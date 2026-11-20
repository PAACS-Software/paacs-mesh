// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/PrimitiveSizeControls.hpp"
#include "PrimitiveMesherInternal.hpp"
#include "paacs/algo/TetToHexMesher.hpp"
#include <chrono>
#include <cmath>
#include <limits>

namespace paacs::algo {
sizing::ProfileSizeControl primitiveFaceSizeControl(const PrimitiveShape &s, int tag, double size,
                                                    double growth) noexcept {
    if ((s.kind != PrimitiveKind::Cylinder && s.kind != PrimitiveKind::Tube) ||
        (tag != 2 && tag != 3) || !std::isfinite(s.height) || s.height <= 0)
        return {0, 0, -1, growth};
    const double z = tag == 2 ? 0 : s.height;
    return {z, z, size, growth};
}
namespace primitive_detail {
PrimitiveMeshingResult profiled(const PrimitiveShape &s, const PrimitiveMeshingOptions &o,
                                mesh::Mesh &destination) {
    // Fallible staging Mesh construction is caught by the only caller,
    // meshPrimitive, which owns the public noexcept/rollback boundary.
    PrimitiveMeshingResult r;
    if (!o.sizeProfile || !o.sizeProfile->valid() || o.sizeProfile->length() != s.height) {
        r.reason = PrimitiveReason::InvalidRequest;
        return r;
    }
    const auto start = std::chrono::steady_clock::now();
    auto limits = o.profileLimits;
    limits.maximumBytes = std::min(limits.maximumBytes, o.maximumWorkingBytes);
    limits.maximumCells = std::min<std::uint64_t>(limits.maximumCells, 1'000'000);
    if (!limits.cancellation)
        limits.cancellation = o.cancellation;
    auto plan = sizing::planProfileDistribution(*o.sizeProfile, limits);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!plan.planned) {
        switch (plan.reason) {
        case sizing::ProfileReason::Cancelled:
            r.reason = PrimitiveReason::Cancelled;
            break;
        case sizing::ProfileReason::BudgetExceeded:
            r.reason = PrimitiveReason::WorkingBudget;
            break;
        case sizing::ProfileReason::AllocationFailed:
            r.reason = PrimitiveReason::AllocationFailed;
            break;
        case sizing::ProfileReason::NumericallyUnresolved:
            r.reason = PrimitiveReason::NumericallyUnresolved;
            break;
        default:
            r.reason = PrimitiveReason::InvalidRequest;
            break;
        }
    } else {
        auto inner = o;
        inner.sizeProfile = nullptr;
        inner.targetSize = o.sizeProfile->minimumSize();
        inner.axialPositions = plan.positions;
        inner.maximumAxialGrowth = limits.maximumAdjacentGrowth == 0
                                       ? o.sizeProfile->maximumGrowth()
                                       : limits.maximumAdjacentGrowth;
        inner.maximumWorkingBytes -= plan.workingBytes;
        if (o.profileStrategy == PrimitiveProfileStrategy::Sweep) {
            r = meshPrimitive(s, inner, destination);
        } else {
            inner.sizeProfile = o.sizeProfile;
            if (o.elements != PrimitiveElements::HexDominant) {
                r = tetrahedra(s, inner, destination);
            } else {
                auto seed = inner;
                seed.elements = PrimitiveElements::Tet4;
                seed.minimumQuality = 0;
                seed.warningQuality = 0;
                seed.qualityPolicy = PrimitiveQualityPolicy::Warn;
                mesh::Mesh linear;
                r = tetrahedra(s, seed, linear);
                if (r.published) {
                    r.published = false;
                    r.destinationUnchanged = true;
                    if (r.outputBytes > inner.maximumWorkingBytes)
                        r.reason = PrimitiveReason::WorkingBudget;
                    else {
                        TetToHexOptions convert;
                        convert.maximumNodes = o.maximumNodes;
                        convert.maximumCells = o.maximumCells;
                        convert.maximumOutputBytes = o.maximumOutputBytes;
                        convert.maximumWorkingBytes = inner.maximumWorkingBytes - r.outputBytes;
                        convert.auditLimits = o.auditLimits;
                        convert.cancellation = o.cancellation;
                        convert.threads = o.threads;
                        convert.minimumQuality = o.minimumQuality;
                        convert.warningQuality = o.warningQuality;
                        convert.qualityPolicy = o.qualityPolicy;
                        mesh::Mesh hex;
                        const auto conversionStart = std::chrono::steady_clock::now();
                        const auto h = convertTet4ToHex8(linear, convert, hex);
                        r.profileConversionSeconds =
                            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                                          conversionStart)
                                .count();
                        r.plannedWorkingBytes =
                            std::max(r.plannedWorkingBytes, r.outputBytes + h.workingBytes);
                        r.reason = h.reason;
                        r.quality = h.quality;
                        r.minimumQuality = h.minimumQuality;
                        r.qualityWarning = h.quality.lowQualityCells > 0;
                        if (h.published) {
                            r.relativeVolumeError =
                                std::abs(h.volume - r.analyticVolume) / r.analyticVolume;
                            if (r.relativeVolumeError > o.maximumRelativeVolumeError)
                                r.reason = PrimitiveReason::GeometryTolerance;
                            else if ((o.cancellation &&
                                      o.cancellation->load(std::memory_order_relaxed)) ||
                                     (o.auditLimits.cancellation &&
                                      o.auditLimits.cancellation->load(std::memory_order_relaxed)))
                                r.reason = PrimitiveReason::Cancelled;
                            else {
                                destination.swap(hex);
                                r.published = true;
                                r.destinationUnchanged = false;
                                r.nodes = h.nodes;
                                r.hexes = h.hexes;
                                r.tet4s = 0;
                                r.boundaryTriangles = 0;
                                r.boundaryQuads = h.boundaryQuads;
                                r.outputBytes = h.outputBytes;
                                r.meshVolume = h.volume;
                                r.minimumCellVolume = h.minimumCellVolume;
                                r.audit = h.audit;
                                r.convertedHex = true;
                            }
                        }
                    }
                }
            }
        }
        r.plannedWorkingBytes += plan.workingBytes;
    }
    r.profiledSizing = true;
    r.profileWorkingBytes = plan.workingBytes;
    r.profilePlanningSeconds = seconds;
    r.transverseTarget = o.sizeProfile->minimumSize();
    r.minimumAxialSpacing = plan.minimumWidth;
    r.maximumAxialSpacing = plan.maximumWidth;
    r.maximumAxialGrowth = plan.maximumAdjacentGrowth;
    r.plateauCellsRemoved = plan.plateauCellsRemoved;
    return r;
}
} // namespace primitive_detail
} // namespace paacs::algo
