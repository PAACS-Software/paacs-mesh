// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-adapter utility. The caller retains native CAD authority.
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"

namespace paacs::geo { class OcctModel; }
namespace paacs::geo::detail {
enum class OcctStarReason : std::uint8_t {
    None, InvalidBoundary, UnsupportedTopology, BoundaryContact,
    KernelNotFound, WorkBudget, OutputBudget, Cancelled, InvalidFan,
    AllocationFailed, InternalError, NativeQueryFailed
};
struct OcctStarOptions {
    std::uint64_t maximumWork = 5'000'000;
    std::uint64_t maximumBoundaryNodes = 100'000;
    std::uint64_t maximumBoundaryTriangles = 100'000;
    unsigned maximumProposalSweeps = 256;
    std::uint64_t maximumNodes = 1'000'000;
    std::uint64_t maximumCells = 5'000'000;
    std::uint64_t maximumOutputBytes = 512ULL * 1024 * 1024;
    // Conservative preflight of requested array/container payload, excluding
    // allocator metadata and native classifier/kernel storage. This is an
    // admission estimate, not a process limit.
    std::uint64_t maximumWorkingBytes = 512ULL * 1024 * 1024;
    mesh::LinearVolumeAuditLimits auditLimits{};
    const std::atomic<bool>* cancellation = nullptr;
    // Explicit candidate admission against an original native body. Both the
    // proposed kernel and every resulting fan Tet4 centroid must classify IN
    // or ON before a proposal is accepted. This supplements the exact PLC
    // certificate; full CAD boundary/volume/source audits remain caller gates.
    // The model must outlive this call. Null preserves the previous search.
    const OcctModel* nativeAdmissionModel = nullptr;
    int nativeAdmissionBodyTag = -1;
};
struct OcctStarOutcome {
    bool published = false;
    bool destinationUnchanged = true;
    bool embedded = false;
    bool kernelFeasible = false;
    OcctStarReason reason = OcctStarReason::InternalError;
    const char* diagnostic = "internal-error";
    Vec3 kernel{};
    double boundaryVolume = 0;
    double minimumQuality = 0;
    std::uint64_t workUsed = 0;
    std::uint64_t exactBoundaryPairs = 0;
    std::uint64_t kernelTrials = 0;
    // Already included in workUsed; failed proposals retain their query count.
    std::uint64_t nativeCandidateQueries = 0;
};
struct OcctBoundaryCertificate {
    bool certified = false;
    bool closedOrientedVertexManifold = false;
    bool embedded = false;
    OcctStarReason reason = OcctStarReason::InternalError;
    const char* diagnostic = "internal-error";
    double signedVolume = 0;
    std::uint64_t nodes = 0, triangles = 0, connectedShells = 0;
    std::int64_t eulerCharacteristic = 0;
    std::uint64_t workUsed = 0, exactBoundaryPairs = 0;
};
// Independently certifies exact topological closure/vertex links and embedded
// geometry of an outward Tri3-only boundary. Disconnected shells, cavities and
// genus are permitted; shell containment and CAD ownership remain caller gates.
// The complete oriented boundary must have finite positive signed volume.
OcctBoundaryCertificate certifyOcctBoundarySurface(const mesh::Mesh& boundary,
                                                   const OcctStarOptions&) noexcept;
// Accepts a single outward Tri3-only shell with native node/face classifications.
// The shell must be an embedded connected orientable closed vertex-manifold
// sphere. A numerical proposal becomes a kernel only when every frozen exact
// orient3d predicate places it strictly inside its outward triangle halfspace.
// The resulting fan retains every input triangle, node and classification.
// Work counts explicit source visits, contact candidates, numerical proposals
// and exact side tests; it does not bound allocations or predicate internals.
// CAD error, CAD inclusion, requested size and body contact checks remain the
// caller's responsibility. KernelNotFound is a bounded-search refusal, not an
// infeasibility proof. Publication uses one nonthrowing swap after all audits.
OcctStarOutcome buildOcctStarVolume(const mesh::Mesh& boundary, int solidTag,
                                   const OcctStarOptions& options,
                                   mesh::Mesh& destination) noexcept;
} // namespace paacs::geo::detail
