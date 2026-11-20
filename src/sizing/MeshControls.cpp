// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/MeshControls.hpp"

#include <algorithm>
#include <cmath>
#include <new>
#include <utility>

namespace paacs::sizing {
namespace {

bool sameSizeKey(const EntitySizeControl& a,
                 const EntitySizeControl& b) noexcept {
    return a.dimension == b.dimension && a.tag == b.tag;
}

bool samePrism(const PrismLayerControl& a,
               const PrismLayerControl& b) noexcept {
    return a.layers == b.layers &&
        a.firstLayerThickness == b.firstLayerThickness &&
        a.cellToCellGrowth == b.cellToCellGrowth &&
        a.maximumTotalThickness == b.maximumTotalThickness;
}

bool validEngine(VolumeEngineKind value) noexcept {
    switch (value) {
        case VolumeEngineKind::Auto:
        case VolumeEngineKind::Tet4:
        case VolumeEngineKind::MappedHex8:
        case VolumeEngineKind::StructuredPolyhedral:
        case VolumeEngineKind::LayeredPrism6:
            return true;
    }
    return false;
}

const EntitySizeControl* findEntity(
        const std::vector<EntitySizeControl>& rows,
        int dimension, int tag) noexcept {
    const auto it = std::lower_bound(
        rows.begin(), rows.end(), EntitySizeControl{dimension, tag, 0.0},
        [](const EntitySizeControl& a, const EntitySizeControl& b) {
            if (a.dimension != b.dimension)
                return a.dimension < b.dimension;
            return a.tag < b.tag;
        });
    return it != rows.end() && it->dimension == dimension && it->tag == tag
        ? &*it : nullptr;
}

}  // namespace

double CompiledMeshControls::effectiveSize(
        int dimension, int tag, double baseSize) const noexcept {
    if (!std::isfinite(baseSize) || !(baseSize > 0.0)) return baseSize;
    double value = baseSize;
    if (const auto* row = findEntity(entitySizes_, dimension, tag))
        value = std::min(value, row->targetSize);
    return std::clamp(value, minimumSize_, maximumSize_);
}

double CompiledMeshControls::effectiveInterfaceSize(
        int faceTag, int firstSolidTag, int secondSolidTag,
        double baseSize) const noexcept {
    if (!std::isfinite(baseSize) || !(baseSize > 0.0)) return baseSize;
    double value = baseSize;
    if (const auto* face = findEntity(entitySizes_, 2, faceTag))
        value = std::min(value, face->targetSize);
    if (const auto* first = findEntity(entitySizes_, 3, firstSolidTag))
        value = std::min(value, first->targetSize);
    if (const auto* second = findEntity(entitySizes_, 3, secondSolidTag))
        value = std::min(value, second->targetSize);
    return std::clamp(value, minimumSize_, maximumSize_);
}

VolumeEngineKind CompiledMeshControls::engineForBody(
        int solidTag) const noexcept {
    const auto it = std::lower_bound(
        bodyMethods_.begin(), bodyMethods_.end(), solidTag,
        [](const BodyMethodControl& row, int tag) {
            return row.solidTag < tag;
        });
    return it != bodyMethods_.end() && it->solidTag == solidTag
        ? it->engine : VolumeEngineKind::Auto;
}

const PrismLayerControl* CompiledMeshControls::prismLayerForFace(
        int faceTag) const noexcept {
    const auto it = std::lower_bound(
        prismLayers_.begin(), prismLayers_.end(), faceTag,
        [](const PrismLayerControl& row, int tag) {
            return row.faceTag < tag;
        });
    return it != prismLayers_.end() && it->faceTag == faceTag ? &*it : nullptr;
}

MeshControlCompileOutcome compileMeshControls(
        const MeshControlSpecification& specification) noexcept {
    MeshControlCompileOutcome out;
    if (!std::isfinite(specification.minimumSize) ||
        !std::isfinite(specification.maximumSize) ||
        specification.minimumSize < 0.0 ||
        !(specification.maximumSize > 0.0) ||
        specification.minimumSize > specification.maximumSize) {
        out.reason = MeshControlCompileReason::InvalidGlobalBounds;
        return out;
    }
    if (!std::isfinite(specification.maximumGrowthRate) ||
        specification.maximumGrowthRate < 1.0) {
        out.reason = MeshControlCompileReason::InvalidGrowthRate;
        return out;
    }

    try {
        CompiledMeshControls compiled;
        compiled.minimumSize_ = specification.minimumSize;
        compiled.maximumSize_ = specification.maximumSize;
        compiled.maximumGrowthRate_ = specification.maximumGrowthRate;
        compiled.entitySizes_ = specification.entitySizes;
        compiled.bodyMethods_ = specification.bodyMethods;
        compiled.prismLayers_ = specification.prismLayers;

        for (const auto& row : compiled.entitySizes_) {
            if (row.dimension < 0 || row.dimension > 3 || row.tag < 0) {
                out.reason = MeshControlCompileReason::InvalidEntity;
                return out;
            }
            if (!std::isfinite(row.targetSize) || !(row.targetSize > 0.0)) {
                out.reason = MeshControlCompileReason::InvalidSize;
                return out;
            }
        }
        std::sort(compiled.entitySizes_.begin(), compiled.entitySizes_.end(),
            [](const EntitySizeControl& a, const EntitySizeControl& b) {
                if (a.dimension != b.dimension)
                    return a.dimension < b.dimension;
                return a.tag < b.tag;
            });
        std::vector<EntitySizeControl> uniqueSizes;
        uniqueSizes.reserve(compiled.entitySizes_.size());
        for (const auto& row : compiled.entitySizes_) {
            if (!uniqueSizes.empty() && sameSizeKey(uniqueSizes.back(), row)) {
                if (uniqueSizes.back().targetSize != row.targetSize) {
                    out.reason =
                        MeshControlCompileReason::ConflictingDuplicateSize;
                    return out;
                }
                continue;
            }
            uniqueSizes.push_back(row);
        }
        compiled.entitySizes_.swap(uniqueSizes);

        for (const auto& row : compiled.bodyMethods_)
            if (row.solidTag < 0 || !validEngine(row.engine)) {
                out.reason = MeshControlCompileReason::InvalidBodyMethod;
                return out;
            }
        std::sort(compiled.bodyMethods_.begin(), compiled.bodyMethods_.end(),
            [](const BodyMethodControl& a, const BodyMethodControl& b) {
                return a.solidTag < b.solidTag;
            });
        std::vector<BodyMethodControl> uniqueMethods;
        uniqueMethods.reserve(compiled.bodyMethods_.size());
        for (const auto& row : compiled.bodyMethods_) {
            if (!uniqueMethods.empty() &&
                uniqueMethods.back().solidTag == row.solidTag) {
                if (uniqueMethods.back().engine != row.engine) {
                    out.reason = MeshControlCompileReason::
                        ConflictingDuplicateBodyMethod;
                    return out;
                }
                continue;
            }
            uniqueMethods.push_back(row);
        }
        compiled.bodyMethods_.swap(uniqueMethods);

        for (const auto& row : compiled.prismLayers_) {
            if (row.faceTag < 0 || row.layers == 0 ||
                !std::isfinite(row.firstLayerThickness) ||
                !(row.firstLayerThickness > 0.0) ||
                !std::isfinite(row.cellToCellGrowth) ||
                row.cellToCellGrowth < 1.0 ||
                row.cellToCellGrowth > specification.maximumGrowthRate ||
                !std::isfinite(row.maximumTotalThickness) ||
                !(row.maximumTotalThickness > 0.0)) {
                out.reason = MeshControlCompileReason::InvalidPrismLayer;
                return out;
            }
        }
        std::sort(compiled.prismLayers_.begin(), compiled.prismLayers_.end(),
            [](const PrismLayerControl& a, const PrismLayerControl& b) {
                return a.faceTag < b.faceTag;
            });
        std::vector<PrismLayerControl> uniqueLayers;
        uniqueLayers.reserve(compiled.prismLayers_.size());
        for (const auto& row : compiled.prismLayers_) {
            if (!uniqueLayers.empty() &&
                uniqueLayers.back().faceTag == row.faceTag) {
                if (!samePrism(uniqueLayers.back(), row)) {
                    out.reason = MeshControlCompileReason::
                        ConflictingDuplicatePrismLayer;
                    return out;
                }
                continue;
            }
            uniqueLayers.push_back(row);
        }
        compiled.prismLayers_.swap(uniqueLayers);

        out.controls = std::move(compiled);
        out.compiled = true;
        out.reason = MeshControlCompileReason::None;
        return out;
    } catch (const std::bad_alloc&) {
        out.reason = MeshControlCompileReason::AllocationFailed;
    } catch (...) {
        out.reason = MeshControlCompileReason::InternalError;
    }
    return out;
}

}  // namespace paacs::sizing
