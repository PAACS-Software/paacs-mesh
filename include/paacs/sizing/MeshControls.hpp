// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/core/Types.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace paacs::sizing {

// Public method selection is stored by solid identity.  Auto is conservative:
// the orchestration layer must prove that a mapped/swept descriptor exists;
// otherwise it routes to the accepted Tet4 engine rather than guessing.
enum class VolumeEngineKind : std::uint8_t {
    Auto,
    Tet4,
    MappedHex8,
    StructuredPolyhedral,
    LayeredPrism6
};

struct EntitySizeControl {
    int dimension = -1;
    int tag = -1;
    double targetSize = 0.0;
};

struct BodyMethodControl {
    int solidTag = -1;
    VolumeEngineKind engine = VolumeEngineKind::Auto;
};

// Geometry-independent layer intent.  The selected face still needs a valid
// sweep/extrusion descriptor before the prism engine can consume this control.
struct PrismLayerControl {
    int faceTag = -1;
    std::uint32_t layers = 0;
    double firstLayerThickness = 0.0;
    double cellToCellGrowth = 1.0;
    double maximumTotalThickness = 0.0;
};

struct MeshControlSpecification {
    double minimumSize = 0.0;
    double maximumSize = 1.0e300;
    double maximumGrowthRate = 1.3;
    std::vector<EntitySizeControl> entitySizes;
    std::vector<BodyMethodControl> bodyMethods;
    std::vector<PrismLayerControl> prismLayers;
};

enum class MeshControlCompileReason : std::uint8_t {
    None,
    InvalidGlobalBounds,
    InvalidGrowthRate,
    InvalidEntity,
    InvalidSize,
    ConflictingDuplicateSize,
    InvalidBodyMethod,
    ConflictingDuplicateBodyMethod,
    InvalidPrismLayer,
    ConflictingDuplicatePrismLayer,
    AllocationFailed,
    InternalError
};

struct MeshControlCompileOutcome;

class CompiledMeshControls {
public:
    double minimumSize() const noexcept { return minimumSize_; }
    double maximumSize() const noexcept { return maximumSize_; }
    double maximumGrowthRate() const noexcept { return maximumGrowthRate_; }

    // Applies one entity-specific cap and the global min/max envelope.  The
    // lookup is deterministic binary search over a compact sorted vector.
    double effectiveSize(int dimension, int tag,
                         double baseSize) const noexcept;

    // Shared-interface planning uses the minimum of the face control and both
    // owning body controls, then applies the global envelope once.
    double effectiveInterfaceSize(int faceTag, int firstSolidTag,
                                  int secondSolidTag,
                                  double baseSize) const noexcept;

    VolumeEngineKind engineForBody(int solidTag) const noexcept;
    const PrismLayerControl* prismLayerForFace(int faceTag) const noexcept;

    std::size_t entitySizeCount() const noexcept { return entitySizes_.size(); }
    std::size_t bodyMethodCount() const noexcept { return bodyMethods_.size(); }
    std::size_t prismLayerCount() const noexcept { return prismLayers_.size(); }

private:
    friend MeshControlCompileOutcome compileMeshControls(
        const MeshControlSpecification&) noexcept;

    double minimumSize_ = 0.0;
    double maximumSize_ = 1.0e300;
    double maximumGrowthRate_ = 1.3;
    std::vector<EntitySizeControl> entitySizes_;
    std::vector<BodyMethodControl> bodyMethods_;
    std::vector<PrismLayerControl> prismLayers_;
};

struct MeshControlCompileOutcome {
    bool compiled = false;
    MeshControlCompileReason reason = MeshControlCompileReason::InternalError;
    CompiledMeshControls controls{};
};

MeshControlCompileOutcome compileMeshControls(
    const MeshControlSpecification& specification) noexcept;

}  // namespace paacs::sizing
