// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctModel.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace paacs::geo {
struct OcctPreviewOptions {
    double linearDeflection = 0, angularDeflection = .4;
    std::uint64_t maximumPoints = 250000, maximumTriangles = 500000;
    std::uint64_t maximumOutputBytes = 64ULL * 1024 * 1024;
};
enum class OcctPreviewReason : std::uint8_t {
    None, InvalidRequest, OutputBudget, KernelFailure, AllocationFailed, SourceChanged
};
struct OcctPreview {
    std::vector<Vec3> points;
    std::vector<std::array<std::uint64_t, 3>> triangles;
    std::vector<int> faceTags;
};
struct OcctPreviewOutcome {
    OcctPreviewReason reason = OcctPreviewReason::KernelFailure;
    bool completed = false;
};
const char* occtPreviewReasonName(OcctPreviewReason reason) noexcept;
// Display tessellation of a deep copy. linearDeflection must be positive and
// uses model coordinates (mm for imported STEP). IDs classify triangles to
// imported face references, rather than a body's outward face uses. Vertices
// may be duplicated on shared edges. This is neither a PAACS
// surface mesh nor a geometric-error certificate. Output limits do not bound
// OCCT tessellator allocations or runtime. maximumOutputBytes counts logical
// point/index/tag storage, excluding container overhead and serialization.
// Failure leaves output unchanged. The caller must serialize access to the
// same model with other native operations.
OcctPreviewOutcome previewOcctCad(const OcctModel& model, const OcctPreviewOptions& options,
                                  OcctPreview& output) noexcept;
} // namespace paacs::geo
