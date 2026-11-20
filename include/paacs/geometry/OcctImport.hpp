// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/PolyhedralModel.hpp"
#include <filesystem>

namespace paacs::geo {
enum class CadFileFormat : std::uint8_t { Step, Brep };
enum class CadImportReason : std::uint8_t {
    None,
    InvalidRequest,
    FileUnavailable,
    FileBudget,
    ParseFailed,
    InvalidShape,
    UnsupportedSolidCount,
    UnsupportedSurface,
    UnsupportedCurve,
    TopologyBudget,
    PolyhedralValidationFailed,
    VolumeMismatch,
    AllocationFailed,
    KernelFailure
};
struct CadImportOptions {
    CadFileFormat format = CadFileFormat::Step;
    std::uint64_t maximumFileBytes = 64ULL * 1024 * 1024;
    PolyhedralLimits polyhedral;
};
struct CadImportOutcome {
    std::unique_ptr<Model> model;
    CadImportReason reason = CadImportReason::KernelFailure;
    PolyhedralReason polyhedralReason = PolyhedralReason::None;
    std::uint64_t fileBytes = 0, solids = 0, faces = 0, edges = 0, vertices = 0;
    std::uint64_t curvedFaces = 0, curvedEdges = 0;
    double cadVolume = 0, polyhedralVolume = 0;
    // STEP coordinates are transferred to millimetres explicitly. BREP has no
    // declared length unit; zero means that its coordinates retain unknown units.
    double millimetresPerModelUnit = 0;
    std::vector<std::string> sourceLengthUnits;
    bool imported() const noexcept { return model != nullptr; }
};
// Optional PAACS::occt target, Open CASCADE 8+. Imports a single valid solid
// with planar faces and straight edges, preserving polygonal holes. Geometry
// is baked into the dependency-free immutable backend. Curved faces, curved
// edges, assemblies, and non-solid extra geometry are explicitly refused.
// STEP always transfers to millimetres, independent of OCCT's global unit
// setting. BREP coordinates have unspecified units. targetSize must use the
// resulting model units. No healing is requested by this API.
// The byte/topology/work bounds do not cap allocations inside OCCT's parser.
CadImportOutcome importPlanarCad(const std::filesystem::path &path,
                                 const CadImportOptions &options = {}) noexcept;
const char *cadImportReasonName(CadImportReason reason) noexcept;
} // namespace paacs::geo
