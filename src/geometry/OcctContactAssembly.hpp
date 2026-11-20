// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Optional adapter implementation; not installed as a public geometry API.
#include "paacs/geometry/OcctMesher.hpp"

namespace paacs::geo::detail {
enum class OcctContactCutPolicy { AllFeatureAxes, PrincipalFeatureAxes };
struct OcctContactAssemblyStrategy {
    OcctContactCutPolicy cuts = OcctContactCutPolicy::AllFeatureAxes;
    unsigned radialCuts = 4, maximumAdaptiveCuts = 4;
    std::uint64_t maximumPieceWork = 5'000'000;
    bool controlDeflectionOnAllSurfaces = false;
};
// Uses one deep-copied native Boolean arrangement and one retained surface
// session for all original bodies. Material overlap, native invalid pieces,
// missing history, unmatched interfaces and exhausted budgets refuse before
// a single destination swap. All original-body CAD checks remain mandatory.
// Logical PAACS admission excludes OCCT transient storage and is not RSS.
NativeCadMeshOutcome buildOcctContactAssembly(const OcctModel&,
    const NativeCadMeshOptions&, mesh::Mesh& destination,
    const std::string* sourceBefore = nullptr) noexcept;
// Additional private strategies keep the original overload and its defaults.
// Every attempt retains the same original geometry/interface publication gates.
NativeCadMeshOutcome buildOcctContactAssembly(const OcctModel&,
    const NativeCadMeshOptions&, const OcctContactAssemblyStrategy&,
    mesh::Mesh& destination, const std::string* sourceBefore = nullptr) noexcept;
// Optional wrapper guard spans the legacy attempt and the assembly fallback.
NativeCadMeshReason snapshotOcctContactSource(const OcctModel&,
    std::uint64_t maximumBytes, std::string& snapshot) noexcept;
}
