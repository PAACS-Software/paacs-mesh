// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional adapter; no public model or mesher layout changes.
#include "paacs/geometry/OcctMesher.hpp"
namespace paacs::geo::detail {
NativeCadMeshOutcome buildOcctIndependentParts(const OcctModel&,const NativeCadMeshOptions&,mesh::Mesh&) noexcept;
// Mesh native face-contact components together and isolated bodies through
// their existing routes. Remap native identities and validate the complete
// original assembly before one destination swap. This route applies only when
// a proper contact component and an isolated component are both present.
NativeCadMeshOutcome buildOcctContactComponents(const OcctModel&,
    const NativeCadMeshOptions&, mesh::Mesh& destination,
    const std::string* sourceBefore = nullptr) noexcept;
}
