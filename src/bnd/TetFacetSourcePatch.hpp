// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/core/Types.hpp"
#include <array>
namespace paacs::bnd::detail {
// One body-local source edit, emitted only after its boundary batch commits.
struct FacetSourcePatch {
    int bodyTag=-1, faceTag=-1;
    Index originalTriangle=-1;
    std::array<Index,3> originalVertices{};
    Index insertedVertex=-1;
};
}
