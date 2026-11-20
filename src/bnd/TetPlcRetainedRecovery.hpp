// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Internal continuation only. Refusal never publishes a complete volume. The
// optional progress output contains checked Tet4 cells, with the source prefix
// intact, so a subsequent bounded strategy can continue completed repairs.
#include "paacs/bnd/TetPlcCavityRecovery.hpp"
namespace paacs::bnd {
TetPlcCavityOutcome recoverTetPlcCavitiesRetaining(const mesh::Mesh&,const mesh::Mesh&,
    const TetPlcCavityOptions&,mesh::Mesh&,mesh::Mesh* retainedProgress) noexcept;
}
