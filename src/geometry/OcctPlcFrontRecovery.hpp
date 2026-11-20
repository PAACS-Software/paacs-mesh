// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "OcctRecoveryMesher.hpp"
namespace paacs::geo::detail {
// Complete one small, immutable native piece PLC. The partition assembler
// retains original-CAD provenance, shared-interface and final geometry gates.
OcctPlcRecoveryOutcome recoverOcctPlcFrontVolume(const mesh::Mesh& boundary,int pieceTag,
    const OcctModel& pieces,const OcctModel& original,int originalBody,
    const OcctStarOptions&,const OcctPartitionRecoveryOptions&,mesh::Mesh& destination) noexcept;
OcctPlcRecoveryOutcome recoverOcctPlcLiveVolume(const mesh::Mesh& boundary,int pieceTag,
    const OcctStarOptions&,const OcctPartitionRecoveryOptions&,mesh::Mesh& destination) noexcept;
}
