// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private producer-to-adapter receipt for CK-2/T2-C. The public candidate
// remains self-contained; this receipt is accepted only after a complete
// value/identity/topology validation and otherwise falls back to the original
// adapter path.

#include "paacs/algo/SurfaceMesherAdapter.hpp"

namespace paacs::algo::detail {

struct PreparedSurfaceCandidate {
    CandidateProviderStamp providerStamp{};
    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    bnd::GeometryProviderStamp geometryStamp{};
    bnd::OwnershipStamp workspaceStamp{};
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::NodeHandle> candidateNodes;
    std::vector<topo::CellHandle> candidateCells;
    bnd::ConstraintSet constraints;
    bnd::SteinerSegmentsReport recovery;
    std::vector<std::size_t> frameBoundary;
    std::size_t exteriorSeedCell = 0;
    bool finalized = false;
};

}  // namespace paacs::algo::detail
