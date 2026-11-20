// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <set>
//
// Topological mesh optimization: 2-3 and 3-2 bistellar flips, interleaved
// with quality-guarded smoothing. This is what lifts the MINIMUM quality -
// vertex smoothing alone cannot fix a sliver whose nodes all sit on the
// boundary; removing it requires changing connectivity.
//
// Every flip is gated on (a) all replacement tets positively oriented,
// (b) total cavity volume conserved to round-off (rejects geometrically
// invalid flips in one check), and (c) the cavity's worst radius-ratio
// strictly improving.
//
#include "paacs/algo/Smoothing.hpp"
#include "paacs/mesh/Mesh.hpp"

namespace paacs::algo {

struct OptimizeOptions {
    int    maxPasses = 5;
    double minGain   = 1e-4;   // required min-quality improvement per flip
    bool   smoothBetweenPasses = true;
    SmoothOptions smooth;

    // Targeted sliver vertex perturbation ("exudation-lite"): for tets below
    // this quality that own at least one interior node, deterministically
    // search small moves of that node that raise the cavity's worst quality.
    // Catches sliver families that flips (connectivity) and boundary peeling
    // (all-boundary tets) cannot reach. 0 disables.
    double perturbBelow = 0.10;

    // Final cleanup: peel remaining sub-threshold slivers that expose a hull
    // face (near-zero-volume caps on curved walls whose base triangle is too
    // skinny for any vertex move to fix). Removes negligible volume; proper
    // fix is frontal surface point placement + boundary recovery (roadmap).
    // 0 disables.
    double peelBelow = 0.05;

    // Peeling is restricted to tets with at least one node classified on
    // one of these (curved) faces - see algo::curvedFaceTags(model). An
    // EMPTY set disables hull peeling entirely: on planar solids the
    // low-quality all-boundary tets can be the legitimate thin domain
    // (PM-ALG-002 guard).
    std::set<int> peelableFaceTags;
};

struct OptimizeResult {
    int    passes  = 0;
    Index  flips23 = 0, flips32 = 0;
    Index  perturbations = 0;
    Index  peeledSlivers = 0;
    double minQualityBefore = 0, minQualityAfter = 0;
    double meanQualityBefore = 0, meanQualityAfter = 0;
};

OptimizeResult optimizeTets(mesh::Mesh& m, const OptimizeOptions& opt = {});

} // namespace paacs::algo
