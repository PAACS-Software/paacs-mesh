// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Quality-guarded ("smart") Laplacian smoothing.
//
// Plain Laplacian smoothing can invert tets near concavities and boundaries;
// here every proposed move is accepted only if the worst quality metric in
// the node's cavity does not decrease and no incident tet inverts. Only
// interior nodes (classification dim == 3) move; boundary fidelity is exact
// by construction.
//
// Metric plumbing is deliberately explicit: the guard metric is the same
// radius-ratio used in reports, so "smoothing improved min quality" is a
// statement about the published number, not a proxy.
//
#include "paacs/mesh/Mesh.hpp"

namespace paacs::algo {

struct SmoothOptions {
    int    maxIterations = 10;
    double relaxation    = 0.6;   // step toward cavity centroid
    double convergenceTol = 1e-4; // stop when max node displacement < tol*h
};

struct SmoothResult {
    int    iterations   = 0;
    Index  acceptedMoves = 0;
    Index  rejectedMoves = 0;
    double minQualityBefore = 0, minQualityAfter = 0;
    double meanQualityBefore = 0, meanQualityAfter = 0;
};

SmoothResult smoothInterior(mesh::Mesh& m,
                            const SmoothOptions& opt = {});

} // namespace paacs::algo
