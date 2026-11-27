// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#pragma once

#include "paacs/geometry/OcctModel.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace paacs::geo {

enum class CadBoundaryEdgeReason : std::uint8_t {
    None, InvalidRequest, InvalidBoundary, MissingAnchor, AmbiguousAnchor,
    ParameterMappingFailed, MissingSurfaceSegment, IncidenceMismatch,
    GeometryError, BudgetExceeded, Cancelled, AllocationFailed, KernelFailure,
    InternalError
};

struct CadBoundaryEdgeOptions {
    std::uint64_t maximumNodes = 1'000'000, maximumTriangles = 2'000'000;
    std::uint64_t maximumSegments = 2'000'000, maximumWork = 20'000'000;
    std::uint64_t maximumWorkingBytes = 256ULL*1024*1024;
    std::uint64_t maximumOutputBytes = 256ULL*1024*1024;
    double maximumSeconds = 60, maximumNodeError = 1e-6, maximumChordError = .05;
    const std::atomic<bool>* cancellation = nullptr;
};

struct CadBoundaryEdgeUse {
    int faceTag = -1;
    std::uint64_t wireIndex = 0, useId = 0;
    CadOrientation orientation = CadOrientation::Forward;
    int sense = 0;
    bool seam = false, degenerated = false;
};

struct CadBoundaryEdgeRecord {
    int edgeTag = -1, startVertexTag = -1, endVertexTag = -1;
    bool periodic = false, closed = false, degenerated = false;
    // Existing boundary node IDs, in increasing native parameter order. A
    // closed curve repeats its first node at the upper parameter endpoint.
    std::vector<Index> naturalNodeOrder, naturalLineOrder;
    std::vector<double> naturalParameters;
    // Source CAD incidence, which can include several owner bodies. Uses are
    // those of the selected boundary body; seams retain both oriented uses.
    std::vector<int> faceTags, bodyTags;
    std::vector<CadBoundaryEdgeUse> uses;
};

struct CadBoundaryEdgeOutcome {
    bool complete = false, destinationUnchanged = true;
    CadBoundaryEdgeReason reason = CadBoundaryEdgeReason::InternalError;
    int solidTag = -1, failureEdgeTag = -1;
    std::uint64_t nodes = 0, triangles = 0, segments = 0, unusedNodes = 0;
    std::uint64_t degeneratedEdges = 0, work = 0;
    double seconds = 0, maximumNodeError = 0, maximumSampledChordError = 0;
    std::vector<CadBoundaryEdgeRecord> edges;
};

// Extract existing CAD-edge segments from an already qualified Tri3 boundary.
// No point is generated, moved, welded or classified by coordinate proximity.
// Declared CAD vertex/edge anchors are checked against their actual source
// curve; consecutive native-parameter anchors must form actual triangle edges
// with matching CAD-face incidence. Inconsistent/ambiguous cases refuse.
//
// The Line2-only result retains ALL boundary nodes/classes and their IDs;
// face-interior points are explicitly unused. Line2 direction is increasing
// native curve parameter, not an individual face's oriented wire use. This is
// an association/export check, not a new surface/volume validity certificate.
// Chord error is sampled at three interior curve parameters, not Hausdorff.
// Limits cover modeled PAACS allocations/work/output and cooperative polling;
// an individual OCCT projection and its workspace are not OS-capped.
// Failure leaves destination unchanged. Aliasing boundary/destination refuses.
CadBoundaryEdgeOutcome exportCadBoundaryEdges(
    const OcctModel&, int solidTag, const mesh::Mesh& qualifiedBoundary,
    const CadBoundaryEdgeOptions&, mesh::Mesh& destination) noexcept;

const char* cadBoundaryEdgeReasonName(CadBoundaryEdgeReason) noexcept;
} // namespace paacs::geo
