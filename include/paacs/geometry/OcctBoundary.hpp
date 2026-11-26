// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/geometry/OcctModel.hpp"
#include "paacs/geometry/CadConstructionBudget.hpp"
#include "paacs/mesh/Mesh.hpp"
#include <atomic>
#include <compare>

namespace paacs::geo {
enum class OcctBoundaryReason : std::uint8_t {
    None, InvalidRequest, UnsupportedTopology, MissingTriangulation,
    NonconformingEdge, GeometryError, InvalidSurface, OpenBoundary,
    Nonmanifold, BudgetExceeded, Cancelled, AllocationFailed, KernelFailure,
    SourceChanged, InternalError
};
struct OcctBoundaryOptions {
    double linearDeflection = .02, angularDeflection = .4;
    // Explicit acceptance bound for sampled physical point-to-trimmed-face
    // distances. Deflection is a kernel request; this bound is independently
    // enforced, with only the owning face's native tolerance added.
    double maximumChordError = .05;
    std::uint64_t maximumNodes = 1'000'000, maximumTriangles = 2'000'000;
    std::uint64_t maximumOutputBytes = 256ULL * 1024 * 1024;
    std::uint64_t maximumWork = 20'000'000;
    const std::atomic<bool>* cancellation = nullptr;
    CadConstructionBudget construction;
};
struct OcctBoundaryNodeKey {
    int dimension = -1, tag = -1;
    // Vertex:0. Edge:increasing native polygon sample ordinal (zero-based).
    // Face:Poly_Triangulation node ordinal (one-based). Keys are valid only
    // within one retained native boundary session, not across imports/meshes.
    std::uint64_t ordinal = 0;
    friend bool operator==(const OcctBoundaryNodeKey&, const OcctBoundaryNodeKey&) = default;
    friend auto operator<=>(const OcctBoundaryNodeKey&, const OcctBoundaryNodeKey&) = default;
};
struct OcctBoundaryOutcome {
    bool published = false, destinationUnchanged = true;
    OcctBoundaryReason reason = OcctBoundaryReason::InternalError;
    int solidTag = -1;
    std::uint64_t nodes = 0, triangles = 0, faces = 0, edgeUses = 0;
    std::uint64_t seamUses = 0, degeneratedUses = 0, collapsedPoleTriangles = 0;
    std::uint64_t closestTrimmedFaceSamples = 0, connectedShells = 0, work = 0;
    double cadVolume = 0, signedSurfaceVolume = 0, relativeVolumeError = 0;
    double maximumStitchError = 0, maximumNodeProjectionError = 0;
    double maximumFacetSampleError = 0, minimumTriangleArea = 0;
    int maximumErrorFace = -1;
    Vec3 maximumErrorSample{}, maximumErrorClosest{};
    bool closed = false, oriented = false, vertexManifold = false;
    bool sourceModelUnmodified = false;
    // Parallel to published boundary nodes. No coordinate-only identity keys.
    std::vector<OcctBoundaryNodeKey> nodeKeys;
};
// Physical closest point on the actual trimmed face, including trim edges and
// vertices. Unlike GFace::projectPoint this accepts an exterior query point.
// Failed queries leave output unchanged; allocation failures propagate.
bool closestPointOnOcctFace(const OcctModel&, int faceTag, const Vec3& point,
                           Vec3& closest);

// Copies native geometry before tessellation, stitches only matching CAD
// vertex/edge identities and curve parameters, preserves native classifications
// and owner face winding, then publishes a Tri3-only boundary by swap.
// Facet centroids and edge midpoints are checked, not a Hausdorff certificate.
// Closed vertex-manifold topology does not certify absence of self-intersection
// or eligibility for volume filling; consumers must audit those separately.
// Output storage and PAACS extraction/check operations are bounded. These
// bounds exclude OCCT copy/tessellation/projection allocations and operations,
// serialization and process RSS. Cancellation is polled by PAACS and forwarded
// through OCCT's cooperative progress interface; it is not a hard time limit.
OcctBoundaryOutcome buildOcctBoundary(const OcctModel&, int solidTag,
                                     const OcctBoundaryOptions&,
                                     mesh::Mesh& destination) noexcept;
const char* occtBoundaryReasonName(OcctBoundaryReason) noexcept;
} // namespace paacs::geo
