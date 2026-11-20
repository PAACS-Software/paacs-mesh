// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::topo::ElementTraits - M4-B.0: the single authoritative semantic
// definition of every element type (directive 7.2/7.3). Topology, audits,
// exports and both volume engines consult THESE tables; no algorithm may
// re-derive orderings locally.
//
// NAMING LOCKS (review requirement - legacy and new names must not become
// competing conventions):
//     mesh::CellType::Line2  <->  ElementType::Segment2
//     mesh::CellType::Line3  <->  ElementType::Segment3
//     "Prism6"               ==   ElementType::Wedge6 (Wedge6 is canonical)
//
// FIRST ROLLOUT includes Tri6 and Segment3 although the directive lists
// them under "future higher-order types": PAACS already GENERATES and
// AUDITS both (order elevation), so excluding them would force the
// element-independent audit to keep hardcoded legacy branches until
// M4-B.9 (review finding). Their traits are frozen to the existing
// MSH/VTK conventions already used by elevation.
//
// HIGH-ORDER NODE INTERPRETATION: node cornerCount+q is the midside of
// edges[q] (MSH/VTK order). Tet10: edges 01,12,20,03,13,23. Tri6/
// Segment3: edges 01,12,20 / 01.
//
// ORIENTATION CONVENTIONS (positive element):
//   Tet4/Tet10: predicates::orient3d(v0,v1,v2,v3) > 0.
//   Tri3/Tri6:  right-hand normal of v0->v1->v2 is the element normal.
//   Quad4:      right-hand normal of the (planar reference) corner loop.
//   Hex8:       corners 0-3 bottom CCW seen from +w, 4-7 above them.
//   Wedge6:     triangle 0-2 bottom, 3-5 above (0->3 etc.).
//   Pyramid5:   quad base 0-3 CCW seen from the apex side, apex 4.
//   All LOCAL FACES are ordered with OUTWARD normals for a positive
//   element (verified by the golden-reference test via centroid tests).
//
// JACOBIAN CONTRACT: for 3D types, the sampled Jacobian determinant at
// every corner (and Gauss points where defined at H0) must be positive
// for a valid element; evaluation routines land with the consumers
// (Quality has Tet10 today; Quad4/Hex8 land at H0). The CONVENTION -
// which map, which sample points, det sign - is fixed here so H0 cannot
// reinvent it.
//
#include <cassert>
#include <stdexcept>
#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/mesh/Mesh.hpp"

namespace paacs::topo {

enum class ElementType : std::uint8_t {
    Segment2 = 0,
    Segment3,
    Tri3,
    Tri6,
    Quad4,
    Tet4,
    Tet10,
    Hex8,
    Wedge6,
    Pyramid5,
    Count_
};

enum class FaceShape : std::uint8_t { None = 0, Tri, Quad };

struct LocalEdge { std::uint8_t a, b; };
struct LocalFace {
    FaceShape shape;
    std::uint8_t n;                   // corner count (3 or 4)
    std::array<std::uint8_t, 4> v;    // corner-local indices, outward order
};
struct RefCoord { double x, y, z; };

struct ElementTraits {
    ElementType type;
    const char* name;
    std::uint8_t dimension;      // 1, 2, or 3
    std::uint8_t cornerCount;
    std::uint8_t nodeCount;      // corners + high-order nodes
    std::uint8_t edgeCount;
    const LocalEdge* edges;
    std::uint8_t faceCount;      // 3D types only; 0 otherwise
    const LocalFace* faces;
    const RefCoord* refCoords;   // nodeCount entries
    std::uint8_t vtkId;          // export mapping (VTK cell type)
    std::uint8_t mshId;          // MSH export element type
};

namespace detail {
// ---- edges -----------------------------------------------------------------
inline constexpr LocalEdge kSegEdges[]  = {{0, 1}};
inline constexpr LocalEdge kTriEdges[]  = {{0, 1}, {1, 2}, {2, 0}};
inline constexpr LocalEdge kQuadEdges[] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
inline constexpr LocalEdge kTetEdges[]  = {{0, 1}, {1, 2}, {2, 0},
                                           {0, 3}, {1, 3}, {2, 3}};
inline constexpr LocalEdge kHexEdges[]  = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                           {4, 5}, {5, 6}, {6, 7}, {7, 4},
                                           {0, 4}, {1, 5}, {2, 6}, {3, 7}};
inline constexpr LocalEdge kWedgeEdges[] = {{0, 1}, {1, 2}, {2, 0},
                                            {3, 4}, {4, 5}, {5, 3},
                                            {0, 3}, {1, 4}, {2, 5}};
inline constexpr LocalEdge kPyrEdges[]  = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                           {0, 4}, {1, 4}, {2, 4}, {3, 4}};
// ---- faces (outward for positive orientation) -------------------------------
inline constexpr LocalFace kTetFaces[] = {
    {FaceShape::Tri, 3, {1, 2, 3, 0}},
    {FaceShape::Tri, 3, {0, 3, 2, 0}},
    {FaceShape::Tri, 3, {0, 1, 3, 0}},
    {FaceShape::Tri, 3, {0, 2, 1, 0}},
};
inline constexpr LocalFace kHexFaces[] = {
    {FaceShape::Quad, 4, {0, 3, 2, 1}},   // bottom (outward -w)
    {FaceShape::Quad, 4, {4, 5, 6, 7}},   // top (outward +w)
    {FaceShape::Quad, 4, {0, 1, 5, 4}},
    {FaceShape::Quad, 4, {1, 2, 6, 5}},
    {FaceShape::Quad, 4, {2, 3, 7, 6}},
    {FaceShape::Quad, 4, {3, 0, 4, 7}},
};
inline constexpr LocalFace kWedgeFaces[] = {
    {FaceShape::Tri, 3, {0, 2, 1, 0}},    // bottom
    {FaceShape::Tri, 3, {3, 4, 5, 0}},    // top
    {FaceShape::Quad, 4, {0, 1, 4, 3}},
    {FaceShape::Quad, 4, {1, 2, 5, 4}},
    {FaceShape::Quad, 4, {2, 0, 3, 5}},
};
inline constexpr LocalFace kPyrFaces[] = {
    {FaceShape::Quad, 4, {0, 3, 2, 1}},   // base
    {FaceShape::Tri, 3, {0, 1, 4, 0}},
    {FaceShape::Tri, 3, {1, 2, 4, 0}},
    {FaceShape::Tri, 3, {2, 3, 4, 0}},
    {FaceShape::Tri, 3, {3, 0, 4, 0}},
};
// ---- reference coordinates ---------------------------------------------------
inline constexpr RefCoord kSeg2Ref[] = {{0,0,0},{1,0,0}};
inline constexpr RefCoord kSeg3Ref[] = {{0,0,0},{1,0,0},{.5,0,0}};
inline constexpr RefCoord kTri3Ref[] = {{0,0,0},{1,0,0},{0,1,0}};
inline constexpr RefCoord kTri6Ref[] = {{0,0,0},{1,0,0},{0,1,0},
                                        {.5,0,0},{.5,.5,0},{0,.5,0}};
inline constexpr RefCoord kQuad4Ref[] = {{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0}};
inline constexpr RefCoord kTet4Ref[] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
inline constexpr RefCoord kTet10Ref[] = {
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},
    {.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
inline constexpr RefCoord kHex8Ref[] = {
    {-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
    {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
inline constexpr RefCoord kWedge6Ref[] = {
    {0,0,-1},{1,0,-1},{0,1,-1},{0,0,1},{1,0,1},{0,1,1}};
inline constexpr RefCoord kPyr5Ref[] = {
    {-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{0,0,1}};

inline constexpr ElementTraits kTraits[] = {
    {ElementType::Segment2, "Segment2", 1, 2, 2, 1, kSegEdges, 0, nullptr,
     kSeg2Ref, 3, 1},
    {ElementType::Segment3, "Segment3", 1, 2, 3, 1, kSegEdges, 0, nullptr,
     kSeg3Ref, 21, 8},
    {ElementType::Tri3, "Tri3", 2, 3, 3, 3, kTriEdges, 0, nullptr,
     kTri3Ref, 5, 2},
    {ElementType::Tri6, "Tri6", 2, 3, 6, 3, kTriEdges, 0, nullptr,
     kTri6Ref, 22, 9},
    {ElementType::Quad4, "Quad4", 2, 4, 4, 4, kQuadEdges, 0, nullptr,
     kQuad4Ref, 9, 3},
    {ElementType::Tet4, "Tet4", 3, 4, 4, 6, kTetEdges, 4, kTetFaces,
     kTet4Ref, 10, 4},
    {ElementType::Tet10, "Tet10", 3, 4, 10, 6, kTetEdges, 4, kTetFaces,
     kTet10Ref, 24, 11},
    {ElementType::Hex8, "Hex8", 3, 8, 8, 12, kHexEdges, 6, kHexFaces,
     kHex8Ref, 12, 5},
    {ElementType::Wedge6, "Wedge6", 3, 6, 6, 9, kWedgeEdges, 5,
     kWedgeFaces, kWedge6Ref, 13, 6},
    {ElementType::Pyramid5, "Pyramid5", 3, 5, 5, 8, kPyrEdges, 5,
     kPyrFaces, kPyr5Ref, 14, 7},
};
} // namespace detail

// Is `t` a valid, in-range element type? (Values can arrive from
// deserialization, corrupted files, C APIs, version mismatches or plugin
// input; an out-of-range enum must never index kTraits out of bounds -
// review found static_cast<ElementType>(255) caused a global-buffer
// overflow.)
inline constexpr bool isValidElementType(ElementType t) {
    return static_cast<std::size_t>(t) <
           static_cast<std::size_t>(ElementType::Count_);
}

// Checked accessor: throws std::out_of_range on an invalid enum. Public
// entry points that receive a possibly-untrusted ElementType must use
// this (or guard with isValidElementType) rather than traits().
inline const ElementTraits& checkedTraits(ElementType t) {
    if (!isValidElementType(t))
        throw std::out_of_range(
            "checkedTraits: ElementType out of range [0, Count_)");
    return detail::kTraits[static_cast<std::size_t>(t)];
}

// Unchecked fast path: callers guarantee t is in range (asserted in debug).
inline const ElementTraits& traits(ElementType t) {
    assert(isValidElementType(t) && "traits(): invalid ElementType");
    return detail::kTraits[static_cast<std::size_t>(t)];
}

// Naming-lock mappings between the legacy public mesh::CellType and
// ElementType. PARTIAL by design: CellType covers only currently
// generated types; generic storage for the rest lands at M4-B.2.
inline bool elementTypeOf(mesh::CellType c, ElementType& out) {
    switch (c) {
        case mesh::CellType::Line2:  out = ElementType::Segment2; return true;
        case mesh::CellType::Line3:  out = ElementType::Segment3; return true;
        case mesh::CellType::Tri3:   out = ElementType::Tri3;     return true;
        case mesh::CellType::Tri6:   out = ElementType::Tri6;     return true;
        case mesh::CellType::Tet4:   out = ElementType::Tet4;     return true;
        case mesh::CellType::Tet10:  out = ElementType::Tet10;    return true;
        case mesh::CellType::Quad4:  out = ElementType::Quad4;    return true;
        case mesh::CellType::Prism6: out = ElementType::Wedge6;   return true;
        case mesh::CellType::Hex8:   out = ElementType::Hex8;     return true;
        case mesh::CellType::Pyramid5: out = ElementType::Pyramid5; return true;
    }
    return false;
}

// Generic connectivity validation per traits (review-corrected: the
// previous form took a bare pointer with no length, so it could not
// detect a short array and could read past it; and it permitted midside
// nodes duplicating corners or each other). Span-based: the caller
// passes the connectivity length; wrong arity is rejected; ALL nodes of
// the element must be pairwise distinct.
// Handle-connectivity validation centralized in traits (review item 9):
// checks arity, that each node satisfies the caller's is-current
// predicate, and pairwise distinctness. Element semantics (arity, and
// later any type-specific rules) stay authoritative in one place. The
// NodeH template parameter avoids a dependency on Storage.hpp; callers
// pass topo::NodeHandle.
template <class NodeH, class IsCurrentFn>
bool validateHandleConnectivity(ElementType t,
                                std::span<const NodeH> nodes,
                                IsCurrentFn&& isCurrent) {
    if (!isValidElementType(t)) return false;
    const ElementTraits& tr = traits(t);
    if (static_cast<int>(nodes.size()) != tr.nodeCount) return false;
    for (const NodeH& nh : nodes)
        if (!isCurrent(nh)) return false;
    for (std::size_t a = 0; a < nodes.size(); ++a)
        for (std::size_t b = a + 1; b < nodes.size(); ++b)
            if (nodes[a] == nodes[b]) return false;
    return true;
}

inline bool validConnectivity(ElementType t, std::span<const Index> conn,
                              Index nodeCount) {
    if (!isValidElementType(t)) return false;
    const ElementTraits& tr = traits(t);
    if (conn.size() != tr.nodeCount) return false;   // arity enforced
    for (const Index v : conn)
        if (v < 0 || v >= nodeCount) return false;
    for (std::size_t a = 0; a < conn.size(); ++a)
        for (std::size_t b = a + 1; b < conn.size(); ++b)
            if (conn[a] == conn[b]) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Orientation-preserving corner automorphisms (the element's rotation
// group), generated ONCE per type from the reference geometry - never
// hand-written (a 24-entry hex table typed by hand is exactly the error
// class the traits test exists to prevent). Used by the canonical
// semantic hash and, later, by face-orientation matching (M4-B.4).
//
// Generation rule: a corner permutation is an automorphism iff it maps
// the local edge set onto itself, the local face set onto itself, and
// preserves orientation - for 3D types the induced linear map on the
// reference corners has positive determinant; for 2D types the corner
// cycle is preserved (cyclic shift); 1D types admit only the identity
// (segments are oriented).
// ---------------------------------------------------------------------------
namespace detail {
std::vector<std::vector<std::uint8_t>> generateRotations(
    const ElementTraits& tr);
} // namespace detail

// Rotation group of the element type (cached; thread-safe via static init).
const std::vector<std::vector<std::uint8_t>>& rotations(ElementType t);

} // namespace paacs::topo
