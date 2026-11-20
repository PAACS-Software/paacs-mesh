// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/algo/ProtectedFeatureRefinement.hpp"

#include "paacs/predicates/Predicates.hpp"
#include "CheckedSizing.hpp"
#include "ProductionSurfaceMesherDetail.hpp"
#include "T2CSurfaceCache.hpp"
#include "TetCandidateDetail.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <numeric>
#include <optional>

namespace paacs::algo {
namespace {

using predicates::Sign;

#if defined(PAACS_TOPO_TESTING)
class BoundaryConformingTestingWallScope {
public:
    explicit BoundaryConformingTestingWallScope(
        boundary_conforming_detail::BoundaryConformingWallPhase phase)
        : phase_(phase),start_(std::chrono::steady_clock::now()) {}
    ~BoundaryConformingTestingWallScope() {
        auto& trace = boundary_conforming_detail::
            boundaryConformingWallTraceForTesting();
        if (!trace.enabled) return;
        const auto elapsed = std::chrono::duration_cast<
            std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-start_).count();
        trace.nanoseconds[static_cast<std::size_t>(phase_)] +=
            static_cast<std::uint64_t>(elapsed);
    }
private:
    boundary_conforming_detail::BoundaryConformingWallPhase phase_;
    std::chrono::steady_clock::time_point start_;
};
#endif

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

bool samePoint(const Vec3& a, const Vec3& b) noexcept {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

std::array<Index, 3> canonicalTriple(std::array<Index, 3> v) noexcept {
    std::sort(v.begin(), v.end());
    return v;
}

bool cyclicEqual(const std::array<Index, 3>& a,
                 const std::array<Index, 3>& b) noexcept {
    for (std::size_t shift = 0; shift < 3; ++shift) {
        if (a[0] == b[shift] && a[1] == b[(shift + 1) % 3] &&
            a[2] == b[(shift + 2) % 3])
            return true;
    }
    return false;
}

bool noncollinear(const Vec3& a, const Vec3& b, const Vec3& c) {
    return predicates::orient2d({a.x, a.y}, {b.x, b.y}, {c.x, c.y}) !=
               Sign::Zero ||
           predicates::orient2d({a.y, a.z}, {b.y, b.z}, {c.y, c.z}) !=
               Sign::Zero ||
           predicates::orient2d({a.z, a.x}, {b.z, b.x}, {c.z, c.x}) !=
               Sign::Zero;
}

BoundaryConformingV3Status statusForSurface(SurfaceSetStatus status) noexcept {
    switch (status) {
        case SurfaceSetStatus::Rejected:
            return BoundaryConformingV3Status::Rejected;
        case SurfaceSetStatus::Refused:
            return BoundaryConformingV3Status::Refused;
        case SurfaceSetStatus::Failed:
            return BoundaryConformingV3Status::Failed;
        case SurfaceSetStatus::Published:
            return BoundaryConformingV3Status::Published;
    }
    return BoundaryConformingV3Status::Failed;
}

BoundaryConformingV3Reason reasonForSurface(SurfaceSetStatus status) noexcept {
    switch (status) {
        case SurfaceSetStatus::Rejected:
            return BoundaryConformingV3Reason::SurfaceConstructionRejected;
        case SurfaceSetStatus::Refused:
            return BoundaryConformingV3Reason::SurfaceConstructionRefused;
        case SurfaceSetStatus::Failed:
            return BoundaryConformingV3Reason::SurfaceConstructionFailed;
        case SurfaceSetStatus::Published:
            return BoundaryConformingV3Reason::None;
    }
    return BoundaryConformingV3Reason::InternalError;
}

ClosedArtifactBuildOutcome artifactFailure(
    BoundaryConformingV3Status status, BoundaryConformingV3Reason reason,
    BoundaryConformingV3Phase phase,
    VolumeProviderCheckpoint checkpoint = VolumeProviderCheckpoint::None) {
    ClosedArtifactBuildOutcome out;
    out.status = status;
    out.failure.reason = reason;
    out.failure.phase = phase;
    out.failure.providerCheckpoint = checkpoint;
    return out;
}

RecoveryEnvelopeBuildOutcome envelopeFailure(
    BoundaryConformingV3Status status, BoundaryConformingV3Reason reason,
    VolumeProviderCheckpoint checkpoint = VolumeProviderCheckpoint::None) {
    RecoveryEnvelopeBuildOutcome out;
    out.status = status;
    out.failure.reason = reason;
    out.failure.phase = BoundaryConformingV3Phase::RecoveryEnvelope;
    out.failure.providerCheckpoint = checkpoint;
    return out;
}

struct PendingFacet {
    FacePersistentId face = 0;
    std::array<Index, 3> directed{};
    mesh::Classification classification{};
};

struct ArtifactAabb {
    Vec3 lo{},hi{};
};

ArtifactAabb mergeArtifactBounds(const ArtifactAabb& a,
                                 const ArtifactAabb& b) noexcept {
    return {{std::min(a.lo.x,b.lo.x),std::min(a.lo.y,b.lo.y),
             std::min(a.lo.z,b.lo.z)},
            {std::max(a.hi.x,b.hi.x),std::max(a.hi.y,b.hi.y),
             std::max(a.hi.z,b.hi.z)}};
}

bool artifactBoundsOverlap(const ArtifactAabb& a,
                           const ArtifactAabb& b) noexcept {
    return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x &&
           a.lo.y <= b.hi.y && b.lo.y <= a.hi.y &&
           a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
}

struct ArtifactLeaf {
    ArtifactAabb bounds{};
    std::size_t facet = 0;
};

struct ArtifactTreeNode {
    ArtifactAabb bounds{};
    std::int32_t left = -1,right = -1;
    std::size_t facet = 0;
    bool leaf = false;
};

std::int32_t buildArtifactTreeNode(
    std::vector<ArtifactTreeNode>& nodes,std::vector<ArtifactLeaf>& leaves,
    std::size_t begin,std::size_t end,std::size_t depth,
    ArtifactPreflightMetrics& metrics) {
    const auto result = static_cast<std::int32_t>(nodes.size());
    nodes.push_back({});
    metrics.treeDepth = std::max(metrics.treeDepth,
        static_cast<std::uint64_t>(depth+1));
    if (end-begin == 1) {
        auto& node = nodes[static_cast<std::size_t>(result)];
        node.bounds = leaves[begin].bounds;
        node.facet = leaves[begin].facet;
        node.leaf = true;
        return result;
    }
    std::array<double,3> centreLo{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double,3> centreHi{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
    for (std::size_t i = begin; i < end; ++i) {
        const auto& bounds = leaves[i].bounds;
        const std::array<double,3> centre{
            std::midpoint(bounds.lo.x,bounds.hi.x),
            std::midpoint(bounds.lo.y,bounds.hi.y),
            std::midpoint(bounds.lo.z,bounds.hi.z)};
        for (std::size_t candidate = 0; candidate < 3; ++candidate) {
            centreLo[candidate] =
                std::min(centreLo[candidate],centre[candidate]);
            centreHi[candidate] =
                std::max(centreHi[candidate],centre[candidate]);
        }
        ++metrics.axisSelectionVisits;
    }
    std::size_t axis = 0;
    for (std::size_t candidate = 1; candidate < 3; ++candidate)
        if (centreHi[candidate]-centreLo[candidate] >
            centreHi[axis]-centreLo[axis])
            axis = candidate;
    const auto coordinate = [axis](const ArtifactAabb& value) noexcept {
        if (axis == 0) return std::midpoint(value.lo.x,value.hi.x);
        if (axis == 1) return std::midpoint(value.lo.y,value.hi.y);
        return std::midpoint(value.lo.z,value.hi.z);
    };
    std::stable_sort(leaves.begin()+static_cast<std::ptrdiff_t>(begin),
                     leaves.begin()+static_cast<std::ptrdiff_t>(end),
                     [&](const ArtifactLeaf& a,const ArtifactLeaf& b) {
        ++metrics.sortComparisons;
        const double av = coordinate(a.bounds),bv = coordinate(b.bounds);
        return av != bv ? av < bv : a.facet < b.facet;
    });
    const std::size_t middle = begin+(end-begin)/2;
    const auto left = buildArtifactTreeNode(
        nodes,leaves,begin,middle,depth+1,metrics);
    const auto right = buildArtifactTreeNode(
        nodes,leaves,middle,end,depth+1,metrics);
    auto& node = nodes[static_cast<std::size_t>(result)];
    node.left = left;
    node.right = right;
    node.bounds = mergeArtifactBounds(
        nodes[static_cast<std::size_t>(left)].bounds,
        nodes[static_cast<std::size_t>(right)].bounds);
    return result;
}

bool auditArtifactContacts(const bnd::ClosedTri3Artifact& artifact,
                           const mesh::Mesh& surface,
                           ArtifactPreflightMetrics& metrics) {
    if (surface.nodes().size() >= bnd::kInvalidTetFacetVertex) return false;
    std::vector<bnd::TetFacetVertex> vertices;
    vertices.reserve(surface.nodes().size());
    for (std::size_t i = 0;i < surface.nodes().size();++i)
        vertices.push_back({surface.nodes()[i],i+1});
    std::vector<ArtifactLeaf> leaves;
    leaves.reserve(artifact.facets.size());
    for (std::size_t fi = 0;fi < artifact.facets.size();++fi) {
        const auto& ids = artifact.facets[fi].directedGlobalNodes;
        ArtifactAabb bounds{surface.node(ids[0]),surface.node(ids[0])};
        for (std::size_t i = 1;i < 3;++i) {
            const Vec3& p = surface.node(ids[i]);
            bounds.lo.x = std::min(bounds.lo.x,p.x);
            bounds.lo.y = std::min(bounds.lo.y,p.y);
            bounds.lo.z = std::min(bounds.lo.z,p.z);
            bounds.hi.x = std::max(bounds.hi.x,p.x);
            bounds.hi.y = std::max(bounds.hi.y,p.y);
            bounds.hi.z = std::max(bounds.hi.z,p.z);
        }
        leaves.push_back({bounds,fi});
    }
    std::vector<ArtifactTreeNode> tree;
    if (!leaves.empty()) {
        tree.reserve(leaves.size()*2-1);
        (void)buildArtifactTreeNode(
            tree,leaves,0,leaves.size(),0,metrics);
    }
    metrics.treeNodes = tree.size();
    std::vector<std::int32_t> stack;
    if (!tree.empty()) stack.reserve(metrics.treeDepth);
    metrics.scratchHighWaterBytes =
        tree.capacity()*sizeof(ArtifactTreeNode)+
        leaves.capacity()*sizeof(ArtifactLeaf)+
        vertices.capacity()*sizeof(bnd::TetFacetVertex)+
        stack.capacity()*sizeof(std::int32_t);
    bnd::TetFacetComplexView view;
    view.vertices = vertices;
    for (std::size_t li = 0;li < leaves.size();++li) {
        const std::size_t fi = leaves[li].facet;
        stack.clear();
        stack.push_back(0);
        while (!stack.empty()) {
            const auto ni = stack.back();
            stack.pop_back();
            const auto& node = tree[static_cast<std::size_t>(ni)];
            ++metrics.nodePairVisits;
            ++metrics.aabbTests;
            if (!artifactBoundsOverlap(leaves[li].bounds,node.bounds)) continue;
            if (!node.leaf) {
                stack.push_back(node.right);
                stack.push_back(node.left);
                continue;
            }
            const std::size_t fj = node.facet;
            if (fj <= fi) continue;
            const auto& a = artifact.facets[fi].directedGlobalNodes;
            const auto& b = artifact.facets[fj].directedGlobalNodes;
            std::array<bnd::TetFacetVertexId,3> ai{},bi{};
            for (std::size_t k = 0;k < 3;++k) {
                ai[k] = static_cast<bnd::TetFacetVertexId>(a[k]);
                bi[k] = static_cast<bnd::TetFacetVertexId>(b[k]);
            }
            ++metrics.facetFacetPairs;
            metrics.edgeEdgePairs += 9;
            metrics.edgeFacetPairs += 6;
            ++metrics.exactFacetFacetTests;
            metrics.exactEdgeEdgeTests += 9;
            metrics.exactEdgeFacetTests += 6;
            const std::array<Vec3,3> ap{surface.node(a[0]),surface.node(a[1]),
                                        surface.node(a[2])};
            const std::array<Vec3,3> bp{surface.node(b[0]),surface.node(b[1]),
                                        surface.node(b[2])};
            const auto contact = bnd::tet_facet_detail::
                auditIdentityAwareClosedContact(
                    view,ai,bnd::tet_facet_detail::kTriangleEdges,
                    bnd::tet_facet_detail::kTriangleFaces,
                    bi,bnd::tet_facet_detail::kTriangleEdges,
                    bnd::tet_facet_detail::kTriangleFaces);
            if (contact.illegal() ||
                bnd::tet_facet_detail::trianglesOpenIntersect(ap,bp)) {
                ++metrics.rejectedContacts;
                return false;
            }
            if (contact.legalIdentityContact) ++metrics.legalSharedContacts;
        }
    }
    metrics.allPairsFallbacks = 0;
    return true;
}

BoundaryConformingV3Reason auditRegionShells(
    const bnd::ClosedTri3Artifact& artifact) {
    struct DirectedUse {
        bnd::VolumeFacetKey facet{};
        Index from = kInvalidIndex;
        Index to = kInvalidIndex;
    };
    for (int region : artifact.regionTags) {
        std::map<bnd::VolumeArtifactEdgeKey, std::vector<DirectedUse>> uses;
        std::map<bnd::VolumeFacetKey, std::set<bnd::VolumeFacetKey>> adjacent;
        std::set<bnd::VolumeFacetKey> regionFacets;
        for (const auto& facet : artifact.facets) {
            auto owner = std::find_if(
                facet.owners.begin(), facet.owners.end(),
                [region](const auto& use) { return use.regionTag == region; });
            if (owner == facet.owners.end()) continue;
            regionFacets.insert(facet.key);
            auto v = facet.directedGlobalNodes;
            if (owner->directedSense < 0) std::swap(v[1], v[2]);
            for (std::size_t i = 0; i < 3; ++i) {
                const Index a = v[i], b = v[(i + 1) % 3];
                uses[bnd::VolumeArtifactEdgeKey::of(a, b)].push_back(
                    {facet.key, a, b});
            }
        }
        if (regionFacets.empty())
            return BoundaryConformingV3Reason::NonClosedShellComponent;
        for (const auto& [key, edgeUses] : uses) {
            (void)key;
            if (edgeUses.size() < 2)
                return BoundaryConformingV3Reason::NonClosedShellComponent;
            if (edgeUses.size() > 2)
                return BoundaryConformingV3Reason::NonManifoldShellComponent;
            if (edgeUses[0].from != edgeUses[1].to ||
                edgeUses[0].to != edgeUses[1].from)
                return BoundaryConformingV3Reason::NonOrientableShellComponent;
            adjacent[edgeUses[0].facet].insert(edgeUses[1].facet);
            adjacent[edgeUses[1].facet].insert(edgeUses[0].facet);
        }
        std::set<bnd::VolumeFacetKey> visited;
        for (const auto& seed : regionFacets) {
            if (!visited.insert(seed).second) continue;
            std::queue<bnd::VolumeFacetKey> q;
            q.push(seed);
            while (!q.empty()) {
                const auto current = q.front();
                q.pop();
                for (const auto& next : adjacent[current])
                    if (visited.insert(next).second) q.push(next);
            }
        }
        if (visited != regionFacets)
            return BoundaryConformingV3Reason::NonClosedShellComponent;

        // Closed two-manifold node-link: incident triangles around each node
        // form one cycle (every opposite edge endpoint has degree two).
        std::map<Index, std::map<Index, std::set<Index>>> links;
        for (const auto& facet : artifact.facets) {
            bool owns = false;
            for (const auto& use : facet.owners)
                owns = owns || use.regionTag == region;
            if (!owns) continue;
            const auto& v = facet.directedGlobalNodes;
            for (std::size_t i = 0; i < 3; ++i) {
                const Index center = v[i];
                const Index a = v[(i + 1) % 3], b = v[(i + 2) % 3];
                links[center][a].insert(b);
                links[center][b].insert(a);
            }
        }
        for (const auto& [center, link] : links) {
            (void)center;
            if (link.empty())
                return BoundaryConformingV3Reason::NonClosedShellComponent;
            for (const auto& [vertex, neighbours] : link) {
                (void)vertex;
                if (neighbours.size() < 2)
                    return BoundaryConformingV3Reason::NonClosedShellComponent;
                if (neighbours.size() > 2)
                    return BoundaryConformingV3Reason::
                        NonManifoldShellComponent;
            }
            std::set<Index> linkVisited;
            std::queue<Index> q;
            q.push(link.begin()->first);
            linkVisited.insert(link.begin()->first);
            while (!q.empty()) {
                const Index current = q.front();
                q.pop();
                for (Index next : link.at(current))
                    if (linkVisited.insert(next).second) q.push(next);
            }
            if (linkVisited.size() != link.size())
                return BoundaryConformingV3Reason::
                    NonManifoldShellComponent;
        }
    }
    return BoundaryConformingV3Reason::None;
}

// A common binary power-of-two rescale leaves every orient3d sign unchanged
// while keeping the determinant away from overflow/underflow in the envelope
// range discriminators.  std::scalbn performs an exact exponent adjustment
// for every representable result; the qualified predicate itself is unchanged.
Sign envelopeOrient3d(const Vec3& a,const Vec3& b,
                      const Vec3& c,const Vec3& d) noexcept {
    const double magnitude = std::max(
        {std::fabs(a.x),std::fabs(a.y),std::fabs(a.z),
         std::fabs(b.x),std::fabs(b.y),std::fabs(b.z),
         std::fabs(c.x),std::fabs(c.y),std::fabs(c.z),
         std::fabs(d.x),std::fabs(d.y),std::fabs(d.z)});
    if (magnitude == 0.0)
        return Sign::Zero;
    int exponent = 0;
    (void)std::frexp(magnitude,&exponent);
    const auto scaled = [exponent](const Vec3& p) noexcept {
        return Vec3{std::scalbn(p.x,-exponent),
                    std::scalbn(p.y,-exponent),
                    std::scalbn(p.z,-exponent)};
    };
    return predicates::orient3d(
        scaled(a),scaled(b),scaled(c),scaled(d));
}

bool validateEnvelopeGeometry(const TetRecoveryEnvelope& envelope,
                              std::span<const Vec3> construction) {
    if (envelope.guardCoordinates.size() < 4 ||
        envelope.outwardHullFaces.size() < 4)
        return false;
    std::set<std::array<std::size_t, 3>> uniqueFaces;
    std::map<std::array<std::size_t, 2>, std::vector<int>> edgeSenses;
    std::map<std::size_t, std::map<std::size_t, std::set<std::size_t>>> links;
    std::vector<std::set<std::size_t>> adjacency(
        envelope.outwardHullFaces.size());
    std::map<std::array<std::size_t, 2>, std::vector<std::size_t>> owners;
    for (const Vec3& p : envelope.guardCoordinates)
        if (!finite(p)) return false;
    for (std::size_t i = 0; i < envelope.guardCoordinates.size(); ++i)
        for (std::size_t j = i + 1; j < envelope.guardCoordinates.size(); ++j)
            if (samePoint(envelope.guardCoordinates[i],
                          envelope.guardCoordinates[j]))
                return false;

    for (std::size_t fi = 0; fi < envelope.outwardHullFaces.size(); ++fi) {
        const auto face = envelope.outwardHullFaces[fi];
        if (face[0] >= envelope.guardCoordinates.size() ||
            face[1] >= envelope.guardCoordinates.size() ||
            face[2] >= envelope.guardCoordinates.size() ||
            face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            return false;
        auto key = face;
        std::sort(key.begin(), key.end());
        if (!uniqueFaces.insert(key).second) return false;
        const Vec3& a = envelope.guardCoordinates[face[0]];
        const Vec3& b = envelope.guardCoordinates[face[1]];
        const Vec3& c = envelope.guardCoordinates[face[2]];
        std::size_t witness = envelope.guardCoordinates.size();
        for (std::size_t q = 0; q < envelope.guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            if (envelopeOrient3d(a, b, c,
                                 envelope.guardCoordinates[q]) !=
                Sign::Zero) {
                witness = q;
                break;
            }
        }
        if (witness == envelope.guardCoordinates.size() ||
            envelopeOrient3d(a, b, c,
                             envelope.guardCoordinates[witness]) !=
                Sign::Negative)
            return false;
        for (std::size_t q = 0; q < envelope.guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            const Sign side = envelopeOrient3d(
                a, b, c, envelope.guardCoordinates[q]);
            if (side == Sign::Positive) return false;
        }
        for (const Vec3& p : construction)
            if (!finite(p) ||
                envelopeOrient3d(a, b, c, p) != Sign::Negative)
                return false;
        for (std::size_t i = 0; i < 3; ++i) {
            const std::size_t from = face[i], to = face[(i + 1) % 3];
            const std::array<std::size_t, 2> edge =
                from < to ? std::array<std::size_t, 2>{from, to}
                          : std::array<std::size_t, 2>{to, from};
            edgeSenses[edge].push_back(from < to ? 1 : -1);
            owners[edge].push_back(fi);
            links[face[i]][face[(i + 1) % 3]].insert(face[(i + 2) % 3]);
            links[face[i]][face[(i + 2) % 3]].insert(face[(i + 1) % 3]);
        }
    }
    for (const auto& [edge, senses] : edgeSenses) {
        (void)edge;
        if (senses.size() != 2 || senses[0] == senses[1]) return false;
    }
    for (const auto& [edge, fs] : owners) {
        (void)edge;
        if (fs.size() != 2) return false;
        adjacency[fs[0]].insert(fs[1]);
        adjacency[fs[1]].insert(fs[0]);
    }
    std::set<std::size_t> visited{0};
    std::queue<std::size_t> q;
    q.push(0);
    while (!q.empty()) {
        const auto f = q.front();
        q.pop();
        for (const auto next : adjacency[f])
            if (visited.insert(next).second) q.push(next);
    }
    if (visited.size() != envelope.outwardHullFaces.size()) return false;
    const std::int64_t chi =
        static_cast<std::int64_t>(envelope.guardCoordinates.size()) -
        static_cast<std::int64_t>(edgeSenses.size()) +
        static_cast<std::int64_t>(envelope.outwardHullFaces.size());
    if (chi != 2) return false;
    for (const auto& [center, link] : links) {
        (void)center;
        for (const auto& [vertex, neighbours] : link) {
            (void)vertex;
            if (neighbours.size() != 2) return false;
        }
        std::set<std::size_t> linkVisited{link.begin()->first};
        std::queue<std::size_t> linkQueue;
        linkQueue.push(link.begin()->first);
        while (!linkQueue.empty()) {
            const auto v = linkQueue.front();
            linkQueue.pop();
            for (const auto next : link.at(v))
                if (linkVisited.insert(next).second) linkQueue.push(next);
        }
        if (linkVisited.size() != link.size()) return false;
    }
    return true;
}

inline constexpr std::array<std::array<std::size_t,3>,12>
    kDerivedBoxFaces{{
        {{0,4,2}},{{4,6,2}},{{1,3,5}},{{3,7,5}},
        {{0,1,4}},{{1,5,4}},{{2,6,3}},{{6,7,3}},
        {{0,2,1}},{{2,3,1}},{{4,5,6}},{{5,7,6}}}};

struct ConstructorFramePoint {
    Vec3 coordinate{};
    std::size_t identity = 0;
};

struct ConstructorFrameScratch {
    std::vector<Vec3> transformedConstruction;
    TetRecoveryEnvelope transformedEnvelope;
    std::vector<ConstructorFramePoint> ordered;

    void prepare(std::size_t constructionCount,
                 std::span<const std::array<std::size_t,3>> faces,
                 std::size_t guardCount) {
        transformedConstruction.resize(constructionCount);
        transformedEnvelope.guardCoordinates.resize(guardCount);
        transformedEnvelope.outwardHullFaces.assign(faces.begin(),faces.end());
        ordered.reserve(constructionCount+guardCount);
    }
};

bool validateEnvelopeSideRules(
    const TetRecoveryEnvelope& envelope,
    std::span<const Vec3> construction) noexcept {
    for (const auto& face : envelope.outwardHullFaces) {
        if (face[0] >= envelope.guardCoordinates.size() ||
            face[1] >= envelope.guardCoordinates.size() ||
            face[2] >= envelope.guardCoordinates.size() ||
            face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            return false;
        const Vec3& a = envelope.guardCoordinates[face[0]];
        const Vec3& b = envelope.guardCoordinates[face[1]];
        const Vec3& c = envelope.guardCoordinates[face[2]];
        std::size_t witness = envelope.guardCoordinates.size();
        for (std::size_t q = 0; q < envelope.guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            if (envelopeOrient3d(
                    a,b,c,envelope.guardCoordinates[q]) != Sign::Zero) {
                witness = q;
                break;
            }
        }
        if (witness == envelope.guardCoordinates.size() ||
            envelopeOrient3d(
                a,b,c,envelope.guardCoordinates[witness]) != Sign::Negative)
            return false;
        for (std::size_t q = 0; q < envelope.guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            if (envelopeOrient3d(
                    a,b,c,envelope.guardCoordinates[q]) == Sign::Positive)
                return false;
        }
        for (const Vec3& point : construction)
            if (envelopeOrient3d(a,b,c,point) != Sign::Negative)
                return false;
    }
    return true;
}

bool validateConstructorFrame(
    const TetRecoveryEnvelope& envelope,
    std::span<const Vec3> construction,
    ConstructorFrameScratch& scratch) noexcept {
    if (construction.empty() || envelope.guardCoordinates.empty() ||
        scratch.transformedConstruction.size() != construction.size() ||
        scratch.transformedEnvelope.guardCoordinates.size() !=
            envelope.guardCoordinates.size())
        return false;
    Vec3 lo = construction.front(), hi = construction.front();
    auto include = [&](const Vec3& point) {
        if (!finite(point)) return false;
        lo.x = std::min(lo.x,point.x); hi.x = std::max(hi.x,point.x);
        lo.y = std::min(lo.y,point.y); hi.y = std::max(hi.y,point.y);
        lo.z = std::min(lo.z,point.z); hi.z = std::max(hi.z,point.z);
        return true;
    };
    for (const Vec3& point : construction)
        if (!include(point)) return false;
    for (const Vec3& point : envelope.guardCoordinates)
        if (!include(point)) return false;

    // This intentionally mirrors frozen delaunay3D preprocessing expression
    // for expression; even algebraically equivalent overflow-avoiding forms
    // would not prove the actual binary-double frame.
    const Vec3 mid = (lo+hi)*0.5;
    const double span = std::max(
        {hi.x-lo.x,hi.y-lo.y,hi.z-lo.z,1e-300});
    const double inv = 1.0/span;
    if (!finite(mid) || !std::isfinite(span) || span <= 0.0 ||
        !std::isfinite(inv))
        return false;

    scratch.ordered.clear();
    std::size_t identity = 0;
    for (std::size_t i = 0; i < construction.size(); ++i) {
        const Vec3 transformed = (construction[i]-mid)*inv;
        if (!finite(transformed)) return false;
        scratch.transformedConstruction[i] = transformed;
        scratch.ordered.push_back({transformed,identity++});
    }
    for (std::size_t i = 0; i < envelope.guardCoordinates.size(); ++i) {
        const Vec3 transformed =
            (envelope.guardCoordinates[i]-mid)*inv;
        if (!finite(transformed)) return false;
        scratch.transformedEnvelope.guardCoordinates[i] = transformed;
        scratch.ordered.push_back({transformed,identity++});
    }
    std::sort(scratch.ordered.begin(),scratch.ordered.end(),
              [](const auto& a,const auto& b) {
        return std::tie(a.coordinate.x,a.coordinate.y,a.coordinate.z,
                        a.identity) <
               std::tie(b.coordinate.x,b.coordinate.y,b.coordinate.z,
                        b.identity);
    });
    for (std::size_t i = 1; i < scratch.ordered.size(); ++i)
        if (samePoint(scratch.ordered[i-1].coordinate,
                      scratch.ordered[i].coordinate))
            return false;
    return validateEnvelopeSideRules(
        scratch.transformedEnvelope,scratch.transformedConstruction);
}

void setDerivedBox(TetRecoveryEnvelope& envelope,
                   const Vec3& guardLo,const Vec3& guardHi) {
    envelope.guardCoordinates.resize(8);
    for (std::size_t i = 0; i < 8; ++i)
        envelope.guardCoordinates[i] =
            {(i&1U) ? guardHi.x : guardLo.x,
             (i&2U) ? guardHi.y : guardLo.y,
             (i&4U) ? guardHi.z : guardLo.z};
}

bool sameDerivedBounds(const Vec3& aLo,const Vec3& aHi,
                       const Vec3& bLo,const Vec3& bHi) noexcept {
    return samePoint(aLo,bLo) && samePoint(aHi,bHi);
}

bool constructionBounds(std::span<const Vec3> construction,
                        Vec3& lo,Vec3& hi,double& referenceSpan) noexcept {
    if (construction.empty() || !finite(construction.front())) return false;
    lo = construction.front();
    hi = construction.front();
    double magnitude = std::max(
        {std::fabs(lo.x),std::fabs(lo.y),std::fabs(lo.z),1.0});
    for (const Vec3& point : construction) {
        if (!finite(point)) return false;
        lo.x = std::min(lo.x,point.x); hi.x = std::max(hi.x,point.x);
        lo.y = std::min(lo.y,point.y); hi.y = std::max(hi.y,point.y);
        lo.z = std::min(lo.z,point.z); hi.z = std::max(hi.z,point.z);
        magnitude = std::max(
            {magnitude,std::fabs(point.x),std::fabs(point.y),
             std::fabs(point.z)});
    }
    referenceSpan = std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
    if (!std::isfinite(referenceSpan)) return false;
    if (referenceSpan == 0.0) referenceSpan = magnitude;
    return std::isfinite(referenceSpan) && referenceSpan > 0.0;
}

}  // namespace

BoundaryConformingV3ReportEnvelope::BoundaryConformingV3ReportEnvelope() =
    default;
BoundaryConformingV3ReportEnvelope::~BoundaryConformingV3ReportEnvelope() =
    default;
BoundaryConformingV3Result::BoundaryConformingV3Result() = default;
BoundaryConformingV3Result::~BoundaryConformingV3Result() = default;
BoundaryConformingV3Outcome::BoundaryConformingV3Outcome() = default;
BoundaryConformingV3Outcome::~BoundaryConformingV3Outcome() = default;
BoundaryConformingV3Outcome::BoundaryConformingV3Outcome(
    BoundaryConformingV3Outcome&&) noexcept = default;
BoundaryConformingV3Outcome& BoundaryConformingV3Outcome::operator=(
    BoundaryConformingV3Outcome&&) noexcept = default;

namespace boundary_conforming_detail {

#if defined(PAACS_TOPO_TESTING)
BoundaryConformingWallTrace& boundaryConformingWallTraceForTesting() noexcept {
    static thread_local BoundaryConformingWallTrace trace;
    return trace;
}
#endif

bool validateTetEnvelopeConstructorFrameForTesting(
    const TetRecoveryEnvelope& envelope,
    std::span<const Vec3> constructionCoordinates) noexcept {
    try {
        ConstructorFrameScratch scratch;
        scratch.prepare(constructionCoordinates.size(),
                        envelope.outwardHullFaces,
                        envelope.guardCoordinates.size());
        return validateConstructorFrame(
            envelope,constructionCoordinates,scratch);
    } catch (...) {
        return false;
    }
}

int derivedTetEnvelopeRungForTesting(
    std::span<const Vec3> constructionCoordinates,
    const TetRecoveryEnvelope& envelope) noexcept {
    if (envelope.guardCoordinates.size() != 8) return -2;
    Vec3 lo{},hi{};
    double referenceSpan = 0.0;
    if (!constructionBounds(
            constructionCoordinates,lo,hi,referenceSpan))
        return -2;
    const double inf = std::numeric_limits<double>::infinity();
    Vec3 candidateLo{std::nextafter(lo.x,-inf),
                     std::nextafter(lo.y,-inf),
                     std::nextafter(lo.z,-inf)};
    Vec3 candidateHi{std::nextafter(hi.x,inf),
                     std::nextafter(hi.y,inf),
                     std::nextafter(hi.z,inf)};
    if (samePoint(envelope.guardCoordinates[0],candidateLo) &&
        samePoint(envelope.guardCoordinates[7],candidateHi))
        return -1;
    Vec3 previousLo = candidateLo,previousHi = candidateHi;
    for (int rung = 0; rung <= 52; ++rung) {
        const double delta = std::ldexp(referenceSpan,rung-52);
        if (!std::isfinite(delta) || delta <= 0.0) continue;
        const Vec3 rawLo{lo.x-delta,lo.y-delta,lo.z-delta};
        const Vec3 rawHi{hi.x+delta,hi.y+delta,hi.z+delta};
        candidateLo = {std::nextafter(rawLo.x,-inf),
                       std::nextafter(rawLo.y,-inf),
                       std::nextafter(rawLo.z,-inf)};
        candidateHi = {std::nextafter(rawHi.x,inf),
                       std::nextafter(rawHi.y,inf),
                       std::nextafter(rawHi.z,inf)};
        if (sameDerivedBounds(candidateLo,candidateHi,
                              previousLo,previousHi))
            continue;
        previousLo = candidateLo;
        previousHi = candidateHi;
        if (samePoint(envelope.guardCoordinates[0],candidateLo) &&
            samePoint(envelope.guardCoordinates[7],candidateHi))
            return rung;
    }
    return -2;
}

bool auditClosedTri3ContactsForTesting(
    const bnd::ClosedTri3Artifact& artifact,const mesh::Mesh& surface,
    ArtifactPreflightMetrics& metrics) noexcept {
    try {
        metrics = {};
        metrics.edgeEntries = artifact.edges.size();
        metrics.facetEntries = artifact.facets.size();
        const bool valid = auditArtifactContacts(artifact,surface,metrics);
        metrics.runComplete = valid;
        return valid;
    } catch (...) {
        metrics.runComplete = false;
        return false;
    }
}

BoundaryConformingV3Reason auditClosedTri3RegionShellsForTesting(
    const bnd::ClosedTri3Artifact& artifact) noexcept {
    try {
        return auditRegionShells(artifact);
    } catch (...) {
        return BoundaryConformingV3Reason::InternalError;
    }
}

}  // namespace boundary_conforming_detail

ClosedArtifactBuildOutcome buildClosedTri3Artifact(
    const SurfaceSetResult& surface, const mesh::Mesh& publishedSurface,
    ArtifactStamp expectedSurfaceSetStamp,
    VolumeBoundaryProviderStamp expectedProviderStamp,
    const VolumeBoundaryOwnershipProvider& provider) noexcept {
    try {
        if (expectedSurfaceSetStamp.persistentId == 0 ||
            surface.setStamp != expectedSurfaceSetStamp)
            return artifactFailure(BoundaryConformingV3Status::Rejected,
                                   BoundaryConformingV3Reason::StaleSurfaceArtifact,
                                   BoundaryConformingV3Phase::Preflight);
        auto artifact = std::make_unique<bnd::ClosedTri3Artifact>();
        artifact->surfaceSetStamp = surface.setStamp;
        std::vector<const SurfaceFaceResult*> faces;
        faces.reserve(surface.faces.size());
        for (const auto& face : surface.faces) {
            if (!face) return artifactFailure(
                BoundaryConformingV3Status::Rejected,
                BoundaryConformingV3Reason::SurfaceCellMappingIncomplete,
                BoundaryConformingV3Phase::SurfaceMapping);
            faces.push_back(face.get());
        }
        std::sort(faces.begin(), faces.end(), [](const auto* a, const auto* b) {
            return a->faceStamp.persistentId < b->faceStamp.persistentId;
        });
        for (std::size_t i = 1; i < faces.size(); ++i)
            if (faces[i - 1]->faceStamp.persistentId ==
                faces[i]->faceStamp.persistentId)
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::DuplicateFacetIdentity,
                    BoundaryConformingV3Phase::SurfaceMapping);

        // Global surface cells have dense public indices.  The previous
        // validation used an ordered set for mapped-cell membership and then
        // rescanned every published cell once per face to prove that no Tri3
        // of the face classification was left unmapped.  That turns a linear
        // artifact audit into O(face_count * surface_cell_count).  Build the
        // classification census once and use a dense byte table for exact
        // alias detection; the proof is equivalent because each mapped cell
        // is already payload-checked below.
        std::vector<std::uint8_t> mappedCells(
            static_cast<std::size_t>(publishedSurface.cellCount()),0);
        std::map<std::pair<int,int>,std::size_t> triCellsByClass;
        for(Index ci=0;ci<publishedSurface.cellCount();++ci) {
            const auto& cell=publishedSurface.cell(ci);
            if(cell.type==mesh::CellType::Tri3)
                ++triCellsByClass[{cell.onEntity.dim,cell.onEntity.tag}];
        }
        std::map<Index, std::set<FacePersistentId>> nodeFaces;
        std::map<FacePersistentId, std::vector<PendingFacet>> byFace;
        std::set<Index> sharedWitness;
        std::map<Index,mesh::Classification> sharedClassByGlobal;
        for (const auto& [key, binding] : surface.sharedLedger) {
            if (binding.existingGlobalIndex < 0 ||
                binding.existingGlobalIndex >= publishedSurface.nodeCount())
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceNodeMappingIncomplete,
                    BoundaryConformingV3Phase::SurfaceMapping);
            if (!(binding.key == key) ||
                !sameClass(binding.classification,
                           publishedSurface.nodeClass(
                               binding.existingGlobalIndex)))
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceNodePayloadMismatch,
                    BoundaryConformingV3Phase::SurfaceMapping);
            const auto [existing,inserted] = sharedClassByGlobal.emplace(
                binding.existingGlobalIndex,binding.classification);
            if (!inserted &&
                (binding.classification.dim != 0 ||
                 !sameClass(existing->second,binding.classification)))
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceNodeMappingAlias,
                    BoundaryConformingV3Phase::SurfaceMapping);
            sharedWitness.insert(binding.existingGlobalIndex);
        }

        for (const SurfaceFaceResult* face : faces) {
            if (face->publishedGlobalByLocalNode.size() !=
                face->localNodes.size())
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceNodeMappingIncomplete,
                    BoundaryConformingV3Phase::SurfaceMapping);
            if (face->publishedGlobalCellByLocalCell.size() !=
                face->localCells.size())
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceCellMappingIncomplete,
                    BoundaryConformingV3Phase::SurfaceMapping);
            std::set<Index> localNodes;
            for (std::size_t i = 0; i < face->localNodes.size(); ++i) {
                const Index global = face->publishedGlobalByLocalNode[i];
                if (global < 0 || global >= publishedSurface.nodeCount())
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceNodeMappingIncomplete,
                        BoundaryConformingV3Phase::SurfaceMapping);
                if (!localNodes.insert(global).second)
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceNodeMappingAlias,
                        BoundaryConformingV3Phase::SurfaceMapping);
                if (!samePoint(face->localNodes[i].coord,
                               publishedSurface.node(global)) ||
                    !sameClass(face->localNodes[i].cls,
                               publishedSurface.nodeClass(global)))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceNodePayloadMismatch,
                        BoundaryConformingV3Phase::SurfaceMapping);
                nodeFaces[global].insert(face->faceStamp.persistentId);
            }
            std::set<std::array<Index, 3>> faceTriangles;
            for (std::size_t ci = 0; ci < face->localCells.size(); ++ci) {
                const auto& local = face->localCells[ci];
                std::array<Index, 3> directed{};
                for (std::size_t n = 0; n < 3; ++n) {
                    if (local.nodes[n] >= face->localNodes.size())
                        return artifactFailure(
                            BoundaryConformingV3Status::Rejected,
                            BoundaryConformingV3Reason::SurfaceCellMappingIncomplete,
                            BoundaryConformingV3Phase::SurfaceMapping);
                    directed[n] =
                        face->publishedGlobalByLocalNode[local.nodes[n]];
                }
                const Index globalCell =
                    face->publishedGlobalCellByLocalCell[ci];
                if (globalCell < 0 || globalCell >= publishedSurface.cellCount() ||
                    mappedCells[static_cast<std::size_t>(globalCell)] != 0)
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceCellMappingAlias,
                        BoundaryConformingV3Phase::SurfaceMapping);
                mappedCells[static_cast<std::size_t>(globalCell)] = 1;
                const auto& cell = publishedSurface.cell(globalCell);
                if (cell.type != mesh::CellType::Tri3 ||
                    !sameClass(cell.onEntity, face->faceClass) ||
                    !sameClass(local.classification, face->faceClass))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceCellPayloadMismatch,
                        BoundaryConformingV3Phase::SurfaceMapping);
                const Index* cn = publishedSurface.cellNodes(globalCell);
                const std::array<Index, 3> actual{cn[0], cn[1], cn[2]};
                if (!cyclicEqual(directed, actual))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::SurfaceCellPayloadMismatch,
                        BoundaryConformingV3Phase::SurfaceMapping);
                const auto canonical = canonicalTriple(directed);
                if (!faceTriangles.insert(canonical).second ||
                    !noncollinear(publishedSurface.node(directed[0]),
                                  publishedSurface.node(directed[1]),
                                  publishedSurface.node(directed[2])))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::DuplicateFacetGeometry,
                        BoundaryConformingV3Phase::SurfaceMapping);
                byFace[face->faceStamp.persistentId].push_back(
                    {face->faceStamp.persistentId, directed,
                     face->faceClass});
            }
            const auto count=triCellsByClass.find(
                {face->faceClass.dim,face->faceClass.tag});
            if(count==triCellsByClass.end()||
               count->second!=face->localCells.size())
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceCellMappingIncomplete,
                    BoundaryConformingV3Phase::SurfaceMapping);
        }
        for (const auto& [global, owners] : nodeFaces)
            if (owners.size() > 1 && sharedWitness.count(global) == 0)
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::SurfaceNodeMappingAlias,
                    BoundaryConformingV3Phase::SurfaceMapping);

        if (expectedProviderStamp.persistentId == 0 ||
            provider.stamp() != expectedProviderStamp)
            return artifactFailure(
                BoundaryConformingV3Status::Rejected,
                BoundaryConformingV3Reason::StaleOwnershipProvider,
                BoundaryConformingV3Phase::Preflight,
                VolumeProviderCheckpoint::Entry);

        if (provider.stamp() != expectedProviderStamp)
            return artifactFailure(
                BoundaryConformingV3Status::Failed,
                BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                BoundaryConformingV3Phase::Ownership,
                VolumeProviderCheckpoint::PreOwnershipTable);
        std::map<FacePersistentId, VolumeFaceOwnership> ownership;
        for (const SurfaceFaceResult* face : faces) {
            VolumeFaceOwnership row;
            try {
                row = provider.ownership(face->faceStamp.persistentId);
            } catch (const std::bad_alloc&) {
                return artifactFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::AllocationFailed,
                    BoundaryConformingV3Phase::Ownership);
            } catch (...) {
                return artifactFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderException,
                    BoundaryConformingV3Phase::Ownership);
            }
            if (row.providerStamp != expectedProviderStamp ||
                row.facePersistentId != face->faceStamp.persistentId)
                return artifactFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderContractViolation,
                    BoundaryConformingV3Phase::Ownership);
            if (row.availability ==
                VolumeProviderAvailability::MissingAuthoritativeOrientation) {
                if (!row.owners.empty()) return artifactFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderContractViolation,
                    BoundaryConformingV3Phase::Ownership);
                return artifactFailure(
                    BoundaryConformingV3Status::Unsupported,
                    BoundaryConformingV3Reason::MissingOwnershipOrientation,
                    BoundaryConformingV3Phase::Ownership);
            }
            if (row.availability != VolumeProviderAvailability::Available)
                return artifactFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderContractViolation,
                    BoundaryConformingV3Phase::Ownership);
            if (row.owners.empty() || row.owners.size() > 2)
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::InvalidOwnerMultiplicity,
                    BoundaryConformingV3Phase::Ownership);
            std::sort(row.owners.begin(), row.owners.end());
            for (const auto& owner : row.owners)
                if (owner.regionTag < 0 ||
                    (owner.directedSense != -1 && owner.directedSense != 1))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::InvalidOwnerSense,
                        BoundaryConformingV3Phase::Ownership);
            if (row.owners.size() == 2 &&
                (row.owners[0].regionTag == row.owners[1].regionTag ||
                 row.owners[0].directedSense == row.owners[1].directedSense))
                return artifactFailure(
                    BoundaryConformingV3Status::Rejected,
                    BoundaryConformingV3Reason::InvalidOwnerSide,
                    BoundaryConformingV3Phase::Ownership);
            ownership.emplace(row.facePersistentId, std::move(row));
        }
        if (provider.stamp() != expectedProviderStamp)
            return artifactFailure(
                BoundaryConformingV3Status::Failed,
                BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                BoundaryConformingV3Phase::Ownership,
                VolumeProviderCheckpoint::PostOwnershipTable);

        std::set<std::array<Index, 3>> completeGeometry;
        std::set<int> regions;
        for (auto& [faceId, pending] : byFace) {
            std::sort(pending.begin(), pending.end(), [](const auto& a,
                                                        const auto& b) {
                return canonicalTriple(a.directed) < canonicalTriple(b.directed);
            });
            for (std::size_t ordinal = 0; ordinal < pending.size(); ++ordinal) {
                if (!completeGeometry.insert(
                        canonicalTriple(pending[ordinal].directed)).second)
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::DuplicateFacetGeometry,
                        BoundaryConformingV3Phase::SurfaceMapping);
                bnd::VolumeFacet facet;
                facet.key = {faceId, ordinal};
                facet.directedGlobalNodes = pending[ordinal].directed;
                facet.classification = pending[ordinal].classification;
                facet.owners = ownership.at(faceId).owners;
                for (const auto& owner : facet.owners)
                    regions.insert(owner.regionTag);
                artifact->facets.push_back(std::move(facet));
            }
        }
        std::sort(artifact->facets.begin(), artifact->facets.end(),
                  [](const auto& a, const auto& b) { return a.key < b.key; });
        artifact->regionTags.assign(regions.begin(), regions.end());

        std::map<bnd::VolumeArtifactEdgeKey, mesh::Classification> curveEdges;
        for (const SurfaceFaceResult* face : faces) {
            for (const auto& edge : face->persistentEdges) {
                if (!edge.publishedGlobalEndpoints) continue;
                const auto key = bnd::VolumeArtifactEdgeKey::of(
                    (*edge.publishedGlobalEndpoints)[0],
                    (*edge.publishedGlobalEndpoints)[1]);
                const auto [it, inserted] =
                    curveEdges.emplace(key, edge.classification);
                if (!inserted && !sameClass(it->second, edge.classification))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::ArtifactEdgeClassificationConflict,
                        BoundaryConformingV3Phase::IndexedPreflight);
            }
        }
        std::map<bnd::VolumeArtifactEdgeKey, bnd::VolumeArtifactEdge> edges;
        for (const auto& facet : artifact->facets) {
            for (std::size_t i = 0; i < 3; ++i) {
                const auto key = bnd::VolumeArtifactEdgeKey::of(
                    facet.directedGlobalNodes[i],
                    facet.directedGlobalNodes[(i + 1) % 3]);
                auto [it, inserted] = edges.emplace(
                    key, bnd::VolumeArtifactEdge{key, facet.classification, {}});
                const auto curve = curveEdges.find(key);
                const mesh::Classification wanted =
                    curve == curveEdges.end() ? facet.classification
                                              : curve->second;
                if (!inserted && !sameClass(it->second.classification, wanted))
                    return artifactFailure(
                        BoundaryConformingV3Status::Rejected,
                        BoundaryConformingV3Reason::ArtifactEdgeClassificationConflict,
                        BoundaryConformingV3Phase::IndexedPreflight);
                it->second.classification = wanted;
                it->second.incidentFacets.push_back(facet.key);
            }
        }
        for (auto& [key, edge] : edges) {
            (void)key;
            std::sort(edge.incidentFacets.begin(), edge.incidentFacets.end());
            artifact->edges.push_back(std::move(edge));
        }
        ArtifactPreflightMetrics preflightMetrics;
        preflightMetrics.edgeEntries = artifact->edges.size();
        preflightMetrics.facetEntries = artifact->facets.size();
        if (!auditArtifactContacts(
                *artifact,publishedSurface,preflightMetrics)) {
            auto failed = artifactFailure(
                BoundaryConformingV3Status::Rejected,
                BoundaryConformingV3Reason::ArtifactContactConflict,
                BoundaryConformingV3Phase::IndexedPreflight);
            failed.metrics = preflightMetrics;
            return failed;
        }
        const auto shellFailure = auditRegionShells(*artifact);
        if (shellFailure != BoundaryConformingV3Reason::None)
            return artifactFailure(
                BoundaryConformingV3Status::Rejected,
                shellFailure,
                BoundaryConformingV3Phase::IndexedPreflight);

        ClosedArtifactBuildOutcome out;
        out.status = BoundaryConformingV3Status::Published;
        out.metrics = preflightMetrics;
        out.metrics.runComplete = true;
        out.artifact = std::move(artifact);
        return out;
    } catch (const std::bad_alloc&) {
        return artifactFailure(BoundaryConformingV3Status::Failed,
                               BoundaryConformingV3Reason::AllocationFailed,
                               BoundaryConformingV3Phase::IndexedPreflight);
    } catch (...) {
        return artifactFailure(BoundaryConformingV3Status::Failed,
                               BoundaryConformingV3Reason::InternalError,
                               BoundaryConformingV3Phase::IndexedPreflight);
    }
}

RecoveryEnvelopeBuildOutcome buildTetRecoveryEnvelope(
    RecoveryEnvelopeMode mode, std::span<const Vec3> constructionCoordinates,
    VolumeBoundaryProviderStamp expectedProviderStamp,
    const VolumeBoundaryOwnershipProvider& provider,
    std::size_t maxProviderNodes, std::size_t maxProviderFaces) noexcept {
    try {
        if (expectedProviderStamp.persistentId == 0 ||
            provider.stamp() != expectedProviderStamp)
            return envelopeFailure(
                BoundaryConformingV3Status::Rejected,
                BoundaryConformingV3Reason::StaleOwnershipProvider,
                VolumeProviderCheckpoint::PreEnvelope);
        if (constructionCoordinates.empty())
            return envelopeFailure(BoundaryConformingV3Status::Rejected,
                                   BoundaryConformingV3Reason::InvalidRecoveryEnvelope);
        auto envelope = std::make_unique<TetRecoveryEnvelope>();
        if (mode == RecoveryEnvelopeMode::DerivedFiniteBox) {
            Vec3 lo{},hi{};
            double referenceSpan = 0.0;
            if (!constructionBounds(
                    constructionCoordinates,lo,hi,referenceSpan))
                return envelopeFailure(
                    BoundaryConformingV3Status::Unsupported,
                    BoundaryConformingV3Reason::EnvelopeRangeUnrepresentable);
            const double inf = std::numeric_limits<double>::infinity();
            envelope->guardCoordinates.resize(8);
            envelope->outwardHullFaces.assign(
                kDerivedBoxFaces.begin(),kDerivedBoxFaces.end());
            ConstructorFrameScratch scratch;
            scratch.prepare(constructionCoordinates.size(),
                            envelope->outwardHullFaces,
                            envelope->guardCoordinates.size());
            auto accept = [&](const Vec3& guardLo,const Vec3& guardHi) {
                if (!finite(guardLo) || !finite(guardHi) ||
                    !(guardLo.x < lo.x && guardLo.y < lo.y &&
                      guardLo.z < lo.z && guardHi.x > hi.x &&
                      guardHi.y > hi.y && guardHi.z > hi.z))
                    return false;
                setDerivedBox(*envelope,guardLo,guardHi);
                return validateEnvelopeGeometry(
                           *envelope,constructionCoordinates) &&
                       validateConstructorFrame(
                           *envelope,constructionCoordinates,scratch);
            };

            Vec3 previousLo{std::nextafter(lo.x,-inf),
                            std::nextafter(lo.y,-inf),
                            std::nextafter(lo.z,-inf)};
            Vec3 previousHi{std::nextafter(hi.x,inf),
                            std::nextafter(hi.y,inf),
                            std::nextafter(hi.z,inf)};
            bool accepted = accept(previousLo,previousHi);
            // Fixed, terminating ladder. Rung zero is approximately one
            // relative binary-double ULP of the construction span; rung 52
            // is exactly one reference span. Outward nextafter guarantees
            // representable strictness even when +/-delta rounds away.
            for (int rung = 0; !accepted && rung <= 52; ++rung) {
                const double delta = std::ldexp(referenceSpan,rung-52);
                if (!std::isfinite(delta) || delta <= 0.0) continue;
                const Vec3 rawLo{lo.x-delta,lo.y-delta,lo.z-delta};
                const Vec3 rawHi{hi.x+delta,hi.y+delta,hi.z+delta};
                const Vec3 candidateLo{std::nextafter(rawLo.x,-inf),
                                       std::nextafter(rawLo.y,-inf),
                                       std::nextafter(rawLo.z,-inf)};
                const Vec3 candidateHi{std::nextafter(rawHi.x,inf),
                                       std::nextafter(rawHi.y,inf),
                                       std::nextafter(rawHi.z,inf)};
                if (sameDerivedBounds(candidateLo,candidateHi,
                                      previousLo,previousHi))
                    continue;
                previousLo = candidateLo;
                previousHi = candidateHi;
                accepted = accept(candidateLo,candidateHi);
            }
            if (!accepted)
                return envelopeFailure(
                    BoundaryConformingV3Status::Unsupported,
                    BoundaryConformingV3Reason::EnvelopeRangeUnrepresentable,
                    VolumeProviderCheckpoint::PostEnvelope);
        } else {
            ProviderRecoveryEnvelope supplied;
            try {
                supplied = provider.recoveryEnvelope(constructionCoordinates);
            } catch (const std::bad_alloc&) {
                return envelopeFailure(BoundaryConformingV3Status::Failed,
                                       BoundaryConformingV3Reason::AllocationFailed);
            } catch (...) {
                return envelopeFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderException);
            }
            if (supplied.providerStamp != expectedProviderStamp)
                return envelopeFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                    VolumeProviderCheckpoint::PostEnvelope);
            if (supplied.availability ==
                VolumeProviderAvailability::RecoveryEnvelopeUnavailable) {
                if (!supplied.guardCoordinates.empty() ||
                    !supplied.outwardHullFaces.empty())
                    return envelopeFailure(
                        BoundaryConformingV3Status::Failed,
                        BoundaryConformingV3Reason::OwnershipProviderContractViolation);
                return envelopeFailure(
                    BoundaryConformingV3Status::Unsupported,
                    BoundaryConformingV3Reason::ProviderEnvelopeUnavailable);
            }
            if (supplied.availability != VolumeProviderAvailability::Available ||
                supplied.guardCoordinates.size() > maxProviderNodes ||
                supplied.outwardHullFaces.size() > maxProviderFaces)
                return envelopeFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::OwnershipProviderContractViolation);
            envelope->guardCoordinates = std::move(supplied.guardCoordinates);
            envelope->outwardHullFaces = std::move(supplied.outwardHullFaces);
        }
        if (provider.stamp() != expectedProviderStamp)
            return envelopeFailure(
                BoundaryConformingV3Status::Failed,
                BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                VolumeProviderCheckpoint::PostEnvelope);
        if (mode == RecoveryEnvelopeMode::ProviderSuppliedConvexHull) {
            ConstructorFrameScratch finalScratch;
            finalScratch.prepare(constructionCoordinates.size(),
                                 envelope->outwardHullFaces,
                                 envelope->guardCoordinates.size());
            if (!validateEnvelopeGeometry(*envelope,constructionCoordinates) ||
                !validateConstructorFrame(
                    *envelope,constructionCoordinates,finalScratch))
                return envelopeFailure(
                    BoundaryConformingV3Status::Failed,
                    BoundaryConformingV3Reason::
                        OwnershipProviderContractViolation,
                    VolumeProviderCheckpoint::PostEnvelope);
        }
        RecoveryEnvelopeBuildOutcome out;
        out.status = BoundaryConformingV3Status::Published;
        out.envelope = std::move(envelope);
        return out;
    } catch (const std::bad_alloc&) {
        return envelopeFailure(BoundaryConformingV3Status::Failed,
                               BoundaryConformingV3Reason::AllocationFailed);
    } catch (...) {
        return envelopeFailure(BoundaryConformingV3Status::Failed,
                               BoundaryConformingV3Reason::InternalError);
    }
}

namespace {

bool validSizing(const sizing::MeshSpec& spec) noexcept {
    return spec.field != nullptr && (spec.order == 1 || spec.order == 2) &&
           std::isfinite(spec.orderFactor) && spec.orderFactor > 0.0 &&
           std::isfinite(spec.gradation) && spec.gradation >= 1.0 &&
           spec.elementsPerCircle >= 0 && spec.threads >= 0;
}

void setFailure(BoundaryConformingV3Outcome& out,
                BoundaryConformingV3Status status,
                BoundaryConformingV3Reason reason,
                BoundaryConformingV3Phase phase,
                VolumeProviderCheckpoint checkpoint =
                    VolumeProviderCheckpoint::None) noexcept {
    out.status = status;
    out.failure.reason = reason;
    out.failure.phase = phase;
    out.failure.providerCheckpoint = checkpoint;
}

BoundaryConformingV3WorkStatus workStatus(
    bnd::TetEdgeRecoveryStatus status) noexcept {
    switch (status) {
        case bnd::TetEdgeRecoveryStatus::AlreadyPresent:
            return BoundaryConformingV3WorkStatus::AlreadyPresent;
        case bnd::TetEdgeRecoveryStatus::Marked:
            return BoundaryConformingV3WorkStatus::Marked;
        case bnd::TetEdgeRecoveryStatus::Recovered:
            return BoundaryConformingV3WorkStatus::Recovered;
        case bnd::TetEdgeRecoveryStatus::Rejected:
            return BoundaryConformingV3WorkStatus::Rejected;
        case bnd::TetEdgeRecoveryStatus::Refused:
            return BoundaryConformingV3WorkStatus::Refused;
        case bnd::TetEdgeRecoveryStatus::Failed:
            return BoundaryConformingV3WorkStatus::Failed;
    }
    return BoundaryConformingV3WorkStatus::Failed;
}

BoundaryConformingV3WorkStatus workStatus(
    bnd::TetFacetRecoveryStatus status) noexcept {
    switch (status) {
        case bnd::TetFacetRecoveryStatus::AlreadyPresent:
            return BoundaryConformingV3WorkStatus::AlreadyPresent;
        case bnd::TetFacetRecoveryStatus::Marked:
            return BoundaryConformingV3WorkStatus::Marked;
        case bnd::TetFacetRecoveryStatus::Recovered:
            return BoundaryConformingV3WorkStatus::Recovered;
        case bnd::TetFacetRecoveryStatus::Rejected:
            return BoundaryConformingV3WorkStatus::Rejected;
        case bnd::TetFacetRecoveryStatus::Refused:
            return BoundaryConformingV3WorkStatus::Refused;
        case bnd::TetFacetRecoveryStatus::Failed:
            return BoundaryConformingV3WorkStatus::Failed;
    }
    return BoundaryConformingV3WorkStatus::Failed;
}

BoundaryConformingV3Phase edgePhase(
    bnd::TetEdgeRecoveryPhase phase) noexcept {
    switch (phase) {
        case bnd::TetEdgeRecoveryPhase::None:
            return BoundaryConformingV3Phase::None;
        case bnd::TetEdgeRecoveryPhase::Preflight:
            return BoundaryConformingV3Phase::IndexedPreflight;
        case bnd::TetEdgeRecoveryPhase::Routing:
            return BoundaryConformingV3Phase::EdgeRouting;
        case bnd::TetEdgeRecoveryPhase::Discovery:
            return BoundaryConformingV3Phase::EdgeDiscovery;
        case bnd::TetEdgeRecoveryPhase::Candidate:
            return BoundaryConformingV3Phase::EdgeCandidate;
        case bnd::TetEdgeRecoveryPhase::Preparation:
            return BoundaryConformingV3Phase::EdgePreparation;
        case bnd::TetEdgeRecoveryPhase::StateCopy:
            return BoundaryConformingV3Phase::EdgeStateCopy;
    }
    return BoundaryConformingV3Phase::EdgeCandidate;
}

BoundaryConformingV3Phase facetPhase(
    bnd::TetFacetRecoveryPhase phase) noexcept {
    switch (phase) {
        case bnd::TetFacetRecoveryPhase::None:
            return BoundaryConformingV3Phase::None;
        case bnd::TetFacetRecoveryPhase::Preflight:
            return BoundaryConformingV3Phase::IndexedPreflight;
        case bnd::TetFacetRecoveryPhase::Routing:
            return BoundaryConformingV3Phase::FacetRouting;
        case bnd::TetFacetRecoveryPhase::Discovery:
            return BoundaryConformingV3Phase::FacetDiscovery;
        case bnd::TetFacetRecoveryPhase::Candidate:
            return BoundaryConformingV3Phase::FacetCandidate;
        case bnd::TetFacetRecoveryPhase::Preparation:
            return BoundaryConformingV3Phase::FacetPreparation;
        case bnd::TetFacetRecoveryPhase::StateCopy:
            return BoundaryConformingV3Phase::FacetStateCopy;
    }
    return BoundaryConformingV3Phase::FacetCandidate;
}

BoundaryConformingV3Reason edgeReason(
    const bnd::TetEdgeRecoveryOutcome& outcome) noexcept {
    using R = bnd::TetEdgeRecoveryReason;
    switch (outcome.reason) {
        case R::None: return BoundaryConformingV3Reason::None;
        case R::StaleWorkspaceStamp:
        case R::StaleRecoveryContext:
            return BoundaryConformingV3Reason::StaleRecoveryContext;
        case R::InvalidIdentity:
            return BoundaryConformingV3Reason::MissingArtifactEdge;
        case R::InvalidClassification:
        case R::StateKeyMismatch:
        case R::StaleBinding:
        case R::WrongBindingEdge:
        case R::WrongBindingClassification:
        case R::IdRebinding:
        case R::EdgeAlias:
            return BoundaryConformingV3Reason::EdgeMetadataConflict;
        case R::ProtectedConstraintConflict:
            return BoundaryConformingV3Reason::ProtectedConstraintConflict;
        case R::BoundaryStar:
        case R::NonManifoldStar:
        case R::DisconnectedStar:
        case R::MixedRegionStar:
        case R::NonCycleLink:
        case R::NonBallClosure:
            return BoundaryConformingV3Reason::NonStarEdgeCavity;
        case R::NonfiniteCoordinate:
        case R::DegenerateTarget:
        case R::UnsupportedCell:
        case R::NonPositiveTet:
        case R::MixedRegion:
        case R::TangentialTouch:
        case R::CoplanarContact:
        case R::IdentityDistinctContact:
            return BoundaryConformingV3Reason::UnsupportedEdgeDegeneracy;
        case R::InvalidCandidate:
            return BoundaryConformingV3Reason::InvalidEdgeCandidate;
        case R::PreparationConflict:
            return BoundaryConformingV3Reason::PreparationConflict;
        case R::AllocationFailed:
            return BoundaryConformingV3Reason::AllocationFailed;
        case R::PendingConstraint:
            return BoundaryConformingV3Reason::InternalError;
        case R::RecoveryConflictCycle:
            return BoundaryConformingV3Reason::RecoveryConflictCycle;
        case R::ResourceBudgetReached:
            return BoundaryConformingV3Reason::ResourceBudgetReached;
    }
    return BoundaryConformingV3Reason::InternalError;
}

BoundaryConformingV3Reason facetReason(
    const bnd::TetFacetRecoveryOutcome& outcome) noexcept {
    using R = bnd::TetFacetRecoveryReason;
    switch (outcome.reason) {
        case R::None: return BoundaryConformingV3Reason::None;
        case R::StaleWorkspaceStamp:
            return BoundaryConformingV3Reason::StaleRecoveryContext;
        case R::InvalidIdentity:
            return BoundaryConformingV3Reason::InvalidRequest;
        case R::InvalidClassification:
        case R::MissingProtectedTargetEdge:
        case R::InvalidEdgeClassification:
        case R::InvalidFaceClassification:
        case R::StateKeyMismatch:
        case R::StaleBinding:
        case R::WrongBindingFace:
        case R::WrongBindingClassification:
        case R::IdRebinding:
        case R::TargetAlias:
        case R::FaceAlias:
            return BoundaryConformingV3Reason::EdgeMetadataConflict;
        case R::ProtectedFacetConflict:
        case R::CrossedProtectedEdge:
        case R::CrossedProtectedFacet:
            return BoundaryConformingV3Reason::ProtectedConstraintConflict;
        case R::BoundaryEdgeStar:
        case R::NonManifoldEdgeStar:
        case R::DisconnectedEdgeStar:
        case R::MixedRegionStar:
        case R::NonBallClosure:
        case R::UnsupportedNonStarCavity:
            return BoundaryConformingV3Reason::UnsupportedFacetDegeneracy;
        case R::NonfiniteCoordinate:
        case R::DegenerateTarget:
        case R::UnsupportedCell:
        case R::NonPositiveTet:
        case R::MixedRegion:
        case R::OneSidedTargetUnsupported:
        case R::InvalidTargetIncidence:
        case R::UnsupportedDegeneracy:
            return BoundaryConformingV3Reason::UnsupportedFacetDegeneracy;
        case R::InvalidCandidate:
            return BoundaryConformingV3Reason::InvalidFacetCandidate;
        case R::PreparationConflict:
            return BoundaryConformingV3Reason::PreparationConflict;
        case R::AllocationFailed:
            return BoundaryConformingV3Reason::AllocationFailed;
        case R::PendingConstraint:
            return BoundaryConformingV3Reason::InternalError;
        case R::RecoveryConflictCycle:
            return BoundaryConformingV3Reason::RecoveryConflictCycle;
    }
    return BoundaryConformingV3Reason::InternalError;
}

std::uint64_t witnessMix(std::uint64_t hash,std::uint64_t value) noexcept {
    for (unsigned byte = 0; byte < 8; ++byte) {
        hash ^= (value >> (byte*8U)) & 0xffU;
        hash *= 1099511628211ULL;
    }
    return hash;
}

void fillWorkFailure(BoundaryConformingV3ReportRow& row) noexcept {
    row.failure.workKind = row.key.kind;
    row.failure.edgeKey = row.key.edgeKey;
    row.failure.facetKey = row.key.facetKey;
    row.failure.scheduleOrdinal = row.scheduleOrdinal;
    std::uint64_t digest = 1469598103934665603ULL;
    digest = witnessMix(digest,static_cast<std::uint64_t>(row.key.kind));
    if (row.key.kind == VolumeWorkKind::ArtifactEdge) {
        row.failure.reason = edgeReason(row.failure.edgeOutcome);
        row.failure.phase = edgePhase(row.failure.edgeOutcome.phase);
        digest = witnessMix(digest,static_cast<std::uint64_t>(
            row.key.edgeKey.canonicalGlobalNodes[0]));
        digest = witnessMix(digest,static_cast<std::uint64_t>(
            row.key.edgeKey.canonicalGlobalNodes[1]));
        digest = witnessMix(digest,static_cast<std::uint64_t>(
            row.failure.edgeOutcome.reason));
    } else if (row.key.kind == VolumeWorkKind::Facet) {
        row.failure.reason = facetReason(row.failure.facetOutcome);
        row.failure.phase = facetPhase(row.failure.facetOutcome.phase);
        digest = witnessMix(digest,row.key.facetKey.facePersistentId);
        digest = witnessMix(digest,row.key.facetKey.normalizedLocalCellOrdinal);
        digest = witnessMix(digest,static_cast<std::uint64_t>(
            row.failure.facetOutcome.reason));
    }
    row.failure.exactWitnessDigest = digest;
}

int permutationParity(const std::array<topo::NodeHandle,3>& directed,
                      const std::array<topo::NodeHandle,3>& canonical) {
    for (std::size_t shift = 0; shift < 3; ++shift)
        if (directed[0] == canonical[shift] &&
            directed[1] == canonical[(shift+1)%3] &&
            directed[2] == canonical[(shift+2)%3])
            return 1;
    for (std::size_t shift = 0; shift < 3; ++shift)
        if (directed[0] == canonical[shift] &&
            directed[1] == canonical[(shift+2)%3] &&
            directed[2] == canonical[(shift+1)%3])
            return -1;
    return 0;
}

struct ExtractedVolume {
    std::unique_ptr<TetAuthoritativeResult> authoritative;
    std::vector<topo::CellHandle> retainedOldCells;
    std::map<topo::CellHandle,int> labels;
    std::map<topo::NodeHandle,Index> constructionGlobalByOldNode;
    VolumeExtractionMetrics metrics{};
};

bool constructionMetricAdd(std::uint64_t& target,
                           std::uint64_t value) noexcept {
    if (value > (std::numeric_limits<std::uint64_t>::max)()-target)
        return false;
    target += value;
    return true;
}

bool aggregateSurfaceConstructionMetrics(
    const std::vector<ProductionSurfaceFaceMetrics>& rows,
    BoundaryConformingV3SurfaceConstructionMetrics& out) noexcept {
    out = {};
    bool allComplete = true;
    const auto add = [](std::uint64_t& target,std::uint64_t value) noexcept {
        return constructionMetricAdd(target,value);
    };
    for (const auto& row : rows) {
        if (!add(out.faceRows,1) ||
            !add(out.boundaryNodes,row.boundaryNodes) ||
            !add(out.interiorCandidates,row.interiorCandidates) ||
            !add(out.spacingQueries,row.spacingQueries) ||
            !add(out.spacingBucketsVisited,row.spacingBucketsVisited) ||
            !add(out.spacingDistanceTests,row.spacingDistanceTests) ||
            !add(out.constraintPointTests,row.constraintPointTests) ||
            !add(out.protectedSegments,row.protectedSegments) ||
            !add(out.recoveredProtectedSegments,
                 row.recoveredProtectedSegments))
            return false;
        const auto& source = row.construction;
        auto& target = out.scalable;
        if (!add(target.inputPoints,source.inputPoints) ||
            !add(target.coordinateSortComparisons,
                 source.coordinateSortComparisons) ||
            !add(target.insertionOrderComparisons,
                 source.insertionOrderComparisons) ||
            !add(target.bootstrapPoints,source.bootstrapPoints) ||
            !add(target.hullOrientationTests,source.hullOrientationTests) ||
            !add(target.bootstrapLegalizationEdgeVisits,
                 source.bootstrapLegalizationEdgeVisits) ||
            !add(target.bootstrapLegalizationOrientationTests,
                 source.bootstrapLegalizationOrientationTests) ||
            !add(target.bootstrapLegalizationIncircleTests,
                 source.bootstrapLegalizationIncircleTests) ||
            !add(target.bootstrapCocircularTieTests,
                 source.bootstrapCocircularTieTests) ||
            !add(target.bootstrapLegalizationFlips,
                 source.bootstrapLegalizationFlips) ||
            !add(target.bootstrapWorklistInsertions,
                 source.bootstrapWorklistInsertions) ||
            !add(target.bootstrapWorklistUpdates,
                 source.bootstrapWorklistUpdates) ||
            !add(target.bootstrapStaleEdgeVisits,
                 source.bootstrapStaleEdgeVisits) ||
            !add(target.pointsInserted,source.pointsInserted) ||
            !add(target.locationTriangleTests,
                 source.locationTriangleTests) ||
            !add(target.walkSteps,source.walkSteps) ||
            !add(target.cavityPredicateTests,
                 source.cavityPredicateTests) ||
            !add(target.cavityMembers,source.cavityMembers) ||
            !add(target.boundaryEdges,source.boundaryEdges) ||
            !add(target.trianglesCommitted,source.trianglesCommitted) ||
            !add(target.fullScanFallbacks,source.fullScanFallbacks) ||
            !add(target.fullScanTriangleTests,
                 source.fullScanTriangleTests))
            return false;
        target.bootstrapWorklistHighWater = (std::max)(
            target.bootstrapWorklistHighWater,
            source.bootstrapWorklistHighWater);
        target.peakLiveTriangles = (std::max)(
            target.peakLiveTriangles,source.peakLiveTriangles);
        target.materializedTriangleSlots = (std::max)(
            target.materializedTriangleSlots,
            source.materializedTriangleSlots);
        target.peakCavitySize = (std::max)(
            target.peakCavitySize,source.peakCavitySize);
        target.peakBoundarySize = (std::max)(
            target.peakBoundarySize,source.peakBoundarySize);
        target.scratchHighWaterBytes = (std::max)(
            target.scratchHighWaterBytes,source.scratchHighWaterBytes);
        target.accountedCapacityHighWaterBytes = (std::max)(
            target.accountedCapacityHighWaterBytes,
            source.accountedCapacityHighWaterBytes);
        allComplete = allComplete && source.adjacencyAuditPassed &&
                      source.runComplete;
    }
    out.scalable.adjacencyAuditPassed = allComplete;
    out.scalable.runComplete = allComplete;
    out.runComplete = allComplete;
    return true;
}

bool multiplyConstructionMetric(std::uint64_t value,std::uint64_t count,
                                std::uint64_t& result) noexcept {
    if (value != 0 && count >
            (std::numeric_limits<std::uint64_t>::max)()/value)
        return false;
    result = value*count;
    return true;
}

bool captureImportConstructionMetrics(
    const TetHandoffMetrics& source,std::uint64_t attempts,
    BoundaryConformingV3ConstructionMetrics& out) noexcept {
    if (!source.workAccountingComplete) return false;
    out.importAttempts = attempts;
    return multiplyConstructionMetric(source.importedNodes,attempts,
                                      out.importNodeVisits) &&
        multiplyConstructionMetric(source.importedCells,attempts,
                                   out.importCellVisits) &&
        multiplyConstructionMetric(source.importEdgeValidationCellVisits,
            attempts,out.importEdgeValidationCellVisits) &&
        multiplyConstructionMetric(source.importEdgeBuildCellVisits,
            attempts,out.importEdgeBuildCellVisits) &&
        multiplyConstructionMetric(source.importEdgeUses,attempts,
                                   out.importEdgeUses) &&
        multiplyConstructionMetric(source.importUniqueEdges,attempts,
                                   out.importUniqueEdges) &&
        multiplyConstructionMetric(source.importFaceValidationCellVisits,
            attempts,out.importFaceValidationCellVisits) &&
        multiplyConstructionMetric(source.importFaceBuildCellVisits,
            attempts,out.importFaceBuildCellVisits) &&
        multiplyConstructionMetric(source.importFaceUses,attempts,
                                   out.importFaceUses) &&
        multiplyConstructionMetric(source.importUniqueFaces,attempts,
                                   out.importUniqueFaces) &&
        multiplyConstructionMetric(source.importEdgeAudits,attempts,
                                   out.importEdgeAudits) &&
        multiplyConstructionMetric(source.importFaceAudits,attempts,
                                   out.importFaceAudits) &&
        multiplyConstructionMetric(source.importPerCellGlobalAudits,attempts,
                                   out.importPerCellGlobalAudits) &&
        multiplyConstructionMetric(
            source.importCandidateValidationCalls,attempts,
            out.importCandidateValidationCalls) &&
        multiplyConstructionMetric(
            source.importCandidateValidationPointRows,attempts,
            out.importCandidateValidationPointRows) &&
        multiplyConstructionMetric(
            source.importCandidateValidationCellRows,attempts,
            out.importCandidateValidationCellRows) &&
        multiplyConstructionMetric(
            source.importCandidateValidationFaceUses,attempts,
            out.importCandidateValidationFaceUses) &&
        multiplyConstructionMetric(
            source.importCandidateValidationFaceRows,attempts,
            out.importCandidateValidationFaceRows) &&
        multiplyConstructionMetric(
            source.importCandidateValidationBoundaryRows,attempts,
            out.importCandidateValidationBoundaryRows) &&
        multiplyConstructionMetric(
            source.importCandidateValidationDualNeighbourVisits,attempts,
            out.importCandidateValidationDualNeighbourVisits) &&
        multiplyConstructionMetric(
            source.importCandidateValidationBoundaryEdgeUses,attempts,
            out.importCandidateValidationBoundaryEdgeUses) &&
        multiplyConstructionMetric(
            source.importCandidateValidationBoundaryNeighbourVisits,attempts,
            out.importCandidateValidationBoundaryNeighbourVisits) &&
        multiplyConstructionMetric(
            source.importCandidateValidationHullTests,attempts,
            out.importCandidateValidationHullTests) &&
        multiplyConstructionMetric(
            source.importCandidateValidationSlotsInitialized,attempts,
            out.importCandidateValidationSlotsInitialized) &&
        multiplyConstructionMetric(source.importSelectionRowsCopied,attempts,
                                   out.importSelectionRowsCopied) &&
        multiplyConstructionMetric(
            source.importSelectionSortComparisons,attempts,
            out.importSelectionSortComparisons) &&
        multiplyConstructionMetric(
            source.importSelectionDuplicateComparisons,attempts,
            out.importSelectionDuplicateComparisons) &&
        multiplyConstructionMetric(source.importSelectionCellVisits,attempts,
                                   out.importSelectionCellVisits) &&
        multiplyConstructionMetric(
            source.importSelectionVertexVisits,attempts,
            out.importSelectionVertexVisits) &&
        multiplyConstructionMetric(
            source.importSelectionPointSlotsInitialized,attempts,
            out.importSelectionPointSlotsInitialized) &&
        multiplyConstructionMetric(source.importGlobalMapRows,attempts,
                                   out.importGlobalMapRows) &&
        multiplyConstructionMetric(source.importGlobalMapComparisons,attempts,
                                   out.importGlobalMapComparisons) &&
        multiplyConstructionMetric(source.importGlobalMapInsertions,attempts,
                                   out.importGlobalMapInsertions) &&
        multiplyConstructionMetric(
            source.importLocalPointSlotsInitialized,attempts,
            out.importLocalPointSlotsInitialized) &&
        multiplyConstructionMetric(
            source.importLocalGlobalSlotsInitialized,attempts,
            out.importLocalGlobalSlotsInitialized) &&
        multiplyConstructionMetric(
            source.importLocalCellSlotsInitialized,attempts,
            out.importLocalCellSlotsInitialized) &&
        multiplyConstructionMetric(
            source.importPersistentNodeSlotsInitialized,attempts,
            out.importPersistentNodeSlotsInitialized) &&
        multiplyConstructionMetric(source.partialImportRowsCopied,attempts,
                                   out.partialImportRowsCopied) &&
        multiplyConstructionMetric(
            source.partialImportGlobalSlotsInitialized,attempts,
            out.partialImportGlobalSlotsInitialized) &&
        multiplyConstructionMetric(
            source.partialImportCandidateMappingVisits,attempts,
            out.partialImportCandidateMappingVisits) &&
        multiplyConstructionMetric(
            source.partialImportEdgeValidationRows,attempts,
            out.partialImportEdgeValidationRows) &&
        multiplyConstructionMetric(
            source.partialImportFaceValidationRows,attempts,
            out.partialImportFaceValidationRows) &&
        multiplyConstructionMetric(
            source.partialImportProtectedEndpointVisits,attempts,
            out.partialImportProtectedEndpointVisits) &&
        multiplyConstructionMetric(
            source.partialImportAuthorityPointValidations,attempts,
            out.partialImportAuthorityPointValidations) &&
        multiplyConstructionMetric(
            source.partialImportDuplicateKeyComparisons,attempts,
            out.partialImportDuplicateKeyComparisons) &&
        multiplyConstructionMetric(
            source.partialImportResolutionNodeVisits,attempts,
            out.partialImportResolutionNodeVisits) &&
        multiplyConstructionMetric(
            source.partialImportEdgeResolutionLookups,attempts,
            out.partialImportEdgeResolutionLookups) &&
        multiplyConstructionMetric(
            source.partialImportFaceResolutionLookups,attempts,
            out.partialImportFaceResolutionLookups) &&
        multiplyConstructionMetric(source.partialImportEdgeDataReads,attempts,
                                   out.partialImportEdgeDataReads) &&
        multiplyConstructionMetric(
            source.partialImportFaceDataCopies,attempts,
            out.partialImportFaceDataCopies) &&
        multiplyConstructionMetric(
            source.partialImportFaceIncidenceRowsCopied,attempts,
            out.partialImportFaceIncidenceRowsCopied) &&
        multiplyConstructionMetric(
            source.partialImportSortComparisons,attempts,
            out.partialImportSortComparisons);
}

bool finalizeConstructionMetrics(
    BoundaryConformingV3ConstructionMetrics& metrics,
    const VolumeExtractionMetrics& extraction) noexcept {
    auto sum = [](std::uint64_t& result,auto... values) noexcept {
        result = 0;
        return (constructionMetricAdd(
            result,static_cast<std::uint64_t>(values)) && ...);
    };
    const auto& surface = metrics.surface;
    const auto& s = surface.scalable;
    if (!sum(metrics.surfaceWorkUnits,
             surface.faceRows,surface.boundaryNodes,
             surface.interiorCandidates,surface.spacingQueries,
             surface.spacingBucketsVisited,surface.spacingDistanceTests,
             surface.constraintPointTests,surface.protectedSegments,
             surface.recoveredProtectedSegments,s.inputPoints,
             s.coordinateSortComparisons,s.insertionOrderComparisons,
             s.bootstrapPoints,s.hullOrientationTests,
             s.bootstrapLegalizationEdgeVisits,
             s.bootstrapLegalizationOrientationTests,
             s.bootstrapLegalizationIncircleTests,
             s.bootstrapCocircularTieTests,s.bootstrapLegalizationFlips,
             s.bootstrapWorklistInsertions,s.bootstrapWorklistUpdates,
             s.bootstrapStaleEdgeVisits,s.pointsInserted,
             s.locationTriangleTests,s.walkSteps,s.cavityPredicateTests,
             s.cavityMembers,s.boundaryEdges,s.trianglesCommitted,
             s.fullScanFallbacks,s.fullScanTriangleTests))
        return false;
    const auto& t = metrics.tetCandidate;
    const auto& d = t.construction;
    if (!sum(metrics.tetCandidateWorkUnits,
             t.inputPoints,t.uniquePoints,t.candidateCells,t.interiorFaces,
             t.boundaryFaces,t.validationFaceUses,t.validationHullTests,t.discardedCandidateWork,
             d.inputPoints,d.pointsInserted,d.pointsDuplicate,d.pointsSkipped,
             d.walkFallbackScans,d.cavityTrimResurrections,
             d.bootstrapDistinctTests,d.bootstrapOrientationTests,
             d.walkCellVisits,d.walkOrientationTests,d.conflictQueries,
             d.conflictOrientationTests,d.conflictInsphereTests,
             d.fallbackCellTests,d.cavityCellVisits,d.cavityNeighbourVisits,
             d.trimCellVisits,d.trimOrientationTests,d.shellFacesEmitted,
             d.adjacencyAuditFaceVisits,d.outputCells,
             d.originalOrientationAuditTests,d.discardedConstructionWork))
        return false;
    if (!sum(metrics.importWorkUnits,
             metrics.importAttempts,metrics.importNodeVisits,
             metrics.importCellVisits,
             metrics.importEdgeValidationCellVisits,
             metrics.importEdgeBuildCellVisits,metrics.importEdgeUses,
             metrics.importUniqueEdges,
             metrics.importFaceValidationCellVisits,
             metrics.importFaceBuildCellVisits,metrics.importFaceUses,
             metrics.importUniqueFaces,metrics.importEdgeAudits,
             metrics.importFaceAudits,metrics.importPerCellGlobalAudits,
             metrics.importCandidateValidationCalls,
             metrics.importCandidateValidationPointRows,
             metrics.importCandidateValidationCellRows,
             metrics.importCandidateValidationFaceUses,
             metrics.importCandidateValidationFaceRows,
             metrics.importCandidateValidationBoundaryRows,
             metrics.importCandidateValidationDualNeighbourVisits,
             metrics.importCandidateValidationBoundaryEdgeUses,
             metrics.importCandidateValidationBoundaryNeighbourVisits,
             metrics.importCandidateValidationHullTests,
             metrics.importCandidateValidationSlotsInitialized,
             metrics.importSelectionRowsCopied,
             metrics.importSelectionSortComparisons,
             metrics.importSelectionDuplicateComparisons,
             metrics.importSelectionCellVisits,
             metrics.importSelectionVertexVisits,
             metrics.importSelectionPointSlotsInitialized,
             metrics.importGlobalMapRows,
             metrics.importGlobalMapComparisons,
             metrics.importGlobalMapInsertions,
             metrics.importLocalPointSlotsInitialized,
             metrics.importLocalGlobalSlotsInitialized,
             metrics.importLocalCellSlotsInitialized,
             metrics.importPersistentNodeSlotsInitialized,
             metrics.partialImportRowsCopied,
             metrics.partialImportGlobalSlotsInitialized,
             metrics.partialImportCandidateMappingVisits,
             metrics.partialImportEdgeValidationRows,
             metrics.partialImportFaceValidationRows,
             metrics.partialImportProtectedEndpointVisits,
             metrics.partialImportAuthorityPointValidations,
             metrics.partialImportDuplicateKeyComparisons,
             metrics.partialImportResolutionNodeVisits,
             metrics.partialImportEdgeResolutionLookups,
             metrics.partialImportFaceResolutionLookups,
             metrics.partialImportEdgeDataReads,
             metrics.partialImportFaceDataCopies,
             metrics.partialImportFaceIncidenceRowsCopied,
             metrics.partialImportSortComparisons))
        return false;
    const auto& edge = extraction.edgeBuild;
    const auto& face = extraction.faceBuild;
    if (!sum(metrics.extractionWorkUnits,
             extraction.barrierFaces,extraction.adjacencyFacesVisited,
             extraction.cellsVisited,extraction.components,
             extraction.exteriorComponents,
             extraction.retainedRegionComponents,extraction.labelSeeds,
             extraction.labelPropagations,extraction.retainedNodes,
             extraction.discardedNodes,extraction.rebuildCellInsertions,
             edge.validationCellVisits,edge.buildCellVisits,edge.edgeUses,
             edge.uniqueEdges,edge.finalEdgeAudits,
             edge.perCellGlobalAudits,face.validationCellVisits,
             face.buildCellVisits,face.faceUses,face.uniqueFaces,
             face.suppliedEdgeAudits,face.finalFaceAudits,
             face.perCellGlobalAudits))
        return false;
    if (!sum(metrics.totalWorkUnits,metrics.surfaceWorkUnits,
             metrics.tetCandidateWorkUnits,metrics.importWorkUnits,
             metrics.extractionWorkUnits))
        return false;
    metrics.runComplete = metrics.surface.runComplete &&
        metrics.tetCandidate.runComplete && extraction.runComplete &&
        metrics.importAttempts != 0 && metrics.importPerCellGlobalAudits == 0 &&
        extraction.edgeBuild.perCellGlobalAudits == 0 &&
        extraction.faceBuild.perCellGlobalAudits == 0;
    return metrics.runComplete;
}

std::size_t workspaceAccountedBytes(
    const topo::TopologyWorkspace& workspace) noexcept {
    return workspace.nodes().stats().totalBytes +
           workspace.cells().stats().totalBytes +
           workspace.edges().stats().totalBytes +
           workspace.faces().stats().totalBytes;
}

std::size_t flatMeshAccountedBytes(const mesh::Mesh& mesh) noexcept {
    std::size_t connectivity = 0;
    for (Index i = 0; i < mesh.cellCount(); ++i)
        connectivity += static_cast<std::size_t>(
            mesh::nodesPerCell(mesh.cell(i).type));
    return static_cast<std::size_t>(mesh.nodeCount()) *
               (sizeof(Vec3)+sizeof(mesh::Classification)) +
           static_cast<std::size_t>(mesh.cellCount()) * sizeof(mesh::Cell) +
           connectivity*sizeof(Index);
}

std::size_t surfaceArtifactAccountedBytes(
    const bnd::ClosedTri3Artifact& artifact) noexcept {
    return artifact.facets.capacity()*sizeof(bnd::VolumeFacet) +
           artifact.edges.capacity()*sizeof(bnd::VolumeArtifactEdge) +
           artifact.regionTags.capacity()*sizeof(int);
}

std::size_t recoveryEnvelopeAccountedBytes(
    const TetRecoveryEnvelope& envelope) noexcept {
    return envelope.guardCoordinates.capacity()*sizeof(Vec3) +
           envelope.outwardHullFaces.capacity()*
               sizeof(std::array<std::size_t,3>);
}

bool auditAuthoritativeBoundaryBijection(
    const bnd::ClosedTri3Artifact& artifact,
    const TetAuthoritativeResult& authoritative,
    BoundaryConformingV3Audit& audit,std::uint64_t& failedGate) {
    failedGate = 0;
    const auto& workspace = *authoritative.workspace;
    std::map<topo::NodeHandle,Index> globalByLocal;
    for (const auto& binding : authoritative.nodeBindings)
        if (!globalByLocal.emplace(binding.local,
                                   binding.existingGlobalIndex).second) {
            failedGate = 1;
            return false;
        }
    std::map<std::array<Index,3>,const bnd::VolumeFacet*> expected;
    for (const auto& facet : artifact.facets) {
        if (!expected.emplace(canonicalTriple(facet.directedGlobalNodes),
                              &facet).second) {
            failedGate = 2;
            return false;
        }
        if (facet.owners.size() == 1) ++audit.prescribedExterior;
        else if (facet.owners.size() == 2) ++audit.prescribedInterface;
        else { failedGate = 3; return false; }
    }
    std::set<std::array<Index,3>> seen;
    for (const auto faceHandle : workspace.faces().activeBySequence()) {
        const auto face = workspace.faces().get(faceHandle);
        if (face.cornerCount != 3 || face.incident.empty() ||
            face.incident.size() > 2) {
            ++audit.nonManifoldFacets;
            failedGate = 4;
            return false;
        }
        std::array<Index,3> key{};
        for (std::size_t i = 0;i < 3;++i) {
            const auto mapped = globalByLocal.find(face.corners[i]);
            if (mapped == globalByLocal.end()) {
                failedGate = 5;
                return false;
            }
            key[i] = mapped->second;
        }
        key = canonicalTriple(key);
        std::vector<int> ownerTags;
        ownerTags.reserve(face.incident.size());
        for (const auto& incidence : face.incident) {
            const topo::CellHandle cell{incidence.cellSlot,incidence.cellGen,
                                         workspace.cells().poolId()};
            if (!workspace.cells().isCurrent(cell)) {
                failedGate = 6;
                return false;
            }
            ownerTags.push_back(workspace.cells().get(cell).cls.tag);
        }
        std::sort(ownerTags.begin(),ownerTags.end());
        const bool boundary = ownerTags.size() == 1;
        const bool interface = ownerTags.size() == 2 &&
                               ownerTags[0] != ownerTags[1];
        if (!boundary && !interface) continue;
        const auto found = expected.find(key);
        if (found == expected.end()) {
            if (boundary) ++audit.unexpectedExterior;
            else ++audit.unexpectedInterface;
            failedGate = 7;
            return false;
        }
        const auto& facet = *found->second;
        std::vector<int> expectedTags;
        for (const auto& owner : facet.owners)
            expectedTags.push_back(owner.regionTag);
        std::sort(expectedTags.begin(),expectedTags.end());
        if (expectedTags != ownerTags ||
            !sameClass(face.cls,facet.classification) ||
            !seen.insert(key).second) {
            ++audit.duplicateFacets;
            failedGate = 8;
            return false;
        }
        if (boundary) ++audit.matchedExterior;
        else ++audit.matchedInterface;
    }
    for (const auto& [key,facet] : expected) {
        if (seen.count(key) != 0) continue;
        if (facet->owners.size() == 1) ++audit.orphanExterior;
        else ++audit.missingInterface;
    }
    const bool passed = audit.orphanExterior == 0 &&
           audit.unexpectedExterior == 0 &&
           audit.missingInterface == 0 && audit.unexpectedInterface == 0 &&
           audit.duplicateFacets == 0 && audit.nonManifoldFacets == 0 &&
           audit.matchedExterior == audit.prescribedExterior &&
           audit.matchedInterface == audit.prescribedInterface;
    if (!passed) failedGate = 9;
    return passed;
}

bool seedLabel(std::map<topo::CellHandle,int>& labels,
               topo::CellHandle cell, int value) {
    auto [it,inserted] = labels.emplace(cell,value);
    if (!inserted &&
        it->second == (std::numeric_limits<int>::min)()) {
        it->second = value;
        return true;
    }
    return inserted || it->second == value;
}

std::unique_ptr<ExtractedVolume> extractVolume(
    const bnd::ClosedTri3Artifact& artifact,
    const bnd::TetFacetSetState& state,
    const std::map<Index,topo::NodeHandle>& localByGlobal,
    const TetAuthoritativeResult& imported, std::size_t guardStart,
    topo::TopologyWorkspace& workspace, const bnd::TetRecoveryContext& context,
    BoundaryConformingV3Reason& failure) {
    auto result = std::make_unique<ExtractedVolume>();
    for (const auto& binding : imported.nodeBindings)
        result->constructionGlobalByOldNode.emplace(
            binding.local,binding.existingGlobalIndex);
    const int unknown = (std::numeric_limits<int>::min)();
    for (const auto& [cell,record] : context.cells) {
        result->labels.emplace(cell,unknown);
        bool guard = false;
        for (const auto dense : record.nodes) {
            const auto old = context.denseToNode[dense];
            const auto global = result->constructionGlobalByOldNode.find(old);
            if (global == result->constructionGlobalByOldNode.end()) {
                failure = BoundaryConformingV3Reason::InvalidTetCandidate;
                return nullptr;
            }
            guard = guard || static_cast<std::size_t>(global->second) >= guardStart;
        }
        if (guard && !seedLabel(result->labels,cell,-1)) {
            failure = BoundaryConformingV3Reason::RegionLabelConflict;
            return nullptr;
        }
    }
    std::set<topo::FaceHandle> barriers;
    for (const auto& facet : artifact.facets) {
        const auto binding = state.bindings.find(facet.key);
        if (binding == state.bindings.end() ||
            !workspace.faces().isCurrent(binding->second.recoveredFace)) {
            failure = BoundaryConformingV3Reason::BoundaryBijectionFailed;
            return nullptr;
        }
        const auto face = workspace.faces().get(binding->second.recoveredFace);
        if (face.incident.size() != 2 ||
            !barriers.insert(binding->second.recoveredFace).second) {
            failure = BoundaryConformingV3Reason::BoundaryBijectionFailed;
            return nullptr;
        }
        std::array<topo::NodeHandle,3> canonical{
            face.corners[0],face.corners[1],face.corners[2]};
        std::sort(canonical.begin(),canonical.end(),[&](const auto a,const auto b) {
            return workspace.nodes().seqOf(a) < workspace.nodes().seqOf(b);
        });
        const std::array<topo::NodeHandle,3> stored{
            face.corners[0],face.corners[1],face.corners[2]};
        std::array<topo::NodeHandle,3> directed{};
        for (std::size_t i = 0; i < 3; ++i) {
            const auto mapped = localByGlobal.find(facet.directedGlobalNodes[i]);
            if (mapped == localByGlobal.end()) {
                failure = BoundaryConformingV3Reason::InvalidOwnerSide;
                return nullptr;
            }
            directed[i] = mapped->second;
        }
        const int storedParity = permutationParity(stored,canonical);
        const int directedParity = permutationParity(directed,canonical);
        if (storedParity == 0 || directedParity == 0) {
            failure = BoundaryConformingV3Reason::InvalidOwnerSide;
            return nullptr;
        }
        std::map<int,topo::CellHandle> ownerCell;
        for (const auto& incidence : face.incident) {
            const topo::CellHandle cell{incidence.cellSlot,incidence.cellGen,
                                        workspace.cells().poolId()};
            const int outward = incidence.sense * storedParity * directedParity;
            for (const auto& owner : facet.owners)
                if (owner.directedSense == outward) {
                    if (!ownerCell.emplace(owner.regionTag,cell).second) {
                        failure = BoundaryConformingV3Reason::InvalidOwnerSide;
                        return nullptr;
                    }
                }
        }
        if (ownerCell.size() != facet.owners.size()) {
            failure = BoundaryConformingV3Reason::InvalidOwnerSide;
            return nullptr;
        }
        for (const auto& owner : facet.owners)
            if (!seedLabel(result->labels,ownerCell.at(owner.regionTag),
                           owner.regionTag)) {
                failure = BoundaryConformingV3Reason::RegionLabelConflict;
                return nullptr;
            }
        if (facet.owners.size() == 1) {
            for (const auto& incidence : face.incident) {
                const topo::CellHandle cell{incidence.cellSlot,incidence.cellGen,
                                            workspace.cells().poolId()};
                if (cell != ownerCell.begin()->second &&
                    !seedLabel(result->labels,cell,-1)) {
                    failure = BoundaryConformingV3Reason::RegionLabelConflict;
                    return nullptr;
                }
            }
        }
    }
    result->metrics.barrierFaces = barriers.size();
    std::queue<topo::CellHandle> queue;
    for (const auto& [cell,label] : result->labels)
        if (label != unknown) {
            queue.push(cell);
            ++result->metrics.labelSeeds;
        }
    while (!queue.empty()) {
        const auto cell = queue.front(); queue.pop();
        const int label = result->labels.at(cell);
        const auto& record = context.cells.at(cell);
        static constexpr std::array<std::array<int,3>,4> faces{{
            {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
        for (const auto& localFace : faces) {
            const auto key = bnd::DenseTetFaceKey::of({
                record.nodes[static_cast<std::size_t>(localFace[0])],
                record.nodes[static_cast<std::size_t>(localFace[1])],
                record.nodes[static_cast<std::size_t>(localFace[2])]});
            const auto row = context.faces.find(key);
            if (row == context.faces.end()) {
                failure = BoundaryConformingV3Reason::RecoveryContextInvalid;
                return nullptr;
            }
            ++result->metrics.adjacencyFacesVisited;
            if (barriers.count(row->second.handle) ||
                row->second.owners.size() != 2)
                continue;
            const auto next = row->second.owners[0].cell == cell
                ? row->second.owners[1].cell : row->second.owners[0].cell;
            int& nextLabel = result->labels.at(next);
            if (nextLabel == unknown) {
                nextLabel = label;
                queue.push(next);
                ++result->metrics.labelPropagations;
            } else if (nextLabel != label) {
                failure = BoundaryConformingV3Reason::RegionLabelConflict;
                return nullptr;
            }
        }
    }
    for (const auto& [cell,label] : result->labels) {
        (void)cell;
        if (label == unknown) {
            failure = BoundaryConformingV3Reason::UnseededDomainComponent;
            return nullptr;
        }
    }

    // A region may own several closed shell components (an outer shell plus
    // sealed cavities), but its retained Tet domain must still be one
    // connected component in this increment.  Reusing one region tag for
    // disconnected solids is the explicit structured unsupported boundary.
    std::set<topo::CellHandle> componentSeen;
    std::map<int,std::size_t> componentsByRegion;
    for (const auto& [seed,label] : result->labels) {
        if (componentSeen.count(seed) != 0) continue;
        ++result->metrics.components;
        if (label < 0) ++result->metrics.exteriorComponents;
        else {
            ++result->metrics.retainedRegionComponents;
            if (++componentsByRegion[label] > 1) {
                failure =
                    BoundaryConformingV3Reason::DisconnectedRegionTagReuse;
                return nullptr;
            }
        }
        std::queue<topo::CellHandle> component;
        component.push(seed);
        componentSeen.insert(seed);
        while (!component.empty()) {
            const auto cell = component.front();
            component.pop();
            ++result->metrics.cellsVisited;
            const auto& record = context.cells.at(cell);
            static constexpr std::array<std::array<int,3>,4> localFaces{{
                {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
            for (const auto& localFace : localFaces) {
                const auto key = bnd::DenseTetFaceKey::of({
                    record.nodes[static_cast<std::size_t>(localFace[0])],
                    record.nodes[static_cast<std::size_t>(localFace[1])],
                    record.nodes[static_cast<std::size_t>(localFace[2])]});
                const auto row = context.faces.find(key);
                if (row == context.faces.end()) {
                    failure =
                        BoundaryConformingV3Reason::RecoveryContextInvalid;
                    return nullptr;
                }
                if (barriers.count(row->second.handle) != 0 ||
                    row->second.owners.size() != 2)
                    continue;
                const auto next = row->second.owners[0].cell == cell
                    ? row->second.owners[1].cell
                    : row->second.owners[0].cell;
                if (result->labels.at(next) != label) continue;
                if (componentSeen.insert(next).second) component.push(next);
            }
        }
    }

    std::set<topo::NodeHandle> used;
    for (const auto& [cell,label] : result->labels) {
        if (label < 0) { ++result->metrics.discardedCells; continue; }
        result->retainedOldCells.push_back(cell);
        ++result->metrics.retainedCells;
        const auto data = workspace.cells().get(cell);
        for (std::size_t i = 0; i < 4; ++i) used.insert(data.nodes[i]);
    }
    if (result->retainedOldCells.empty()) {
        failure = BoundaryConformingV3Reason::UnseededDomainComponent;
        return nullptr;
    }
    std::vector<topo::NodeHandle> oldNodes(used.begin(),used.end());
    std::sort(oldNodes.begin(),oldNodes.end(),[&](const auto a,const auto b) {
        return workspace.nodes().seqOf(a) < workspace.nodes().seqOf(b);
    });
    topo::NodePool nodes;
    std::map<topo::NodeHandle,topo::NodeHandle> newByOld;
    auto authoritative = std::make_unique<TetAuthoritativeResult>();
    for (const auto old : oldNodes) {
        const auto data = workspace.nodes().get(old);
        const auto local = nodes.add(data.coord,data.cls);
        newByOld.emplace(old,local);
        const Index global = result->constructionGlobalByOldNode.at(old);
        authoritative->nodeBindings.push_back({kInvalidIndex,global,local});
        if (authoritative->existingGlobalByLocalNodeSlot.size() <=
            static_cast<std::size_t>(local.slot))
            authoritative->existingGlobalByLocalNodeSlot.resize(
                static_cast<std::size_t>(local.slot)+1,kInvalidIndex);
        authoritative->existingGlobalByLocalNodeSlot[
            static_cast<std::size_t>(local.slot)] = global;
    }
    topo::CellPool cells;
    std::sort(result->retainedOldCells.begin(),result->retainedOldCells.end(),
              [&](const auto a,const auto b) {
        return workspace.cells().seqOf(a) < workspace.cells().seqOf(b);
    });
    for (const auto old : result->retainedOldCells) {
        const auto data = workspace.cells().get(old);
        std::array<topo::NodeHandle,4> conn{};
        for (std::size_t i = 0; i < 4; ++i) conn[i] = newByOld.at(data.nodes[i]);
        const auto local = cells.add(topo::ElementType::Tet4,conn,
                                     {3,result->labels.at(old)},nodes);
        authoritative->selectedCellBindings.push_back({kInvalidIndex,local});
        ++result->metrics.rebuildCellInsertions;
    }
    auto edges = topo::EdgePool::buildComplete(
        nodes,cells,&result->metrics.edgeBuild);
    auto faces = topo::FacePool::buildComplete(
        nodes,cells,edges,&result->metrics.faceBuild);
    for (const auto& [key,binding] : state.edgeBindings) {
        const auto ga = localByGlobal.find(key.canonicalGlobalNodes[0]);
        const auto gb = localByGlobal.find(key.canonicalGlobalNodes[1]);
        if (ga == localByGlobal.end() || gb == localByGlobal.end() ||
            newByOld.count(ga->second) == 0 || newByOld.count(gb->second) == 0)
            continue;
        const auto edge = edges.find(nodes,newByOld.at(ga->second),
                                     newByOld.at(gb->second));
        if (!edge) { failure = BoundaryConformingV3Reason::AuthoritativeRebuildFailed; return nullptr; }
        edges.setClassification(*edge,binding.classification);
        edges.setRequiredConstraint(*edge,true);
    }
    for (const auto& [key,binding] : state.bindings) {
        (void)key;
        std::array<topo::NodeHandle,3> mapped{};
        bool present = true;
        for (std::size_t i = 0; i < 3; ++i) {
            present = present && newByOld.count(binding.canonicalVertices[i]);
            if (present) mapped[i] = newByOld.at(binding.canonicalVertices[i]);
        }
        if (!present) continue;
        const auto face = faces.find(nodes,mapped);
        if (!face) { failure = BoundaryConformingV3Reason::AuthoritativeRebuildFailed; return nullptr; }
        faces.setClassification(*face,binding.classification);
    }
    authoritative->workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    result->metrics.retainedNodes = oldNodes.size();
    result->metrics.discardedNodes = context.denseToNode.size()-oldNodes.size();
    result->metrics.centroidClassifierCalls = 0;
    result->metrics.boundaryPeelingCalls = 0;
    result->metrics.hullCapPeelingCalls = 0;
    result->metrics.runComplete = true;
    result->authoritative = std::move(authoritative);
    failure = BoundaryConformingV3Reason::None;
    return result;
}

}  // namespace

BoundaryConformingV3Outcome
t2c_detail::meshModelBoundaryConformingWithProtectedSplits(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const BoundaryConformingV3Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& provider,
    std::span<const t2c_detail::DyadicEdgeSplit> protectedSplits,
    std::span<const t2c_detail::PrivateFacePoint> protectedFacePoints,
    t2c_detail::PrivateBuildTrace* protectedTrace) noexcept {
    BoundaryConformingV3Outcome out;
#if defined(PAACS_TOPO_TESTING)
    auto& testingWall = boundary_conforming_detail::
        boundaryConformingWallTraceForTesting();
    if (testingWall.enabled) testingWall.clear();
    struct ProtectedTraceWallCopy {
        t2c_detail::PrivateBuildTrace* trace = nullptr;
        boundary_conforming_detail::BoundaryConformingWallTrace* wall = nullptr;
        ~ProtectedTraceWallCopy() noexcept {
            if (trace && wall && wall->enabled)
                for (std::size_t i=0;i<wall->nanoseconds.size();++i)
                    trace->phaseNanoseconds[i]=wall->nanoseconds[i];
        }
    } protectedTraceWallCopy{protectedTrace,&testingWall};
    std::optional<BoundaryConformingTestingWallScope> testingWallScope;
#endif
    auto markAllocationPhase = [](BoundaryConformingV3Phase phase) noexcept {
#if defined(PAACS_TOPO_TESTING)
        boundary_conforming_detail::testingAllocationPhaseMarker() =
            static_cast<std::uint8_t>(phase);
#else
        (void)phase;
#endif
    };
    markAllocationPhase(BoundaryConformingV3Phase::ResultAssembly);
    BoundaryConformingV3MemoryMetrics memory;
    BoundaryConformingV3ConstructionMetrics constructionMetrics;
    auto recordOverlap = [&](std::size_t bytes,
                             BoundaryConformingV3Phase phase) noexcept {
        if (bytes > memory.peakAccountedOverlapBytes) {
            memory.peakAccountedOverlapBytes = bytes;
            memory.peakOverlapStage = phase;
        }
    };
    try {
        out.report = std::make_unique<BoundaryConformingV3ReportEnvelope>();
        out.report->reportComplete = true;
        out.reportComplete = true;
    } catch (const std::bad_alloc&) {
        out.status = BoundaryConformingV3Status::Failed;
        out.failure.reason = BoundaryConformingV3Reason::AllocationFailed;
        out.failure.phase = BoundaryConformingV3Phase::ResultAssembly;
        return out;
    }
    try {
        out.surfaceOutcome = std::make_unique<SurfaceSetOutcome>();
    } catch (const std::bad_alloc&) {
        out.status = BoundaryConformingV3Status::Failed;
        out.failure.reason = BoundaryConformingV3Reason::AllocationFailed;
        out.failure.phase = BoundaryConformingV3Phase::ResultAssembly;
        return out;
    }
    markAllocationPhase(BoundaryConformingV3Phase::Preflight);
    if (!validSizing(spec)) {
        setFailure(out,BoundaryConformingV3Status::Rejected,
                   BoundaryConformingV3Reason::InvalidSizingPolicy,
                   BoundaryConformingV3Phase::Preflight);
        return out;
    }
    if (request.expectedSurfaceSetStamp.persistentId == 0 ||
        request.expectedOwnershipProviderStamp.persistentId == 0) {
        setFailure(out,BoundaryConformingV3Status::Rejected,
                   BoundaryConformingV3Reason::InvalidRequest,
                   BoundaryConformingV3Phase::Preflight);
        return out;
    }
    markAllocationPhase(BoundaryConformingV3Phase::SurfaceConstruction);
    try {
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::Surface);
#endif
        const sizing::MeshSpec checkedSpec = detail::makeCheckedSpec(spec);
        auto surfaceMesh = std::make_unique<mesh::Mesh>();
        core::setThreadCount(spec.threads);
        if (protectedTrace) {
            // The v4 coordinator needs the private trace even when no edge
            // split is active. Build the solver-facing edge mesh and its
            // trace in one pass instead of first meshing into `surfaceMesh`
            // and then repeating the complete edge mesher into `traceOnly`.
            t2c_detail::meshEdgesWithProtectedSplits(
                model,checkedSpec,protectedSplits,*surfaceMesh,*protectedTrace,
                protectedTrace->proposedEdgeSplits);
        } else {
            meshEdges(model,checkedSpec,*surfaceMesh);
        }
        auto surfaceStage = detail::buildProductionSurfaceSetForV3(
            model,checkedSpec,request.expectedSurfaceSetStamp,*surfaceMesh,
            protectedTrace
                ? std::span<const t2c_detail::PrivateCurveNodeIdentity>(
                      protectedTrace->curveNodeIdentities)
                : std::span<const t2c_detail::PrivateCurveNodeIdentity>{},
            protectedFacePoints,
            protectedTrace ? protectedTrace->surfaceCacheSession : nullptr,
            protectedTrace
                ? protectedTrace->proposedEdgeSplits
                : std::span<const t2c_detail::DyadicEdgeSplit>{},
            protectedTrace
                ? protectedTrace->proposedFacePoints
                : std::span<const t2c_detail::PrivateFacePoint>{},
            protectedTrace&&protectedTrace->preferDirectConstraintInitialization
                ? detail::SurfaceCandidateConstructionMode::
                    ScalableDelaunayDirectConstraints
                : request.usePlanarSurfaceRecovery
                    ? detail::SurfaceCandidateConstructionMode::
                        ScalablePlanarRecovery
                    : detail::SurfaceCandidateConstructionMode::
                        ScalableDelaunay);
        if(protectedTrace&&protectedTrace->surfaceCacheSession) {
            protectedTrace->surfaceCacheDiagnostic=
                protectedTrace->surfaceCacheSession->diagnostic();
            protectedTrace->surfaceCacheMetrics=
                protectedTrace->surfaceCacheSession->metrics();
            if(protectedTrace->surfaceCacheChannel) {
                auto& channel=*protectedTrace->surfaceCacheChannel;
                channel.status=protectedTrace->surfaceCacheDiagnostic.status;
                channel.reason=protectedTrace->surfaceCacheDiagnostic.reason;
                channel.phase=protectedTrace->surfaceCacheDiagnostic.phase;
                channel.cacheWorkUnits=
                    protectedTrace->surfaceCacheDiagnostic.workUnits;
                channel.cacheAccountedBytes=
                    protectedTrace->surfaceCacheDiagnostic.accountedBytes;
            }
        }
        *out.surfaceOutcome = std::move(surfaceStage.surface);
        if (protectedTrace) {
            protectedTrace->surfaceNodes =
                static_cast<std::uint64_t>(surfaceMesh->nodeCount());
            protectedTrace->surfaceSegments = static_cast<std::uint64_t>(
                surfaceMesh->countCells(mesh::CellType::Line2));
        }
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        if (out.surfaceOutcome->setStatus != SurfaceSetStatus::Published ||
            !out.surfaceOutcome->retainedArtifact) {
            const auto status = statusForSurface(out.surfaceOutcome->setStatus);
            setFailure(out,status,reasonForSurface(out.surfaceOutcome->setStatus),
                       BoundaryConformingV3Phase::SurfaceConstruction);
            out.report->reportComplete = true;
            out.reportComplete = true;
            return out;
        }
        if (!aggregateSurfaceConstructionMetrics(
                surfaceStage.metrics,constructionMetrics.surface)) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::ConstructionMetricsOverflow,
                       BoundaryConformingV3Phase::SurfaceConstruction);
            return out;
        }
        markAllocationPhase(BoundaryConformingV3Phase::SurfaceMapping);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::Artifact);
#endif
        auto artifactBuilt = buildClosedTri3Artifact(
            *out.surfaceOutcome->retainedArtifact,*surfaceMesh,
            request.expectedSurfaceSetStamp,
            request.expectedOwnershipProviderStamp,provider);
        if (!artifactBuilt.artifact) {
            out.status = artifactBuilt.status;
            out.failure = artifactBuilt.failure;
            return out;
        }
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        auto artifact = std::move(artifactBuilt.artifact);
        memory.surfaceArtifactBytes = surfaceArtifactAccountedBytes(*artifact);
        memory.flatMeshBytes = flatMeshAccountedBytes(*surfaceMesh);
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes,
                      BoundaryConformingV3Phase::ArtifactIndexBuild);
        const std::size_t edgeCount = artifact->edges.size();
        const std::size_t facetCount = artifact->facets.size();
        if (edgeCount > (std::numeric_limits<std::size_t>::max)()-facetCount) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::InternalError,
                       BoundaryConformingV3Phase::ResultAssembly);
            out.report->reportComplete = false;
            out.reportComplete = false;
            return out;
        }
        try {
            const std::size_t total = edgeCount+facetCount;
            out.report->processed.reserve(total);
            out.report->unattempted.reserve(total);
            for (const auto& edge : artifact->edges)
                out.report->unattempted.push_back(
                    {VolumeWorkKind::ArtifactEdge,edge.key,{}});
            for (const auto& facet : artifact->facets)
                out.report->unattempted.push_back(
                    {VolumeWorkKind::Facet,{},facet.key});
        } catch (const std::bad_alloc&) {
            out.report->processed.clear();
            out.report->unattempted.clear();
            out.report->reportComplete = false;
            out.reportComplete = false;
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::AllocationFailed,
                       BoundaryConformingV3Phase::ResultAssembly);
            return out;
        }
        markAllocationPhase(BoundaryConformingV3Phase::RecoveryEnvelope);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::Envelope);
#endif
        auto envelopeBuilt = buildTetRecoveryEnvelope(
            request.envelopeMode,surfaceMesh->nodes(),
            request.expectedOwnershipProviderStamp,provider,
            request.maxProviderEnvelopeNodes,request.maxProviderEnvelopeFaces);
        if (!envelopeBuilt.envelope) {
            out.status = envelopeBuilt.status;
            out.failure = envelopeBuilt.failure;
            return out;
        }
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        auto envelope = std::move(envelopeBuilt.envelope);
        memory.recoveryEnvelopeBytes =
            recoveryEnvelopeAccountedBytes(*envelope);
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes+
                          memory.recoveryEnvelopeBytes,
                      BoundaryConformingV3Phase::RecoveryEnvelope);
        const std::size_t guardStart =
            static_cast<std::size_t>(surfaceMesh->nodeCount());
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::
                TetConstruction);
#endif
        mesh::Mesh construction = *surfaceMesh;
        for (const Vec3& guard : envelope->guardCoordinates)
            construction.addNode(guard,{3,0});
        markAllocationPhase(BoundaryConformingV3Phase::TetConstruction);
        auto candidateBuilt = buildTetCandidate(construction.nodes());
        if (protectedTrace && candidateBuilt.candidate)
            protectedTrace->candidateTetCells =
                candidateBuilt.candidate->cells.size();
#if defined(PAACS_TOPO_TESTING)
        {
            auto& trace = boundary_conforming_detail::
                boundaryConformingWallTraceForTesting();
            if (trace.enabled)
                trace.tetCandidateMetrics = candidateBuilt.metrics;
        }
#endif
        if (!candidateBuilt.built()) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       candidateBuilt.reason ==
                               TetCandidateReason::AllocationFailed
                           ? BoundaryConformingV3Reason::AllocationFailed
                           : BoundaryConformingV3Reason::TetConstructionFailed,
                       BoundaryConformingV3Phase::TetConstruction);
            return out;
        }
        constructionMetrics.tetCandidate = candidateBuilt.metrics;
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        auto candidate = std::move(candidateBuilt.candidate);
        memory.tetCandidateBytes = static_cast<std::size_t>(
            candidate->metrics.candidateAccountedBytes);
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes+
                          memory.recoveryEnvelopeBytes+
                          memory.tetCandidateBytes+
                          flatMeshAccountedBytes(construction),
                      BoundaryConformingV3Phase::TetConstruction);
        std::vector<Index> selected(candidate->cells.size());
        std::iota(selected.begin(),selected.end(),Index{0});
        std::vector<Index> globalByPoint(candidate->points.size(),kInvalidIndex);
        std::vector<mesh::Classification> classByPoint(candidate->points.size());
        for (std::size_t point = 0; point < candidate->points.size(); ++point) {
            const Index global = candidate->pointToInput[point];
            globalByPoint[point] = global;
            classByPoint[point] = construction.nodeClass(global);
        }
        std::vector<TetImportProtectedEdge> protectedImportEdges;
        std::vector<TetImportProtectedFace> protectedImportFaces;
        protectedImportEdges.reserve(artifact->edges.size());
        protectedImportFaces.reserve(artifact->facets.size());
        for (const auto& edge : artifact->edges)
            protectedImportEdges.push_back(
                {edge.key.canonicalGlobalNodes,edge.classification});
        for (const auto& facet : artifact->facets)
            protectedImportFaces.push_back(
                {facet.directedGlobalNodes,facet.classification});
        markAllocationPhase(BoundaryConformingV3Phase::AuthoritativeImport);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::
                AuthoritativeImport);
#endif
        const bool mixedExactPresentImport=protectedTrace&&
            protectedTrace->mixedPresentBatchSink;
        auto imported = mixedExactPresentImport
            ? tet_candidate_detail::
                  importTetCandidateWithPartialProtectedBoundary(
                *candidate,selected,globalByPoint,classByPoint,construction,
                {3,0},protectedImportEdges,protectedImportFaces)
            : importTetCandidateWithProtectedBoundary(
                *candidate,selected,globalByPoint,classByPoint,construction,
                {3,0},protectedImportEdges,protectedImportFaces);
        std::uint64_t protectedImportFallbacks = 0;
        if (!mixedExactPresentImport&&!imported.imported() &&
            imported.reason == TetCandidateReason::InvalidImportSelection) {
            // A missing prescribed edge/facet is the normal T2-B recovery
            // case, not a handoff failure. Rebuild the still-private generic
            // authoritative candidate and let indexed recovery handle it.
            imported = importTetCandidate(
                *candidate,selected,globalByPoint,classByPoint,construction,
                {3,0});
            protectedImportFallbacks = 1;
        }
        if (!imported.imported()) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       imported.reason == TetCandidateReason::AllocationFailed
                           ? BoundaryConformingV3Reason::AllocationFailed
                           : BoundaryConformingV3Reason::AuthoritativeImportFailed,
                       BoundaryConformingV3Phase::AuthoritativeImport);
            return out;
        }
        if (!captureImportConstructionMetrics(
                imported.result->metrics,protectedImportFallbacks+1,
                constructionMetrics)) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::ConstructionMetricsOverflow,
                       BoundaryConformingV3Phase::AuthoritativeImport);
            return out;
        }
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        auto authoritative = std::move(imported.result);
        memory.authoritativeWorkspaceBytes =
            workspaceAccountedBytes(*authoritative->workspace);
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes+
                          memory.recoveryEnvelopeBytes+
                          memory.tetCandidateBytes+
                          memory.authoritativeWorkspaceBytes,
                      BoundaryConformingV3Phase::AuthoritativeImport);
        candidate.reset();
        std::map<Index,topo::NodeHandle> localByGlobal;
        for (const auto& binding : authoritative->nodeBindings)
            localByGlobal.emplace(binding.existingGlobalIndex,binding.local);
        markAllocationPhase(BoundaryConformingV3Phase::RecoveryContextBuild);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::
                RecoveryContext);
#endif
        auto contextBuilt = bnd::buildTetRecoveryContext(
            *authoritative->workspace,{3,0});
        if (!contextBuilt.built()) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       contextBuilt.reason ==
                               bnd::TetRecoveryContextReason::AllocationFailed
                           ? BoundaryConformingV3Reason::AllocationFailed
                           : BoundaryConformingV3Reason::RecoveryContextInvalid,
                       BoundaryConformingV3Phase::RecoveryContextBuild);
            return out;
        }
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        auto context = std::move(contextBuilt.context);
        memory.recoveryIndexBytes = context->metrics.indexBytes;
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes+
                          memory.recoveryEnvelopeBytes+
                          memory.authoritativeWorkspaceBytes+
                          memory.recoveryIndexBytes,
                      BoundaryConformingV3Phase::RecoveryContextBuild);
        bnd::TetFacetSetState state;
        state.workspaceStamp = topo::WorkspaceStamp::of(*authoritative->workspace);
        state.artifactStamp = artifact->surfaceSetStamp;
        if (provider.stamp() != request.expectedOwnershipProviderStamp) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                       BoundaryConformingV3Phase::EdgeRouting,
                       VolumeProviderCheckpoint::PreEdgeRecovery);
            return out;
        }
        markAllocationPhase(BoundaryConformingV3Phase::IndexedPreflight);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::Recovery);
#endif
        auto recovered = request.edgeGeneralRecovery.enabled
            ? bnd::recoverTetFacetSetIndexed(
                *artifact,localByGlobal,{3,0},state,*authoritative->workspace,
                *context,request.edgeGeneralRecovery,
                protectedTrace?protectedTrace->frontierSink:nullptr,
                protectedTrace?protectedTrace->facetConflictCycleSink:nullptr,
                protectedTrace?protectedTrace->mixedPresentBatchSink:nullptr)
            : bnd::recoverTetFacetSetIndexed(
                *artifact,localByGlobal,{3,0},state,*authoritative->workspace,
                *context,protectedTrace?protectedTrace->frontierSink:nullptr,
                protectedTrace?protectedTrace->facetConflictCycleSink:nullptr,
                protectedTrace?protectedTrace->mixedPresentBatchSink:nullptr);
#if defined(PAACS_TOPO_TESTING)
        bnd::T2CEdgeLineageRunForTesting lineageRun;
        if(protectedTrace&&protectedTrace->frontierSink&&
           bnd::t2cEdgeLineageSinkForTesting())
            lineageRun=bnd::traceT2CEdgeLineageForTesting(
                *artifact,localByGlobal,{3,0},state,
                *authoritative->workspace,*context,recovered,
                *protectedTrace->frontierSink);
        (void)lineageRun;
#endif
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        if (protectedTrace) {
            auto partialConstruction = constructionMetrics;
            const VolumeExtractionMetrics noExtraction{};
            (void)finalizeConstructionMetrics(
                partialConstruction,noExtraction);
            protectedTrace->constructionWorkUnits =
                partialConstruction.totalWorkUnits;
            if(protectedTrace->surfaceCacheSession&&
               !constructionMetricAdd(
                    protectedTrace->constructionWorkUnits,
                    protectedTrace->surfaceCacheSession->attemptWorkUnits()))
                protectedTrace->workComplete=false;
            std::uint64_t recoveryWork = 0;
            const auto addRecovery = [&](std::uint64_t value) noexcept {
                if (!constructionMetricAdd(recoveryWork,value)) {
                    protectedTrace->workComplete = false;
                    return false;
                }
                return true;
            };
            const auto& scheduler = recovered.schedulerMetrics;
            bool recoveryComplete =
                addRecovery(scheduler.routeAssessments) &&
                addRecovery(scheduler.dependencyInvalidations) &&
                addRecovery(scheduler.selectedRebuilds) &&
                addRecovery(scheduler.batchPreflightAssessments) &&
                addRecovery(scheduler.alreadyPresentBatchEntities) &&
                addRecovery(scheduler.protectedConflicts) &&
                addRecovery(scheduler.conflictSccVisits) &&
                addRecovery(scheduler.fullStageRescans);
            if(protectedTrace->mixedPresentBatchSink) {
                const auto& mixed=
                    protectedTrace->mixedPresentBatchSink->metrics;
                recoveryComplete=recoveryComplete&&mixed.complete&&
                    addRecovery(mixed.copiedRows)&&
                    addRecovery(mixed.schedulerBindingAuditRows)&&
                    addRecovery(mixed.pendingRowsVisited)&&
                    addRecovery(mixed.dependencyPairsTested)&&
                    addRecovery(mixed.presentBindingEntities)&&
                    addRecovery(mixed.nontrivialRecoveries);
                const auto addBindingWork=[&](const auto& work) {
                    recoveryComplete=recoveryComplete&&work.complete&&
                        addRecovery(work.edgeBindingAuditCalls)&&
                        addRecovery(work.edgeBindingRowsVisited)&&
                        addRecovery(work.edgeAliasRowsVisited)&&
                        addRecovery(work.facetBindingAuditCalls)&&
                        addRecovery(work.facetBindingRowsVisited)&&
                        addRecovery(work.facetBindingRequestRowsVisited)&&
                        addRecovery(work.facetBindingDuplicateRowsCompared)&&
                        addRecovery(work.facetBindingReceiptBuilds)&&
                        addRecovery(work.facetBindingReceiptReuses)&&
                        addRecovery(work.facetBindingReceiptRowsCopied)&&
                        addRecovery(work.facetBindingReceiptLookups)&&
                        addRecovery(work.facetBindingReceiptReserveCalls)&&
                        addRecovery(work.edgeDataCopyRows)&&
                        addRecovery(work.faceDataCopyRows)&&
                        addRecovery(work.faceIncidenceRowsCopied)&&
                        addRecovery(work.stateCopies)&&
                        addRecovery(work.stateCopyRows);
                };
                addBindingWork(mixed.assessmentBindingWork);
                addBindingWork(mixed.recoveryBindingWork);
            } else
                recoveryComplete=recoveryComplete&&
                    addRecovery(scheduler.dependencyEdges);
            const auto addCounters=[&](const auto& counters) {
                recoveryComplete = recoveryComplete &&
                    addRecovery(counters.cellsScanned) &&
                    addRecovery(counters.crossedEdgeStars) &&
                    addRecovery(counters.closureInsertions) &&
                    addRecovery(counters.fixedPointPasses) &&
                    addRecovery(counters.apexAttempts) &&
                    addRecovery(counters.tetTetIntersectionTests);
            };
            if(protectedTrace->mixedPresentBatchSink) {
                const auto& mixed=
                    protectedTrace->mixedPresentBatchSink->metrics;
                addCounters(mixed.assessmentCounters);
                addCounters(mixed.recoveryCounters);
            } else if(protectedTrace->frontierSink&&
                      protectedTrace->frontierSink->rowsUsed!=0) {
                const auto& sink=*protectedTrace->frontierSink;
                protectedTrace->workComplete=
                    protectedTrace->workComplete&&sink.complete;
                addCounters(sink.assessmentCounters);
            } else for (const auto& row : recovered.processed) {
                addCounters(row.kind == bnd::TetFacetSetWorkKind::Edge
                    ? row.edgeOutcome.counters : row.facetOutcome.counters);
            }
            protectedTrace->recoveryWorkUnits = recoveryWork;
            protectedTrace->workComplete =
                protectedTrace->workComplete && recoveryComplete;
        }
        if (recovered.reportComplete) {
            // The exact buffers were reserved and prefilled before any Tet
            // work.  Rewriting them here cannot allocate.
            out.report->processed.clear();
            out.report->unattempted.clear();
            std::size_t ordinal = 0;
            for (const auto& source : recovered.processed) {
                BoundaryConformingV3ReportRow row;
                row.scheduleOrdinal = ordinal++;
                row.key.kind = source.kind == bnd::TetFacetSetWorkKind::Edge
                    ? VolumeWorkKind::ArtifactEdge : VolumeWorkKind::Facet;
                row.key.edgeKey = source.edgeKey;
                row.key.facetKey = source.facetKey;
                if (source.kind == bnd::TetFacetSetWorkKind::Edge) {
                    row.status = workStatus(source.edgeOutcome.status);
                    row.failure.edgeOutcome = source.edgeOutcome;
                } else {
                    row.status = workStatus(source.facetOutcome.status);
                    row.failure.facetOutcome = source.facetOutcome;
                }
                fillWorkFailure(row);
                out.report->processed.push_back(row);
            }
            for (const auto& key : recovered.unattemptedEdges)
                out.report->unattempted.push_back(
                    {VolumeWorkKind::ArtifactEdge,key,{}});
            for (const auto& key : recovered.unattemptedFacets)
                out.report->unattempted.push_back(
                    {VolumeWorkKind::Facet,{},key});
        }
        // An internal set-report allocation failure cannot erase the already
        // complete top-level multiplicity envelope.  In that case every work
        // item remains in its canonical unattempted row.
        out.report->reportComplete = true;
        out.reportComplete = true;
        if (recovered.status != bnd::TetFacetSetRecoveryStatus::Complete) {
            const auto status = recovered.status ==
                    bnd::TetFacetSetRecoveryStatus::Rejected
                ? BoundaryConformingV3Status::Rejected
                : recovered.status == bnd::TetFacetSetRecoveryStatus::Refused
                    ? BoundaryConformingV3Status::Refused
                    : BoundaryConformingV3Status::Failed;
            if (recovered.reportComplete && !out.report->processed.empty()) {
                out.status = status;
                out.failure = out.report->processed.back().failure;
            } else {
                setFailure(out,status,
                           BoundaryConformingV3Reason::AllocationFailed,
                           BoundaryConformingV3Phase::ResultAssembly);
            }
            return out;
        }
        if (provider.stamp() != request.expectedOwnershipProviderStamp) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                       BoundaryConformingV3Phase::DomainExtraction,
                       VolumeProviderCheckpoint::PreDomainExtraction);
            return out;
        }
        BoundaryConformingV3Reason extractionFailure{};
        markAllocationPhase(BoundaryConformingV3Phase::DomainExtraction);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::Extraction);
#endif
        auto extracted = extractVolume(
            *artifact,state,localByGlobal,*authoritative,guardStart,
            *authoritative->workspace,*context,extractionFailure);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        if (!extracted) {
            setFailure(out,
                       extractionFailure ==
                               BoundaryConformingV3Reason::DisconnectedRegionTagReuse
                           ? BoundaryConformingV3Status::Unsupported
                           : BoundaryConformingV3Status::Failed,
                       extractionFailure,
                       BoundaryConformingV3Phase::DomainExtraction);
            return out;
        }
        memory.extractedWorkspaceBytes =
            workspaceAccountedBytes(*extracted->authoritative->workspace);
        recordOverlap(memory.surfaceArtifactBytes+memory.flatMeshBytes+
                          memory.recoveryEnvelopeBytes+
                          memory.authoritativeWorkspaceBytes+
                          memory.recoveryIndexBytes+
                          memory.extractedWorkspaceBytes,
                      BoundaryConformingV3Phase::DomainExtraction);
        BoundaryConformingV3Audit topologyAudit;
        std::uint64_t topologyAuditGate = 0;
        markAllocationPhase(BoundaryConformingV3Phase::BoundaryAudit);
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.emplace(
            boundary_conforming_detail::BoundaryConformingWallPhase::
                AuditAndPublication);
#endif
        if (!auditAuthoritativeBoundaryBijection(
                *artifact,*extracted->authoritative,topologyAudit,
                topologyAuditGate)) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::BoundaryBijectionFailed,
                       BoundaryConformingV3Phase::BoundaryAudit);
            // Compact deterministic category witness retained on failure;
            // each low byte saturates rather than truncating a large count.
            const auto byte = [](std::uint64_t value) noexcept {
                return std::min<std::uint64_t>(value,255U);
            };
            out.failure.exactWitnessDigest =
                byte(topologyAudit.prescribedExterior) |
                (byte(topologyAudit.matchedExterior)<<8U) |
                (byte(topologyAudit.orphanExterior)<<16U) |
                (byte(topologyAudit.unexpectedExterior)<<24U) |
                (byte(topologyAudit.prescribedInterface)<<32U) |
                (byte(topologyAudit.matchedInterface)<<40U) |
                (byte(topologyAudit.missingInterface)<<48U) |
                (byte(topologyAudit.unexpectedInterface)<<56U);
            if (out.failure.exactWitnessDigest == 0)
                out.failure.exactWitnessDigest = topologyAuditGate;
            return out;
        }
        markAllocationPhase(BoundaryConformingV3Phase::MeshBuild);
        mesh::Mesh published = *surfaceMesh;
        published.reserveAppend(0,extracted->retainedOldCells.size(),
                                extracted->retainedOldCells.size()*4);
        std::vector<Index> publishedCells;
        publishedCells.reserve(extracted->retainedOldCells.size());
        for (const auto old : extracted->retainedOldCells) {
            const auto data = authoritative->workspace->cells().get(old);
            std::vector<Index> conn(4);
            for (std::size_t i = 0; i < 4; ++i) {
                const Index global =
                    extracted->constructionGlobalByOldNode.at(data.nodes[i]);
                if (global < 0 || static_cast<std::size_t>(global) >= guardStart)
                    throw std::logic_error("guard node retained");
                conn[i] = global;
            }
            publishedCells.push_back(published.addCell(
                mesh::CellType::Tet4,conn,{3,extracted->labels.at(old)}));
        }
        if (provider.stamp() != request.expectedOwnershipProviderStamp) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::OwnershipProviderStampChanged,
                       BoundaryConformingV3Phase::BoundaryAudit,
                       VolumeProviderCheckpoint::PreBoundaryAudit);
            return out;
        }
        mesh::AuditOptions auditOptions;
        auditOptions.boundaryMode = mesh::BoundaryComparisonMode::Required;
        auditOptions.boundaryOrphanBudget = 0;
        auditOptions.boundaryUnexpectedBudget = 0;
        auditOptions.duplicateFacetsAreErrors = true;
        auditOptions.unusedNodesAreErrors = true;
        const auto audit = mesh::auditMesh(published,auditOptions);
        memory.flatMeshBytes = flatMeshAccountedBytes(published);
        recordOverlap(memory.surfaceArtifactBytes+
                          memory.extractedWorkspaceBytes+
                          memory.flatMeshBytes,
                      BoundaryConformingV3Phase::BoundaryAudit);
        // MeshAudit's Required mode intentionally treats every surface Tri3
        // as an exterior prescription. In v3, an exactly audited interface
        // Tri3 has two opposite-region Tet owners and is therefore reported
        // by that legacy exterior-only comparison as one orphan. The exact
        // authoritative bijection above has already proved every such key,
        // both owner tags, classification and incidence. No exterior budget
        // is relaxed: interface-free fixtures still require literal zero.
        const auto expectedLegacyInterfaceOrphans = static_cast<Index>(
            topologyAudit.prescribedInterface);
        if (audit.boundaryOrphans != expectedLegacyInterfaceOrphans ||
            audit.boundaryUnexpected != 0 ||
            audit.boundaryDuplicates != 0 || audit.nonManifoldFaces != 0 ||
            audit.nInverted != 0 || audit.nNonFinite != 0 ||
            audit.unusedNodes != 0) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::BoundaryBijectionFailed,
                       BoundaryConformingV3Phase::BoundaryAudit);
            const auto byte = [](Index value) noexcept {
                return static_cast<std::uint64_t>(
                    std::min<Index>(value,255));
            };
            out.failure.exactWitnessDigest =
                byte(audit.boundaryOrphans) |
                (byte(audit.boundaryUnexpected)<<8U) |
                (byte(audit.boundaryDuplicates)<<16U) |
                (byte(audit.nonManifoldFaces)<<24U) |
                (byte(audit.nInverted)<<32U) |
                (byte(audit.nNonFinite)<<40U) |
                (byte(audit.unusedNodes)<<48U);
            return out;
        }
        if (!finalizeConstructionMetrics(
                constructionMetrics,extracted->metrics)) {
            setFailure(out,BoundaryConformingV3Status::Failed,
                       BoundaryConformingV3Reason::ConstructionMetricsOverflow,
                       BoundaryConformingV3Phase::ResultAssembly);
            return out;
        }
        if (protectedTrace)
            protectedTrace->constructionWorkUnits =
                constructionMetrics.totalWorkUnits;
        if(protectedTrace&&protectedTrace->surfaceCacheSession&&
           !constructionMetricAdd(
                protectedTrace->constructionWorkUnits,
                protectedTrace->surfaceCacheSession->attemptWorkUnits()))
            protectedTrace->workComplete=false;
        markAllocationPhase(BoundaryConformingV3Phase::ResultAssembly);
        auto retained = std::make_unique<BoundaryConformingV3Result>();
        retained->finalAuthoritative = std::move(extracted->authoritative);
        retained->constraintState = std::move(state);
        retained->artifactPreflightMetrics = artifactBuilt.metrics;
        retained->constructionMetrics = constructionMetrics;
        retained->recoveryMetrics.contextNodes = context->metrics.nodes;
        retained->recoveryMetrics.contextCells = context->metrics.cells;
        retained->recoveryMetrics.contextEdges = context->metrics.edges;
        retained->recoveryMetrics.contextFaces = context->metrics.faces;
        retained->recoveryMetrics.contextAuditVisits =
            context->metrics.auditVisits;
        retained->recoveryMetrics.indexBytes = context->metrics.indexBytes;
        retained->recoveryMetrics.edgeAabbQueries =
            context->metrics.segmentQueries;
        retained->recoveryMetrics.facetAabbQueries =
            context->metrics.facetQueries + context->metrics.boundsQueries;
        retained->recoveryMetrics.indexNodesVisited =
            context->metrics.indexNodesVisited;
        retained->recoveryMetrics.returnedTetCandidates =
            context->metrics.returnedCandidates;
        retained->recoveryMetrics.indexRemovals =
            context->metrics.indexRemovals;
        retained->recoveryMetrics.indexInsertions =
            context->metrics.indexInsertions;
        retained->recoveryMetrics.indexRebalances =
            context->metrics.indexRebuilds;
        retained->recoveryMetrics.indexHighWater =
            context->metrics.deltaHighWater;
        retained->recoveryMetrics.routeAssessments =
            recovered.schedulerMetrics.routeAssessments;
        retained->recoveryMetrics.dependencyEdges =
            recovered.schedulerMetrics.dependencyEdges;
        retained->recoveryMetrics.dependencyInvalidations =
            recovered.schedulerMetrics.dependencyInvalidations;
        retained->recoveryMetrics.selectedRebuilds =
            recovered.schedulerMetrics.selectedRebuilds;
        retained->recoveryMetrics.batchPreflightAssessments =
            recovered.schedulerMetrics.batchPreflightAssessments;
        retained->recoveryMetrics.alreadyPresentBatchEntities =
            recovered.schedulerMetrics.alreadyPresentBatchEntities;
        retained->recoveryMetrics.protectedImportFallbacks =
            protectedImportFallbacks;
        retained->recoveryMetrics.protectedConflicts =
            recovered.schedulerMetrics.protectedConflicts;
        retained->recoveryMetrics.conflictSccVisits =
            recovered.schedulerMetrics.conflictSccVisits;
        retained->recoveryMetrics.fullStageRescans =
            recovered.schedulerMetrics.fullStageRescans;
        for (const auto& row : recovered.processed) {
            const auto& counters = row.kind == bnd::TetFacetSetWorkKind::Edge
                ? row.edgeOutcome.counters : row.facetOutcome.counters;
            retained->recoveryMetrics.targetTetTests += counters.cellsScanned;
            retained->recoveryMetrics.completeStarQueries +=
                counters.crossedEdgeStars;
            retained->recoveryMetrics.closureInsertions +=
                counters.closureInsertions;
            retained->recoveryMetrics.apexPairs += counters.apexAttempts;
            retained->recoveryMetrics.candidateSurvivorTests +=
                counters.tetTetIntersectionTests;
            if (row.kind == bnd::TetFacetSetWorkKind::Edge &&
                row.edgeOutcome.status == bnd::TetEdgeRecoveryStatus::Recovered) {
                ++retained->recoveryMetrics.edgeCavityCommits;
                ++retained->recoveryMetrics.predictedEdgeMatches;
            }
            if (row.kind == bnd::TetFacetSetWorkKind::Facet &&
                row.facetOutcome.status ==
                    bnd::TetFacetRecoveryStatus::Recovered)
                ++retained->recoveryMetrics.facetCommits;
        }
        retained->recoveryMetrics.productionGlobalScanFallbacks = 0;
        retained->recoveryMetrics.runComplete = true;
        retained->extractionMetrics = extracted->metrics;
        memory.finalRetainedOverlapBytes =
            memory.extractedWorkspaceBytes+memory.flatMeshBytes;
        memory.runComplete = true;
        retained->memoryMetrics = memory;
        retained->boundaryAudit = topologyAudit;
        retained->boundaryAudit.invertedTets = audit.nInverted;
        retained->boundaryAudit.unusedNodes = audit.unusedNodes;
        retained->boundaryAudit.passed = true;
        retained->globalByAuthoritativeNode =
            retained->finalAuthoritative->existingGlobalByLocalNodeSlot;
        retained->globalByAuthoritativeCell = std::move(publishedCells);
        markAllocationPhase(BoundaryConformingV3Phase::Publication);
        boundary_conforming_detail::publicationTailActive() = true;
        destination.swap(published);
        out.retainedResult = std::move(retained);
        out.status = BoundaryConformingV3Status::Published;
        out.failure = {};
        out.destinationUnchanged = false;
        boundary_conforming_detail::publicationTailActive() = false;
#if defined(PAACS_TOPO_TESTING)
        testingWallScope.reset();
#endif
        return out;
    } catch (const detail::SizingEvaluationError&) {
        if (out.surfaceOutcome) {
            out.surfaceOutcome->setStatus = SurfaceSetStatus::Failed;
            out.surfaceOutcome->setFailure.reason = SurfaceReason::InternalError;
            out.surfaceOutcome->reportComplete = true;
        }
        setFailure(out,BoundaryConformingV3Status::Failed,
                   BoundaryConformingV3Reason::SizingEvaluationFailed,
                   BoundaryConformingV3Phase::SurfaceConstruction);
        return out;
    } catch (const std::bad_alloc&) {
        setFailure(out,BoundaryConformingV3Status::Failed,
                   BoundaryConformingV3Reason::AllocationFailed,
                   out.surfaceOutcome &&
                           out.surfaceOutcome->setStatus == SurfaceSetStatus::Published
                       ? BoundaryConformingV3Phase::ResultAssembly
                       : BoundaryConformingV3Phase::SurfaceConstruction);
        return out;
    } catch (...) {
        setFailure(out,BoundaryConformingV3Status::Failed,
                   BoundaryConformingV3Reason::InternalError,
                   BoundaryConformingV3Phase::ResultAssembly);
        return out;
    }
}

BoundaryConformingV3Outcome meshModelBoundaryConforming(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const BoundaryConformingV3Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& provider) noexcept {
    return t2c_detail::meshModelBoundaryConformingWithProtectedSplits(
        model,spec,request,destination,provider,{},{},nullptr);
}

}  // namespace paacs::algo
