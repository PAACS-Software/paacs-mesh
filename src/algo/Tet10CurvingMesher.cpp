// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Tet10CurvingMesher.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <numeric>
#include <set>
#include <utility>

namespace paacs::algo {
namespace {

bool finite(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool sameClass(const mesh::Classification &a, const mesh::Classification &b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

T4BStatus sourceStatus(T4AStatus status) noexcept {
    switch (status) {
    case T4AStatus::Published:
        return T4BStatus::Published;
    case T4AStatus::Incomplete:
        return T4BStatus::Incomplete;
    case T4AStatus::Rejected:
        return T4BStatus::Rejected;
    case T4AStatus::Unsupported:
        return T4BStatus::Unsupported;
    case T4AStatus::Failed:
        return T4BStatus::Failed;
    }
    return T4BStatus::Failed;
}

bool checkedAdd(std::uint64_t a, std::uint64_t b, std::uint64_t &out) noexcept {
    if (a > std::numeric_limits<std::uint64_t>::max() - b)
        return false;
    out = a + b;
    return true;
}
bool checkedMul(std::uint64_t a, std::uint64_t b, std::uint64_t &out) noexcept {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a)
        return false;
    out = a * b;
    return true;
}
bool addAccounted(std::uint64_t &total, std::uint64_t count, std::uint64_t width) noexcept {
    std::uint64_t bytes = 0, next = 0;
    if (!checkedMul(count, width, bytes) || !checkedAdd(total, bytes, next))
        return false;
    total = next;
    return true;
}

Tet10SourceEdgeKey makeKey(RefinementVertexKey a, RefinementVertexKey b) noexcept {
    if (b < a)
        std::swap(a, b);
    return {{{a, b}}};
}

struct EdgeWork {
    enum class Kind : std::uint8_t { Interior, Curve, Surface };
    Tet10SourceEdgeKey key{};
    std::array<Index, 2> source{{kInvalidIndex, kInvalidIndex}};
    mesh::Classification classification{3, 0};
    Kind kind = Kind::Interior;
    bool boundary = false;
    Index midside = kInvalidIndex;
    std::vector<bnd::VolumeFacetKey> incidentFacets;
};

constexpr std::array<std::array<unsigned, 2>, 6> tetEdges{
    {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
constexpr std::array<std::array<unsigned, 2>, 3> triEdges{{{{0, 1}}, {{1, 2}}, {{2, 0}}}};

bool appendCellEdges(const mesh::Mesh &source, Index cellIndex,
                     const std::vector<RefinementVertexKey> &keys,
                     std::map<Tet10SourceEdgeKey, EdgeWork> &edges, T4BMetrics &metrics) {
    const auto cell = source.cell(cellIndex);
    const Index *nodes = source.cellNodes(cellIndex);
    auto add = [&](unsigned a, unsigned b, bool boundary, EdgeWork::Kind kind) -> bool {
        Index ia = nodes[a], ib = nodes[b];
        if (ia < 0 || ib < 0 || ia >= source.nodeCount() || ib >= source.nodeCount() || ia == ib)
            return false;
        const auto &keyA = keys[static_cast<std::size_t>(ia)];
        const auto &keyB = keys[static_cast<std::size_t>(ib)];
        if (keyB < keyA)
            std::swap(ia, ib);
        const auto key =
            makeKey(keys[static_cast<std::size_t>(ia)], keys[static_cast<std::size_t>(ib)]);
        auto [it, inserted] = edges.emplace(
            key, EdgeWork{key, {{ia, ib}}, cell.onEntity, kind, boundary, kInvalidIndex, {}});
        if (!inserted && boundary) {
            it->second.boundary = true;
            if (kind == EdgeWork::Kind::Curve || it->second.kind == EdgeWork::Kind::Curve)
                it->second.kind = EdgeWork::Kind::Curve;
            else
                it->second.kind = EdgeWork::Kind::Surface;
        }
        ++metrics.sourceEdgeUses;
        return true;
    };
    switch (cell.type) {
    case mesh::CellType::Line2:
        return add(0, 1, true, EdgeWork::Kind::Curve);
    case mesh::CellType::Tri3:
        for (const auto &e : triEdges)
            if (!add(e[0], e[1], true, EdgeWork::Kind::Surface))
                return false;
        return true;
    case mesh::CellType::Tet4:
        for (const auto &e : tetEdges)
            if (!add(e[0], e[1], false, EdgeWork::Kind::Interior))
                return false;
        return true;
    default:
        return false;
    }
}

bool buildSourceKeys(const mesh::Mesh &source, const T4AMutableArtifact &artifact,
                     std::vector<RefinementVertexKey> &keys) {
    if (!artifact.authoritative || !artifact.authoritative->workspace)
        return false;
    keys.assign(static_cast<std::size_t>(source.nodeCount()), {});
    std::vector<unsigned char> seen(static_cast<std::size_t>(source.nodeCount()), 0);
    const auto &ws = *artifact.authoritative->workspace;
    for (const auto &[handle, key] : artifact.vertexKeyByHandle) {
        if (!ws.nodes().isCurrent(handle) || handle.slot < 0)
            return false;
        const auto slot = static_cast<std::size_t>(handle.slot);
        if (slot >= artifact.publishedGlobalByAuthoritativeNodeSlot.size())
            return false;
        const Index output = artifact.publishedGlobalByAuthoritativeNodeSlot[slot];
        if (output < 0 || output >= source.nodeCount() || seen[static_cast<std::size_t>(output)])
            return false;
        const auto data = ws.nodes().get(handle);
        if (data.coord.x != source.node(output).x || data.coord.y != source.node(output).y ||
            data.coord.z != source.node(output).z || !sameClass(data.cls, source.nodeClass(output)))
            return false;
        keys[static_cast<std::size_t>(output)] = key;
        seen[static_cast<std::size_t>(output)] = 1;
    }
    return std::all_of(seen.begin(), seen.end(), [](unsigned char v) { return v != 0; });
}

using OutputFaceKey = std::array<Index, 3>;

OutputFaceKey outputFaceKey(OutputFaceKey value) noexcept {
    std::sort(value.begin(), value.end());
    return value;
}

bool sameDirectedCycle(const OutputFaceKey &a, const OutputFaceKey &b) noexcept {
    for (unsigned shift = 0; shift < 3; ++shift)
        if (a[0] == b[shift] && a[1] == b[(shift + 1) % 3] && a[2] == b[(shift + 2) % 3])
            return true;
    return false;
}

bool validateBoundaryBindingsImpl(std::map<Tet10SourceEdgeKey, EdgeWork> &edges,
                                  const T4AMutableArtifact &artifact, const mesh::Mesh &source,
                                  const std::vector<RefinementVertexKey> &sourceKeys) {
    const auto &state = artifact.constraintState;
    if (!artifact.authoritative || !artifact.authoritative->workspace)
        return false;
    const auto &ws = *artifact.authoritative->workspace;
    if (!(state.workspaceStamp == topo::WorkspaceStamp::of(ws)))
        return false;
    if (state.bindings.size() != artifact.immutableBoundary.size() ||
        state.bindings.size() != state.canonicalSchedule.size() ||
        state.edgeBindings.size() != state.canonicalEdgeSchedule.size() ||
        !std::is_sorted(state.canonicalSchedule.begin(), state.canonicalSchedule.end()) ||
        !std::is_sorted(state.canonicalEdgeSchedule.begin(), state.canonicalEdgeSchedule.end()))
        return false;
    std::map<Index, topo::NodeHandle> handleByGlobal;
    for (const auto &[handle, key] : artifact.vertexKeyByHandle)
        if (key.kind == RefinementVertexKind::ExistingGlobal)
            if (!handleByGlobal.emplace(static_cast<Index>(key.value), handle).second)
                return false;
    std::map<Index, Index> outputByGlobal;
    for (Index output = 0; output < source.nodeCount(); ++output) {
        const auto &key = sourceKeys[static_cast<std::size_t>(output)];
        if (key.kind == RefinementVertexKind::ExistingGlobal &&
            !outputByGlobal.emplace(static_cast<Index>(key.value), output).second)
            return false;
    }
    static constexpr std::array<std::array<unsigned, 3>, 4> tetFaces{
        {{{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};
    std::map<OutputFaceKey, std::vector<OutputFaceKey>> tetOwners;
    std::map<OutputFaceKey, std::pair<OutputFaceKey, mesh::Classification>> surface;
    for (Index ci = 0; ci < source.cellCount(); ++ci) {
        const auto cell = source.cell(ci);
        const Index *nodes = source.cellNodes(ci);
        if (cell.type == mesh::CellType::Tet4) {
            for (const auto &local : tetFaces) {
                const OutputFaceKey directed{nodes[local[0]], nodes[local[1]], nodes[local[2]]};
                auto &owners = tetOwners[outputFaceKey(directed)];
                owners.push_back(directed);
                if (owners.size() > 2)
                    return false;
            }
        } else if (cell.type == mesh::CellType::Tri3) {
            const OutputFaceKey directed{nodes[0], nodes[1], nodes[2]};
            if (!surface.emplace(outputFaceKey(directed), std::make_pair(directed, cell.onEntity))
                     .second)
                return false;
        }
    }
    std::set<bnd::VolumeArtifactEdgeKey> scheduled;
    for (const auto &key : state.canonicalEdgeSchedule)
        if (!scheduled.insert(key).second || state.edgeBindings.count(key) != 1)
            return false;
    for (const auto &boundary : artifact.immutableBoundary) {
        const auto binding = state.bindings.find(boundary.facetKey);
        if (binding == state.bindings.end() ||
            !sameClass(binding->second.classification, boundary.classification) ||
            binding->second.recoveredFace != boundary.authoritativeFace ||
            !ws.faces().isCurrent(binding->second.recoveredFace))
            return false;
        const auto face = ws.faces().get(binding->second.recoveredFace);
        if (face.cornerCount != 3 || face.incident.size() != 1 ||
            !sameClass(face.cls, boundary.classification))
            return false;
        std::array<topo::NodeHandle, 3> directed{};
        for (unsigned i = 0; i < 3; ++i) {
            const auto found = handleByGlobal.find(boundary.directedGlobalNodes[i]);
            if (found == handleByGlobal.end())
                return false;
            directed[i] = found->second;
        }
        std::array<topo::NodeHandle, 3> canonical = directed;
        std::sort(canonical.begin(), canonical.end(), [&](const auto a, const auto b) {
            return ws.nodes().seqOf(a) < ws.nodes().seqOf(b);
        });
        if (canonical != binding->second.canonicalVertices)
            return false;
        auto poolKey =
            std::array<topo::NodeHandle, 3>{face.corners[0], face.corners[1], face.corners[2]};
        std::sort(poolKey.begin(), poolKey.end(), [&](const auto a, const auto b) {
            return ws.nodes().seqOf(a) < ws.nodes().seqOf(b);
        });
        if (canonical != poolKey)
            return false;
        OutputFaceKey prescribed{};
        for (unsigned i = 0; i < 3; ++i) {
            const auto found = outputByGlobal.find(boundary.directedGlobalNodes[i]);
            if (found == outputByGlobal.end())
                return false;
            prescribed[i] = found->second;
        }
        const auto key = outputFaceKey(prescribed);
        const auto surfaceRow = surface.find(key);
        const auto tetRow = tetOwners.find(key);
        if (surfaceRow == surface.end() || tetRow == tetOwners.end() ||
            tetRow->second.size() != 1 ||
            !sameDirectedCycle(surfaceRow->second.first, prescribed) ||
            !sameClass(surfaceRow->second.second, boundary.classification))
            return false;
        OutputFaceKey protectedDirected{};
        for (unsigned i = 0; i < 3; ++i) {
            const auto slot = static_cast<std::size_t>(face.corners[i].slot);
            if (slot >= artifact.publishedGlobalByAuthoritativeNodeSlot.size())
                return false;
            protectedDirected[i] = artifact.publishedGlobalByAuthoritativeNodeSlot[slot];
            if (protectedDirected[i] < 0 || protectedDirected[i] >= source.nodeCount())
                return false;
        }
        if (face.incident.front().sense < 0)
            std::swap(protectedDirected[1], protectedDirected[2]);
        if (!sameDirectedCycle(tetRow->second.front(), protectedDirected))
            return false;
        for (unsigned e = 0; e < 3; ++e) {
            const auto persistent = bnd::VolumeArtifactEdgeKey::of(
                boundary.directedGlobalNodes[e], boundary.directedGlobalNodes[(e + 1) % 3]);
            const auto row = edges.find(makeKey(
                RefinementVertexKey{RefinementVertexKind::ExistingGlobal,
                                    static_cast<std::uint64_t>(persistent.canonicalGlobalNodes[0])},
                RefinementVertexKey{
                    RefinementVertexKind::ExistingGlobal,
                    static_cast<std::uint64_t>(persistent.canonicalGlobalNodes[1])}));
            if (row == edges.end())
                return false;
            row->second.incidentFacets.push_back(boundary.facetKey);
        }
    }
    for (auto &[key, row] : edges) {
        if (!row.boundary)
            continue;
        if (key.vertices[0].kind != RefinementVertexKind::ExistingGlobal ||
            key.vertices[1].kind != RefinementVertexKind::ExistingGlobal)
            return false;
        const auto persistent = bnd::VolumeArtifactEdgeKey::of(
            static_cast<Index>(key.vertices[0].value), static_cast<Index>(key.vertices[1].value));
        const auto it = state.edgeBindings.find(persistent);
        if (it == state.edgeBindings.end() || !ws.edges().isCurrent(it->second.edge))
            return false;
        const auto data = ws.edges().get(it->second.edge);
        if (!it->second.required || !data.requiredConstraint ||
            !sameClass(it->second.classification, data.cls))
            return false;
        const auto ha = handleByGlobal.find(persistent.canonicalGlobalNodes[0]);
        const auto hb = handleByGlobal.find(persistent.canonicalGlobalNodes[1]);
        if (ha == handleByGlobal.end() || hb == handleByGlobal.end() ||
            !((data.a == ha->second && data.b == hb->second) ||
              (data.a == hb->second && data.b == ha->second)))
            return false;
        row.classification = it->second.classification;
        std::sort(row.incidentFacets.begin(), row.incidentFacets.end());
        row.incidentFacets.erase(std::unique(row.incidentFacets.begin(), row.incidentFacets.end()),
                                 row.incidentFacets.end());
        if (row.incidentFacets.empty())
            return false;
        if (row.kind == EdgeWork::Kind::Curve && row.classification.dim != 1)
            return false;
        if (row.kind == EdgeWork::Kind::Surface && row.classification.dim != 2)
            return false;
        if (row.kind == EdgeWork::Kind::Surface)
            for (const auto &facet : row.incidentFacets) {
                const auto owner = state.bindings.find(facet);
                if (owner == state.bindings.end() ||
                    !sameClass(owner->second.classification, row.classification))
                    return false;
            }
    }
    const auto boundaryCount = static_cast<std::size_t>(std::count_if(
        edges.begin(), edges.end(), [](const auto &value) { return value.second.boundary; }));
    if (boundaryCount != state.edgeBindings.size())
        return false;
    return true;
}

T4BReason validateBoundaryBindings(std::map<Tet10SourceEdgeKey, EdgeWork> &edges,
                                   const T4AMutableArtifact &artifact, const mesh::Mesh &source,
                                   const std::vector<RefinementVertexKey> &sourceKeys,
                                   const ArtifactStamp &expectedArtifact) {
    if (!artifact.authoritative || !artifact.authoritative->workspace)
        return T4BReason::InvalidAuthoritativeSource;
    const auto &state = artifact.constraintState;
    const auto &ws = *artifact.authoritative->workspace;
    if (state.artifactStamp.persistentId != expectedArtifact.persistentId)
        return T4BReason::ForeignEdgeBinding;
    if (state.artifactStamp.revision != expectedArtifact.revision ||
        !(state.workspaceStamp == topo::WorkspaceStamp::of(ws)))
        return T4BReason::StaleEdgeBinding;
    if (!std::is_sorted(state.canonicalSchedule.begin(), state.canonicalSchedule.end()) ||
        !std::is_sorted(state.canonicalEdgeSchedule.begin(), state.canonicalEdgeSchedule.end()) ||
        std::adjacent_find(state.canonicalSchedule.begin(), state.canonicalSchedule.end()) !=
            state.canonicalSchedule.end() ||
        std::adjacent_find(state.canonicalEdgeSchedule.begin(),
                           state.canonicalEdgeSchedule.end()) != state.canonicalEdgeSchedule.end())
        return T4BReason::DuplicateEdgeBinding;
    if (state.bindings.size() < artifact.immutableBoundary.size() ||
        state.canonicalSchedule.size() < artifact.immutableBoundary.size() ||
        state.edgeBindings.size() <
            static_cast<std::size_t>(std::count_if(
                edges.begin(), edges.end(), [](const auto &row) { return row.second.boundary; })))
        return T4BReason::MissingEdgeBinding;
    if (state.bindings.size() > artifact.immutableBoundary.size() ||
        state.canonicalSchedule.size() > artifact.immutableBoundary.size() ||
        state.edgeBindings.size() >
            static_cast<std::size_t>(std::count_if(
                edges.begin(), edges.end(), [](const auto &row) { return row.second.boundary; })) ||
        state.canonicalEdgeSchedule.size() != state.edgeBindings.size())
        return T4BReason::ForeignEdgeBinding;
    for (const auto &row : artifact.immutableBoundary) {
        const auto binding = state.bindings.find(row.facetKey);
        if (binding == state.bindings.end())
            return T4BReason::MissingEdgeBinding;
        if (binding->second.recoveredFace != row.authoritativeFace)
            return T4BReason::ForeignEdgeBinding;
        if (!ws.faces().isCurrent(binding->second.recoveredFace))
            return T4BReason::StaleEdgeBinding;
        if (!sameClass(binding->second.classification, row.classification))
            return T4BReason::ConflictingEdgeClassification;
    }
    for (const auto &key : state.canonicalEdgeSchedule) {
        const auto binding = state.edgeBindings.find(key);
        if (binding == state.edgeBindings.end())
            return T4BReason::MissingEdgeBinding;
        if (!ws.edges().isCurrent(binding->second.edge))
            return T4BReason::StaleEdgeBinding;
        const auto data = ws.edges().get(binding->second.edge);
        if (!sameClass(data.cls, binding->second.classification))
            return T4BReason::ConflictingEdgeClassification;
        if (!binding->second.required || !data.requiredConstraint)
            return T4BReason::ForeignEdgeBinding;
        if (binding->second.classification.dim != 1 && binding->second.classification.dim != 2)
            return T4BReason::AmbiguousBoundaryOwnership;
    }
    return validateBoundaryBindingsImpl(edges, artifact, source, sourceKeys)
               ? T4BReason::None
               : T4BReason::ForeignEdgeBinding;
}

T4BReason certificateReason(Tet10CertificateStatus status) noexcept {
    switch (status) {
    case Tet10CertificateStatus::PositiveBelowThreshold:
        return T4BReason::JacobianBelowThreshold;
    case Tet10CertificateStatus::Nonpositive:
        return T4BReason::NonpositiveJacobian;
    case Tet10CertificateStatus::RangeUnrepresentable:
        return T4BReason::CertificateRangeUnrepresentable;
    case Tet10CertificateStatus::IncompleteDepth:
    case Tet10CertificateStatus::IncompleteSubcells:
    case Tet10CertificateStatus::IncompleteCoefficients:
    case Tet10CertificateStatus::IncompleteWork:
        return T4BReason::JacobianCertificationBudgetReached;
    case Tet10CertificateStatus::AllocationFailed:
        return T4BReason::AllocationFailed;
    case Tet10CertificateStatus::InvalidInput:
        return T4BReason::InvalidElevatedTopology;
    case Tet10CertificateStatus::Certified:
        return T4BReason::None;
    }
    return T4BReason::InternalError;
}

T4BStatus certificateStatus(Tet10CertificateStatus status) noexcept {
    switch (status) {
    case Tet10CertificateStatus::IncompleteDepth:
    case Tet10CertificateStatus::IncompleteSubcells:
    case Tet10CertificateStatus::IncompleteCoefficients:
    case Tet10CertificateStatus::IncompleteWork:
        return T4BStatus::Incomplete;
    case Tet10CertificateStatus::RangeUnrepresentable:
        return T4BStatus::Unsupported;
    default:
        return T4BStatus::Failed;
    }
}

bool buildElevatedCells(const mesh::Mesh &source,
                        const std::map<Tet10SourceEdgeKey, EdgeWork> &edges,
                        const std::vector<RefinementVertexKey> &keys, mesh::Mesh &result,
                        std::vector<Index> &scratch) {
    auto mid = [&](Index a, Index b) -> Index {
        const auto it = edges.find(
            makeKey(keys[static_cast<std::size_t>(a)], keys[static_cast<std::size_t>(b)]));
        return it == edges.end() ? kInvalidIndex : it->second.midside;
    };
    for (Index ci = 0; ci < source.cellCount(); ++ci) {
        const auto c = source.cell(ci);
        const Index *n = source.cellNodes(ci);
        scratch.clear();
        if (c.type == mesh::CellType::Line2) {
            const Index m = mid(n[0], n[1]);
            if (m < 0)
                return false;
            scratch.push_back(n[0]);
            scratch.push_back(n[1]);
            scratch.push_back(m);
            result.addCell(mesh::CellType::Line3, scratch, c.onEntity);
        } else if (c.type == mesh::CellType::Tri3) {
            const Index m01 = mid(n[0], n[1]), m12 = mid(n[1], n[2]), m20 = mid(n[2], n[0]);
            if (m01 < 0 || m12 < 0 || m20 < 0)
                return false;
            scratch.insert(scratch.end(), {n[0], n[1], n[2], m01, m12, m20});
            result.addCell(mesh::CellType::Tri6, scratch, c.onEntity);
        } else if (c.type == mesh::CellType::Tet4) {
            scratch.insert(scratch.end(), {n[0], n[1], n[2], n[3]});
            for (const auto &e : tetEdges) {
                const Index m = mid(n[e[0]], n[e[1]]);
                if (m < 0)
                    return false;
                scratch.push_back(m);
            }
            result.addCell(mesh::CellType::Tet10, scratch, c.onEntity);
        } else
            return false;
    }
    return true;
}

void cloneMesh(const mesh::Mesh &source, mesh::Mesh &destination) {
    std::uint64_t connectivity = 0;
    for (Index i = 0; i < source.cellCount(); ++i)
        connectivity += static_cast<std::uint64_t>(mesh::nodesPerCell(source.cell(i).type));
    destination.reserveAppend(static_cast<std::size_t>(source.nodeCount()),
                              static_cast<std::size_t>(source.cellCount()),
                              static_cast<std::size_t>(connectivity));
    for (Index i = 0; i < source.nodeCount(); ++i)
        destination.addNode(source.node(i), source.nodeClass(i));
    for (Index i = 0; i < source.cellCount(); ++i) {
        const auto cell = source.cell(i);
        const int count = mesh::nodesPerCell(cell.type);
        std::vector<Index> nodes(static_cast<std::size_t>(count));
        const Index *input = source.cellNodes(i);
        for (int j = 0; j < count; ++j)
            nodes[static_cast<std::size_t>(j)] = input[j];
        destination.addCell(cell.type, nodes, cell.onEntity);
    }
}

} // namespace

#if defined(PAACS_TOPO_TESTING)
T4BReason
t4b_detail::validateSourceAuthorityForTesting(const mesh::Mesh &source,
                                              T4AMutableArtifact &artifact,
                                              const ArtifactStamp &expectedArtifact) noexcept {
    try {
        std::vector<RefinementVertexKey> sourceKeys;
        if (!buildSourceKeys(source, artifact, sourceKeys))
            return T4BReason::InvalidAuthoritativeSource;
        std::map<Tet10SourceEdgeKey, EdgeWork> edges;
        T4BMetrics metrics{};
        metrics.providerStampCalls = 1;
        for (Index ci = 0; ci < source.cellCount(); ++ci)
            if (!appendCellEdges(source, ci, sourceKeys, edges, metrics))
                return T4BReason::UnsupportedSourceCell;
        return validateBoundaryBindings(edges, artifact, source, sourceKeys, expectedArtifact);
    } catch (const std::bad_alloc &) {
        return T4BReason::AllocationFailed;
    } catch (...) {
        return T4BReason::InternalError;
    }
}
#endif

T4BOutcome meshModelCurvedTet10(const geo::Model &model, const sizing::MeshSpec &spec,
                                const T4BRequest &request, mesh::Mesh &destination,
                                const VolumeBoundaryOwnershipProvider &ownershipProvider,
                                const Tet10PlacementProvider &placementProvider) noexcept {
#if defined(PAACS_TOPO_TESTING)
    struct RouteTrackingScope {
        bool previous = false;
        RouteTrackingScope() noexcept : previous(t4b_detail::routeActive()) {
            t4b_detail::routeActive() = true;
        }
        ~RouteTrackingScope() { t4b_detail::routeActive() = previous; }
    } routeTrackingScope;
#endif
    T4BOutcome out;
    T4BPhase phase = T4BPhase::ReportShell;
    T4BMetrics metrics{};
    try {
        out.report = std::make_unique<T4BReport>();
        phase = T4BPhase::RequestPreflight;
        if (!std::isfinite(request.policy.minimumCertifiedJacobianRatio) ||
            request.policy.minimumCertifiedJacobianRatio < 0.0 ||
            request.policy.minimumCertifiedJacobianRatio > 1.0 ||
            request.policy.maxCertificateDepth > 20) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::InvalidRequest, phase};
            return out;
        }
        auto consumeWork = [&](std::uint64_t units, T4BPhase checkpoint,
                               const Tet10SourceEdgeKey &edge = {}) noexcept {
            std::uint64_t next = 0;
            if (!checkedAdd(metrics.workUnits, units, next) || next > request.policy.maxWorkUnits) {
                out.status = T4BStatus::Incomplete;
                out.failure = {T4BReason::JacobianCertificationBudgetReached, checkpoint, edge};
                return false;
            }
            metrics.workUnits = next;
            return true;
        };
        Tet10PlacementStamp initialStamp{};
        if (!consumeWork(1, phase)) {
            out.reportComplete = true;
            return out;
        }
        ++metrics.providerStampCalls;
        try {
            initialStamp = placementProvider.stamp();
        } catch (const std::bad_alloc &) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::AllocationFailed, phase};
            return out;
        } catch (...) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::PlacementProviderThrew, phase};
            return out;
        }
        if (initialStamp.persistentId != request.expectedPlacementStamp.persistentId) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::ForeignPlacementProvider, phase};
            return out;
        }
        if (initialStamp.revision != request.expectedPlacementStamp.revision) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::StalePlacementProvider, phase};
            return out;
        }
        auto checkedStamp = [&](T4BPhase checkpoint, const Tet10SourceEdgeKey &edge = {}) noexcept {
            if (!consumeWork(1, checkpoint, edge))
                return false;
            ++metrics.providerStampCalls;
            Tet10PlacementStamp observed{};
            try {
                observed = placementProvider.stamp();
            } catch (const std::bad_alloc &) {
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::AllocationFailed, checkpoint, edge};
                return false;
            } catch (...) {
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::PlacementProviderThrew, checkpoint, edge};
                return false;
            }
            if (observed.persistentId != request.expectedPlacementStamp.persistentId) {
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::PlacementProviderIdentityChanged, checkpoint, edge};
                return false;
            }
            if (observed.revision != request.expectedPlacementStamp.revision) {
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::StalePlacementProvider, checkpoint, edge};
                return false;
            }
            return true;
        };
        phase = T4BPhase::SourceT4A;
        mesh::Mesh sourceMesh;
#if defined(PAACS_TOPO_TESTING)
        struct NestedSourceTrackingPause {
            bool previous = false;
            NestedSourceTrackingPause() noexcept : previous(t4b_detail::routeActive()) {
                t4b_detail::routeActive() = false;
            }
            ~NestedSourceTrackingPause() { t4b_detail::routeActive() = previous; }
        } sourceTrackingPause;
#endif
        auto source =
            meshModelQualityTet4(model, spec, request.source, sourceMesh, ownershipProvider);
#if defined(PAACS_TOPO_TESTING)
        t4b_detail::routeActive() = true;
#endif
        out.sourceOutcome = std::make_unique<T4AOutcome>(std::move(source));
        if (!out.sourceOutcome->published()) {
            out.status = sourceStatus(out.sourceOutcome->status);
            out.failure = {T4BReason::SourceT4AStopped, phase};
            return out;
        }
        phase = T4BPhase::SourceValidation;
        if (!out.sourceOutcome->retainedResult ||
            !out.sourceOutcome->retainedResult->publicationReady ||
            !out.sourceOutcome->retainedResult->artifact ||
            sourceMesh.countCells(mesh::CellType::Tet4) <= 0) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::InvalidAuthoritativeSource, phase};
            return out;
        }
        auto &sourceResult = *out.sourceOutcome->retainedResult;
        auto &artifact = *sourceResult.artifact;
        std::vector<RefinementVertexKey> sourceKeys;
        if (!buildSourceKeys(sourceMesh, artifact, sourceKeys)) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::InvalidAuthoritativeSource, phase};
            return out;
        }
        std::map<Tet10SourceEdgeKey, EdgeWork> edges;
        phase = T4BPhase::EdgeExtraction;
        for (Index ci = 0; ci < sourceMesh.cellCount(); ++ci) {
            ++metrics.sourceCellsScanned;
            if (!appendCellEdges(sourceMesh, ci, sourceKeys, edges, metrics)) {
                out.status = T4BStatus::Unsupported;
                out.failure = {T4BReason::UnsupportedSourceCell, phase};
                return out;
            }
        }
        const auto bindingStatus = validateBoundaryBindings(
            edges, artifact, sourceMesh, sourceKeys,
            request.source.source.conformingRequest.expectedSurfaceSetStamp);
        if (bindingStatus != T4BReason::None) {
            out.status = T4BStatus::Rejected;
            out.failure = {bindingStatus, phase};
            return out;
        }
        metrics.uniqueEdges = edges.size();
        for (const auto &[key, row] : edges) {
            (void)key;
            if (!row.boundary)
                ++metrics.interiorEdges;
            else if (row.kind == EdgeWork::Kind::Curve)
                ++metrics.curveEdges;
            else
                ++metrics.surfaceEdges;
        }
        const std::uint64_t sourceNodes = static_cast<std::uint64_t>(sourceMesh.nodeCount());
        const std::uint64_t sourceCells = static_cast<std::uint64_t>(sourceMesh.cellCount());
        std::uint64_t requiredNodes = 0, requiredConn = 0, unattempted = 0;
        if (!checkedAdd(sourceNodes, metrics.uniqueEdges, requiredNodes)) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::ConnectivityOverflow, phase};
            return out;
        }
        for (Index ci = 0; ci < sourceMesh.cellCount(); ++ci) {
            std::uint64_t addCount = 0;
            switch (sourceMesh.cell(ci).type) {
            case mesh::CellType::Line2:
                addCount = 3;
                break;
            case mesh::CellType::Tri3:
                addCount = 6;
                break;
            case mesh::CellType::Tet4:
                addCount = 10;
                break;
            default:
                break;
            }
            if (!checkedAdd(requiredConn, addCount, requiredConn)) {
                out.status = T4BStatus::Rejected;
                out.failure = {T4BReason::ConnectivityOverflow, phase};
                return out;
            }
        }
        const std::uint64_t tetCount =
            static_cast<std::uint64_t>(sourceMesh.countCells(mesh::CellType::Tet4));
        std::vector<std::uint64_t> tetCellOrdinals;
        tetCellOrdinals.reserve(static_cast<std::size_t>(tetCount));
        for (Index ci = 0; ci < sourceMesh.cellCount(); ++ci)
            if (sourceMesh.cell(ci).type == mesh::CellType::Tet4)
                tetCellOrdinals.push_back(static_cast<std::uint64_t>(ci));
        if (tetCellOrdinals.size() != static_cast<std::size_t>(tetCount)) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::InvalidAuthoritativeSource, phase};
            return out;
        }
        std::uint64_t depthPower = 1;
        for (std::uint32_t i = 0; i <= request.policy.maxCertificateDepth; ++i)
            if (!checkedMul(depthPower, 8, depthPower)) {
                out.status = T4BStatus::Rejected;
                out.failure = {T4BReason::CapacityOverflow, phase};
                return out;
            }
        std::uint64_t maxVisited = 0, pending = 0;
        if (!checkedAdd(metrics.uniqueEdges, tetCount, unattempted) ||
            !checkedAdd(1, 7ull * request.policy.maxCertificateDepth, pending)) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::CapacityOverflow, phase};
            return out;
        }
        maxVisited = (depthPower - 1) / 7;
        (void)maxVisited;
        (void)pending;
        const auto indexMaximum = static_cast<std::uint64_t>((std::numeric_limits<Index>::max)());
        if (requiredNodes > indexMaximum || sourceCells > indexMaximum ||
            requiredConn > indexMaximum || requiredNodes > std::vector<Vec3>{}.max_size() ||
            sourceCells > std::vector<mesh::Cell>{}.max_size() ||
            requiredConn > std::vector<Index>{}.max_size() ||
            metrics.uniqueEdges > std::vector<T4BPlacementRow>{}.max_size() ||
            tetCount > std::vector<T4BCertificateRow>{}.max_size() ||
            tetCount > std::vector<std::uint64_t>{}.max_size() ||
            unattempted > std::vector<Tet10SourceEdgeKey>{}.max_size()) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::CapacityOverflow, phase};
            return out;
        }
        if (requiredNodes > request.policy.maxResultNodes ||
            sourceCells > request.policy.maxResultCells ||
            requiredConn > request.policy.maxResultConnectivity) {
            out.status = T4BStatus::Incomplete;
            out.failure = {T4BReason::ResultBudgetReached, phase};
            return out;
        }
        if (metrics.uniqueEdges > request.policy.maxPlacementRows ||
            tetCount > request.policy.maxCertificateRows ||
            unattempted > request.policy.maxUnattemptedRows) {
            out.status = T4BStatus::Incomplete;
            out.failure = {T4BReason::ReportBudgetReached, T4BPhase::CertificatePreparation};
            return out;
        }
        // Conservative accounted high-water for the additive T4-B stage.
        // It includes the carried T4-A bound, source/candidate/publication flat
        // buffers at their maximum overlap, edge-map nodes and vector payload,
        // both report/result rows, and fixed certificate DFS scratch.
        std::uint64_t accounted = sourceResult.metrics.transientAccountedHighWaterBytes;
        const std::uint64_t flatNodeWidth = sizeof(Vec3) + sizeof(mesh::Classification);
        const std::uint64_t edgeNodeWidth =
            sizeof(std::pair<const Tet10SourceEdgeKey, EdgeWork>) + 4 * sizeof(void *);
        const std::uint64_t pendingBytes = 141ull * (40ull * sizeof(Tet10Interval) +
                                                     sizeof(std::uint64_t) + sizeof(std::uint32_t));
        if (!addAccounted(accounted, requiredNodes, 3 * flatNodeWidth) ||
            !addAccounted(accounted, sourceCells, 3 * sizeof(mesh::Cell)) ||
            !addAccounted(accounted, requiredConn, 3 * sizeof(Index)) ||
            !addAccounted(accounted, metrics.uniqueEdges, edgeNodeWidth) ||
            !addAccounted(accounted, metrics.uniqueEdges,
                          2 * sizeof(T4BPlacementRow) + sizeof(Tet10SourceEdgeKey)) ||
            !addAccounted(accounted, tetCount,
                          2 * sizeof(T4BCertificateRow) + 2 * sizeof(std::uint64_t)) ||
            !addAccounted(accounted, sourceNodes, sizeof(RefinementVertexKey)) ||
            !addAccounted(accounted, sourceNodes, sizeof(Index)) ||
            !addAccounted(accounted, sourceCells, sizeof(Index)) ||
            !addAccounted(accounted, 1,
                          pendingBytes + 10 * sizeof(Index) + sizeof(T4BReport) +
                              sizeof(T4BResult))) {
            out.status = T4BStatus::Rejected;
            out.failure = {T4BReason::CapacityOverflow, T4BPhase::CertificatePreparation};
            return out;
        }
        for (const auto &[key, row] : edges) {
            (void)key;
            if (!addAccounted(accounted, row.incidentFacets.capacity(),
                              sizeof(bnd::VolumeFacetKey))) {
                out.status = T4BStatus::Rejected;
                out.failure = {T4BReason::CapacityOverflow, T4BPhase::CertificatePreparation};
                return out;
            }
        }
        metrics.accountedHighWaterBytes = accounted;
        phase = T4BPhase::CertificatePreparation;
        out.report->placements.reserve(static_cast<std::size_t>(metrics.uniqueEdges));
        out.report->certificates.reserve(static_cast<std::size_t>(tetCount));
        out.report->unattemptedEdges.reserve(static_cast<std::size_t>(metrics.uniqueEdges));
        out.report->unattemptedCells.reserve(static_cast<std::size_t>(tetCount));
        auto receipt = std::make_unique<T4BSourceReceipt>();
        auto working = std::make_unique<T4BResult>();
        auto candidate = std::make_unique<mesh::Mesh>();
        candidate->reserveAppend(static_cast<std::size_t>(requiredNodes),
                                 static_cast<std::size_t>(sourceCells),
                                 static_cast<std::size_t>(requiredConn));
        working->edgeToMidside.reserve(static_cast<std::size_t>(metrics.uniqueEdges));
        working->certificates.reserve(static_cast<std::size_t>(tetCount));
        working->sourceNodeToOutput.reserve(static_cast<std::size_t>(sourceNodes));
        working->sourceCellToOutput.reserve(static_cast<std::size_t>(sourceCells));
        for (Index i = 0; i < sourceMesh.nodeCount(); ++i)
            working->sourceNodeToOutput.push_back(i);
        for (Index i = 0; i < sourceMesh.cellCount(); ++i)
            working->sourceCellToOutput.push_back(i);
        std::vector<Index> cellScratch;
        cellScratch.reserve(10);
        auto allUnattempted = [&]() noexcept {
            for (const auto &[key, row] : edges) {
                (void)row;
                out.report->unattemptedEdges.push_back(key);
            }
            out.report->unattemptedCells.insert(out.report->unattemptedCells.end(),
                                                tetCellOrdinals.begin(), tetCellOrdinals.end());
            out.reportComplete = true;
        };
        auto allCellsUnattempted = [&]() noexcept {
            out.report->unattemptedCells.insert(out.report->unattemptedCells.end(),
                                                tetCellOrdinals.begin(), tetCellOrdinals.end());
            out.reportComplete = true;
        };
        std::uint64_t extractionWork = 0;
        if (!checkedAdd(metrics.sourceCellsScanned, metrics.sourceEdgeUses, extractionWork) ||
            !consumeWork(extractionWork, T4BPhase::EdgeExtraction)) {
            allUnattempted();
            return out;
        }
        if (!checkedStamp(phase)) {
            allUnattempted();
            return out;
        }
        // All allocations and immutable source checks have completed. Heavy
        // ownership now crosses through unique_ptr moves only.
        receipt->sourceOutcome = std::move(out.sourceOutcome->sourceOutcome);
        receipt->sourceReceipt = std::move(out.sourceOutcome->sourceReceipt);
        receipt->sourceReport = std::move(out.sourceOutcome->report);
        receipt->before = sourceResult.before;
        receipt->after = sourceResult.after;
        receipt->metrics = sourceResult.metrics;
        receipt->boundaryAudit = sourceResult.boundaryAudit;
        receipt->strictHash = sourceResult.finalStrictHash;
        receipt->semanticHash = sourceResult.finalSemanticHash;
        working->sourceArtifact = std::move(sourceResult.artifact);
        out.sourceOutcome.reset();
        out.sourceReceipt = std::move(receipt);
        for (Index i = 0; i < sourceMesh.nodeCount(); ++i)
            candidate->addNode(sourceMesh.node(i), sourceMesh.nodeClass(i));
        phase = T4BPhase::Placement;
        for (auto &[key, row] : edges) {
            T4BPlacementRow reportRow{};
            reportRow.edge = key;
            auto finalizePlacementStop = [&](T4BPlacementDisposition disposition) noexcept {
                reportRow.disposition = disposition;
                out.report->placements.push_back(reportRow);
                auto remaining = edges.find(key);
                if (remaining != edges.end())
                    ++remaining;
                for (; remaining != edges.end(); ++remaining)
                    out.report->unattemptedEdges.push_back(remaining->first);
                out.report->unattemptedCells.insert(out.report->unattemptedCells.end(),
                                                    tetCellOrdinals.begin(), tetCellOrdinals.end());
                out.reportComplete = true;
            };
            const Vec3 a = sourceMesh.node(row.source[0]), b = sourceMesh.node(row.source[1]);
            Vec3 coordinate{};
            if (!row.boundary) {
                coordinate = {std::midpoint(a.x, b.x), std::midpoint(a.y, b.y),
                              std::midpoint(a.z, b.z)};
                reportRow.disposition = T4BPlacementDisposition::InteriorMidpoint;
            } else {
                if (!checkedStamp(phase, key)) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    return out;
                }
                Tet10PlacementRequest placementRequest{
                    key.vertices,
                    {{a, b}},
                    row.classification,
                    row.kind == EdgeWork::Kind::Curve ? Tet10PlacementKind::Curve
                                                      : Tet10PlacementKind::Surface,
                    std::span<const bnd::VolumeFacetKey>(row.incidentFacets)};
                Tet10PlacementResult placed{};
                if (!consumeWork(1, phase, key)) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    return out;
                }
                try {
                    ++metrics.providerCalls;
                    placed = placementProvider.place(placementRequest);
                } catch (const std::bad_alloc &) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    out.status = T4BStatus::Failed;
                    out.failure = {T4BReason::AllocationFailed, phase, key};
                    return out;
                } catch (...) {
                    finalizePlacementStop(T4BPlacementDisposition::Failed);
                    out.status = T4BStatus::Failed;
                    out.failure = {T4BReason::PlacementProviderThrew, phase, key};
                    return out;
                }
                reportRow.providerStatus = placed.status;
                reportRow.providerReason = placed.providerReason;
                if (placed.observedStamp.persistentId !=
                    request.expectedPlacementStamp.persistentId) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    out.status = T4BStatus::Failed;
                    out.failure = {T4BReason::PlacementProviderIdentityChanged, phase, key};
                    return out;
                }
                if (placed.observedStamp.revision != request.expectedPlacementStamp.revision) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    out.status = T4BStatus::Failed;
                    out.failure = {T4BReason::StalePlacementProvider, phase, key};
                    return out;
                }
                if (!checkedStamp(phase, key)) {
                    finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                    return out;
                }
                if (placed.status != Tet10PlacementStatus::Placed) {
                    out.status =
                        placed.status == Tet10PlacementStatus::Unsupported
                            ? T4BStatus::Unsupported
                            : (placed.status == Tet10PlacementStatus::Rejected ? T4BStatus::Rejected
                                                                               : T4BStatus::Failed);
                    out.failure = {placed.status == Tet10PlacementStatus::Unsupported
                                       ? T4BReason::PlacementUnsupported
                                       : (placed.status == Tet10PlacementStatus::Rejected
                                              ? T4BReason::PlacementRejected
                                              : T4BReason::PlacementFailed),
                                   phase, key};
                    finalizePlacementStop(placed.status == Tet10PlacementStatus::Unsupported
                                              ? T4BPlacementDisposition::Unsupported
                                              : (placed.status == Tet10PlacementStatus::Rejected
                                                     ? T4BPlacementDisposition::Rejected
                                                     : T4BPlacementDisposition::Failed));
                    return out;
                }
                coordinate = placed.coordinate;
                reportRow.disposition = T4BPlacementDisposition::Placed;
            }
            if (!finite(coordinate)) {
                reportRow.coordinate = coordinate;
                finalizePlacementStop(T4BPlacementDisposition::Failed);
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::NonfinitePlacement, phase, key};
                return out;
            }
            reportRow.coordinate = coordinate;
            reportRow.classification = row.classification;
            if (!consumeWork(1, phase, key)) {
                finalizePlacementStop(T4BPlacementDisposition::TerminalFailure);
                return out;
            }
            row.midside = candidate->addNode(coordinate, row.classification);
            reportRow.midsideIndex = row.midside;
            out.report->placements.push_back(reportRow);
            working->edgeToMidside.push_back(reportRow);
        }
        if (!checkedStamp(phase)) {
            allCellsUnattempted();
            return out;
        }
        phase = T4BPhase::CertificatePreparation;
        try {
            if (!buildElevatedCells(sourceMesh, edges, sourceKeys, *candidate, cellScratch)) {
                allCellsUnattempted();
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::InvalidElevatedTopology, phase};
                return out;
            }
        } catch (const std::bad_alloc &) {
            allCellsUnattempted();
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::AllocationFailed, phase};
            return out;
        } catch (...) {
            allCellsUnattempted();
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::InvalidElevatedTopology, phase};
            return out;
        }
        working->privateMesh = std::move(candidate);
        working->metrics = metrics;
        out.retainedResult = std::move(working);
        phase = T4BPhase::JacobianCertification;
        Tet10CertificatePolicy certificatePolicy;
        certificatePolicy.minimumJacobianRatio = request.policy.minimumCertifiedJacobianRatio;
        certificatePolicy.maxDepth = request.policy.maxCertificateDepth;
        for (std::size_t tetPosition = 0; tetPosition < tetCellOrdinals.size(); ++tetPosition) {
            const auto sourceOrdinal = tetCellOrdinals[tetPosition];
            const Index ci = static_cast<Index>(sourceOrdinal);
            const auto cell = out.retainedResult->privateMesh->cell(ci);
            if (cell.type != mesh::CellType::Tet10) {
                out.retainedResult->metrics = metrics;
                out.status = T4BStatus::Failed;
                out.failure = {T4BReason::InvalidElevatedTopology, phase, {}, sourceOrdinal};
                return out;
            }
            std::array<Vec3, 10> coordinates{};
            const Index *n = out.retainedResult->privateMesh->cellNodes(ci);
            for (unsigned j = 0; j < 10; ++j)
                coordinates[j] = out.retainedResult->privateMesh->node(n[j]);
            certificatePolicy.maxVisitedSubcells =
                metrics.certificateSubcells >= request.policy.maxCertificateCells
                    ? 0
                    : request.policy.maxCertificateCells - metrics.certificateSubcells;
            certificatePolicy.maxCoefficientTests =
                metrics.coefficientTests >= request.policy.maxCoefficientTests
                    ? 0
                    : request.policy.maxCoefficientTests - metrics.coefficientTests;
            certificatePolicy.maxWorkUnits = metrics.workUnits >= request.policy.maxWorkUnits
                                                 ? 0
                                                 : request.policy.maxWorkUnits - metrics.workUnits;
            const auto cert = certifyTet10Jacobian(coordinates, certificatePolicy);
            T4BCertificateRow row{sourceOrdinal, cert.status, cert.witness, cert.metrics,
                                  cert.conservativeMinimumLowerBound};
            out.report->certificates.push_back(row);
            out.retainedResult->certificates.push_back(row);
            metrics.coefficientTests += cert.metrics.coefficientTests;
            metrics.certificateSubcells += cert.metrics.visitedSubcells;
            metrics.workUnits += cert.metrics.workUnits;
            metrics.maximumCertificateDepth =
                std::max(metrics.maximumCertificateDepth, cert.metrics.maximumDepth);
            if (cert.status != Tet10CertificateStatus::Certified) {
                out.retainedResult->metrics = metrics;
                out.report->unattemptedCells.insert(
                    out.report->unattemptedCells.end(),
                    std::next(tetCellOrdinals.begin(),
                              static_cast<std::ptrdiff_t>(tetPosition + 1)),
                    tetCellOrdinals.end());
                out.reportComplete = true;
                out.status = certificateStatus(cert.status);
                out.failure = {certificateReason(cert.status), phase, {}, sourceOrdinal};
                return out;
            }
        }
        out.retainedResult->certificationComplete = true;
        out.reportComplete = true;
        phase = T4BPhase::FinalAudit;
        if (!checkedAdd(static_cast<std::uint64_t>(out.retainedResult->privateMesh->nodeCount()),
                        static_cast<std::uint64_t>(out.retainedResult->privateMesh->cellCount()),
                        metrics.auditRows) ||
            metrics.auditRows > request.policy.maxWorkUnits ||
            metrics.workUnits > request.policy.maxWorkUnits - metrics.auditRows) {
            out.retainedResult->metrics = metrics;
            out.status = T4BStatus::Incomplete;
            out.failure = {T4BReason::JacobianCertificationBudgetReached, phase};
            return out;
        }
        metrics.workUnits += metrics.auditRows;
        mesh::AuditOptions options;
        options.boundaryMode = mesh::BoundaryComparisonMode::Required;
        options.boundaryOrphanBudget = 0;
        options.boundaryUnexpectedBudget = 0;
        out.retainedResult->audit = mesh::auditMesh(*out.retainedResult->privateMesh, options);
        if (out.retainedResult->audit.boundaryOrphans != 0 ||
            out.retainedResult->audit.boundaryUnexpected != 0) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::BoundaryBijectionFailed, phase};
            return out;
        }
        if (!out.retainedResult->audit.passed() || out.retainedResult->audit.nInverted != 0 ||
            out.retainedResult->audit.nInvalidTet10 != 0) {
            out.status = T4BStatus::Failed;
            out.failure = {T4BReason::InvalidElevatedTopology, phase};
            return out;
        }
        out.retainedResult->strictHash = mesh::meshHash(*out.retainedResult->privateMesh);
        out.retainedResult->semanticHash = mesh::semanticHashEx(*out.retainedResult->privateMesh);
        out.retainedResult->metrics = metrics;
        out.retainedResult->publicationReady = true;
        if (!checkedStamp(T4BPhase::PublicationPreparation)) {
            out.retainedResult->publicationReady = false;
            out.retainedResult->metrics = metrics;
            return out;
        }
        phase = T4BPhase::PublicationPreparation;
        std::uint64_t publicationWork = 0;
        if (!checkedAdd(requiredNodes, sourceCells, publicationWork) ||
            !checkedAdd(publicationWork, requiredConn, publicationWork) ||
            publicationWork > request.policy.maxWorkUnits ||
            metrics.workUnits > request.policy.maxWorkUnits - publicationWork) {
            out.retainedResult->publicationReady = false;
            out.retainedResult->metrics = metrics;
            out.status = T4BStatus::Incomplete;
            out.failure = {T4BReason::JacobianCertificationBudgetReached, phase};
            return out;
        }
        metrics.workUnits += publicationWork;
#if defined(PAACS_TOPO_TESTING)
        struct PublicationPreparationTrackingScope {
            bool previous = false;
            bool active = true;
            PublicationPreparationTrackingScope() noexcept
                : previous(t4b_detail::publicationPreparationActive()) {
                t4b_detail::publicationPreparationActive() = true;
            }
            void finish() noexcept {
                if (active) {
                    t4b_detail::publicationPreparationActive() = previous;
                    active = false;
                }
            }
            ~PublicationPreparationTrackingScope() { finish(); }
        } publicationPreparationTracking;
#endif
        // Tracking is active before the target object's constructor.  Its
        // destructor restores the marker on every exceptional exit.
        mesh::Mesh publication;
        cloneMesh(*out.retainedResult->privateMesh, publication);
#if defined(PAACS_TOPO_TESTING)
        publicationPreparationTracking.finish();
#endif
        phase = T4BPhase::Publication;
#if defined(PAACS_TOPO_TESTING)
        t4b_detail::publicationTailActive() = true;
#endif
        destination.swap(publication);
#if defined(PAACS_TOPO_TESTING)
        t4b_detail::publicationTailActive() = false;
#endif
        out.retainedResult->metrics = metrics;
        out.retainedResult->metrics.runComplete = true;
        out.status = T4BStatus::Published;
        out.destinationUnchanged = false;
        out.reportComplete = true;
        return out;
    } catch (const std::bad_alloc &) {
        out.status = T4BStatus::Failed;
        out.failure = {T4BReason::AllocationFailed, phase};
        return out;
    } catch (...) {
        out.status = T4BStatus::Failed;
        out.failure = {phase == T4BPhase::PublicationPreparation
                           ? T4BReason::PublicationPreparationFailed
                           : T4BReason::InternalError,
                       phase};
        return out;
    }
}

} // namespace paacs::algo
