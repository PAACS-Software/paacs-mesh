// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/bnd/TetRecoveryIndex.hpp"

#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <queue>
#include <set>

namespace paacs::bnd {
namespace {

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

TetAabb pointBounds(const Vec3& p) noexcept { return {p, p}; }

TetAabb merge(const TetAabb& a, const TetAabb& b) noexcept {
    return {{std::min(a.lo.x,b.lo.x), std::min(a.lo.y,b.lo.y),
             std::min(a.lo.z,b.lo.z)},
            {std::max(a.hi.x,b.hi.x), std::max(a.hi.y,b.hi.y),
             std::max(a.hi.z,b.hi.z)}};
}

bool overlaps(const TetAabb& a, const TetAabb& b) noexcept {
    return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x &&
           a.lo.y <= b.hi.y && b.lo.y <= a.hi.y &&
           a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
}

double axisLo(const TetAabb& a, std::size_t axis) noexcept {
    return axis == 0 ? a.lo.x : axis == 1 ? a.lo.y : a.lo.z;
}
double axisHi(const TetAabb& a, std::size_t axis) noexcept {
    return axis == 0 ? a.hi.x : axis == 1 ? a.hi.y : a.hi.z;
}

struct Leaf {
    topo::CellHandle cell{};
    TetAabb bounds{};
    std::uint64_t sequence = 0;
};

std::int32_t buildNode(std::vector<DynamicTetAabbIndex::Node>& nodes,
                       std::vector<Leaf>& leaves, std::size_t begin,
                       std::size_t end, std::size_t depth) {
    const std::int32_t result = static_cast<std::int32_t>(nodes.size());
    nodes.push_back({});
    if (end - begin == 1) {
        auto& node = nodes[static_cast<std::size_t>(result)];
        node.bounds = leaves[begin].bounds;
        node.cell = leaves[begin].cell;
        node.minimumSequence = leaves[begin].sequence;
        node.leaf = true;
        return result;
    }
    const std::size_t axis = depth % 3;
    std::stable_sort(leaves.begin() + static_cast<std::ptrdiff_t>(begin),
                     leaves.begin() + static_cast<std::ptrdiff_t>(end),
                     [axis](const Leaf& a, const Leaf& b) {
        if (axisLo(a.bounds,axis) != axisLo(b.bounds,axis))
            return axisLo(a.bounds,axis) < axisLo(b.bounds,axis);
        if (axisHi(a.bounds,axis) != axisHi(b.bounds,axis))
            return axisHi(a.bounds,axis) < axisHi(b.bounds,axis);
        return a.sequence < b.sequence;
    });
    const std::size_t middle = begin + (end - begin) / 2;
    const std::int32_t left = buildNode(nodes, leaves, begin, middle, depth+1);
    const std::int32_t right = buildNode(nodes, leaves, middle, end, depth+1);
    auto& node = nodes[static_cast<std::size_t>(result)];
    node.left = left;
    node.right = right;
    const auto& l = nodes[static_cast<std::size_t>(left)];
    const auto& r = nodes[static_cast<std::size_t>(right)];
    node.bounds = merge(l.bounds, r.bounds);
    node.minimumSequence = std::min(l.minimumSequence, r.minimumSequence);
    return result;
}

}  // namespace

void DynamicTetAabbIndex::build(
    const std::map<topo::CellHandle, CurrentTetRecord>& cells) {
    std::vector<Leaf> leaves;
    leaves.reserve(cells.size());
    for (const auto& [handle, cell] : cells)
        leaves.push_back({handle, cell.bounds, cell.sequence});
    std::sort(leaves.begin(), leaves.end(), [](const Leaf& a, const Leaf& b) {
        return a.sequence < b.sequence;
    });
    nodes_.clear();
    sortedLeafHandles_.clear();
    if (leaves.empty()) return;
    nodes_.reserve(leaves.size() * 2 - 1);
    (void)buildNode(nodes_, leaves, 0, leaves.size(), 0);
    sortedLeafHandles_.reserve(leaves.size());
    for (const auto& node : nodes_)
        if (node.leaf) sortedLeafHandles_.push_back(node.cell);
    std::sort(sortedLeafHandles_.begin(),sortedLeafHandles_.end());
}

std::vector<topo::CellHandle> DynamicTetAabbIndex::query(
    const TetAabb& bounds, std::uint64_t* nodesVisited) const {
    std::vector<topo::CellHandle> result;
    if (nodes_.empty()) return result;
    std::vector<std::int32_t> stack{0};
    while (!stack.empty()) {
        const std::int32_t index = stack.back();
        stack.pop_back();
        if (nodesVisited) ++*nodesVisited;
        const auto& node = nodes_[static_cast<std::size_t>(index)];
        if (!overlaps(node.bounds, bounds)) continue;
        if (node.leaf) {
            result.push_back(node.cell);
        } else {
            const auto& left = nodes_[static_cast<std::size_t>(node.left)];
            const auto& right = nodes_[static_cast<std::size_t>(node.right)];
            if (left.minimumSequence <= right.minimumSequence) {
                stack.push_back(node.right);
                stack.push_back(node.left);
            } else {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
    }
    return result;
}

bool DynamicTetAabbIndex::covers(
    const std::map<topo::CellHandle, CurrentTetRecord>& cells) const noexcept {
    if (cells.empty())
        return nodes_.empty() && sortedLeafHandles_.empty();
    std::size_t leaves = 0;
    for (const auto& node : nodes_)
        if (node.leaf) {
            ++leaves;
        }
    if (leaves != cells.size() ||
        sortedLeafHandles_.size() != cells.size()) return false;
    auto expected = cells.begin();
    for (const auto handle : sortedLeafHandles_) {
        if (expected == cells.end() || !(handle == expected->first))
            return false;
        ++expected;
    }
    if (expected != cells.end()) return false;
    return true;
}

std::vector<topo::CellHandle> TetRecoveryContext::querySegment(
    const Vec3& a, const Vec3& b) const {
    ++metrics.segmentQueries;
    const TetAabb bounds = merge(pointBounds(a), pointBounds(b));
    auto result = index.query(bounds,&metrics.indexNodesVisited);
    result.erase(std::remove_if(result.begin(), result.end(),
        [&](const auto cell) { return cells.count(cell) == 0; }), result.end());
    for (const auto cell : deltaCells) {
        const auto found = cells.find(cell);
        if (found != cells.end() && overlaps(found->second.bounds, bounds))
            result.push_back(cell);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    std::sort(result.begin(), result.end(), [&](const auto& x, const auto& y) {
        return cells.at(x).sequence < cells.at(y).sequence;
    });
    metrics.returnedCandidates += result.size();
    return result;
}

std::vector<topo::CellHandle> TetRecoveryContext::queryTriangle(
    const std::array<Vec3, 3>& triangle) const {
    ++metrics.facetQueries;
    TetAabb bounds = pointBounds(triangle[0]);
    bounds = merge(bounds, pointBounds(triangle[1]));
    bounds = merge(bounds, pointBounds(triangle[2]));
    auto result = index.query(bounds,&metrics.indexNodesVisited);
    result.erase(std::remove_if(result.begin(), result.end(),
        [&](const auto cell) { return cells.count(cell) == 0; }), result.end());
    for (const auto cell : deltaCells) {
        const auto found = cells.find(cell);
        if (found != cells.end() && overlaps(found->second.bounds, bounds))
            result.push_back(cell);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    std::sort(result.begin(), result.end(), [&](const auto& x, const auto& y) {
        return cells.at(x).sequence < cells.at(y).sequence;
    });
    metrics.returnedCandidates += result.size();
    return result;
}

std::vector<topo::CellHandle> TetRecoveryContext::queryBounds(
    const TetAabb& bounds) const {
    ++metrics.boundsQueries;
    auto result = index.query(bounds,&metrics.indexNodesVisited);
    result.erase(std::remove_if(result.begin(), result.end(),
        [&](const auto cell) { return cells.count(cell) == 0; }), result.end());
    for (const auto cell : deltaCells) {
        const auto found = cells.find(cell);
        if (found != cells.end() && overlaps(found->second.bounds,bounds))
            result.push_back(cell);
    }
    std::sort(result.begin(),result.end());
    result.erase(std::unique(result.begin(),result.end()),result.end());
    std::sort(result.begin(),result.end(),[&](const auto& a,const auto& b) {
        return cells.at(a).sequence < cells.at(b).sequence;
    });
    metrics.returnedCandidates += result.size();
    return result;
}

bool TetRecoveryContext::currentFor(
    const topo::TopologyWorkspace& workspace) const noexcept {
    if (poisoned) return false;
    const auto current = topo::WorkspaceStamp::of(workspace);
    return stamp == current && contextRevision == current.revision;
}

TetRecoveryContextBuildOutcome buildTetRecoveryContext(
    const topo::TopologyWorkspace& workspace,
    mesh::Classification temporaryRegion) noexcept {
    TetRecoveryContextBuildOutcome outcome;
    try {
        if (temporaryRegion.dim != 3 || temporaryRegion.tag < 0) {
            outcome.reason = TetRecoveryContextReason::InvalidClassification;
            return outcome;
        }
        auto context = std::make_unique<TetRecoveryContext>();
        context->stamp = topo::WorkspaceStamp::of(workspace);
        context->contextRevision = context->stamp.revision;
        context->region = temporaryRegion;
        context->denseToNode = workspace.nodes().activeBySequence();
        if (context->denseToNode.size() >= kInvalidDenseTetNode) {
            outcome.reason = TetRecoveryContextReason::InvalidConnectivity;
            return outcome;
        }
        context->nodeUseCount.assign(context->denseToNode.size(), 0);
        context->nodeStars.resize(context->denseToNode.size());
        context->plainVertices.reserve(context->denseToNode.size());
        std::set<std::uint64_t> nodeSequences;
        for (std::size_t i = 0; i < context->denseToNode.size(); ++i) {
            const auto handle = context->denseToNode[i];
            const auto data = workspace.nodes().get(handle);
            if (!finite(data.coord)) {
                outcome.reason = TetRecoveryContextReason::NonfiniteNode;
                return outcome;
            }
            if (!nodeSequences.insert(workspace.nodes().seqOf(handle)).second) {
                outcome.reason = TetRecoveryContextReason::DuplicateNodeSequence;
                return outcome;
            }
            context->nodeToDense.emplace(handle,
                static_cast<DenseTetNodeId>(i));
            context->plainVertices.push_back(
                {data.coord,workspace.nodes().seqOf(handle)});
            ++context->metrics.auditVisits;
        }
        std::set<std::array<DenseTetNodeId,4>> cellKeys;
        for (const auto handle : workspace.cells().activeBySequence()) {
            const auto data = workspace.cells().get(handle);
            if (data.type != topo::ElementType::Tet4) {
                outcome.reason = TetRecoveryContextReason::UnsupportedCell;
                return outcome;
            }
            if (!sameClass(data.cls, temporaryRegion)) {
                outcome.reason = TetRecoveryContextReason::InvalidClassification;
                return outcome;
            }
            CurrentTetRecord record;
            record.handle = handle;
            record.sequence = workspace.cells().seqOf(handle);
            std::array<Vec3,4> points{};
            record.bounds.lo = {std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()};
            record.bounds.hi = {-std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()};
            for (std::size_t i = 0; i < 4; ++i) {
                const auto found = context->nodeToDense.find(data.nodes[i]);
                if (found == context->nodeToDense.end()) {
                    outcome.reason = TetRecoveryContextReason::InvalidConnectivity;
                    return outcome;
                }
                record.nodes[i] = found->second;
                ++context->nodeUseCount[found->second];
                context->nodeStars[found->second].push_back(handle);
                points[i] = workspace.nodes().get(data.nodes[i]).coord;
                record.bounds = i == 0 ? pointBounds(points[i])
                                       : merge(record.bounds, pointBounds(points[i]));
            }
            auto key = record.nodes;
            std::sort(key.begin(), key.end());
            if (std::adjacent_find(key.begin(), key.end()) != key.end() ||
                !cellKeys.insert(key).second) {
                outcome.reason = TetRecoveryContextReason::DuplicateTet;
                return outcome;
            }
            if (predicates::orient3d(points[0], points[1], points[2], points[3]) !=
                predicates::Sign::Positive) {
                outcome.reason = TetRecoveryContextReason::NonPositiveTet;
                return outcome;
            }
            record.classification = data.cls;
            context->cells.emplace(handle, record);
            ++context->metrics.auditVisits;
        }
        for (std::uint64_t uses : context->nodeUseCount)
            if (uses == 0) {
                outcome.reason = TetRecoveryContextReason::OrphanNode;
                return outcome;
            }
        for (const auto edge : workspace.edges().activeBySequence()) {
            const auto data = workspace.edges().get(edge);
            const auto a = context->nodeToDense.find(data.a);
            const auto b = context->nodeToDense.find(data.b);
            if (a == context->nodeToDense.end() ||
                b == context->nodeToDense.end()) {
                outcome.reason = TetRecoveryContextReason::MissingEdgeIncidence;
                return outcome;
            }
            context->edges.emplace(DenseTetEdgeKey::of(a->second,b->second),
                EdgeResolverRow{edge, {}, data.cls, data.requiredConstraint});
        }
        for (const auto& [cell, record] : context->cells) {
            (void)record;
            if (!workspace.edges().hasCellEdges(workspace.cells(), cell)) {
                outcome.reason = TetRecoveryContextReason::MissingEdgeIncidence;
                return outcome;
            }
            for (const auto edge : workspace.edges().cellEdges(
                     workspace.cells(), cell)) {
                const auto data = workspace.edges().get(edge);
                const auto key = DenseTetEdgeKey::of(
                    context->nodeToDense.at(data.a),
                    context->nodeToDense.at(data.b));
                auto found = context->edges.find(key);
                if (found == context->edges.end()) {
                    outcome.reason = TetRecoveryContextReason::MissingEdgeIncidence;
                    return outcome;
                }
                found->second.completeStar.push_back(cell);
            }
        }
        for (auto& [key,row] : context->edges) {
            (void)key;
            std::sort(row.completeStar.begin(), row.completeStar.end(),
                [&](const auto& a, const auto& b) {
                    return context->cells.at(a).sequence <
                           context->cells.at(b).sequence;
                });
        }
        for (const auto face : workspace.faces().activeBySequence()) {
            const auto data = workspace.faces().get(face);
            if (data.cornerCount != 3 || data.incident.empty() ||
                data.incident.size() > 2) {
                outcome.reason = TetRecoveryContextReason::NonManifoldFace;
                return outcome;
            }
            std::array<DenseTetNodeId,3> dense{};
            for (std::size_t i = 0; i < 3; ++i)
                dense[i] = context->nodeToDense.at(data.corners[i]);
            FaceResolverRow row;
            row.handle = face;
            row.classification = data.cls;
            for (const auto& owner : data.incident) {
                topo::CellHandle cell{owner.cellSlot, owner.cellGen,
                                      workspace.cells().poolId()};
                if (context->cells.count(cell) == 0 ||
                    (owner.sense != -1 && owner.sense != 1)) {
                    outcome.reason = TetRecoveryContextReason::InvalidFaceSense;
                    return outcome;
                }
                row.owners.push_back({cell, owner.sense});
            }
            if (row.owners.size() == 2 &&
                row.owners[0].sense * row.owners[1].sense != -1) {
                outcome.reason = TetRecoveryContextReason::InvalidFaceSense;
                return outcome;
            }
            context->faces.emplace(DenseTetFaceKey::of(dense), std::move(row));
        }
        for (const auto& [cell, record] : context->cells) {
            (void)record;
            if (!workspace.faces().hasCellFaces(workspace.cells(), cell)) {
                outcome.reason = TetRecoveryContextReason::MissingFaceIncidence;
                return outcome;
            }
        }
        for (const auto& [key,row] : context->edges)
            if (row.required)
                context->protectedEdges.push_back(
                    canonicalTetFacetEdge(key.nodes[0],key.nodes[1]));
        for (const auto& [key,row] : context->faces)
            if (row.classification.dim == 2 && row.classification.tag >= 0)
                context->protectedFaces.push_back({
                    TetFacetFaceKey{key.nodes},
                    static_cast<std::uint64_t>(row.handle.slot)+1U,
                    row.classification});
        context->index.build(context->cells);
        if (!context->index.covers(context->cells)) {
            outcome.reason = TetRecoveryContextReason::IndexCoverageMismatch;
            return outcome;
        }
        context->metrics.nodes = context->denseToNode.size();
        context->metrics.cells = context->cells.size();
        context->metrics.edges = context->edges.size();
        context->metrics.faces = context->faces.size();
        context->metrics.indexNodes = context->index.nodes().size();
        context->metrics.indexBytes = context->index.accountedBytes();
        context->metrics.runComplete = true;
        outcome.reason = TetRecoveryContextReason::None;
        outcome.context = std::move(context);
        return outcome;
    } catch (const std::bad_alloc&) {
        outcome.reason = TetRecoveryContextReason::AllocationFailed;
        return outcome;
    } catch (...) {
        outcome.reason = TetRecoveryContextReason::InternalError;
        return outcome;
    }
}

TetRecoveryContextUpdateOutcome updateTetRecoveryContext(
    TetRecoveryContext& context, const topo::TopologyWorkspace& workspace,
    std::span<const topo::CellHandle> retired,
    std::span<const topo::CellHandle> created) noexcept {
    TetRecoveryContextUpdateOutcome outcome;
    try {
        if (context.poisoned) {
            outcome.reason = TetRecoveryContextReason::StaleWorkspace;
            return outcome;
        }
        std::set<topo::CellHandle> retiredSet;
        std::set<DenseTetEdgeKey> changedEdges;
        std::set<DenseTetFaceKey> changedFaces;
        std::set<DenseTetNodeId> changedNodes;
        for (const auto cell : retired) {
            ++context.metrics.indexRemovals;
            const auto found = context.cells.find(cell);
            if (found == context.cells.end() || !retiredSet.insert(cell).second) {
                outcome.reason = TetRecoveryContextReason::InvalidConnectivity;
                context.poisoned = true;
                return outcome;
            }
            const auto& nodes = found->second.nodes;
            static constexpr std::array<std::array<int,2>,6> edges{{
                {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
            static constexpr std::array<std::array<int,3>,4> faces{{
                {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
            for (const auto& edge : edges) {
                const auto key = DenseTetEdgeKey::of(
                    nodes[static_cast<std::size_t>(edge[0])],
                    nodes[static_cast<std::size_t>(edge[1])]);
                changedEdges.insert(key);
                auto row = context.edges.find(key);
                if (row == context.edges.end()) throw std::logic_error("edge delta");
                auto& star = row->second.completeStar;
                star.erase(std::remove(star.begin(), star.end(), cell), star.end());
            }
            for (const auto& face : faces) {
                const auto key = DenseTetFaceKey::of({
                    nodes[static_cast<std::size_t>(face[0])],
                    nodes[static_cast<std::size_t>(face[1])],
                    nodes[static_cast<std::size_t>(face[2])]});
                changedFaces.insert(key);
                auto row = context.faces.find(key);
                if (row == context.faces.end()) throw std::logic_error("face delta");
                auto& owners = row->second.owners;
                owners.erase(std::remove_if(owners.begin(), owners.end(),
                    [&](const auto& owner) { return owner.cell == cell; }),
                    owners.end());
            }
            for (const auto node : nodes) {
                changedNodes.insert(node);
                if (context.nodeUseCount[node] == 0)
                    throw std::logic_error("node-use delta");
                --context.nodeUseCount[node];
                auto& star = context.nodeStars[node];
                star.erase(std::remove(star.begin(),star.end(),cell),star.end());
            }
        }
        for (const auto cell : retired) context.cells.erase(cell);

        for (const auto cell : created) {
            ++context.metrics.indexInsertions;
            if (!workspace.cells().isCurrent(cell) || context.cells.count(cell))
                throw std::logic_error("created cell identity");
            const auto data = workspace.cells().get(cell);
            if (data.type != topo::ElementType::Tet4 ||
                !sameClass(data.cls, context.region))
                throw std::logic_error("created cell payload");
            CurrentTetRecord record;
            record.handle = cell;
            record.sequence = workspace.cells().seqOf(cell);
            std::array<Vec3,4> points{};
            for (std::size_t i = 0; i < 4; ++i) {
                const auto found = context.nodeToDense.find(data.nodes[i]);
                if (found == context.nodeToDense.end())
                    throw std::logic_error("created node identity");
                record.nodes[i] = found->second;
                points[i] = workspace.nodes().get(data.nodes[i]).coord;
                record.bounds = i == 0 ? pointBounds(points[i])
                                       : merge(record.bounds, pointBounds(points[i]));
                ++context.nodeUseCount[found->second];
                context.nodeStars[found->second].push_back(cell);
            }
            auto unique = record.nodes;
            std::sort(unique.begin(), unique.end());
            if (std::adjacent_find(unique.begin(), unique.end()) != unique.end() ||
                predicates::orient3d(points[0],points[1],points[2],points[3]) !=
                    predicates::Sign::Positive)
                throw std::logic_error("created cell geometry");
            record.classification = data.cls;
            context.cells.emplace(cell, record);
            context.deltaCells.push_back(cell);
            context.metrics.deltaHighWater = std::max<std::uint64_t>(
                context.metrics.deltaHighWater,context.deltaCells.size());

            if (!workspace.edges().hasCellEdges(workspace.cells(), cell) ||
                !workspace.faces().hasCellFaces(workspace.cells(), cell))
                throw std::logic_error("created incidence");
            for (const auto edge : workspace.edges().cellEdges(
                     workspace.cells(), cell)) {
                const auto edgeData = workspace.edges().get(edge);
                const auto key = DenseTetEdgeKey::of(
                    context.nodeToDense.at(edgeData.a),
                    context.nodeToDense.at(edgeData.b));
                changedEdges.insert(key);
                auto& row = context.edges[key];
                row.handle = edge;
                row.classification = edgeData.cls;
                row.required = edgeData.requiredConstraint;
                row.completeStar.push_back(cell);
            }
            for (const auto face : workspace.faces().cellFaces(
                     workspace.cells(), cell)) {
                const auto faceData = workspace.faces().get(face);
                std::array<DenseTetNodeId,3> dense{};
                for (std::size_t i = 0; i < 3; ++i)
                    dense[i] = context.nodeToDense.at(faceData.corners[i]);
                const auto key = DenseTetFaceKey::of(dense);
                changedFaces.insert(key);
                auto& row = context.faces[key];
                row.handle = face;
                row.classification = faceData.cls;
                auto owner = std::find_if(
                    faceData.incident.begin(), faceData.incident.end(),
                    [&](const auto& incidence) {
                        return incidence.cellSlot == cell.slot &&
                               incidence.cellGen == cell.gen;
                    });
                if (owner == faceData.incident.end())
                    throw std::logic_error("created face owner");
                row.owners.push_back({cell, owner->sense});
            }
            for (const auto node : record.nodes) changedNodes.insert(node);
        }

        for (const auto& key : changedEdges) {
            auto it = context.edges.find(key);
            if (it == context.edges.end())
                throw std::logic_error("changed edge missing");
            auto& star = it->second.completeStar;
            std::sort(star.begin(), star.end(), [&](const auto& a, const auto& b) {
                return context.cells.at(a).sequence < context.cells.at(b).sequence;
            });
            star.erase(std::unique(star.begin(), star.end()), star.end());
            if (star.empty()) {
                context.edges.erase(it);
                continue;
            }
            const auto a = context.denseToNode[it->first.nodes[0]];
            const auto b = context.denseToNode[it->first.nodes[1]];
            const auto current = workspace.edges().find(workspace.nodes(), a, b);
            if (!current) throw std::logic_error("updated edge missing");
            const auto data = workspace.edges().get(*current);
            it->second.handle = *current;
            it->second.classification = data.cls;
            it->second.required = data.requiredConstraint;
        }
        for (const auto& key : changedFaces) {
            auto it = context.faces.find(key);
            if (it == context.faces.end())
                throw std::logic_error("changed face missing");
            auto& owners = it->second.owners;
            std::sort(owners.begin(), owners.end(), [&](const auto& a, const auto& b) {
                return context.cells.at(a.cell).sequence <
                       context.cells.at(b.cell).sequence;
            });
            if (owners.empty()) {
                context.faces.erase(it);
                continue;
            }
            if (owners.size() > 2 ||
                (owners.size() == 2 && owners[0].sense * owners[1].sense != -1))
                throw std::logic_error("updated face manifold");
            const std::array<topo::NodeHandle,3> nodes{
                context.denseToNode[it->first.nodes[0]],
                context.denseToNode[it->first.nodes[1]],
                context.denseToNode[it->first.nodes[2]]};
            const auto current = workspace.faces().find(workspace.nodes(), nodes);
            if (!current) throw std::logic_error("updated face missing");
            it->second.handle = *current;
            it->second.classification = workspace.faces().get(*current).cls;
        }
        for (const auto node : changedNodes) {
            if (context.nodeUseCount[node] == 0)
                throw std::logic_error("updated orphan");
            auto& star = context.nodeStars[node];
            std::sort(star.begin(),star.end(),[&](const auto& a,const auto& b) {
                return context.cells.at(a).sequence < context.cells.at(b).sequence;
            });
            star.erase(std::unique(star.begin(),star.end()),star.end());
        }

        const std::size_t rebuildThreshold =
            std::max<std::size_t>(64, context.cells.size() / 4);
        if (context.deltaCells.size() > rebuildThreshold) {
            context.index.build(context.cells);
            context.deltaCells.clear();
            outcome.rebuilt = true;
            ++context.metrics.indexRebuilds;
        }
        context.stamp = topo::WorkspaceStamp::of(workspace);
        context.contextRevision = context.stamp.revision;
        context.metrics.cells = context.cells.size();
        context.metrics.edges = context.edges.size();
        context.metrics.faces = context.faces.size();
        outcome.reason = TetRecoveryContextReason::None;
        return outcome;
    } catch (const std::bad_alloc&) {
        context.poisoned = true;
        outcome.reason = TetRecoveryContextReason::AllocationFailed;
        return outcome;
    } catch (...) {
        context.poisoned = true;
        outcome.reason = TetRecoveryContextReason::InternalError;
        return outcome;
    }
}

}  // namespace paacs::bnd
