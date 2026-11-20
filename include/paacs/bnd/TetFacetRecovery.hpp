// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T2-A -- protected-edge recovery of one prescribed interior Tri3 facet.
// Design authority: docs/T2A_DESIGN.md revision 2 (accepted 2026-08-03).

#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/topo/TetCavity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace paacs::bnd {

struct TetFacetRequest {
    topo::WorkspaceStamp stamp{};
    std::uint64_t entityId = 0;
    std::array<topo::NodeHandle, 3> vertices{};
    mesh::Classification classification{};
};

struct PersistentTetFacet {
    std::uint64_t entityId = 0;
    std::array<topo::NodeHandle, 3> canonicalVertices{};
    mesh::Classification classification{};
    topo::FaceHandle recoveredFace{};
};

struct TetFacetConstraintState {
    topo::WorkspaceStamp stamp{};
    std::map<std::uint64_t, PersistentTetFacet> facets;

    static TetFacetConstraintState forWorkspace(
        const topo::TopologyWorkspace& workspace) {
        TetFacetConstraintState state;
        state.stamp = topo::WorkspaceStamp::of(workspace);
        return state;
    }
    void swap(TetFacetConstraintState& other) noexcept {
        using std::swap;
        swap(stamp, other.stamp);
        facets.swap(other.facets);
    }
};

static_assert(std::is_nothrow_swappable_v<topo::WorkspaceStamp>);
static_assert(std::is_nothrow_swappable_v<
              decltype(TetFacetConstraintState::facets)>);

enum class TetFacetRecoveryStatus : std::uint8_t {
    AlreadyPresent,
    Marked,
    Recovered,
    Refused,
    Rejected,
    Failed
};

enum class TetFacetRecoveryPhase : std::uint8_t {
    None,
    Preflight,
    Routing,
    Discovery,
    Candidate,
    Preparation,
    StateCopy
};

enum class TetFacetRecoveryReason : std::uint8_t {
    None,
    StaleWorkspaceStamp,
    InvalidIdentity,
    InvalidClassification,
    NonfiniteCoordinate,
    DegenerateTarget,
    UnsupportedCell,
    NonPositiveTet,
    MixedRegion,
    MissingProtectedTargetEdge,
    InvalidEdgeClassification,
    InvalidFaceClassification,
    StateKeyMismatch,
    StaleBinding,
    WrongBindingFace,
    WrongBindingClassification,
    IdRebinding,
    TargetAlias,
    FaceAlias,
    OneSidedTargetUnsupported,
    InvalidTargetIncidence,
    ProtectedFacetConflict,
    CrossedProtectedEdge,
    CrossedProtectedFacet,
    BoundaryEdgeStar,
    NonManifoldEdgeStar,
    DisconnectedEdgeStar,
    MixedRegionStar,
    UnsupportedDegeneracy,
    NonBallClosure,
    UnsupportedNonStarCavity,
    InvalidCandidate,
    PreparationConflict,
    AllocationFailed,
    // T2-B set-scheduler outcomes. The accepted single-facet entry never
    // emits either value.
    PendingConstraint,
    RecoveryConflictCycle
};

struct TetFacetRecoveryOutcome {
    TetFacetRecoveryStatus status = TetFacetRecoveryStatus::Failed;
    TetFacetRecoveryReason reason = TetFacetRecoveryReason::None;
    TetFacetRecoveryPhase phase = TetFacetRecoveryPhase::None;
    topo::FaceHandle face{};
    TetFacetGeometryReason geometryReason = TetFacetGeometryReason::None;
    TetFacetGeometryCounters counters{};
    bool capacityResiduePermitted = false;

    bool succeeded() const noexcept {
        return status == TetFacetRecoveryStatus::AlreadyPresent ||
               status == TetFacetRecoveryStatus::Marked ||
               status == TetFacetRecoveryStatus::Recovered;
    }
};

static_assert(std::is_standard_layout_v<TetFacetRecoveryOutcome>);
static_assert(std::is_trivially_copyable_v<TetFacetRecoveryOutcome>);

namespace tet_facet_recovery_detail {

inline bool sameClass(const mesh::Classification& a,
                      const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}

inline bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

inline bool defaultClass(const mesh::Classification& value) noexcept {
    const mesh::Classification unset{};
    return sameClass(value, unset);
}

inline topo::WorkspaceStamp currentStamp(
    const topo::TopologyWorkspace& workspace) {
    return topo::WorkspaceStamp::of(workspace);
}

inline topo::WorkspaceStamp nextStamp(const topo::WorkspaceStamp& current) {
    topo::WorkspaceStamp next = current;
    ++next.revision;
    return next;
}

inline std::array<topo::NodeHandle, 3> canonicalTriple(
    const topo::NodePool& nodes,
    std::array<topo::NodeHandle, 3> vertices) {
    std::sort(vertices.begin(), vertices.end(),
              [&](topo::NodeHandle a, topo::NodeHandle b) {
                  return nodes.seqOf(a) < nodes.seqOf(b);
              });
    return vertices;
}

inline bool sameTriple(const std::array<topo::NodeHandle, 3>& a,
                       const std::array<topo::NodeHandle, 3>& b) noexcept {
    return a == b;
}

inline TetFacetRecoveryOutcome rejected(TetFacetRecoveryReason reason) {
    TetFacetRecoveryOutcome result;
    result.status = TetFacetRecoveryStatus::Rejected;
    result.reason = reason;
    result.phase = TetFacetRecoveryPhase::Preflight;
    return result;
}

inline TetFacetRecoveryOutcome refused(
    TetFacetRecoveryReason reason, TetFacetRecoveryPhase phase,
    TetFacetGeometryReason geometry = TetFacetGeometryReason::None,
    const TetFacetGeometryCounters& counters = {}) {
    TetFacetRecoveryOutcome result;
    result.status = TetFacetRecoveryStatus::Refused;
    result.reason = reason;
    result.phase = phase;
    result.geometryReason = geometry;
    result.counters = counters;
    return result;
}

inline TetFacetRecoveryReason geometryRefusalReason(
    TetFacetGeometryReason reason) noexcept {
    switch (reason) {
        case TetFacetGeometryReason::CrossedProtectedEdge:
            return TetFacetRecoveryReason::CrossedProtectedEdge;
        case TetFacetGeometryReason::CrossedProtectedFace:
        case TetFacetGeometryReason::ProtectedFaceLost:
            return TetFacetRecoveryReason::CrossedProtectedFacet;
        case TetFacetGeometryReason::BoundaryEdgeStar:
            return TetFacetRecoveryReason::BoundaryEdgeStar;
        case TetFacetGeometryReason::NonManifoldEdgeStar:
            return TetFacetRecoveryReason::NonManifoldEdgeStar;
        case TetFacetGeometryReason::DisconnectedEdgeStar:
            return TetFacetRecoveryReason::DisconnectedEdgeStar;
        case TetFacetGeometryReason::MixedRegionStar:
            return TetFacetRecoveryReason::MixedRegionStar;
        case TetFacetGeometryReason::NonBallClosure:
            return TetFacetRecoveryReason::NonBallClosure;
        case TetFacetGeometryReason::UnsupportedNonStarCavity:
            return TetFacetRecoveryReason::UnsupportedNonStarCavity;
        default:
            return TetFacetRecoveryReason::UnsupportedDegeneracy;
    }
}

struct ComplexStorage {
    std::vector<topo::NodeHandle> nodeHandles;
    std::vector<topo::CellHandle> cellHandles;
    std::vector<TetFacetVertex> vertices;
    std::vector<TetFacetCell> cells;
    std::vector<TetFacetEdgeKey> protectedEdges;
    std::vector<TetFacetProtectedFace> protectedFaces;
    TetFacetTriangle target{};
    mesh::Classification region{};

    TetFacetComplexView view() const {
        return {vertices, cells, target, protectedEdges, protectedFaces,
                region};
    }
};

inline TetFacetRecoveryReason auditPersistentState(
    const topo::TopologyWorkspace& workspace,
    const TetFacetConstraintState& state,
    const TetFacetRequest& request,
    const std::array<topo::NodeHandle, 3>& requestKey) {
    const auto current = currentStamp(workspace);
    if (!(state.stamp == current) || !(request.stamp == current))
        return TetFacetRecoveryReason::StaleWorkspaceStamp;

    std::set<std::array<topo::NodeHandle, 3>> triples;
    std::set<topo::FaceHandle> faces;
    for (const auto& [mapId, binding] : state.facets) {
        if (mapId != binding.entityId)
            return TetFacetRecoveryReason::StateKeyMismatch;
        for (const auto vertex : binding.canonicalVertices)
            if (!workspace.nodes().isCurrent(vertex))
                return TetFacetRecoveryReason::StaleBinding;
        if (canonicalTriple(workspace.nodes(), binding.canonicalVertices) !=
            binding.canonicalVertices)
            return TetFacetRecoveryReason::WrongBindingFace;
        if (!workspace.faces().isCurrent(binding.recoveredFace))
            return TetFacetRecoveryReason::StaleBinding;
        const auto face = workspace.faces().get(binding.recoveredFace);
        if (face.cornerCount != 3)
            return TetFacetRecoveryReason::WrongBindingFace;
        std::array<topo::NodeHandle, 3> faceNodes{
            face.corners[0], face.corners[1], face.corners[2]};
        if (canonicalTriple(workspace.nodes(), faceNodes) !=
            binding.canonicalVertices)
            return TetFacetRecoveryReason::WrongBindingFace;
        if (!sameClass(face.cls, binding.classification))
            return TetFacetRecoveryReason::WrongBindingClassification;
        if (!faces.insert(binding.recoveredFace).second)
            return TetFacetRecoveryReason::FaceAlias;
        if (!triples.insert(binding.canonicalVertices).second)
            return TetFacetRecoveryReason::TargetAlias;
    }

    const auto id = state.facets.find(request.entityId);
    if (id != state.facets.end()) {
        if (!sameTriple(id->second.canonicalVertices, requestKey) ||
            !sameClass(id->second.classification, request.classification))
            return TetFacetRecoveryReason::IdRebinding;
    }
    for (const auto& [entity, binding] : state.facets)
        if (entity != request.entityId &&
            sameTriple(binding.canonicalVertices, requestKey))
            return TetFacetRecoveryReason::TargetAlias;
    return TetFacetRecoveryReason::None;
}

inline TetFacetRecoveryReason buildComplex(
    const topo::TopologyWorkspace& workspace, const TetFacetRequest& request,
    ComplexStorage& output) {
    output.nodeHandles = workspace.nodes().activeBySequence();
    // A dense ID may use [0, kInvalidTetFacetVertex); the sentinel itself
    // must never be produced by the size_t-to-ID conversion below.
    if (output.nodeHandles.size() >
        static_cast<std::size_t>(kInvalidTetFacetVertex))
        return TetFacetRecoveryReason::InvalidIdentity;
    output.vertices.reserve(output.nodeHandles.size());
    std::map<topo::NodeHandle, TetFacetVertexId> nodeIndex;
    for (std::size_t i = 0; i < output.nodeHandles.size(); ++i) {
        const auto handle = output.nodeHandles[i];
        const auto node = workspace.nodes().get(handle);
        if (!finite(node.coord))
            return TetFacetRecoveryReason::NonfiniteCoordinate;
        const auto id = static_cast<TetFacetVertexId>(i);
        nodeIndex.emplace(handle, id);
        output.vertices.push_back(
            TetFacetVertex{node.coord, workspace.nodes().seqOf(handle)});
    }
    for (std::size_t i = 0; i < 3; ++i) {
        const auto found = nodeIndex.find(request.vertices[i]);
        if (found == nodeIndex.end())
            return TetFacetRecoveryReason::InvalidIdentity;
        output.target.vertices[i] = found->second;
    }

    output.cellHandles = workspace.cells().activeBySequence();
    output.cells.reserve(output.cellHandles.size());
    bool haveRegion = false;
    for (const auto handle : output.cellHandles) {
        const auto cell = workspace.cells().get(handle);
        if (cell.type != topo::ElementType::Tet4 || cell.count != 4)
            return TetFacetRecoveryReason::UnsupportedCell;
        if (cell.cls.dim != 3 || cell.cls.tag < 0)
            return TetFacetRecoveryReason::MixedRegion;
        if (!haveRegion) {
            output.region = cell.cls;
            haveRegion = true;
        } else if (!sameClass(output.region, cell.cls)) {
            return TetFacetRecoveryReason::MixedRegion;
        }
        TetFacetCell plain;
        plain.sequence = workspace.cells().seqOf(handle);
        plain.classification = cell.cls;
        for (std::size_t i = 0; i < 4; ++i) {
            const auto found = nodeIndex.find(cell.nodes[i]);
            if (found == nodeIndex.end())
                return TetFacetRecoveryReason::InvalidIdentity;
            plain.vertices[i] = found->second;
        }
        const auto& a = output.vertices[plain.vertices[0]].position;
        const auto& b = output.vertices[plain.vertices[1]].position;
        const auto& c = output.vertices[plain.vertices[2]].position;
        const auto& d = output.vertices[plain.vertices[3]].position;
        if (predicates::orient3d(a, b, c, d) !=
            predicates::Sign::Positive)
            return TetFacetRecoveryReason::NonPositiveTet;
        output.cells.push_back(plain);
    }
    if (!haveRegion) return TetFacetRecoveryReason::UnsupportedCell;

    for (const auto edge : workspace.edges().activeBySequence()) {
        const auto data = workspace.edges().get(edge);
        const auto ia = nodeIndex.find(data.a);
        const auto ib = nodeIndex.find(data.b);
        if (ia == nodeIndex.end() || ib == nodeIndex.end())
            return TetFacetRecoveryReason::InvalidIdentity;
        if (data.requiredConstraint) {
            if (data.cls.dim != 1 || data.cls.tag < 0)
                return TetFacetRecoveryReason::InvalidEdgeClassification;
            output.protectedEdges.push_back(
                canonicalTetFacetEdge(ia->second, ib->second));
        }
    }
    for (std::size_t i = 0; i < 3; ++i) {
        const auto a = request.vertices[i];
        const auto b = request.vertices[(i + 1) % 3];
        const auto edge = workspace.edges().find(workspace.nodes(), a, b);
        if (!edge || !workspace.edges().get(*edge).requiredConstraint)
            return TetFacetRecoveryReason::MissingProtectedTargetEdge;
        const auto cls = workspace.edges().get(*edge).cls;
        if (cls.dim != 1 || cls.tag < 0)
            return TetFacetRecoveryReason::InvalidEdgeClassification;
    }

    for (const auto faceHandle : workspace.faces().activeBySequence()) {
        const auto face = workspace.faces().get(faceHandle);
        if (!defaultClass(face.cls) &&
            (face.cls.dim != 2 || face.cls.tag < 0))
            return TetFacetRecoveryReason::InvalidFaceClassification;
        if (face.cls.dim != 2) continue;
        if (face.cornerCount != 3)
            return TetFacetRecoveryReason::InvalidFaceClassification;
        std::array<TetFacetVertexId, 3> ids{};
        for (std::size_t i = 0; i < 3; ++i) {
            const auto found = nodeIndex.find(face.corners[i]);
            if (found == nodeIndex.end())
                return TetFacetRecoveryReason::InvalidIdentity;
            ids[i] = found->second;
        }
        output.protectedFaces.push_back(
            {canonicalTetFacetFace(ids), workspace.faces().seqOf(faceHandle),
             face.cls});
    }
    return TetFacetRecoveryReason::None;
}

inline TetFacetRecoveryReason targetShapePreflight(
    const topo::TopologyWorkspace& workspace,
    const TetFacetRequest& request) {
    if (request.classification.dim != 2 || request.classification.tag < 0)
        return TetFacetRecoveryReason::InvalidClassification;
    for (const auto vertex : request.vertices) {
        if (!workspace.nodes().isCurrent(vertex))
            return TetFacetRecoveryReason::InvalidIdentity;
        if (!finite(workspace.nodes().get(vertex).coord))
            return TetFacetRecoveryReason::NonfiniteCoordinate;
    }
    if (request.vertices[0] == request.vertices[1] ||
        request.vertices[1] == request.vertices[2] ||
        request.vertices[0] == request.vertices[2])
        return TetFacetRecoveryReason::InvalidIdentity;
    const auto a = workspace.nodes().get(request.vertices[0]).coord;
    const auto b = workspace.nodes().get(request.vertices[1]).coord;
    const auto c = workspace.nodes().get(request.vertices[2]).coord;
    if (predicates::orient2d({a.x, a.y}, {b.x, b.y}, {c.x, c.y}) ==
            predicates::Sign::Zero &&
        predicates::orient2d({a.y, a.z}, {b.y, b.z}, {c.y, c.z}) ==
            predicates::Sign::Zero &&
        predicates::orient2d({a.z, a.x}, {b.z, b.x}, {c.z, c.x}) ==
            predicates::Sign::Zero)
        return TetFacetRecoveryReason::DegenerateTarget;
    return TetFacetRecoveryReason::None;
}

inline void publishBindingScratch(
    TetFacetConstraintState& scratch,
    const topo::WorkspaceStamp& stamp, std::uint64_t entityId,
    const std::array<topo::NodeHandle, 3>& key,
    const mesh::Classification& classification, topo::FaceHandle face) {
    scratch.stamp = stamp;
    const PersistentTetFacet value{entityId, key, classification, face};
    const auto existing = scratch.facets.find(entityId);
    if (existing == scratch.facets.end())
        scratch.facets.emplace(entityId, value);
    else
        existing->second = value;
}

}  // namespace tet_facet_recovery_detail

inline TetFacetRecoveryOutcome recoverTetFacet(
    topo::TopologyWorkspace& workspace, TetFacetConstraintState& state,
    const TetFacetRequest& request) noexcept {
    using namespace tet_facet_recovery_detail;
    TetFacetRecoveryPhase phase = TetFacetRecoveryPhase::Preflight;
    const auto setPhase = [&](TetFacetRecoveryPhase next) {
        phase = next;
    };
    bool preparationStarted = false;
    try {
        const auto preflightStamp = currentStamp(workspace);
        if (!(request.stamp == preflightStamp) ||
            !(state.stamp == preflightStamp))
            return rejected(TetFacetRecoveryReason::StaleWorkspaceStamp);
        const auto shape = targetShapePreflight(workspace, request);
        if (shape != TetFacetRecoveryReason::None) return rejected(shape);
        const auto requestKey = canonicalTriple(workspace.nodes(),
                                                request.vertices);
        const auto stateAudit =
            auditPersistentState(workspace, state, request, requestKey);
        if (stateAudit != TetFacetRecoveryReason::None)
            return rejected(stateAudit);

        ComplexStorage storage;
        const auto complexAudit = buildComplex(workspace, request, storage);
        if (complexAudit != TetFacetRecoveryReason::None)
            return rejected(complexAudit);

        setPhase(TetFacetRecoveryPhase::Routing);
        const auto existingFace =
            workspace.faces().find(workspace.nodes(), request.vertices);
        if (existingFace) {
            const auto face = workspace.faces().get(*existingFace);
            if (face.incident.size() == 1)
                return refused(
                    TetFacetRecoveryReason::OneSidedTargetUnsupported, phase);
            if (face.incident.size() != 2 ||
                face.incident[0].sense * face.incident[1].sense != -1)
                return refused(TetFacetRecoveryReason::InvalidTargetIncidence,
                               phase);
            for (const auto& incidence : face.incident) {
                const topo::CellHandle owner{incidence.cellSlot,
                                             incidence.cellGen,
                                             workspace.cells().poolId()};
                if (!workspace.cells().isCurrent(owner) ||
                    !sameClass(workspace.cells().get(owner).cls,
                               storage.region))
                    return refused(
                        TetFacetRecoveryReason::InvalidTargetIncidence, phase);
            }

            const auto current = currentStamp(workspace);
            if (sameClass(face.cls, request.classification)) {
                const auto exact = state.facets.find(request.entityId);
                if (exact == state.facets.end()) {
                    setPhase(TetFacetRecoveryPhase::StateCopy);
                    TetFacetConstraintState scratch = state;
                    publishBindingScratch(scratch, current, request.entityId,
                                          requestKey, request.classification,
                                          *existingFace);
                    state.swap(scratch);
                } else if (exact->second.recoveredFace != *existingFace) {
                    return rejected(TetFacetRecoveryReason::FaceAlias);
                }
                TetFacetRecoveryOutcome result;
                result.status = TetFacetRecoveryStatus::AlreadyPresent;
                result.face = *existingFace;
                return result;
            }
            if (!defaultClass(face.cls))
                return refused(TetFacetRecoveryReason::ProtectedFacetConflict,
                               phase);

            auto transaction = workspace.beginTransaction();
            transaction.setFaceClassification(*existingFace,
                                              request.classification);
            setPhase(TetFacetRecoveryPhase::Preparation);
            preparationStarted = true;
            transaction.prepare();
            setPhase(TetFacetRecoveryPhase::StateCopy);
            TetFacetConstraintState scratch = state;
            publishBindingScratch(scratch, nextStamp(current),
                                  request.entityId, requestKey,
                                  request.classification, *existingFace);
            const auto committed = transaction.commit();
            (void)committed;
            state.swap(scratch);
            TetFacetRecoveryOutcome result;
            result.status = TetFacetRecoveryStatus::Marked;
            result.face = *existingFace;
            return result;
        }

        setPhase(TetFacetRecoveryPhase::Discovery);
        const auto plan = planTetFacetCavity(storage.view());
        if (!plan.planned()) {
            if (plan.reason == TetFacetGeometryReason::AllocationFailed)
                throw std::bad_alloc{};
            return refused(geometryRefusalReason(plan.reason), phase,
                           plan.reason, plan.plan.counters);
        }
        setPhase(TetFacetRecoveryPhase::Candidate);
        const auto candidate =
            reconstructTetFacetCandidate(storage.view(), plan.plan);
        if (!candidate.built()) {
            if (candidate.reason == TetFacetGeometryReason::AllocationFailed)
                throw std::bad_alloc{};
            TetFacetRecoveryOutcome result;
            result.status = candidate.status == TetFacetPlanStatus::Failed
                                ? TetFacetRecoveryStatus::Failed
                                : TetFacetRecoveryStatus::Refused;
            result.reason = candidate.status == TetFacetPlanStatus::Failed
                                ? TetFacetRecoveryReason::InvalidCandidate
                                : geometryRefusalReason(candidate.reason);
            result.phase = phase;
            result.geometryReason = candidate.reason;
            result.counters = candidate.candidate.counters;
            return result;
        }

        topo::TetCavitySpec cavity;
        cavity.stamp = request.stamp;
        cavity.target = request.vertices;
        cavity.targetClassification = request.classification;
        cavity.region = storage.region;
        for (const auto index : plan.plan.selectedCells)
            cavity.removed.push_back(storage.cellHandles[index]);
        for (const auto& plain : candidate.candidate.cells) {
            topo::TetCavityCell cell;
            for (std::size_t i = 0; i < 4; ++i)
                cell.nodes[i] = storage.nodeHandles[plain.vertices[i]];
            cell.cls = plain.classification;
            cavity.added.push_back(cell);
        }

        auto transaction = workspace.beginTetCavityTransaction();
        try {
            (void)transaction.stage(cavity);
        } catch (const std::bad_alloc&) {
            throw;
        } catch (...) {
            TetFacetRecoveryOutcome result;
            result.status = TetFacetRecoveryStatus::Failed;
            result.reason = TetFacetRecoveryReason::PreparationConflict;
            result.phase = TetFacetRecoveryPhase::Preparation;
            result.counters = candidate.candidate.counters;
            return result;
        }
        setPhase(TetFacetRecoveryPhase::Preparation);
        preparationStarted = true;
        transaction.prepare();
        const auto predictedFace = transaction.predictedTargetFace();
        setPhase(TetFacetRecoveryPhase::StateCopy);
        TetFacetConstraintState scratch = state;
        publishBindingScratch(scratch, nextStamp(request.stamp),
                              request.entityId, requestKey,
                              request.classification, predictedFace);
        const auto commit = transaction.commit();
        state.swap(scratch);
        TetFacetRecoveryOutcome result;
        result.status = TetFacetRecoveryStatus::Recovered;
        result.face = commit.targetFace;
        result.counters = candidate.candidate.counters;
        return result;
    } catch (const std::bad_alloc&) {
        TetFacetRecoveryOutcome result;
        result.status = TetFacetRecoveryStatus::Failed;
        result.reason = TetFacetRecoveryReason::AllocationFailed;
        result.phase = phase;
        result.capacityResiduePermitted = preparationStarted;
        return result;
    } catch (...) {
        if (phase == TetFacetRecoveryPhase::Preflight)
            return rejected(TetFacetRecoveryReason::InvalidIdentity);
        TetFacetRecoveryOutcome result;
        result.status = TetFacetRecoveryStatus::Failed;
        result.reason = TetFacetRecoveryReason::PreparationConflict;
        result.phase = phase == TetFacetRecoveryPhase::StateCopy
                           ? TetFacetRecoveryPhase::Preparation
                           : phase;
        result.capacityResiduePermitted = preparationStarted;
        return result;
    }
}

}  // namespace paacs::bnd
