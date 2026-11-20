// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "test_harness.hpp"
#include "sealed_cavity_fixture.hpp"
#include "workspace_digest.hpp"

#include "paacs/algo/TetRefinementMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "paacs/topo/TetRefinement.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string_view>
#include <vector>
#if defined(_WIN32)
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

namespace {
using namespace paacs;

class ReversedFaceModel final : public geo::Model {
public:
    explicit ReversedFaceModel(const geo::Model& source) : source_(source) {}
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        return source_.edges();
    }
    std::vector<const geo::GFace*> faces() const override {
        auto result=source_.faces();
        std::reverse(result.begin(),result.end());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override { return source_.edge(tag); }
    const geo::GFace* face(int tag) const override { return source_.face(tag); }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }
private:
    const geo::Model& source_;
};

void checkedPowerPredicates() {
    TEST("T3-A exact checked power predicates");
    using predicates::PowerPredicateStatus;
    using predicates::Sign;
    const Vec3 a{0,0,0}, b{2,0,0};
    CHECK(predicates::edgeDiametralPower(a,b,{1,0,0}).sign == Sign::Negative);
    CHECK(predicates::edgeDiametralPower(a,b,{1,1,0}).sign == Sign::Zero);
    CHECK(predicates::edgeDiametralPower(a,b,{3,0,0}).sign == Sign::Positive);
    CHECK(predicates::edgeDiametralPower(b,a,{1,0,0}).sign == Sign::Negative);
    CHECK(predicates::edgeDiametralPower(a,a,{1,0,0}).status ==
          PowerPredicateStatus::DegenerateEdge);
    CHECK(predicates::edgeDiametralPower(
              a,b,{(std::numeric_limits<double>::infinity)(),0,0}).status ==
          PowerPredicateStatus::NonFiniteInput);

    const Vec3 c{0,2,0};
    CHECK(predicates::triangleCircumballPower(a,b,c,{1,1,0}).sign ==
          Sign::Negative);
    CHECK(predicates::triangleCircumballPower(a,b,c,{2,2,0}).sign ==
          Sign::Zero);
    CHECK(predicates::triangleCircumballPower(a,c,b,{1,1,0}).sign ==
          Sign::Negative);
    const std::array<Vec3,3> triangle{a,b,c};
    const std::array<std::array<std::size_t,3>,6> permutations{{
        {{0,1,2}},{{0,2,1}},{{1,0,2}},{{1,2,0}},{{2,0,1}},{{2,1,0}}}};
    for (const auto& permutation:permutations)
        CHECK(predicates::triangleCircumballPower(
            triangle[permutation[0]],triangle[permutation[1]],
            triangle[permutation[2]],{1,1,0}).sign==Sign::Negative);
    CHECK(predicates::triangleCircumballPower(a,b,{4,0,0},{1,0,0}).status ==
          PowerPredicateStatus::DegenerateTriangle);

    const double tiny = (std::numeric_limits<double>::denorm_min)();
    const auto extreme = predicates::edgeDiametralPower(
        {-(std::numeric_limits<double>::max)(),tiny,0},
        { (std::numeric_limits<double>::max)(),tiny,0},{0,tiny,0});
    CHECK(extreme.status == PowerPredicateStatus::Valid);
    CHECK(extreme.sign == Sign::Negative);
    const double translated=0x1p50;
    const auto translatedResult=predicates::triangleCircumballPower(
        {translated,translated,translated},
        {translated+2.0,translated,translated},
        {translated,translated+2.0,translated},
        {translated+1.0,translated+1.0,translated});
    CHECK(translatedResult.status==PowerPredicateStatus::Valid);
    CHECK(translatedResult.sign==Sign::Negative);
    CHECK(predicates::triangleCircumballPower(
        a,b,c,{0,0,(std::numeric_limits<double>::quiet_NaN)()}).status==
        PowerPredicateStatus::NonFiniteInput);
}

void certifiedProtectedBounds() {
    TEST("T3-A certified protected AABB contains witness missed by one ULP");
    const std::array<Vec3,2> edge{{
        {-0x1.e57f18445ad90p-47,-0x1.8ad78b8b9d320p-31,
          0x1.a4501ba705cf4p-56},
        { 0x1.bddb1b106edcep-41,-0x1.baf75d8cfe134p-47,
         -0x1.a7b30c43bcda4p-62}}};
    const Vec3 centre=(edge[0]+edge[1])*0.5;
    const double radius=(edge[1]-edge[0]).norm()*0.5;
    const double naive=std::nextafter(
        centre.y+radius,(std::numeric_limits<double>::infinity)());
    const Vec3 witness{
        centre.x,std::nextafter(
            naive,(std::numeric_limits<double>::infinity)()),centre.z};
    CHECK(witness.y>naive);
    const auto exact=predicates::edgeDiametralPower(
        edge[0],edge[1],witness);
    CHECK(exact.status==predicates::PowerPredicateStatus::Valid);
    CHECK(exact.sign==predicates::Sign::Negative);
    Vec3 lo{},hi{};
    CHECK(algo::t3a_detail::certifiedProtectedEdgeBoundsForTesting(
        edge,lo,hi));
    CHECK(witness.x>=lo.x&&witness.x<=hi.x);
    CHECK(witness.y>=lo.y&&witness.y<=hi.y);
    CHECK(witness.z>=lo.z&&witness.z<=hi.z);

    const double maximum=(std::numeric_limits<double>::max)();
    const std::array<Vec3,2> unrepresentable{{
        {maximum,0,0},{std::nextafter(maximum,0.0),0,0}}};
    CHECK(!algo::t3a_detail::certifiedProtectedEdgeBoundsForTesting(
        unrepresentable,lo,hi));

    const std::array<Vec3,3> facet{{
        {-0x1.1661b1faf92d9p-207,0x1.e2146a6980a80p-208,0},
        { 0x1.55f509f581541p-207,0x1.e2146a6980a80p-208,0},
        { 0x1.fc9abfd441340p-211,0x1.139ac9967eca6p-206,0}}};
    const Vec3 u=facet[1]-facet[0],v=facet[2]-facet[0];
    const Vec3 n=u.cross(v);
    const Vec3 relative=(u.norm2()*v.cross(n)+v.norm2()*n.cross(u)) /
                        (2.0*n.norm2());
    const Vec3 computedCentre=facet[0]+relative;
    const double naiveFacetHi=std::nextafter(
        computedCentre.y+relative.norm(),
        (std::numeric_limits<double>::infinity)());
    CHECK(facet[2].y>naiveFacetHi);
    const auto onCircumball=predicates::triangleCircumballPower(
        facet[0],facet[1],facet[2],facet[2]);
    CHECK(onCircumball.status==predicates::PowerPredicateStatus::Valid);
    CHECK(onCircumball.sign==predicates::Sign::Zero);
    CHECK(algo::t3a_detail::certifiedProtectedFacetBoundsForTesting(
        facet,lo,hi));
    CHECK(facet[2].x>=lo.x&&facet[2].x<=hi.x);
    CHECK(facet[2].y>=lo.y&&facet[2].y<=hi.y);
    CHECK(facet[2].z>=lo.z&&facet[2].z<=hi.z);
}

void violationIdentityAndVolumeFloor() {
    TEST("T3-A exact size witness and inherited created-volume floor");
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(1.0);
    algo::TetRefinementPolicy policy;
    policy.maxEdgeToTarget=0.5;
    policy.maxRadiusEdge=100.0;
    const std::array<Vec3,4> points{{
        {0,0,0},{4,0,0},{0,3,0},{0,0,2}}};
    const std::array<Index,4> global{{10,11,12,13}};
    algo::TetViolationRow row;
    CHECK(algo::t3a_detail::evaluateTetViolationForTesting(
        points,global,9,spec,policy,row));
    CHECK(row.key.kind==algo::TetViolationKind::Size);
    CHECK(row.key.hasSizeEdge);
    CHECK(row.key.sizeEdge.vertices[0].value==11);
    CHECK(row.key.sizeEdge.vertices[1].value==12);
    CHECK(row.normalizedSeverity==100.0);
    CHECK(row.hasSizeSamples);
    CHECK(row.sizeSamples[0].value==1.0&&row.sizeSamples[1].value==1.0&&
          row.sizeSamples[2].value==1.0);

    sizing::MeshSpec relaxed;
    relaxed.field=std::make_shared<sizing::ConstantField>(1000.0);
    algo::TetRefinementPolicy radiusPolicy;
    radiusPolicy.maxEdgeToTarget=100.0;
    radiusPolicy.maxRadiusEdge=1.0;
    const std::array<Vec3,4> skinny{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,0.01}}};
    CHECK(algo::t3a_detail::evaluateTetViolationForTesting(
        skinny,global,9,relaxed,radiusPolicy,row));
    CHECK(row.key.kind==algo::TetViolationKind::RadiusEdge);
    CHECK(!row.key.hasSizeEdge);

    auto sliverPolicy=radiusPolicy;
    sliverPolicy.maxRadiusEdge=1.0e9;
    sliverPolicy.sliverPolicy=algo::TetSliverPolicy::BoundedInsertion;
    sliverPolicy.sliverRadiusRatio=0.5;
    sliverPolicy.minimumSliverImprovement=0.01;
    sliverPolicy.maxSliverCandidateAttempts=5;
    CHECK(algo::t3a_detail::evaluateTetViolationForTesting(
        skinny,global,9,relaxed,sliverPolicy,row));
    CHECK(row.key.kind==algo::TetViolationKind::Sliver);

    sliverPolicy.sliverPolicy=algo::TetSliverPolicy::ReportOnly;
    CHECK(!algo::t3a_detail::evaluateTetViolationForTesting(
        skinny,global,9,relaxed,sliverPolicy,row));
}

struct OneTetFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::array<topo::NodeHandle,4> nodes{};
    topo::CellHandle cell{};
    topo::EdgeHandle protectedEdge{};
    topo::FaceHandle protectedFace{};

    OneTetFixture() {
        topo::NodePool np;
        std::array<topo::NodeHandle,4> churnNodes{
            np.add({10,0,0},{3,8}),np.add({11,0,0},{3,8}),
            np.add({10,1,0},{3,8}),np.add({10,0,1},{3,8})};
        topo::CellPool cp;
        const auto churnCell=cp.add(
            topo::ElementType::Tet4,churnNodes,{3,8},np);
        topo::EdgePool ep;
        (void)ep.insertCellEdges(np,cp,churnCell);
        topo::FacePool fp;
        (void)fp.insertCellFaces(np,cp,ep,churnCell);
        fp.removeCellFaces(cp,churnCell);
        ep.removeCellEdges(cp,churnCell);
        cp.remove(churnCell);
        for (const auto node:churnNodes) np.remove(node);

        nodes[0] = np.add({0,0,0},{3,9});
        nodes[1] = np.add({1,0,0},{3,9});
        nodes[2] = np.add({0,1,0},{3,9});
        nodes[3] = np.add({0,0,1},{3,9});
        cell = cp.add(topo::ElementType::Tet4,nodes,{3,9},np);
        (void)ep.insertCellEdges(np,cp,cell);
        const auto edge=ep.find(np,nodes[0],nodes[1]);
        if (!edge) throw std::logic_error("fixture edge missing");
        protectedEdge=*edge;
        ep.setRequiredConstraint(protectedEdge,true);
        ep.setClassification(protectedEdge,{1,66});
        (void)fp.insertCellFaces(np,cp,ep,cell);
        const std::array<topo::NodeHandle,3> boundary{
            nodes[0],nodes[2],nodes[1]};
        const auto found = fp.find(np,boundary);
        if (!found) throw std::logic_error("fixture face missing");
        protectedFace = *found;
        fp.setClassification(protectedFace,{2,77});
        workspace = std::make_unique<topo::TopologyWorkspace>(
            std::move(np),std::move(cp),std::move(ep),std::move(fp));
    }
};

topo::TetRefinementCavitySpec splitSpec(const OneTetFixture& f) {
    topo::TetRefinementCavitySpec spec;
    spec.stamp = topo::WorkspaceStamp::of(*f.workspace);
    spec.stagedNode = {0,{0.2,0.2,0.2},{3,9}};
    spec.removed.push_back(f.cell);
    static constexpr std::array<std::array<int,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    for (const auto& face : faces) {
        topo::TetRefinementCellSpec cell;
        cell.vertices[0] = topo::TetRefinementNodeRef::existingNode(
            f.nodes[static_cast<std::size_t>(face[0])]);
        cell.vertices[1] = topo::TetRefinementNodeRef::existingNode(
            f.nodes[static_cast<std::size_t>(face[1])]);
        cell.vertices[2] = topo::TetRefinementNodeRef::existingNode(
            f.nodes[static_cast<std::size_t>(face[2])]);
        cell.vertices[3] = topo::TetRefinementNodeRef::stagedNode();
        cell.classification = {3,9};
        spec.added.push_back(cell);
    }
    return spec;
}

void stagedTetRefinementTransaction() {
    TEST("T3-A staged Tet refinement transaction and boundary owner change");
    OneTetFixture f;
    const auto faceBefore = f.workspace->faces().get(f.protectedFace);
    const auto edgeBefore = f.workspace->edges().get(f.protectedEdge);
    CHECK(edgeBefore.requiredConstraint);
    CHECK(edgeBefore.useCount==1);
    CHECK(faceBefore.incident.size() == 1);
    const topo::CellHandle ownerBefore{
        faceBefore.incident[0].cellSlot,faceBefore.incident[0].cellGen,
        f.workspace->cells().poolId()};
    const auto revision = f.workspace->topologyRevision();
    const auto spec = splitSpec(f);
    auto tx = f.workspace->beginTetRefinementTransaction();
    const auto& preview = tx.stage(spec);
    CHECK(preview.valid);
    CHECK(preview.removedCellCount == 1);
    CHECK(preview.addedCellCount == 4);
    const auto& prediction = tx.prepare();
    CHECK(prediction.postCommitStamp.revision == revision + 1);
    CHECK(prediction.stagedNode.isSpecified());
    CHECK(prediction.cells.size() == 5);
    CHECK(!prediction.edges.empty());
    CHECK(!prediction.faces.empty());
    const auto predictedEdge=std::find_if(
        prediction.edges.begin(),prediction.edges.end(),[&](const auto& row) {
            return row.predicted==f.protectedEdge;
        });
    CHECK(predictedEdge!=prediction.edges.end());
    if (predictedEdge!=prediction.edges.end()) {
        CHECK(predictedEdge->transition==topo::TetRefinementTransition::Retained);
        CHECK(predictedEdge->initialUseCount==1);
        CHECK(predictedEdge->finalUseCount==2);
        CHECK(predictedEdge->required);
        CHECK(predictedEdge->classification.dim==1&&
              predictedEdge->classification.tag==66);
    }
    const auto commit = tx.commit();
    CHECK(commit.actualNode == commit.predictedNode);
    CHECK(commit.stamp.revision == revision + 1);
    CHECK(f.workspace->nodes().isCurrent(commit.actualNode));
    CHECK(!f.workspace->cells().isCurrent(f.cell));
    CHECK(f.workspace->faces().isCurrent(f.protectedFace));
    CHECK(f.workspace->edges().isCurrent(f.protectedEdge));
    const auto edgeAfter=f.workspace->edges().get(f.protectedEdge);
    CHECK(edgeAfter.requiredConstraint);
    CHECK(edgeAfter.cls.dim==1&&edgeAfter.cls.tag==66);
    CHECK(edgeAfter.useCount==2);
    const auto faceAfter = f.workspace->faces().get(f.protectedFace);
    CHECK(faceAfter.cls.dim == 2 && faceAfter.cls.tag == 77);
    CHECK(faceAfter.incident.size() == 1);
    const topo::CellHandle ownerAfter{
        faceAfter.incident[0].cellSlot,faceAfter.incident[0].cellGen,
        f.workspace->cells().poolId()};
    CHECK(ownerAfter != ownerBefore);
    CHECK(f.workspace->cells().isCurrent(ownerAfter));
    for (const auto& row : commit.cells) CHECK(row.actual == row.predicted);
    for (const auto& row : commit.edges) CHECK(row.actual == row.predicted);
    for (const auto& row : commit.faces) CHECK(row.actual == row.predicted);
}

struct PlainTetFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::CellHandle> cells;

    PlainTetFixture(const std::vector<Vec3>& points,
                    const std::vector<std::array<std::size_t,4>>& connectivity) {
        topo::NodePool nodes;
        std::vector<topo::NodeHandle> handles;
        handles.reserve(points.size());
        for (const auto point:points) handles.push_back(nodes.add(point,{3,17}));
        topo::CellPool cellPool;
        for (const auto& tet:connectivity) {
            std::array<topo::NodeHandle,4> corners{
                handles[tet[0]],handles[tet[1]],handles[tet[2]],handles[tet[3]]};
            cells.push_back(cellPool.add(
                topo::ElementType::Tet4,corners,{3,17},nodes));
        }
        topo::EdgePool edges;
        for (const auto cell:cells) (void)edges.insertCellEdges(nodes,cellPool,cell);
        topo::FacePool faces;
        for (const auto cell:cells)
            (void)faces.insertCellFaces(nodes,cellPool,edges,cell);
        workspace=std::make_unique<topo::TopologyWorkspace>(
            std::move(nodes),std::move(cellPool),std::move(edges),std::move(faces));
    }
};

void indexedLocationAndCandidateContracts() {
    TEST("T3-A indexed location and complete plain candidate contracts");
    PlainTetFixture separated(
        {{0,0,0},{1,0,0},{0,1,0},{0,0,1},
         {3,0,0},{4,0,0},{3,1,0},{3,0,1}},
        {{{0,1,2,3}},{{4,5,6,7}}});
    topo::CellHandle located{};
    CHECK(algo::t3a_detail::exactLocationForTesting(
        *separated.workspace,{3.1,0.1,0.1},located));
    CHECK(located==separated.cells[1]);
    CHECK(!algo::t3a_detail::validateRefinementCavityForTesting(
        *separated.workspace,separated.cells,{0.2,0.2,0.2}));

    OneTetFixture valid;
    const std::array<topo::CellHandle,1> one{{valid.cell}};
    CHECK(algo::t3a_detail::validateRefinementCavityForTesting(
        *valid.workspace,one,{0.2,0.2,0.2}));

    // The survivor is geometrically coincident but identity-distinct.  A
    // topology-only transaction could stage the split; the complete T3-A
    // validator must reject its illegal closed contact first.
    PlainTetFixture coincident(
        {{0,0,0},{1,0,0},{0,1,0},{0,0,1},
         {0,0,0},{1,0,0},{0,1,0},{0,0,1}},
        {{{0,1,2,3}},{{4,5,6,7}}});
    const std::array<topo::CellHandle,1> first{{coincident.cells[0]}};
    CHECK(!algo::t3a_detail::validateRefinementCavityForTesting(
        *coincident.workspace,first,{0.2,0.2,0.2}));
}

class Provider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{701,4};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;

    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        return current;
    }
    algo::VolumeFaceOwnership ownership(
        algo::FacePersistentId face) const override {
        algo::VolumeFaceOwnership result;
        result.providerStamp=current;
        result.facePersistentId=face;
        const auto found=rows.find(face);
        if (found!=rows.end()) result.owners=found->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
        std::span<const Vec3>) const override {
        algo::ProviderRecoveryEnvelope result;
        result.providerStamp=current;
        result.availability=
            algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

void configureProvider(const geo::Model& model,Provider& provider) {
    const auto solids=model.solids();
    CHECK(solids.size()==1);
    if (solids.size()!=1) return;
    Vec3 lo{},hi{};
    solids.front()->bbox(lo,hi);
    const Vec3 centre=(lo+hi)*0.5;
    for (const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{};
        face->paramBounds(uvLo,uvHi);
        const Vec2 uv=(uvLo+uvHi)*0.5;
        const int sense=face->normalAt(uv).dot(
            face->pointAt(uv)-centre)>0.0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,sense}};
    }
}

algo::T3ARequest requestFor(const Provider& provider) {
    algo::T3ARequest request;
    request.conformingRequest.expectedSurfaceSetStamp={9801,1};
    request.conformingRequest.expectedOwnershipProviderStamp=provider.current;
    request.conformingRequest.envelopeMode=
        algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    request.refinement.maxEdgeToTarget=2.0;
    request.refinement.maxRadiusEdge=100.0;
    request.refinement.maxInsertions=8;
    request.refinement.maxWorkUnits=10000;
    request.refinement.maxOutputNodes=1000;
    request.refinement.maxOutputTets=10000;
    request.refinement.maxReportRows=45;
    request.refinement.maxReportWitnesses=10000;
    request.refinement.maxReportTopologyHandles=10008;
    return request;
}

void requireProcessedXorUnattempted(const algo::T3AReport& report);
const algo::t3a_detail::T3ABlockerStateDigest* findBlockerDigest(
    const algo::t3a_detail::T3ABlockerDigestTrace& trace,
    algo::t3a_detail::T3ABlockerDigestCheckpoint checkpoint,
    std::size_t occurrence=0) noexcept;
void requireSameBlockerPreparationState(
    const algo::t3a_detail::T3ABlockerStateDigest& before,
    const algo::t3a_detail::T3ABlockerStateDigest& after,
    bool terminalFinalizationAllowed);

void boundedPublicRoute() {
    TEST("T3-A bounded public route preserves v3 boundary and publishes atomically");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(100.0);
    auto request=requestFor(provider);
    mesh::Mesh destination;
    destination.addNode({-8,-7,-6},{0,441});
    const auto before=mesh::meshHash(destination);
    const auto outcome=algo::meshModelRefinedTet4(
        *model,spec,request,destination,provider);
    if (!outcome.published())
        std::printf("  t3a diagnostic status=%d reason=%d phase=%d source=%d/%d\n",
            static_cast<int>(outcome.status),static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase),
            outcome.sourceReceipt?static_cast<int>(outcome.sourceReceipt->status):-1,
            outcome.sourceReceipt?static_cast<int>(outcome.sourceReceipt->failure.reason):-1);
    CHECK(outcome.published());
    CHECK(!outcome.destinationUnchanged);
    CHECK(mesh::meshHash(destination)!=before);
    CHECK(outcome.retainedResult!=nullptr);
    CHECK(outcome.report!=nullptr&&outcome.report->reportComplete);
    CHECK(outcome.sourceReceipt!=nullptr&&
          outcome.sourceReceipt->status==algo::BoundaryConformingV3Status::Published);
    if (outcome.retainedResult) {
        CHECK(std::string_view(outcome.retainedResult->routeId)==
              "interior-refined-tet4-t3a");
        CHECK(outcome.retainedResult->routeVersion==algo::kT3ARouteVersion);
        CHECK(outcome.retainedResult->publicationReady);
        CHECK(outcome.retainedResult->metrics.runComplete);
        CHECK(outcome.retainedResult->metrics.committedInsertions==0);
        CHECK(outcome.retainedResult->metrics.authoritativeAccountedBytes>0);
        CHECK(outcome.retainedResult->metrics.transientAccountedHighWaterBytes>0);
        CHECK(outcome.retainedResult->metrics.publicationAccountedBytes>0);
        CHECK(!outcome.retainedResult->immutableBoundary.empty());
        CHECK(outcome.retainedResult->boundaryAudit.orphanExterior==0);
        CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
        CHECK(outcome.retainedResult->quality.inverted==0);
    }

    auto floorRequest=request;
    floorRequest.refinement.minimumCreatedTetVolume=1000.0;
    mesh::Mesh floorDestination;
    const auto floorOutcome=algo::meshModelRefinedTet4(
        *model,spec,floorRequest,floorDestination,provider);
    CHECK(floorOutcome.published());
    CHECK(floorOutcome.retainedResult!=nullptr);
    CHECK(floorOutcome.retainedResult&&
          floorOutcome.retainedResult->quality.inheritedBelowCreatedFloor>0);
    CHECK(floorOutcome.retainedResult&&
          floorOutcome.retainedResult->quality.createdBelowFloor==0);

    auto invalid=request;
    invalid.refinement.maxEdgeToTarget=0.0;
    mesh::Mesh unchanged;
    unchanged.addNode({-8,-7,-6},{0,441});
    const auto invalidBefore=mesh::meshHash(unchanged);
    const auto rejected=algo::meshModelRefinedTet4(
        *model,spec,invalid,unchanged,provider);
    CHECK(rejected.status==algo::T3AStatus::Rejected);
    CHECK(rejected.failure.reason==algo::T3AReason::InvalidRefinementPolicy);
    CHECK(rejected.destinationUnchanged);
    CHECK(mesh::meshHash(unchanged)==invalidBefore);

    sizing::MeshSpec refiningSpec;
    refiningSpec.field=std::make_shared<sizing::ConstantField>(0.55);
    auto refiningRequest=requestFor(provider);
    refiningRequest.refinement.maxInsertions=16;
    refiningRequest.refinement.maxEdgeToTarget=0.5;
    refiningRequest.refinement.maxWorkUnits=100000;
    refiningRequest.refinement.maxReportRows=85;
    refiningRequest.refinement.maxReportWitnesses=100000;
    refiningRequest.refinement.maxReportTopologyHandles=100016;
    mesh::Mesh refinementDestination;
    refinementDestination.addNode({-4,-3,-2},{0,991});
    const auto refinementBefore=mesh::meshHash(refinementDestination);
    const auto refined=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,refinementDestination,provider);
    CHECK(refined.status==algo::T3AStatus::Incomplete);
    CHECK(refined.failure.reason==algo::T3AReason::BoundaryRefinementRequired);
    CHECK(refined.failure.phase==algo::T3APhase::Encroachment);
    CHECK(refined.destinationUnchanged);
    CHECK(mesh::meshHash(refinementDestination)==refinementBefore);
    CHECK(refined.retainedResult!=nullptr);
    CHECK(refined.report!=nullptr&&refined.report->reportComplete);
    if (refined.retainedResult) {
        CHECK(refined.retainedResult->metrics.committedInsertions==2);
        CHECK(refined.retainedResult->metrics.perInsertionGlobalTetScans==0);
        CHECK(refined.retainedResult->metrics.coordinatorRowsPrepared>0);
        CHECK(refined.retainedResult->metrics.locationIndexUpdateVisits>0);
        CHECK(refined.retainedResult->metrics.locationIndexUpdateVisits<=
              refined.retainedResult->metrics.locationIndexUpdateWorkCharged);
        CHECK(refined.retainedResult->metrics.coordinatorTailWorkCharged>0);
        CHECK(!refined.retainedResult->publicationReady);
        CHECK(refined.retainedResult->quality.inverted==0);
        CHECK(refined.retainedResult->boundaryAudit.orphanExterior==0);
        CHECK(refined.retainedResult->boundaryAudit.unexpectedExterior==0);
    }
    CHECK(refined.report&&refined.report->committedPrefixUsed==2);
    CHECK(refined.report&&refined.report->protectedWitnessesUsed>0);
    if (refined.report) {
        const auto attempts=refined.report->attempts();
        const auto encroached=std::find_if(
            attempts.begin(),attempts.end(),[](const auto& row) {
                return row.disposition==algo::T3AAttemptDisposition::Encroached;
            });
        CHECK(encroached!=attempts.end());
        if (encroached!=attempts.end()) {
            CHECK(std::any_of(encroached+1,attempts.end(),[](const auto& row) {
                return row.disposition==algo::T3AAttemptDisposition::Committed;
            }));
            CHECK(encroached->witnessCount>1);
            for (std::uint64_t i=1;i<encroached->witnessCount;++i) {
                const auto& a=refined.report->protectedWitnessStorage[
                    static_cast<std::size_t>(encroached->witnessOffset+i-1)];
                const auto& b=refined.report->protectedWitnessStorage[
                    static_cast<std::size_t>(encroached->witnessOffset+i)];
                CHECK(a.kind<b.kind||
                      (a.kind==b.kind&&
                       (a.kind==algo::T3AProtectedFeatureKind::Edge
                            ? !(b.edge<a.edge):!(b.facet<a.facet))));
            }
        }
    }
    if (refined.report)
        std::printf("  t3a prefix attempts=%llu commits=%llu witnesses=%llu work=%llu\n",
            static_cast<unsigned long long>(refined.report->attemptsUsed),
            static_cast<unsigned long long>(refined.report->committedPrefixUsed),
            static_cast<unsigned long long>(refined.report->protectedWitnessesUsed),
            static_cast<unsigned long long>(
                refined.retainedResult?refined.retainedResult->metrics.workUnitsConsumed:0));

    mesh::Mesh witnessFaultDestination;
    witnessFaultDestination.addNode({-4,-3,-2},{0,993});
    const auto witnessFaultBefore=mesh::meshHash(witnessFaultDestination);
    auto& blockerTrace=algo::t3a_detail::blockerDigestTraceForTesting();
    const auto blockerTransitionMask=
        (std::uint32_t{1}<<static_cast<std::uint32_t>(
            algo::t3a_detail::T3ABlockerDigestCheckpoint::
                BlockTransitionBefore))|
        (std::uint32_t{1}<<static_cast<std::uint32_t>(
            algo::t3a_detail::T3ABlockerDigestCheckpoint::
                BlockTransitionTerminalAfter));
    blockerTrace.reset(true,blockerTransitionMask);
    algo::t3a_detail::reportFaultCountdownForTesting()=1;
    const auto witnessFault=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,witnessFaultDestination,provider);
    algo::t3a_detail::reportFaultCountdownForTesting()=0;
    CHECK(witnessFault.status==algo::T3AStatus::Failed);
    CHECK(witnessFault.failure.reason==algo::T3AReason::AllocationFailed);
    CHECK(witnessFault.failure.phase==algo::T3APhase::Encroachment);
    CHECK(witnessFault.destinationUnchanged);
    CHECK(mesh::meshHash(witnessFaultDestination)==witnessFaultBefore);
    CHECK(witnessFault.report&&witnessFault.report->reportComplete);
    CHECK(witnessFault.report&&witnessFault.report->attemptsUsed>0);
    if (witnessFault.report&&witnessFault.report->attemptsUsed>0)
        CHECK(witnessFault.report->attemptStorage[
            static_cast<std::size_t>(witnessFault.report->attemptsUsed-1)].
                disposition==algo::T3AAttemptDisposition::TerminalFailure);
    const auto* blockerBefore=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionBefore);
    const auto* blockerAfter=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionTerminalAfter);
    CHECK(blockerBefore!=nullptr);
    CHECK(blockerAfter!=nullptr);
    if (blockerBefore&&blockerAfter)
        requireSameBlockerPreparationState(
            *blockerBefore,*blockerAfter,true);
    blockerTrace.reset(false);
    mesh::Mesh witnessRetryDestination;
    const auto witnessRetry=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,witnessRetryDestination,provider);
    CHECK(witnessRetry.status==algo::T3AStatus::Incomplete);
    CHECK(witnessRetry.retainedResult&&
          witnessRetry.retainedResult->metrics.committedInsertions==2);

    const auto checkBlockedBudgetStop=[&](const algo::T3AOutcome& stopped,
                                          algo::T3AReason reason,
                                          std::uint64_t destinationDigest) {
        CHECK(stopped.status==algo::T3AStatus::Incomplete);
        CHECK(stopped.failure.reason==reason);
        CHECK(stopped.failure.phase==algo::T3APhase::Encroachment);
        CHECK(stopped.destinationUnchanged);
        CHECK(stopped.retainedResult!=nullptr);
        CHECK(stopped.report&&stopped.report->reportComplete);
        CHECK(stopped.report&&stopped.report->blockedUsed==0);
        CHECK(stopped.report&&stopped.report->attemptsUsed>0);
        CHECK(stopped.retainedResult&&
              stopped.retainedResult->metrics.phase3Allocations==0);
        CHECK(stopped.retainedResult&&
              stopped.retainedResult->metrics.publicationTailAllocations==0);
        if (!stopped.report||stopped.report->attemptsUsed==0) return;
        const auto& terminal=stopped.report->attemptStorage[
            static_cast<std::size_t>(stopped.report->attemptsUsed-1)];
        CHECK(terminal.disposition==
              algo::T3AAttemptDisposition::TerminalFailure);
        CHECK(terminal.reason==reason);
        CHECK(stopped.report->hasTerminalViolation);
        CHECK(terminal.violation==stopped.report->terminalViolation);
        for (std::uint64_t i=0;i<stopped.report->unattemptedUsed;++i)
            CHECK(!(stopped.report->unattemptedStorage[
                        static_cast<std::size_t>(i)]==terminal.violation));
        CHECK(destinationDigest!=0);
        if (stopped.report) requireProcessedXorUnattempted(*stopped.report);
    };

    mesh::Mesh reportBudgetDestination;
    reportBudgetDestination.addNode({-4,-3,-2},{0,994});
    const auto reportBudgetBefore=mesh::meshHash(reportBudgetDestination);
    algo::t3a_detail::blockerReportBudgetCountdownForTesting()=1;
    blockerTrace.reset(true,blockerTransitionMask);
    const auto reportBudget=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,reportBudgetDestination,provider);
    algo::t3a_detail::blockerReportBudgetCountdownForTesting()=0;
    CHECK(mesh::meshHash(reportBudgetDestination)==reportBudgetBefore);
    checkBlockedBudgetStop(
        reportBudget,algo::T3AReason::ReportBudgetReached,reportBudgetBefore);
    blockerBefore=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionBefore);
    blockerAfter=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionTerminalAfter);
    CHECK(blockerBefore!=nullptr&&blockerAfter!=nullptr);
    if (blockerBefore&&blockerAfter)
        requireSameBlockerPreparationState(
            *blockerBefore,*blockerAfter,true);
    blockerTrace.reset(false);

    mesh::Mesh workBudgetDestination;
    workBudgetDestination.addNode({-4,-3,-2},{0,995});
    const auto workBudgetBefore=mesh::meshHash(workBudgetDestination);
    algo::t3a_detail::blockerWorkBudgetCountdownForTesting()=1;
    blockerTrace.reset(true,blockerTransitionMask);
    const auto workBudget=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,workBudgetDestination,provider);
    algo::t3a_detail::blockerWorkBudgetCountdownForTesting()=0;
    CHECK(mesh::meshHash(workBudgetDestination)==workBudgetBefore);
    checkBlockedBudgetStop(
        workBudget,algo::T3AReason::WorkBudgetReached,workBudgetBefore);
    blockerBefore=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionBefore);
    blockerAfter=findBlockerDigest(
        blockerTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            BlockTransitionTerminalAfter);
    CHECK(blockerBefore!=nullptr&&blockerAfter!=nullptr);
    if (blockerBefore&&blockerAfter)
        requireSameBlockerPreparationState(
            *blockerBefore,*blockerAfter,true);
    blockerTrace.reset(false);
    if (reportBudget.retainedResult&&workBudget.retainedResult&&
        reportBudget.retainedResult->authoritative&&
        workBudget.retainedResult->authoritative&&
        reportBudget.retainedResult->authoritative->workspace&&
        workBudget.retainedResult->authoritative->workspace) {
        CHECK(paacs_tests::observableDigest(
                  *reportBudget.retainedResult->authoritative->workspace,
                  false,false)==
              paacs_tests::observableDigest(
                  *workBudget.retainedResult->authoritative->workspace,
                  false,false));
    }

    auto volumeBlockedRequest=refiningRequest;
    volumeBlockedRequest.refinement.minimumCreatedTetVolume=1.0;
    volumeBlockedRequest.refinement.maxWorkUnits=100000;
    volumeBlockedRequest.refinement.maxReportTopologyHandles=100016;
    mesh::Mesh volumeBlockedDestination;
    volumeBlockedDestination.addNode({-4,-3,-2},{0,992});
    const auto volumeBlockedBefore=mesh::meshHash(volumeBlockedDestination);
    const auto volumeBlocked=algo::meshModelRefinedTet4(
        *model,refiningSpec,volumeBlockedRequest,volumeBlockedDestination,
        provider);
    std::printf("  t3a volume block status=%d reason=%d phase=%d commits=%llu "
                "attempts=%llu work=%llu\n",
        static_cast<int>(volumeBlocked.status),
        static_cast<int>(volumeBlocked.failure.reason),
        static_cast<int>(volumeBlocked.failure.phase),
        static_cast<unsigned long long>(volumeBlocked.retainedResult?
            volumeBlocked.retainedResult->metrics.committedInsertions:0),
        static_cast<unsigned long long>(volumeBlocked.retainedResult?
            volumeBlocked.retainedResult->metrics.candidateAttempts:0),
        static_cast<unsigned long long>(volumeBlocked.retainedResult?
            volumeBlocked.retainedResult->metrics.workUnitsConsumed:0));
    CHECK(volumeBlocked.status==algo::T3AStatus::Incomplete);
    CHECK(volumeBlocked.failure.reason==algo::T3AReason::CreatedTetVolumeBlocked);
    CHECK(volumeBlocked.destinationUnchanged);
    CHECK(mesh::meshHash(volumeBlockedDestination)==volumeBlockedBefore);
    CHECK(volumeBlocked.retainedResult&&
          volumeBlocked.retainedResult->metrics.committedInsertions==0);

    mesh::Mesh quarantinedDestination;
    quarantinedDestination.addNode({-4,-3,-2},{0,991});
    const auto quarantineBefore=mesh::meshHash(quarantinedDestination);
    algo::t3a_detail::injectPostCommitAuditFailureForTesting()=true;
    const auto quarantined=algo::meshModelRefinedTet4(
        *model,refiningSpec,refiningRequest,quarantinedDestination,provider);
    algo::t3a_detail::injectPostCommitAuditFailureForTesting()=false;
    CHECK(quarantined.status==algo::T3AStatus::Failed);
    CHECK(quarantined.failure.reason==algo::T3AReason::InternalError);
    CHECK(quarantined.failure.phase==algo::T3APhase::LocalAudit);
    CHECK(quarantined.destinationUnchanged);
    CHECK(mesh::meshHash(quarantinedDestination)==quarantineBefore);
    CHECK(quarantined.retainedResult==nullptr);
}

class BallField final : public sizing::SizeField {
public:
    BallField(Vec3 centre,double radius,double inside,double outside)
        : centre_(centre),radius2_(radius*radius),inside_(inside),
          outside_(outside) {}
    double eval(const Vec3& point) const override {
        return (point-centre_).norm2()<=radius2_?inside_:outside_;
    }
private:
    Vec3 centre_{};
    double radius2_=0.0,inside_=1.0,outside_=1.0;
};

class RefinementPhaseBallField final : public sizing::SizeField {
public:
    RefinementPhaseBallField(Vec3 centre,double radius,double construction,
                             double inside,double outside)
        : centre_(centre),radius2_(radius*radius),
          construction_(construction),inside_(inside),outside_(outside) {}
    double eval(const Vec3& point) const override {
        if (!algo::t3a_detail::refinementSizingPhaseForTesting())
            return construction_;
        return (point-centre_).norm2()<=radius2_?inside_:outside_;
    }
private:
    Vec3 centre_{};
    double radius2_=0.0,construction_=1.0,inside_=1.0,outside_=1.0;
};

class ExactPointSetField final : public sizing::SizeField {
public:
    ExactPointSetField(std::vector<Vec3> points,double construction,
                       double selected,double outside)
        : construction_(construction),selected_(selected),outside_(outside) {
        for (const auto point:points) selectedBits_.insert(bits(point));
    }
    double eval(const Vec3& point) const override {
        if (!algo::t3a_detail::refinementSizingPhaseForTesting())
            return construction_;
        return selectedBits_.count(bits(point))?selected_:outside_;
    }
private:
    static std::array<std::uint64_t,3> bits(Vec3 point) noexcept {
        const auto canonical=[](double value) {
            if (value==0.0) value=0.0;
            return std::bit_cast<std::uint64_t>(value);
        };
        return {canonical(point.x),canonical(point.y),canonical(point.z)};
    }
    std::set<std::array<std::uint64_t,3>> selectedBits_;
    double construction_=1.0,selected_=1.0,outside_=1.0;
};

bool testViolationKeyLess(const algo::TetViolationKey& a,
                          const algo::TetViolationKey& b) noexcept {
    if (a.kind!=b.kind) return a.kind<b.kind;
    if (a.tet<b.tet) return true;
    if (b.tet<a.tet) return false;
    if (a.hasSizeEdge!=b.hasSizeEdge) return a.hasSizeEdge>b.hasSizeEdge;
    return a.sizeEdge<b.sizeEdge;
}

struct ObservableGeneration {
    algo::TetViolationKey key{};
    std::uint64_t generation=0;
};

bool observableGenerationLess(const ObservableGeneration& a,
                              const ObservableGeneration& b) noexcept {
    if (testViolationKeyLess(a.key,b.key)) return true;
    if (testViolationKeyLess(b.key,a.key)) return false;
    return a.generation<b.generation;
}

void requireProcessedXorUnattempted(const algo::T3AReport& report) {
    std::vector<ObservableGeneration> processed;
    processed.reserve(report.attemptsUsed);
    std::set<std::uint64_t> attemptOrdinals;
    for (const auto& attempt:report.attempts()) {
        CHECK(attemptOrdinals.insert(attempt.attemptOrdinal).second);
        processed.push_back({attempt.violation,attempt.violationGeneration});
    }
    std::sort(processed.begin(),processed.end(),observableGenerationLess);
    processed.erase(std::unique(processed.begin(),processed.end(),
        [](const auto& a,const auto& b) {
            return a.key==b.key&&a.generation==b.generation;
        }),processed.end());
    for (std::size_t i=1;i<processed.size();++i)
        CHECK(observableGenerationLess(processed[i-1],processed[i]));

    const auto unattempted=std::span<const algo::TetViolationKey>(
        report.unattemptedStorage.data(),
        static_cast<std::size_t>(report.unattemptedUsed));
    for (std::size_t i=1;i<unattempted.size();++i)
        CHECK(testViolationKeyLess(unattempted[i-1],unattempted[i]));
    for (const auto& blocked:report.blocked()) {
        CHECK(std::none_of(unattempted.begin(),unattempted.end(),
            [&](const auto& pending) { return pending==blocked.violation; }));
        CHECK(std::count_if(processed.begin(),processed.end(),
            [&](const auto& selected) {
                return selected.key==blocked.violation&&
                       selected.generation==blocked.violationGeneration;
            })==1);
    }
    if (report.hasTerminalViolation)
        CHECK(std::none_of(unattempted.begin(),unattempted.end(),
            [&](const auto& pending) {
                return pending==report.terminalViolation;
            }));
}

void requireLiveGenerationXor(
    const algo::T3AReport& report,
    const algo::t3a_detail::T3ABlockerDigestTrace& trace) {
    requireProcessedXorUnattempted(report);
    CHECK(trace.finalLiveComplete);
    CHECK(trace.finalLiveUsed==report.unattemptedUsed);
    if (!trace.finalLiveComplete||
        trace.finalLiveUsed!=report.unattemptedUsed) return;
    std::vector<ObservableGeneration> processed;
    processed.reserve(report.attemptsUsed);
    for (const auto& attempt:report.attempts())
        processed.push_back({attempt.violation,attempt.violationGeneration});
    std::sort(processed.begin(),processed.end(),observableGenerationLess);
    processed.erase(std::unique(processed.begin(),processed.end(),
        [](const auto& a,const auto& b) {
            return a.key==b.key&&a.generation==b.generation;
        }),processed.end());
    for (std::size_t i=0;i<trace.finalLiveUsed;++i) {
        const auto& live=trace.finalLive[i];
        CHECK(live.key==report.unattemptedStorage[i]);
        CHECK(std::none_of(processed.begin(),processed.end(),
            [&](const auto& selected) {
                return selected.key==live.key&&
                       selected.generation==live.generation;
            }));
        if (i!=0) {
            const auto& prior=trace.finalLive[i-1];
            CHECK(testViolationKeyLess(prior.key,live.key)||
                  (!(testViolationKeyLess(live.key,prior.key))&&
                   prior.generation<live.generation));
        }
    }
}

const algo::t3a_detail::T3ABlockerStateDigest* findBlockerDigest(
    const algo::t3a_detail::T3ABlockerDigestTrace& trace,
    algo::t3a_detail::T3ABlockerDigestCheckpoint checkpoint,
    std::size_t occurrence) noexcept {
    for (std::size_t i=0;i<trace.used;++i) {
        if (trace.rows[i].checkpoint!=checkpoint) continue;
        if (occurrence--==0) return &trace.rows[i];
    }
    return nullptr;
}

void requireSameBlockerPreparationState(
    const algo::t3a_detail::T3ABlockerStateDigest& before,
    const algo::t3a_detail::T3ABlockerStateDigest& after,
    bool terminalFinalizationAllowed) {
    CHECK(before.topologySemantic==after.topologySemantic);
    CHECK(before.topologyStrict==after.topologyStrict);
    CHECK(before.topologyRevision==after.topologyRevision);
    CHECK(before.poolState==after.poolState);
    CHECK(before.queueStructural==after.queueStructural);
    CHECK(before.queueLogical==after.queueLogical);
    CHECK(before.heapStructural==after.heapStructural);
    CHECK(before.heapLogical==after.heapLogical);
    CHECK(before.positionStructural==after.positionStructural);
    CHECK(before.history==after.history);
    CHECK(before.blockerMetrics==after.blockerMetrics);
    CHECK(before.attemptsUsed==after.attemptsUsed);
    CHECK(before.protectedWitnessesUsed==after.protectedWitnessesUsed);
    CHECK(before.blockedUsed==after.blockedUsed);
    CHECK(before.currentBlocked==after.currentBlocked);
    CHECK(before.committedPrefixUsed==after.committedPrefixUsed);
    if (!terminalFinalizationAllowed) {
        CHECK(before.reportState==after.reportState);
        CHECK(before.structural==after.structural);
        CHECK(before.logical==after.logical);
    }
}

void appendViolationKey(std::vector<std::uint64_t>& words,
                        const algo::TetViolationKey& key) {
    words.push_back(static_cast<std::uint64_t>(key.kind));
    for (const auto vertex:key.tet.vertices) {
        words.push_back(static_cast<std::uint64_t>(vertex.kind));
        words.push_back(vertex.value);
    }
    words.push_back(key.hasSizeEdge?1:0);
    for (const auto vertex:key.sizeEdge.vertices) {
        words.push_back(static_cast<std::uint64_t>(vertex.kind));
        words.push_back(vertex.value);
    }
}

template<class Tag>
void appendLocalHandle(std::vector<std::uint64_t>& words,
                       topo::Handle<Tag> handle) {
    words.push_back(static_cast<std::uint32_t>(handle.slot));
    words.push_back(handle.gen);
}

void appendCoordinate(std::vector<std::uint64_t>& words,Vec3 point) {
    words.push_back(std::bit_cast<std::uint64_t>(point.x));
    words.push_back(std::bit_cast<std::uint64_t>(point.y));
    words.push_back(std::bit_cast<std::uint64_t>(point.z));
}

std::vector<std::uint64_t> normalizedReportWords(
    const algo::T3AReport& report) {
    std::vector<std::uint64_t> words;
    words.reserve(static_cast<std::size_t>(
        32+report.attemptsUsed*40+report.blockedUsed*28+
        report.topologyHandlesUsed*20));
    words.insert(words.end(),{
        report.envelopePrepared?1ull:0ull,report.reportComplete?1ull:0ull,
        report.attemptStorage.size(),report.protectedWitnessStorage.size(),
        report.blockedStorage.size(),report.unattemptedStorage.size(),
        report.committedPrefixStorage.size(),
        report.topologyHandleStorage.size(),report.attemptsUsed,
        report.protectedWitnessesUsed,report.blockedUsed,
        report.unattemptedUsed,report.committedPrefixUsed,
        report.topologyHandlesUsed,report.hasTerminalViolation?1ull:0ull});
    if (report.hasTerminalViolation)
        appendViolationKey(words,report.terminalViolation);
    for (const auto& row:report.attempts()) {
        words.push_back(row.attemptOrdinal);
        appendViolationKey(words,row.violation);
        words.push_back(row.violationGeneration);
        words.push_back(static_cast<std::uint64_t>(row.candidateKind));
        words.push_back(row.hasFace?1:0);
        for (const auto vertex:row.face.vertices) {
            words.push_back(static_cast<std::uint64_t>(vertex.kind));
            words.push_back(vertex.value);
        }
        appendCoordinate(words,row.proposedCoordinate);
        words.push_back(static_cast<std::uint64_t>(row.disposition));
        words.push_back(static_cast<std::uint64_t>(row.reason));
        words.push_back(row.witnessOffset);
        words.push_back(row.witnessCount);
        words.push_back(row.cavityDigest);
        words.push_back(row.workBefore);
        words.push_back(row.workAfter);
    }
    for (std::uint64_t i=0;i<report.protectedWitnessesUsed;++i) {
        const auto& row=report.protectedWitnessStorage[
            static_cast<std::size_t>(i)];
        words.push_back(static_cast<std::uint64_t>(row.kind));
        words.push_back(row.edge.canonicalGlobalNodes[0]);
        words.push_back(row.edge.canonicalGlobalNodes[1]);
        words.push_back(row.facet.facePersistentId);
        words.push_back(row.facet.normalizedLocalCellOrdinal);
    }
    for (const auto& row:report.blocked()) {
        appendViolationKey(words,row.violation);
        words.push_back(row.violationGeneration);
        words.push_back(row.attemptOrdinal);
        words.push_back(row.witnessOffset);
        words.push_back(row.witnessCount);
        words.push_back(row.invalidatedByInsertionOrdinal);
    }
    for (std::uint64_t i=0;i<report.unattemptedUsed;++i)
        appendViolationKey(words,report.unattemptedStorage[
            static_cast<std::size_t>(i)]);
    for (std::uint64_t i=0;i<report.committedPrefixUsed;++i) {
        const auto& row=report.committedPrefixStorage[
            static_cast<std::size_t>(i)];
        words.push_back(row.insertionOrdinal);
        appendViolationKey(words,row.violation);
        words.push_back(row.stagedLocalId);
        appendLocalHandle(words,row.predictedNode);
        appendLocalHandle(words,row.actualNode);
        words.push_back(row.topologyOffset);
        words.push_back(row.topologyCount);
        words.push_back(row.postCommitRevision);
    }
    for (std::uint64_t i=0;i<report.topologyHandlesUsed;++i) {
        const auto& row=report.topologyHandleStorage[
            static_cast<std::size_t>(i)];
        words.push_back(static_cast<std::uint64_t>(row.kind));
        words.push_back(static_cast<std::uint64_t>(row.transition));
        appendLocalHandle(words,row.predictedNode);
        appendLocalHandle(words,row.actualNode);
        appendLocalHandle(words,row.predictedCell);
        appendLocalHandle(words,row.actualCell);
        appendLocalHandle(words,row.predictedEdge);
        appendLocalHandle(words,row.actualEdge);
        appendLocalHandle(words,row.predictedFace);
        appendLocalHandle(words,row.actualFace);
    }
    return words;
}

std::vector<std::uint64_t> metricWords(
    const algo::TetRefinementMetrics& value) {
    return {
        value.initialCellsScanned,value.finalCellsAudited,
        value.sizingEvaluations,value.violationEvaluations,
        value.queuePushes,value.queuePops,value.staleQueueRows,
        value.queueHighWater,value.candidateAttempts,
        value.rejectedCandidates,value.protectedIndexVisits,
        value.protectedEdgeTests,value.protectedFacetTests,
        value.locationIndexVisits,value.locationIndexUpdateVisits,
        value.locationIndexUpdateWorkCharged,
        value.coordinatorRowsPrepared,value.coordinatorTailWorkCharged,
        value.locationTetTests,value.cavityCellsVisited,
        value.cavityFaceVisits,value.cavityHighWater,
        value.orient3dTests,value.insphereTests,value.survivorIndexVisits,
        value.candidateSurvivorTests,value.preparedTransactions,
        value.committedInsertions,value.blockedRowsRecorded,
        value.blockedRowsInvalidated,value.blockedRowsCurrent,
        value.blockedRowsCurrentHighWater,value.contextDeltaUpdates,
        value.contextRebuilds,value.locallyReevaluatedCells,
        value.perInsertionGlobalTetScans,value.phase3Allocations,
        value.publicationTailAllocations,value.workUnitsConsumed,
        value.transientAccountedHighWaterBytes,
        value.authoritativeAccountedBytes,value.publicationAccountedBytes,
        value.runComplete?1ull:0ull,value.qualityComplete?1ull:0ull};
}

std::vector<std::uint64_t> immutableBoundaryWords(
    const algo::T3AResult& result) {
    std::vector<std::uint64_t> words;
    words.reserve(result.immutableBoundary.size()*10);
    for (const auto& row:result.immutableBoundary) {
        words.push_back(row.facetKey.facePersistentId);
        words.push_back(row.facetKey.normalizedLocalCellOrdinal);
        for (const auto node:row.directedGlobalNodes) words.push_back(node);
        words.push_back(static_cast<std::uint64_t>(row.classification.dim));
        words.push_back(static_cast<std::uint64_t>(row.classification.tag));
        appendLocalHandle(words,row.authoritativeFace);
    }
    return words;
}

std::vector<std::uint64_t> failureWords(const algo::T3AFailure& failure) {
    std::vector<std::uint64_t> words{
        static_cast<std::uint64_t>(failure.reason),
        static_cast<std::uint64_t>(failure.phase),
        failure.hasViolation?1ull:0ull};
    if (failure.hasViolation) appendViolationKey(words,failure.violation);
    words.push_back(failure.hasProtectedWitness?1ull:0ull);
    if (failure.hasProtectedWitness) {
        words.push_back(static_cast<std::uint64_t>(
            failure.protectedWitness.kind));
        words.push_back(
            failure.protectedWitness.edge.canonicalGlobalNodes[0]);
        words.push_back(
            failure.protectedWitness.edge.canonicalGlobalNodes[1]);
        words.push_back(failure.protectedWitness.facet.facePersistentId);
        words.push_back(
            failure.protectedWitness.facet.normalizedLocalCellOrdinal);
    }
    words.push_back(failure.insertionOrdinal);
    words.push_back(failure.completedWorkUnits);
    return words;
}

void multiInsertionPublicRoute() {
    TEST("T3-A public route commits four independent exact-edge violations");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);

    sizing::MeshSpec baselineSpec;
    baselineSpec.field=std::make_shared<sizing::ConstantField>(0.2);
    auto baselineRequest=requestFor(provider);
    baselineRequest.conformingRequest.expectedSurfaceSetStamp={9903,1};
    baselineRequest.refinement.maxEdgeToTarget=1.0e9;
    baselineRequest.refinement.maxRadiusEdge=1.0e9;
    baselineRequest.refinement.maxInsertions=0;
    baselineRequest.refinement.maxWorkUnits=1000000;
    baselineRequest.refinement.maxOutputNodes=100000;
    baselineRequest.refinement.maxOutputTets=100000;
    baselineRequest.refinement.maxReportRows=5;
    baselineRequest.refinement.maxReportWitnesses=200000;
    baselineRequest.refinement.maxReportTopologyHandles=1000000;
    mesh::Mesh baselineDestination;
    const auto baseline=algo::meshModelRefinedTet4(
        *model,baselineSpec,baselineRequest,baselineDestination,provider);
    CHECK(baseline.published());
    CHECK(baseline.retainedResult&&baseline.retainedResult->authoritative&&
          baseline.retainedResult->authoritative->workspace);
    if (!baseline.retainedResult||!baseline.retainedResult->authoritative||
        !baseline.retainedResult->authoritative->workspace) return;

    const auto& workspace=
        *baseline.retainedResult->authoritative->workspace;
    static constexpr std::array<std::array<std::size_t,2>,6> edges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    struct CandidateEdge {
        std::array<topo::NodeHandle,2> nodes{};
        Vec3 midpoint{};
    };
    std::set<std::array<topo::NodeHandle,2>> seen;
    std::vector<CandidateEdge> candidates;
    for (const auto cell:workspace.cells().activeBySequence()) {
        const auto data=workspace.cells().get(cell);
        if (data.type!=topo::ElementType::Tet4) continue;
        const std::array<topo::NodeHandle,4> nodes{
            data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
        for (const auto& local:edges) {
            std::array<topo::NodeHandle,2> edge{
                nodes[local[0]],nodes[local[1]]};
            if (edge[1]<edge[0]) std::swap(edge[0],edge[1]);
            if (!seen.insert(edge).second) continue;
            const auto a=workspace.nodes().get(edge[0]);
            const auto b=workspace.nodes().get(edge[1]);
            if ((b.coord-a.coord).norm2()<=0.16) continue;
            const Vec3 midpoint=(a.coord+b.coord)*0.5;
            const double boundaryDistance=(std::min)({
                midpoint.x,1.0-midpoint.x,midpoint.y,1.0-midpoint.y,
                midpoint.z,1.0-midpoint.z});
            if (boundaryDistance<=0.15) continue;
            candidates.push_back({edge,midpoint});
        }
    }
    std::set<topo::NodeHandle> usedNodes;
    std::vector<Vec3> screened;
    for (const auto& candidate:candidates) {
        if (usedNodes.count(candidate.nodes[0])||
            usedNodes.count(candidate.nodes[1])) continue;
        bool separated=true;
        for (const auto prior:screened)
            separated=separated&&
                (candidate.midpoint-prior).norm2()>0.01;
        if (!separated) continue;
        screened.push_back(candidate.midpoint);
        usedNodes.insert(candidate.nodes[0]);
        usedNodes.insert(candidate.nodes[1]);
        if (screened.size()==8) break;
    }
    std::vector<Vec3> selectedMidpoints;
    if (screened.size()==8)
        for (const auto index:std::array<std::size_t,4>{{0,2,4,7}})
            selectedMidpoints.push_back(screened[index]);
    std::printf("  t3a exact-edge selected=%zu seedNodes=%zu seedCells=%zu\n",
        selectedMidpoints.size(),workspace.nodes().activeBySequence().size(),
        workspace.cells().activeBySequence().size());
    CHECK(selectedMidpoints.size()==4);
    if (selectedMidpoints.size()!=4) return;

    sizing::MeshSpec refinedSpec;
    refinedSpec.field=std::make_shared<ExactPointSetField>(
        selectedMidpoints,0.2,0.15,1.0);
    auto request=requestFor(provider);
    request.conformingRequest.expectedSurfaceSetStamp={9904,1};
    request.refinement.maxEdgeToTarget=2.5;
    request.refinement.maxRadiusEdge=1.0e9;
    request.refinement.maxInsertions=8;
    request.refinement.maxWorkUnits=1000000;
    request.refinement.maxOutputNodes=100000;
    request.refinement.maxOutputTets=100000;
    request.refinement.maxReportRows=45;
    request.refinement.maxReportWitnesses=200000;
    request.refinement.maxReportTopologyHandles=1000008;
    mesh::Mesh destination;
    const auto refined=algo::meshModelRefinedTet4(
        *model,refinedSpec,request,destination,provider);
    CHECK(refined.status==algo::T3AStatus::Incomplete||refined.published());
    CHECK(refined.retainedResult!=nullptr);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.committedInsertions>=4);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.perInsertionGlobalTetScans==0);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.contextRebuilds==0);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.phase3Allocations==0);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.publicationTailAllocations==0);
    CHECK(refined.retainedResult&&
          refined.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(refined.retainedResult&&
          refined.retainedResult->boundaryAudit.unexpectedExterior==0);
    CHECK(refined.report!=nullptr);
    CHECK(refined.report&&refined.report->blockedStorage.size()==
          refined.report->attemptStorage.size());
    CHECK(refined.report&&refined.report->blockedStorage.size()!=
          request.refinement.maxOutputTets);
    bool sawBlock=false,sawCommitAfterBlock=false;
    if (refined.report) {
        for (const auto& attempt:refined.report->attempts()) {
            if (attempt.disposition==algo::T3AAttemptDisposition::Encroached)
                sawBlock=true;
            if (sawBlock&&attempt.disposition==
                    algo::T3AAttemptDisposition::Committed)
                sawCommitAfterBlock=true;
        }
        std::uint64_t current=0,invalidated=0;
        for (const auto& blocked:refined.report->blocked()) {
            CHECK(blocked.attemptOrdinal<refined.report->attemptsUsed);
            if (blocked.attemptOrdinal>=refined.report->attemptsUsed) continue;
            const auto& attempt=refined.report->attemptStorage[
                static_cast<std::size_t>(blocked.attemptOrdinal)];
            CHECK(attempt.disposition==
                  algo::T3AAttemptDisposition::Encroached);
            CHECK(attempt.violation==blocked.violation);
            CHECK(attempt.violationGeneration==blocked.violationGeneration);
            CHECK(attempt.witnessOffset==blocked.witnessOffset);
            CHECK(attempt.witnessCount==blocked.witnessCount);
            if (blocked.invalidatedByInsertionOrdinal==
                    algo::kCurrentBlockedInsertion) ++current;
            else {
                ++invalidated;
                CHECK(blocked.invalidatedByInsertionOrdinal<
                      refined.retainedResult->metrics.committedInsertions);
            }
        }
        CHECK(refined.retainedResult->metrics.blockedRowsRecorded==
              refined.report->blockedUsed);
        CHECK(refined.retainedResult->metrics.blockedRowsCurrent==current);
        CHECK(refined.retainedResult->metrics.blockedRowsInvalidated==invalidated);
    }
    CHECK(sawBlock);
    CHECK(sawCommitAfterBlock);
    CHECK(refined.retainedResult&&
          refined.retainedResult->metrics.blockedRowsInvalidated>0);
    std::printf("  t3a exact-edge commits=%llu blocked=%llu invalidated=%llu "
                "attempts=%llu work=%llu\n",
        static_cast<unsigned long long>(refined.retainedResult?
            refined.retainedResult->metrics.committedInsertions:0),
        static_cast<unsigned long long>(refined.report?
            refined.report->blockedUsed:0),
        static_cast<unsigned long long>(refined.retainedResult?
            refined.retainedResult->metrics.blockedRowsInvalidated:0),
        static_cast<unsigned long long>(refined.report?
            refined.report->attemptsUsed:0),
        static_cast<unsigned long long>(refined.retainedResult?
            refined.retainedResult->metrics.workUnitsConsumed:0));

    mesh::Mesh invalidationFaultDestination;
    invalidationFaultDestination.addNode({-9,-8,-7},{0,9906});
    const auto invalidationFaultBefore=
        mesh::meshHash(invalidationFaultDestination);
    auto& invalidationTrace=
        algo::t3a_detail::blockerDigestTraceForTesting();
    const auto invalidationMask=
        (std::uint32_t{1}<<static_cast<std::uint32_t>(
            algo::t3a_detail::T3ABlockerDigestCheckpoint::
                InvalidationBeforeTransaction))|
        (std::uint32_t{1}<<static_cast<std::uint32_t>(
            algo::t3a_detail::T3ABlockerDigestCheckpoint::
                InvalidationFailureAfterAbandon));
    invalidationTrace.reset(true,invalidationMask);
    algo::t3a_detail::blockerInvalidationFaultCountdownForTesting()=1;
    const auto invalidationFault=algo::meshModelRefinedTet4(
        *model,refinedSpec,request,invalidationFaultDestination,provider);
    algo::t3a_detail::blockerInvalidationFaultCountdownForTesting()=0;
    CHECK(invalidationFault.status==algo::T3AStatus::Failed);
    CHECK(invalidationFault.failure.reason==algo::T3AReason::AllocationFailed);
    CHECK(invalidationFault.failure.phase==
          algo::T3APhase::CoordinatorPreparation);
    CHECK(invalidationFault.destinationUnchanged);
    CHECK(mesh::meshHash(invalidationFaultDestination)==
          invalidationFaultBefore);
    CHECK(invalidationFault.report&&invalidationFault.report->reportComplete);
    CHECK(invalidationFault.report&&invalidationFault.report->blockedUsed>0);
    CHECK(invalidationFault.report&&invalidationFault.report->attemptsUsed>0);
    if (invalidationFault.report&&
        invalidationFault.report->attemptsUsed>0) {
        const auto& terminal=invalidationFault.report->attemptStorage[
            static_cast<std::size_t>(
                invalidationFault.report->attemptsUsed-1)];
        CHECK(terminal.disposition==
              algo::T3AAttemptDisposition::TerminalFailure);
        CHECK(terminal.reason==algo::T3AReason::AllocationFailed);
    }
    const auto* invalidationBefore=findBlockerDigest(
        invalidationTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            InvalidationBeforeTransaction);
    const auto* invalidationAfter=findBlockerDigest(
        invalidationTrace,algo::t3a_detail::T3ABlockerDigestCheckpoint::
            InvalidationFailureAfterAbandon);
    CHECK(invalidationBefore!=nullptr);
    CHECK(invalidationAfter!=nullptr);
    if (invalidationBefore&&invalidationAfter)
        requireSameBlockerPreparationState(
            *invalidationBefore,*invalidationAfter,false);
    invalidationTrace.reset(false);

    mesh::Mesh invalidationRetryDestination;
    const auto invalidationRetry=algo::meshModelRefinedTet4(
        *model,refinedSpec,request,invalidationRetryDestination,provider);
    CHECK(invalidationRetry.status==refined.status);
    CHECK(invalidationRetry.failure.reason==refined.failure.reason);
    CHECK(invalidationRetry.retainedResult&&refined.retainedResult&&
          invalidationRetry.retainedResult->metrics.committedInsertions==
              refined.retainedResult->metrics.committedInsertions);
    CHECK(invalidationRetry.report&&refined.report&&
          invalidationRetry.report->attemptsUsed==refined.report->attemptsUsed);
    CHECK(invalidationRetry.report&&refined.report&&
          invalidationRetry.report->blockedUsed==refined.report->blockedUsed);
    if (invalidationRetry.retainedResult&&refined.retainedResult&&
        invalidationRetry.retainedResult->authoritative&&
        refined.retainedResult->authoritative&&
        invalidationRetry.retainedResult->authoritative->workspace&&
        refined.retainedResult->authoritative->workspace) {
        CHECK(paacs_tests::observableDigest(
                  *invalidationRetry.retainedResult->authoritative->workspace,
                  false,false)==
              paacs_tests::observableDigest(
                  *refined.retainedResult->authoritative->workspace,
                  false,false));
    }
    if (refined.report) requireProcessedXorUnattempted(*refined.report);
    if (invalidationFault.report)
        requireProcessedXorUnattempted(*invalidationFault.report);
    if (invalidationRetry.report)
        requireProcessedXorUnattempted(*invalidationRetry.report);

    std::vector<Vec3> mixedMidpoints(screened.begin(),screened.begin()+4);
    sizing::MeshSpec mixedSpec;
    mixedSpec.field=std::make_shared<ExactPointSetField>(
        mixedMidpoints,0.2,0.15,1.0);
    auto mixedRequest=request;
    mixedRequest.conformingRequest.expectedSurfaceSetStamp={9905,1};
    mesh::Mesh mixedDestination;
    const auto mixed=algo::meshModelRefinedTet4(
        *model,mixedSpec,mixedRequest,mixedDestination,provider);
    CHECK(mixed.status==algo::T3AStatus::Incomplete);
    CHECK(mixed.failure.reason==algo::T3AReason::InvalidCavity);
    CHECK(mixed.report&&mixed.report->blockedUsed>0);
    CHECK(mixed.report&&mixed.report->hasTerminalViolation);
    CHECK(mixed.report&&mixed.report->attemptsUsed>0);
    if (mixed.report&&mixed.report->attemptsUsed>0) {
        const auto& terminal=mixed.report->attemptStorage[
            static_cast<std::size_t>(mixed.report->attemptsUsed-1)];
        CHECK(terminal.reason==algo::T3AReason::InvalidCavity);
        CHECK(terminal.violation==mixed.report->terminalViolation);
    }
}

void multiInsertionScaling() {
    TEST("T3-A four-row public multi-insertion and local-delta work scaling");
    const std::array<std::size_t,4> contexts{{16,64,256,1024}};
    std::array<std::uint64_t,4> affected{},coordinatorWork{},
                                removedBytes{},addedBytes{},
                                indexUpdates{},indexQueries{};
    for (std::size_t row=0;row<contexts.size();++row) {
        CHECK(algo::t3a_detail::coordinatorDeltaScalingForTesting(
            contexts[row],affected[row],coordinatorWork[row],
            removedBytes[row],addedBytes[row]));
        CHECK(algo::t3a_detail::dynamicLocationScalingForTesting(
            2*contexts[row],indexUpdates[row],indexQueries[row]));
        CHECK(affected[row]>0);
        CHECK(coordinatorWork[row]>0);
        CHECK(removedBytes[row]>0);
        CHECK(addedBytes[row]>0);
        CHECK(indexUpdates[row]>0);
        CHECK(indexQueries[row]>0);
    }
    for (std::size_t row=1;row<contexts.size();++row) {
        CHECK(affected[row]==affected[0]);
        CHECK(removedBytes[row]==removedBytes[0]);
        CHECK(addedBytes[row]==addedBytes[0]);
        CHECK(coordinatorWork[row]<12*coordinatorWork[row-1]);
        CHECK(indexUpdates[row]<12*indexUpdates[row-1]);
        CHECK(indexQueries[row]<12*indexQueries[row-1]);
    }
    std::printf("  t3a delta scaling contexts=%zu/%zu/%zu/%zu "
                "work=%llu/%llu/%llu/%llu updates=%llu/%llu/%llu/%llu "
                "queries=%llu/%llu/%llu/%llu\n",
        contexts[0],contexts[1],contexts[2],contexts[3],
        static_cast<unsigned long long>(coordinatorWork[0]),
        static_cast<unsigned long long>(coordinatorWork[1]),
        static_cast<unsigned long long>(coordinatorWork[2]),
        static_cast<unsigned long long>(coordinatorWork[3]),
        static_cast<unsigned long long>(indexUpdates[0]),
        static_cast<unsigned long long>(indexUpdates[1]),
        static_cast<unsigned long long>(indexUpdates[2]),
        static_cast<unsigned long long>(indexUpdates[3]),
        static_cast<unsigned long long>(indexQueries[0]),
        static_cast<unsigned long long>(indexQueries[1]),
        static_cast<unsigned long long>(indexQueries[2]),
        static_cast<unsigned long long>(indexQueries[3]));

    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec baselineSpec;
    baselineSpec.field=std::make_shared<sizing::ConstantField>(0.2);
    auto baselineRequest=requestFor(provider);
    baselineRequest.conformingRequest.expectedSurfaceSetStamp={9910,1};
    baselineRequest.refinement.maxEdgeToTarget=1.0e9;
    baselineRequest.refinement.maxRadiusEdge=1.0e9;
    baselineRequest.refinement.maxInsertions=0;
    baselineRequest.refinement.maxWorkUnits=1000000;
    baselineRequest.refinement.maxOutputNodes=100000;
    baselineRequest.refinement.maxOutputTets=100000;
    baselineRequest.refinement.maxReportRows=5;
    baselineRequest.refinement.maxReportWitnesses=200000;
    baselineRequest.refinement.maxReportTopologyHandles=1000000;
    mesh::Mesh baselineDestination;
    const auto baseline=algo::meshModelRefinedTet4(
        *model,baselineSpec,baselineRequest,baselineDestination,provider);
    CHECK(baseline.published());
    CHECK(baseline.retainedResult&&baseline.retainedResult->authoritative&&
          baseline.retainedResult->authoritative->workspace);
    if (!baseline.retainedResult||!baseline.retainedResult->authoritative||
        !baseline.retainedResult->authoritative->workspace) return;
    const auto& workspace=*baseline.retainedResult->authoritative->workspace;
    static constexpr std::array<std::array<std::size_t,2>,6> tetEdges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    struct ScalingEdge {
        std::array<topo::NodeHandle,2> nodes{};
        Vec3 midpoint{};
    };
    std::set<std::array<topo::NodeHandle,2>> seenEdges;
    std::vector<ScalingEdge> candidateEdges;
    for (const auto cell:workspace.cells().activeBySequence()) {
        const auto data=workspace.cells().get(cell);
        if (data.type!=topo::ElementType::Tet4) continue;
        const std::array<topo::NodeHandle,4> nodes{
            data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
        for (const auto& local:tetEdges) {
            std::array<topo::NodeHandle,2> edge{
                nodes[local[0]],nodes[local[1]]};
            if (edge[1]<edge[0]) std::swap(edge[0],edge[1]);
            if (!seenEdges.insert(edge).second) continue;
            const auto a=workspace.nodes().get(edge[0]).coord;
            const auto b=workspace.nodes().get(edge[1]).coord;
            if ((b-a).norm2()<=0.16) continue;
            const auto midpoint=(a+b)*0.5;
            const double boundaryDistance=(std::min)({
                midpoint.x,1.0-midpoint.x,midpoint.y,1.0-midpoint.y,
                midpoint.z,1.0-midpoint.z});
            if (boundaryDistance>0.15)
                candidateEdges.push_back({edge,midpoint});
        }
    }
    std::set<topo::NodeHandle> usedNodes;
    std::vector<Vec3> screened;
    for (const auto& candidate:candidateEdges) {
        if (usedNodes.count(candidate.nodes[0])||
            usedNodes.count(candidate.nodes[1])) continue;
        bool separated=true;
        for (const auto prior:screened)
            separated=separated&&
                (candidate.midpoint-prior).norm2()>0.01;
        if (!separated) continue;
        screened.push_back(candidate.midpoint);
        usedNodes.insert(candidate.nodes[0]);
        usedNodes.insert(candidate.nodes[1]);
    }
    std::printf("  t3a scaling candidate midpoints=%zu\n",
                candidateEdges.size());
    const std::vector<Vec3> basePoints{
        screened[0],screened[2],screened[4],screened[7]};
    const std::array<std::vector<std::size_t>,4> extraSets{{
        {8},{9,24},{38,46,53,37},{38,46,53,52,67}}};
    std::array<std::uint64_t,4> publicCommits{},publicWork{};
    ReversedFaceModel reversed(*model);
    const std::array<int,3> threadCounts{{1,4,0}};
    for (std::size_t scalingRow=0;scalingRow<extraSets.size();++scalingRow) {
        auto points=basePoints;
        for (const auto extra:extraSets[scalingRow]) {
            CHECK(extra<candidateEdges.size());
            if (extra<candidateEdges.size())
                points.push_back(candidateEdges[extra].midpoint);
        }
        sizing::MeshSpec spec;
        spec.field=std::make_shared<ExactPointSetField>(
            points,0.2,0.15,1.0);
        std::vector<std::uint64_t> referenceReport,referenceMetrics,
                                   referenceBoundary,referenceFailure;
        std::uint64_t referenceTopology=0,referenceBlockerLogical=0;
        algo::T3AReason referenceReason=algo::T3AReason::None;
        for (std::size_t run=0;run<threadCounts.size();++run) {
            topo::NodePool churnNodes;
            topo::CellPool churnCells;
            for (std::size_t i=0;i<run*11;++i)
                (void)churnNodes.add(
                    {static_cast<double>(i),-17.0,23.0},{0,700});
            CHECK(churnCells.stats().live==0);
            const geo::Model& presented=run==1
                ? static_cast<const geo::Model&>(reversed)
                : static_cast<const geo::Model&>(*model);
            Provider runProvider;
            configureProvider(presented,runProvider);
            spec.threads=threadCounts[run];
            core::setThreadCount(threadCounts[run]);
            auto request=requestFor(runProvider);
            request.conformingRequest.expectedSurfaceSetStamp=
                {9920+scalingRow,1};
            request.refinement.maxEdgeToTarget=2.5;
            request.refinement.maxRadiusEdge=1.0e9;
            request.refinement.maxInsertions=20;
            request.refinement.maxWorkUnits=200000;
            request.refinement.maxOutputNodes=100000;
            request.refinement.maxOutputTets=100000;
            request.refinement.maxReportRows=105;
            request.refinement.maxReportWitnesses=200000;
            request.refinement.maxReportTopologyHandles=200020;
            auto& trace=algo::t3a_detail::blockerDigestTraceForTesting();
            const auto finalMask=std::uint32_t{1}<<
                static_cast<std::uint32_t>(
                    algo::t3a_detail::T3ABlockerDigestCheckpoint::Final);
            trace.reset(true,finalMask);
            auto& wall=algo::t3a_detail::wallTimingTraceForTesting();
            wall.reset(run==0);
            auto& sourceWall=algo::boundary_conforming_detail::
                boundaryConformingWallTraceForTesting();
            sourceWall.reset(run==0);
            mesh::Mesh destination;
            const auto outcome=algo::meshModelRefinedTet4(
                presented,spec,request,destination,runProvider);
            const auto* finalDigest=findBlockerDigest(
                trace,algo::t3a_detail::T3ABlockerDigestCheckpoint::Final);
            CHECK(finalDigest!=nullptr);
            CHECK(outcome.status==algo::T3AStatus::Incomplete);
            CHECK(outcome.retainedResult!=nullptr);
            CHECK(outcome.report&&outcome.report->reportComplete);
            if (!outcome.retainedResult||!outcome.report||!finalDigest)
                continue;
            CHECK(outcome.retainedResult->metrics.perInsertionGlobalTetScans==0);
            CHECK(outcome.retainedResult->metrics.contextRebuilds==0);
            CHECK(outcome.retainedResult->metrics.phase3Allocations==0);
            CHECK(outcome.retainedResult->metrics.publicationTailAllocations==0);
            CHECK(outcome.retainedResult->boundaryAudit.orphanExterior==0);
            CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
            CHECK(outcome.sourceReceipt!=nullptr);
            CHECK(outcome.sourceReceipt&&
                  outcome.sourceReceipt->recoveryMetrics.
                      alreadyPresentBatchEntities>0);
            CHECK(outcome.sourceReceipt&&
                  outcome.sourceReceipt->recoveryMetrics.selectedRebuilds==0);
            requireLiveGenerationXor(*outcome.report,trace);
            const auto report=normalizedReportWords(*outcome.report);
            const auto metrics=metricWords(outcome.retainedResult->metrics);
            const auto boundary=immutableBoundaryWords(*outcome.retainedResult);
            const auto failure=failureWords(outcome.failure);
            const auto topology=
                outcome.retainedResult->authoritative->workspace->
                    semanticTopologyHash();
            if (run==0) {
                publicCommits[scalingRow]=
                    outcome.retainedResult->metrics.committedInsertions;
                publicWork[scalingRow]=
                    outcome.retainedResult->metrics.workUnitsConsumed;
                referenceReason=outcome.failure.reason;
                referenceReport=report;
                referenceMetrics=metrics;
                referenceBoundary=boundary;
                referenceFailure=failure;
                referenceTopology=topology;
                referenceBlockerLogical=finalDigest->logical;
                std::printf("  t3a public row=%zu targets=%zu commits=%llu "
                            "blocked=%llu attempts=%llu work=%llu\n",
                    scalingRow,points.size(),
                    static_cast<unsigned long long>(publicCommits[scalingRow]),
                    static_cast<unsigned long long>(outcome.report->blockedUsed),
                    static_cast<unsigned long long>(outcome.report->attemptsUsed),
                    static_cast<unsigned long long>(publicWork[scalingRow]));
                std::printf("    wall_ms source/setup/protect/location/cavity/"
                            "candidate/topology/coordinator/local/final="
                            "%.2f/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f\n",
                    wall.nanoseconds[0]/1.0e6,wall.nanoseconds[1]/1.0e6,
                    wall.nanoseconds[2]/1.0e6,wall.nanoseconds[3]/1.0e6,
                    wall.nanoseconds[4]/1.0e6,wall.nanoseconds[5]/1.0e6,
                    wall.nanoseconds[6]/1.0e6,wall.nanoseconds[7]/1.0e6,
                    wall.nanoseconds[8]/1.0e6,wall.nanoseconds[9]/1.0e6);
                std::printf("    source_wall_ms surface/artifact/envelope/tet/"
                            "import/context/recovery=%.2f/%.2f/%.2f/%.2f/"
                            "%.2f/%.2f/%.2f\n",
                    sourceWall.nanoseconds[0]/1.0e6,
                    sourceWall.nanoseconds[1]/1.0e6,
                    sourceWall.nanoseconds[2]/1.0e6,
                    sourceWall.nanoseconds[3]/1.0e6,
                    sourceWall.nanoseconds[4]/1.0e6,
                    sourceWall.nanoseconds[5]/1.0e6,
                    sourceWall.nanoseconds[6]/1.0e6);
                if (outcome.sourceReceipt)
                    std::printf("    source work artifactPairs=%llu exact=%llu "
                                "recoveryVisits=%llu segTet=%llu targetTet=%llu "
                                "survivor=%llu rebuilds=%llu batch=%llu "
                                "cells=%llu\n",
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                edgeEdgePairs+
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                edgeFacetPairs+
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                facetFacetPairs),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                exactEdgeEdgeTests+
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                exactEdgeFacetTests+
                            outcome.sourceReceipt->artifactPreflightMetrics.
                                exactFacetFacetTests),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.
                                indexNodesVisited),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.segmentTetTests),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.targetTetTests),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.
                                candidateSurvivorTests),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.
                                selectedRebuilds),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.
                                alreadyPresentBatchEntities),
                        static_cast<unsigned long long>(
                            outcome.sourceReceipt->recoveryMetrics.contextCells));
            } else {
                CHECK(outcome.failure.reason==referenceReason);
                CHECK(report==referenceReport);
                CHECK(metrics==referenceMetrics);
                CHECK(boundary==referenceBoundary);
                CHECK(failure==referenceFailure);
                CHECK(topology==referenceTopology);
                CHECK(finalDigest->logical==referenceBlockerLogical);
            }
            trace.reset(false);
        }
        algo::t3a_detail::wallTimingTraceForTesting().reset(false);
        algo::boundary_conforming_detail::
            boundaryConformingWallTraceForTesting().reset(false);
    }
    core::setThreadCount(0);
    CHECK(publicCommits[0]>=4);
    for (std::size_t row=1;row<publicCommits.size();++row) {
        CHECK(publicCommits[row]>publicCommits[row-1]);
        CHECK(publicWork[row]<12*publicWork[row-1]);
    }
}

void publicFixtureScaling() {
    TEST("T3-A four increasing public fixtures satisfy W(4n) scaling");
    auto kernel=geo::makeBuiltinKernel();
    const std::array<double,4> spacings{{0.2,0.16,0.105,0.09}};
    constexpr std::uint64_t insertionCap=20;
    const std::array<int,3> threadCounts{{1,4,0}};
    std::array<std::uint64_t,4> seedCells{},commits{},work{};
    static constexpr std::array<std::array<std::size_t,2>,6> tetEdges{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for (std::size_t row=0;row<spacings.size();++row) {
        const double spacing=spacings[row];
        auto model=kernel->makeBox({0,0,0},{1,1,1});
        Provider provider;
        configureProvider(*model,provider);
        sizing::MeshSpec baselineSpec;
        baselineSpec.field=std::make_shared<sizing::ConstantField>(spacing);
        baselineSpec.threads=1;
        auto baselineRequest=requestFor(provider);
        baselineRequest.conformingRequest.expectedSurfaceSetStamp=
            {9940+row,1};
        baselineRequest.refinement.maxEdgeToTarget=1.0e9;
        baselineRequest.refinement.maxRadiusEdge=1.0e9;
        baselineRequest.refinement.maxInsertions=0;
        baselineRequest.refinement.maxWorkUnits=2000000;
        baselineRequest.refinement.maxOutputNodes=200000;
        baselineRequest.refinement.maxOutputTets=200000;
        baselineRequest.refinement.maxReportRows=5;
        baselineRequest.refinement.maxReportWitnesses=1000000;
        baselineRequest.refinement.maxReportTopologyHandles=2000000;
        mesh::Mesh baselineDestination;
        const auto baseline=algo::meshModelRefinedTet4(
            *model,baselineSpec,baselineRequest,baselineDestination,provider);
        if (!baseline.published())
            std::printf("  t3a W4n baseline row=%zu status=%d reason=%d "
                        "phase=%d\n",row,static_cast<int>(baseline.status),
                static_cast<int>(baseline.failure.reason),
                static_cast<int>(baseline.failure.phase));
        CHECK(baseline.published());
        CHECK(baseline.retainedResult&&baseline.retainedResult->authoritative&&
              baseline.retainedResult->authoritative->workspace);
        if (!baseline.retainedResult||!baseline.retainedResult->authoritative||
            !baseline.retainedResult->authoritative->workspace) return;
        const auto& workspace=
            *baseline.retainedResult->authoritative->workspace;
        seedCells[row]=workspace.cells().activeBySequence().size();
        struct CandidateEdge {
            std::array<topo::NodeHandle,2> nodes{};
            Vec3 midpoint{};
        };
        std::set<std::array<topo::NodeHandle,2>> seen;
        std::vector<CandidateEdge> candidates;
        for (const auto cell:workspace.cells().activeBySequence()) {
            const auto data=workspace.cells().get(cell);
            if (data.type!=topo::ElementType::Tet4) continue;
            const std::array<topo::NodeHandle,4> nodes{
                data.nodes[0],data.nodes[1],data.nodes[2],data.nodes[3]};
            for (const auto& local:tetEdges) {
                std::array<topo::NodeHandle,2> edge{
                    nodes[local[0]],nodes[local[1]]};
                if (edge[1]<edge[0]) std::swap(edge[0],edge[1]);
                if (!seen.insert(edge).second) continue;
                const auto a=workspace.nodes().get(edge[0]).coord;
                const auto b=workspace.nodes().get(edge[1]).coord;
                const double threshold=row==0?4.0:2.25;
                if ((b-a).norm2()<=threshold*spacing*spacing) continue;
                const auto midpoint=(a+b)*0.5;
                const double boundaryDistance=(std::min)({
                    midpoint.x,1.0-midpoint.x,
                    midpoint.y,1.0-midpoint.y,
                    midpoint.z,1.0-midpoint.z});
                if (boundaryDistance>0.75*spacing)
                    candidates.push_back({edge,midpoint});
            }
        }
        std::vector<Vec3> selected;
        if (row==0) {
            std::set<topo::NodeHandle> usedNodes;
            std::vector<Vec3> screened;
            for (const auto& candidate:candidates) {
                if (usedNodes.count(candidate.nodes[0])||
                    usedNodes.count(candidate.nodes[1])) continue;
                bool separated=true;
                for (const auto prior:screened)
                    separated=separated&&
                        (candidate.midpoint-prior).norm2()>0.01;
                if (!separated) continue;
                screened.push_back(candidate.midpoint);
                usedNodes.insert(candidate.nodes[0]);
                usedNodes.insert(candidate.nodes[1]);
            }
            CHECK(screened.size()>7&&candidates.size()>8);
            if (screened.size()<=7||candidates.size()<=8) return;
            selected={screened[0],screened[2],screened[4],screened[7],
                      candidates[8].midpoint};
        } else {
            selected.reserve(candidates.size());
            for (const auto& candidate:candidates) {
                selected.push_back(candidate.midpoint);
                if (selected.size()==1024) break;
            }
        }
        const std::size_t minimumTargets=row==0?5:8;
        CHECK(selected.size()>=minimumTargets);
        if (selected.size()<minimumTargets) return;

        ReversedFaceModel reversed(*model);
        std::vector<std::uint64_t> referenceReport,referenceMetrics,
                                   referenceBoundary,referenceFailure;
        std::uint64_t referenceTopology=0,referenceBlockerLogical=0;
        algo::T3AReason referenceReason=algo::T3AReason::None;
        for (std::size_t run=0;run<threadCounts.size();++run) {
            topo::NodePool churnNodes;
            topo::CellPool churnCells;
            for (std::size_t i=0;i<run*13;++i)
                (void)churnNodes.add(
                    {static_cast<double>(i),-31.0,41.0},{0,703});
            CHECK(churnCells.stats().live==0);
            const geo::Model& presented=run==1
                ? static_cast<const geo::Model&>(reversed)
                : static_cast<const geo::Model&>(*model);
            Provider runProvider;
            configureProvider(presented,runProvider);
            sizing::MeshSpec spec;
            const double selectedSize=(row==0?0.75:0.5)*spacing;
            spec.field=std::make_shared<ExactPointSetField>(
                selected,spacing,selectedSize,1.0);
            spec.threads=threadCounts[run];
            auto request=requestFor(runProvider);
            request.conformingRequest.expectedSurfaceSetStamp={9950+row,1};
            request.refinement.maxEdgeToTarget=2.5;
            request.refinement.maxRadiusEdge=1.0e9;
            // One common nonbinding cap: increasing commit counts must come
            // from the selected violations, never from row-specific budgets.
            request.refinement.maxInsertions=insertionCap;
            request.refinement.maxWorkUnits=2000000;
            request.refinement.maxOutputNodes=200000;
            request.refinement.maxOutputTets=200000;
            request.refinement.maxReportRows=4096;
            request.refinement.maxReportWitnesses=1000000;
            request.refinement.maxReportTopologyHandles=2000020;
            auto& trace=algo::t3a_detail::blockerDigestTraceForTesting();
            const auto finalMask=std::uint32_t{1}<<
                static_cast<std::uint32_t>(
                    algo::t3a_detail::T3ABlockerDigestCheckpoint::Final);
            trace.reset(true,finalMask);
            mesh::Mesh destination;
            const auto outcome=algo::meshModelRefinedTet4(
                presented,spec,request,destination,runProvider);
            if (outcome.status!=algo::T3AStatus::Incomplete&&
                !outcome.published())
                std::printf("  t3a W4n route row=%zu run=%zu status=%d "
                            "reason=%d phase=%d\n",row,run,
                    static_cast<int>(outcome.status),
                    static_cast<int>(outcome.failure.reason),
                    static_cast<int>(outcome.failure.phase));
            const auto* finalDigest=findBlockerDigest(
                trace,algo::t3a_detail::T3ABlockerDigestCheckpoint::Final);
            CHECK(finalDigest!=nullptr);
            CHECK(outcome.status==algo::T3AStatus::Incomplete||
                  outcome.published());
            CHECK(outcome.retainedResult!=nullptr);
            CHECK(outcome.report&&outcome.report->reportComplete);
            if (!outcome.retainedResult||!outcome.report||!finalDigest)
                continue;
            CHECK(outcome.retainedResult->metrics.perInsertionGlobalTetScans==0);
            CHECK(outcome.retainedResult->metrics.contextRebuilds==0);
            CHECK(outcome.retainedResult->metrics.phase3Allocations==0);
            CHECK(outcome.retainedResult->metrics.publicationTailAllocations==0);
            CHECK(outcome.retainedResult->boundaryAudit.orphanExterior==0);
            CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
            requireLiveGenerationXor(*outcome.report,trace);
            const auto report=normalizedReportWords(*outcome.report);
            const auto metrics=metricWords(outcome.retainedResult->metrics);
            const auto boundary=immutableBoundaryWords(*outcome.retainedResult);
            const auto failure=failureWords(outcome.failure);
            const auto topology=outcome.retainedResult->authoritative->
                workspace->semanticTopologyHash();
            if (run==0) {
                commits[row]=
                    outcome.retainedResult->metrics.committedInsertions;
                work[row]=outcome.retainedResult->metrics.workUnitsConsumed;
                referenceReason=outcome.failure.reason;
                referenceReport=report;
                referenceMetrics=metrics;
                referenceBoundary=boundary;
                referenceFailure=failure;
                referenceTopology=topology;
                referenceBlockerLogical=finalDigest->logical;
            } else {
                CHECK(outcome.failure.reason==referenceReason);
                CHECK(report==referenceReport);
                CHECK(metrics==referenceMetrics);
                CHECK(boundary==referenceBoundary);
                CHECK(failure==referenceFailure);
                CHECK(topology==referenceTopology);
                CHECK(finalDigest->logical==referenceBlockerLogical);
            }
            trace.reset(false);
        }
        CHECK(commits[row]>=4);
        std::printf("  t3a W4n row=%zu spacing=%.6f targets=%zu seed=%llu "
                    "commits=%llu work=%llu\n",row,spacing,selected.size(),
            static_cast<unsigned long long>(seedCells[row]),
            static_cast<unsigned long long>(commits[row]),
            static_cast<unsigned long long>(work[row]));
    }
    core::setThreadCount(0);
    for (std::size_t row=1;row<spacings.size();++row) {
        CHECK(seedCells[row]>seedCells[row-1]);
        CHECK(commits[row]>commits[row-1]);
        CHECK(work[row]<12*work[row-1]);
    }
    CHECK(seedCells.back()>=4*seedCells.front());
    CHECK(work.back()<12*work.front());
}
void shapeAndLocalFieldCoverage() {
    TEST("T3-A cylinder thin-box and locally concentrated field coverage");
    auto kernel=geo::makeBuiltinKernel();
    struct Fixture {
        std::unique_ptr<geo::Model> model;
        double h=1.0;
        std::uint64_t stamp=0;
    };
    std::array<Fixture,2> fixtures{{
        {kernel->makeCylinder({0,0,0},0.5,1.0),0.15,9901},
        {kernel->makeBox({0,0,0},{1,1,1e-3}),0.2,9902}}};
    for (auto& fixture:fixtures) {
        Provider provider;
        configureProvider(*fixture.model,provider);
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(fixture.h);
        auto request=requestFor(provider);
        request.conformingRequest.expectedSurfaceSetStamp={fixture.stamp,1};
        request.refinement.maxEdgeToTarget=100.0;
        request.refinement.maxRadiusEdge=1.0e9;
        request.refinement.maxInsertions=0;
        request.refinement.maxWorkUnits=100000;
        request.refinement.maxReportRows=5;
        request.refinement.maxReportWitnesses=200000;
        request.refinement.maxReportTopologyHandles=100000;
        mesh::Mesh destination;
        const auto outcome=algo::meshModelRefinedTet4(
            *fixture.model,spec,request,destination,provider);
        CHECK(outcome.published());
        CHECK(outcome.retainedResult!=nullptr);
        if (!outcome.retainedResult) continue;
        CHECK(outcome.retainedResult->boundaryAudit.orphanExterior==0);
        CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior==0);
        CHECK(outcome.retainedResult->quality.inverted==0);
        CHECK(outcome.retainedResult->quality.nonfinite==0);
        CHECK(outcome.retainedResult->metrics.perInsertionGlobalTetScans==0);
        CHECK(outcome.retainedResult->metrics.contextRebuilds==0);
    }

    auto box=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*box,provider);
    sizing::MeshSpec fieldSpec;
    fieldSpec.field=std::make_shared<BallField>(
        Vec3{0.5,0.5,0.5},0.3,0.2,0.55);
    auto request=requestFor(provider);
    request.refinement.maxEdgeToTarget=0.5;
    request.refinement.maxInsertions=16;
    request.refinement.maxWorkUnits=100000;
    request.refinement.maxReportRows=85;
    request.refinement.maxReportWitnesses=200000;
    request.refinement.maxReportTopologyHandles=100016;
    mesh::Mesh destination;
    const auto local=algo::meshModelRefinedTet4(
        *box,fieldSpec,request,destination,provider);
    CHECK(local.status==algo::T3AStatus::Incomplete);
    CHECK(local.failure.reason==algo::T3AReason::BoundaryRefinementRequired);
    CHECK(local.retainedResult!=nullptr);
    CHECK(local.retainedResult&&
          local.retainedResult->metrics.committedInsertions==1);
    CHECK(local.report&&local.report->blockedUsed>0);
    CHECK(local.retainedResult&&local.report&&
          local.retainedResult->metrics.blockedRowsCurrent==
              local.report->blockedUsed);
    CHECK(local.retainedResult&&
          local.retainedResult->metrics.perInsertionGlobalTetScans==0);
    CHECK(local.retainedResult&&local.retainedResult->metrics.contextRebuilds==0);
    CHECK(local.retainedResult&&local.retainedResult->boundaryAudit.orphanExterior==0);
    CHECK(local.retainedResult&&local.retainedResult->boundaryAudit.unexpectedExterior==0);
    if (local.retainedResult)
        std::printf("  local field status=%d reason=%d phase=%d "
                    "commits=%llu blocked=%llu attempts=%llu "
                    "unattempted=%llu work=%llu\n",
            static_cast<int>(local.status),
            static_cast<int>(local.failure.reason),
            static_cast<int>(local.failure.phase),
            static_cast<unsigned long long>(
                local.retainedResult->metrics.committedInsertions),
            static_cast<unsigned long long>(local.report?
                local.report->blockedUsed:0),
            static_cast<unsigned long long>(local.report?
                local.report->attemptsUsed:0),
            static_cast<unsigned long long>(local.report?
                local.report->unattemptedUsed:0),
            static_cast<unsigned long long>(
                local.retainedResult->metrics.workUnitsConsumed));
}

void sealedCavityBoundary() {
    TEST("T3-A outer shell plus sealed cavity preserves immutable components");
    paacs_tests::fixtures::SealedCavityModel model;
    Provider provider;
    for (const auto* face:model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,1}};
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(10.0);
    auto request=requestFor(provider);
    request.conformingRequest.expectedSurfaceSetStamp={9940,1};
    request.refinement.maxEdgeToTarget=1.0e9;
    request.refinement.maxRadiusEdge=1.0e9;
    request.refinement.maxInsertions=0;
    request.refinement.maxWorkUnits=100000;
    request.refinement.maxReportRows=5;
    request.refinement.maxReportWitnesses=200000;
    request.refinement.maxReportTopologyHandles=100000;
    mesh::Mesh destination;
    const auto outcome=algo::meshModelRefinedTet4(
        model,spec,request,destination,provider);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult!=nullptr);
    if (!outcome.retainedResult||!outcome.retainedResult->authoritative||
        !outcome.retainedResult->authoritative->workspace) return;
    const auto& result=*outcome.retainedResult;
    const auto& workspace=*result.authoritative->workspace;
    CHECK(result.boundaryAudit.orphanExterior==0);
    CHECK(result.boundaryAudit.unexpectedExterior==0);
    CHECK(result.immutableBoundary.size()>=24);
    for (const auto& boundary:result.immutableBoundary) {
        CHECK(workspace.faces().isCurrent(boundary.authoritativeFace));
        if (workspace.faces().isCurrent(boundary.authoritativeFace)) {
            const auto face=workspace.faces().get(boundary.authoritativeFace);
            CHECK(face.incident.size()==1);
            CHECK(face.cls.dim==boundary.classification.dim&&
                  face.cls.tag==boundary.classification.tag);
        }
    }
}

void dynamicLocationScaling() {
    TEST("T3-A evolving location index applies many deltas subquadratically");
    CHECK(algo::t3a_detail::blockedSchedulerReuseForTesting());
    const std::array<std::size_t,3> contextSizes{{16,64,256}};
    std::array<std::uint64_t,3> affectedRows{},coordinatorWork{},
                                accountedRemoved{},accountedAdded{};
    for (std::size_t i=0;i<contextSizes.size();++i)
        CHECK(algo::t3a_detail::coordinatorDeltaScalingForTesting(
            contextSizes[i],affectedRows[i],coordinatorWork[i],
            accountedRemoved[i],accountedAdded[i]));
    CHECK(affectedRows[0]==affectedRows[1]);
    CHECK(affectedRows[1]==affectedRows[2]);
    CHECK(coordinatorWork[1]>coordinatorWork[0]);
    CHECK(coordinatorWork[2]>coordinatorWork[1]);
    CHECK(coordinatorWork[1]<2*coordinatorWork[0]);
    CHECK(coordinatorWork[2]<2*coordinatorWork[1]);
    CHECK(accountedRemoved[0]==accountedRemoved[1]);
    CHECK(accountedRemoved[1]==accountedRemoved[2]);
    CHECK(accountedAdded[0]==accountedAdded[1]);
    CHECK(accountedAdded[1]==accountedAdded[2]);

    const std::array<std::size_t,3> sizes{{128,512,2048}};
    std::array<std::uint64_t,3> updates{},queries{};
    for (std::size_t i=0;i<sizes.size();++i) {
        CHECK(algo::t3a_detail::dynamicLocationScalingForTesting(
            sizes[i],updates[i],queries[i]));
        CHECK(updates[i]>0);
        CHECK(queries[i]>0);
        CHECK(updates[i]<sizes[i]*64);
        CHECK(queries[i]<sizes[i]*128);
    }
    for (std::size_t i=1;i<sizes.size();++i) {
        CHECK(updates[i]<8*updates[i-1]);
        CHECK(queries[i]<8*queries[i-1]);
    }
    std::printf("  t3a local coordinator rows=%llu work=%llu/%llu/%llu "
                "accounting=%llu->%llu; "
                "dynamic index updates=%llu/%llu/%llu "
                "queries=%llu/%llu/%llu\n",
        static_cast<unsigned long long>(affectedRows[0]),
        static_cast<unsigned long long>(coordinatorWork[0]),
        static_cast<unsigned long long>(coordinatorWork[1]),
        static_cast<unsigned long long>(coordinatorWork[2]),
        static_cast<unsigned long long>(accountedRemoved[0]),
        static_cast<unsigned long long>(accountedAdded[0]),
        static_cast<unsigned long long>(updates[0]),
        static_cast<unsigned long long>(updates[1]),
        static_cast<unsigned long long>(updates[2]),
        static_cast<unsigned long long>(queries[0]),
        static_cast<unsigned long long>(queries[1]),
        static_cast<unsigned long long>(queries[2]));

}

void operationScaling() {
    TEST("T3-A uniform and local-field operation counters reject quadratic growth");
    auto kernel=geo::makeBuiltinKernel();
    const std::array<double,4> sizes{{0.8,0.5,0.315,0.2}};
    std::array<std::uint64_t,4> uniformWork{},uniformCells{},localWork{};
    for (std::size_t i=0;i<sizes.size();++i) {
        auto model=kernel->makeBox({0,0,0},{1,1,1});
        Provider provider;
        configureProvider(*model,provider);
        auto request=requestFor(provider);
        request.conformingRequest.expectedSurfaceSetStamp={9950+i,1};
        request.refinement.maxEdgeToTarget=1.0e9;
        request.refinement.maxRadiusEdge=1.0e9;
        request.refinement.maxInsertions=0;
        request.refinement.maxWorkUnits=1000000;
        request.refinement.maxReportRows=5;
        request.refinement.maxReportWitnesses=200000;
        request.refinement.maxReportTopologyHandles=1000000;
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(sizes[i]);
        mesh::Mesh destination;
        const auto outcome=algo::meshModelRefinedTet4(
            *model,spec,request,destination,provider);
        CHECK(outcome.published());
        CHECK(outcome.retainedResult!=nullptr);
        if (!outcome.retainedResult) return;
        const auto& metrics=outcome.retainedResult->metrics;
        uniformWork[i]=metrics.workUnitsConsumed;
        uniformCells[i]=metrics.initialCellsScanned;
        CHECK(uniformWork[i]>=uniformCells[i]);
        CHECK(metrics.perInsertionGlobalTetScans==0);
        CHECK(metrics.contextRebuilds==0);

        auto localRequest=request;
        localRequest.conformingRequest.expectedSurfaceSetStamp={9960+i,1};
        localRequest.refinement.maxEdgeToTarget=0.5;
        localRequest.refinement.maxInsertions=1;
        localRequest.refinement.maxReportRows=10;
        localRequest.refinement.maxReportTopologyHandles=1000001;
        sizing::MeshSpec localSpec;
        localSpec.field=std::make_shared<BallField>(
            Vec3{0.5,0.5,0.5},0.3,sizes[i]*0.4,sizes[i]);
        mesh::Mesh localDestination;
        const auto local=algo::meshModelRefinedTet4(
            *model,localSpec,localRequest,localDestination,provider);
        CHECK(local.status==algo::T3AStatus::Published||
              local.status==algo::T3AStatus::Incomplete);
        CHECK(local.retainedResult!=nullptr);
        if (!local.retainedResult) return;
        localWork[i]=local.retainedResult->metrics.workUnitsConsumed;
        CHECK(localWork[i]>0);
        CHECK(localWork[i]<20*uniformWork[i]);
        CHECK(local.retainedResult->metrics.perInsertionGlobalTetScans==0);
        CHECK(local.retainedResult->metrics.contextRebuilds==0);
    }
    for (std::size_t i=1;i<sizes.size();++i) {
        CHECK(uniformCells[i]>uniformCells[i-1]);
        CHECK(uniformWork[i]<12*uniformWork[i-1]);
    }
    std::printf("  t3a scaling uniform cells=%llu/%llu/%llu/%llu "
                "work=%llu/%llu/%llu/%llu local=%llu/%llu/%llu/%llu\n",
        static_cast<unsigned long long>(uniformCells[0]),
        static_cast<unsigned long long>(uniformCells[1]),
        static_cast<unsigned long long>(uniformCells[2]),
        static_cast<unsigned long long>(uniformCells[3]),
        static_cast<unsigned long long>(uniformWork[0]),
        static_cast<unsigned long long>(uniformWork[1]),
        static_cast<unsigned long long>(uniformWork[2]),
        static_cast<unsigned long long>(uniformWork[3]),
        static_cast<unsigned long long>(localWork[0]),
        static_cast<unsigned long long>(localWork[1]),
        static_cast<unsigned long long>(localWork[2]),
        static_cast<unsigned long long>(localWork[3]));
}

class PointFaultField final : public sizing::SizeField {
public:
    enum class Fault { Zero, Nonfinite, BadAlloc, Standard, Nonstandard };
    PointFaultField(Vec3 point,Fault fault) : bits_(bits(point)),fault_(fault) {}
    double eval(const Vec3& point) const override {
        if (bits(point)!=bits_) return 100.0;
        switch (fault_) {
            case Fault::Zero: return 0.0;
            case Fault::Nonfinite:
                return (std::numeric_limits<double>::infinity)();
            case Fault::BadAlloc: throw std::bad_alloc();
            case Fault::Standard: throw std::runtime_error("sizing fault");
            case Fault::Nonstandard: throw 17;
        }
        return 100.0;
    }
private:
    static std::array<std::uint64_t,3> bits(Vec3 p) noexcept {
        return {std::bit_cast<std::uint64_t>(p.x),
                std::bit_cast<std::uint64_t>(p.y),
                std::bit_cast<std::uint64_t>(p.z)};
    }
    std::array<std::uint64_t,3> bits_{};
    Fault fault_;
};

class ChangingPointField final : public sizing::SizeField {
public:
    ChangingPointField(Vec3 point,std::uint64_t stableCalls)
        : bits_(bits(point)),stableCalls_(stableCalls) {}
    double eval(const Vec3& point) const override {
        if (bits(point)!=bits_) return 100.0;
        ++hits_;
        return hits_<=stableCalls_?100.0:101.0;
    }
    std::uint64_t hits() const noexcept { return hits_; }
private:
    static std::array<std::uint64_t,3> bits(Vec3 p) noexcept {
        return {std::bit_cast<std::uint64_t>(p.x),
                std::bit_cast<std::uint64_t>(p.y),
                std::bit_cast<std::uint64_t>(p.z)};
    }
    std::array<std::uint64_t,3> bits_{};
    std::uint64_t stableCalls_=0;
    mutable std::uint64_t hits_=0;
};

void sizingTaxonomyAndBudgets() {
    TEST("T3-A sizing taxonomy and every hard envelope/budget boundary");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec baselineSpec;
    baselineSpec.field=std::make_shared<sizing::ConstantField>(100.0);
    mesh::Mesh baselineMesh;
    const auto baseline=algo::meshModelRefinedTet4(
        *model,baselineSpec,requestFor(provider),baselineMesh,provider);
    CHECK(baseline.published());
    CHECK(baseline.retainedResult&&baseline.retainedResult->authoritative&&
          baseline.retainedResult->authoritative->workspace);
    if (!baseline.retainedResult||!baseline.retainedResult->authoritative||
        !baseline.retainedResult->authoritative->workspace) return;
    const auto& ws=*baseline.retainedResult->authoritative->workspace;
    const auto cells=ws.cells().activeBySequence();
    CHECK(!cells.empty());
    if (cells.empty()) return;
    const auto cell=ws.cells().get(cells.front());
    const Vec3 target=(ws.nodes().get(cell.nodes[0]).coord+
                       ws.nodes().get(cell.nodes[1]).coord)*0.5;

    const std::array<std::pair<PointFaultField::Fault,algo::T3AReason>,5> faults{{
        {PointFaultField::Fault::Zero,algo::T3AReason::InvalidSizingEvaluation},
        {PointFaultField::Fault::Nonfinite,algo::T3AReason::InvalidSizingEvaluation},
        {PointFaultField::Fault::BadAlloc,algo::T3AReason::AllocationFailed},
        {PointFaultField::Fault::Standard,algo::T3AReason::SizingFieldThrew},
        {PointFaultField::Fault::Nonstandard,algo::T3AReason::SizingFieldThrew}}};
    for (const auto& [fault,reason]:faults) {
        sizing::MeshSpec spec;
        spec.field=std::make_shared<PointFaultField>(target,fault);
        mesh::Mesh destination;
        destination.addNode({-3,-2,-1},{0,808});
        const auto before=mesh::meshHash(destination);
        const auto outcome=algo::meshModelRefinedTet4(
            *model,spec,requestFor(provider),destination,provider);
        CHECK(outcome.status==algo::T3AStatus::Failed);
        CHECK(outcome.failure.reason==reason);
        CHECK(outcome.failure.phase==algo::T3APhase::SizingEvaluation);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
    }

    auto v3Counter=std::make_shared<ChangingPointField>(
        target,(std::numeric_limits<std::uint64_t>::max)());
    sizing::MeshSpec v3CountSpec; v3CountSpec.field=v3Counter;
    mesh::Mesh v3Mesh;
    const auto v3=algo::meshModelBoundaryConforming(
        *model,v3CountSpec,requestFor(provider).conformingRequest,
        v3Mesh,provider);
    CHECK(v3.published());
    const auto v3Hits=v3Counter->hits();

    auto preparationChange=std::make_shared<ChangingPointField>(target,v3Hits+1);
    sizing::MeshSpec preparationChangeSpec;
    preparationChangeSpec.field=preparationChange;
    mesh::Mesh changedDestination;
    const auto changed=algo::meshModelRefinedTet4(
        *model,preparationChangeSpec,requestFor(provider),changedDestination,
        provider);
    CHECK(changed.status==algo::T3AStatus::Failed);
    CHECK(changed.failure.reason==algo::T3AReason::SizingFieldChanged);
    CHECK(changed.failure.phase==algo::T3APhase::SizingEvaluation);
    CHECK(changed.destinationUnchanged);

    auto fullCounter=std::make_shared<ChangingPointField>(
        target,(std::numeric_limits<std::uint64_t>::max)());
    sizing::MeshSpec fullCountSpec; fullCountSpec.field=fullCounter;
    mesh::Mesh fullCountMesh;
    const auto fullCount=algo::meshModelRefinedTet4(
        *model,fullCountSpec,requestFor(provider),fullCountMesh,provider);
    CHECK(fullCount.published());
    CHECK(fullCounter->hits()>v3Hits+1);
    auto finalChange=std::make_shared<ChangingPointField>(
        target,fullCounter->hits()-1);
    sizing::MeshSpec finalChangeSpec; finalChangeSpec.field=finalChange;
    mesh::Mesh finalDestination;
    const auto finalChanged=algo::meshModelRefinedTet4(
        *model,finalChangeSpec,requestFor(provider),finalDestination,provider);
    CHECK(finalChanged.status==algo::T3AStatus::Failed);
    CHECK(finalChanged.failure.reason==algo::T3AReason::SizingFieldChanged);
    CHECK(finalChanged.failure.phase==algo::T3APhase::FinalQualityAudit);
    CHECK(finalChanged.destinationUnchanged);

    sizing::MeshSpec refiningSpec;
    refiningSpec.field=std::make_shared<sizing::ConstantField>(0.55);
    const auto requireIncomplete=[&](algo::T3ARequest request,
                                     algo::T3AReason reason) {
        mesh::Mesh destination;
        destination.addNode({-6,-5,-4},{0,809});
        const auto before=mesh::meshHash(destination);
        const auto outcome=algo::meshModelRefinedTet4(
            *model,refiningSpec,request,destination,provider);
        CHECK(outcome.status==algo::T3AStatus::Incomplete);
        CHECK(outcome.failure.reason==reason);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
    };
    auto insertion=requestFor(provider);
    insertion.refinement.maxInsertions=0;
    insertion.refinement.maxEdgeToTarget=0.5;
    requireIncomplete(insertion,algo::T3AReason::InsertionBudgetReached);
    auto work=requestFor(provider);
    work.refinement.maxInsertions=2;
    work.refinement.maxWorkUnits=1;
    work.refinement.maxEdgeToTarget=0.5;
    {
        mesh::Mesh destination;
        const auto outcome=algo::meshModelRefinedTet4(
            *model,refiningSpec,work,destination,provider);
        if (outcome.status!=algo::T3AStatus::Incomplete)
            std::printf("  hard-work status=%d reason=%d phase=%d complete=%d\n",
                static_cast<int>(outcome.status),
                static_cast<int>(outcome.failure.reason),
                static_cast<int>(outcome.failure.phase),
                outcome.reportComplete?1:0);
        CHECK(outcome.status==algo::T3AStatus::Incomplete);
        CHECK(outcome.failure.reason==algo::T3AReason::WorkBudgetReached);
        CHECK(outcome.destinationUnchanged);
        CHECK(outcome.retainedResult!=nullptr);
        CHECK(outcome.retainedResult&&
              outcome.retainedResult->metrics.workUnitsConsumed<=
                  work.refinement.maxWorkUnits);
    }

    const auto requireRejected=[&](algo::T3ARequest request) {
        mesh::Mesh destination;
        destination.addNode({-6,-5,-4},{0,810});
        const auto before=mesh::meshHash(destination);
        const auto outcome=algo::meshModelRefinedTet4(
            *model,refiningSpec,request,destination,provider);
        CHECK(outcome.status==algo::T3AStatus::Rejected);
        CHECK(outcome.failure.reason==algo::T3AReason::InvalidRefinementPolicy);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==before);
    };
    auto nodes=requestFor(provider); nodes.refinement.maxOutputNodes=1;
    requireRejected(nodes);
    auto tets=requestFor(provider); tets.refinement.maxOutputTets=1;
    requireRejected(tets);
    auto rows=requestFor(provider); rows.refinement.maxReportRows=1;
    requireRejected(rows);
    auto witnesses=requestFor(provider); witnesses.refinement.maxReportWitnesses=1;
    requireRejected(witnesses);
    auto topology=requestFor(provider); topology.refinement.maxReportTopologyHandles=1;
    requireRejected(topology);
    auto overflow=requestFor(provider);
    overflow.refinement.maxInsertions=(std::numeric_limits<std::uint64_t>::max)();
    requireRejected(overflow);

    auto reportOnly=requestFor(provider);
    reportOnly.refinement.maxEdgeToTarget=100.0;
    reportOnly.refinement.maxRadiusEdge=1.0e9;
    reportOnly.refinement.sliverRadiusRatio=0.999;
    mesh::Mesh reportedMesh;
    const auto reported=algo::meshModelRefinedTet4(
        *model,baselineSpec,reportOnly,reportedMesh,provider);
    CHECK(reported.published());
    CHECK(reported.retainedResult&&
          reported.retainedResult->quality.sliversRemaining>0);
    CHECK(reported.retainedResult&&
          reported.retainedResult->quality.sliversAttempted==0);

    auto boundedSliver=reportOnly;
    boundedSliver.refinement.sliverPolicy=algo::TetSliverPolicy::BoundedInsertion;
    boundedSliver.refinement.minimumSliverImprovement=0.01;
    boundedSliver.refinement.maxSliverCandidateAttempts=1;
    mesh::Mesh boundedDestination;
    const auto bounded=algo::meshModelRefinedTet4(
        *model,baselineSpec,boundedSliver,boundedDestination,provider);
    CHECK(bounded.status==algo::T3AStatus::Incomplete);
    CHECK(bounded.destinationUnchanged);
    CHECK(bounded.retainedResult&&
          bounded.retainedResult->quality.sliversAttempted>0);
    CHECK(bounded.retainedResult&&
          bounded.retainedResult->quality.sliversRemaining>0);
    CHECK(bounded.retainedResult&&
          !bounded.retainedResult->quality.sliverGoalComplete);
}

void threadDeterminism() {
    TEST("T3-A deterministic publication and private prefix at 1/4/auto");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    ReversedFaceModel reversed(*model);
    sizing::MeshSpec publishSpec;
    publishSpec.field=std::make_shared<sizing::ConstantField>(100.0);
    sizing::MeshSpec refineSpec;
    refineSpec.field=std::make_shared<sizing::ConstantField>(0.55);
    std::array<std::uint64_t,3> hashes{};
    std::array<std::vector<algo::T3AAttemptRow>,3> attempts;
    const std::array<int,3> threadCounts{{1,4,0}};
    for (std::size_t run=0;run<threadCounts.size();++run) {
        topo::NodePool churnNodes;
        topo::CellPool churnCells;
        for (std::size_t i=0;i<run*13;++i)
            (void)churnNodes.add({static_cast<double>(i),-3,7},
                                 {0,static_cast<int>(i)});
        CHECK(churnCells.stats().live==0);
        core::setThreadCount(threadCounts[run]);
        const geo::Model& presented=run==1
            ? static_cast<const geo::Model&>(reversed)
            : static_cast<const geo::Model&>(*model);
        Provider runProvider;
        configureProvider(presented,runProvider);
        auto runPublishSpec=publishSpec;
        runPublishSpec.threads=threadCounts[run];
        auto runRefineSpec=refineSpec;
        runRefineSpec.threads=threadCounts[run];
        mesh::Mesh published;
        const auto completed=algo::meshModelRefinedTet4(
            presented,runPublishSpec,requestFor(runProvider),published,
            runProvider);
        CHECK(completed.published());
        hashes[run]=mesh::meshHash(published);

        auto request=requestFor(runProvider);
        request.refinement.maxInsertions=2;
        request.refinement.maxEdgeToTarget=0.5;
        request.refinement.maxWorkUnits=100000;
        request.refinement.maxReportTopologyHandles=100002;
        mesh::Mesh privateDestination;
        const auto incomplete=algo::meshModelRefinedTet4(
            presented,runRefineSpec,request,privateDestination,runProvider);
        CHECK(incomplete.status==algo::T3AStatus::Incomplete);
        CHECK(incomplete.retainedResult&&
              incomplete.retainedResult->metrics.committedInsertions==1);
        if (incomplete.report)
            attempts[run].assign(incomplete.report->attempts().begin(),
                                 incomplete.report->attempts().end());
    }
    core::setThreadCount(0);
    CHECK(hashes[0]==hashes[1]&&hashes[1]==hashes[2]);
    CHECK(attempts[0].size()==attempts[1].size()&&
          attempts[1].size()==attempts[2].size());
    for (std::size_t run=1;run<attempts.size();++run)
        for (std::size_t i=0;i<attempts[0].size();++i) {
            CHECK(attempts[run][i].attemptOrdinal==
                  attempts[0][i].attemptOrdinal);
            CHECK(attempts[run][i].violation==attempts[0][i].violation);
            CHECK(attempts[run][i].candidateKind==
                  attempts[0][i].candidateKind);
            CHECK(attempts[run][i].disposition==
                  attempts[0][i].disposition);
            CHECK(attempts[run][i].reason==attempts[0][i].reason);
            CHECK(attempts[run][i].proposedCoordinate.x==
                  attempts[0][i].proposedCoordinate.x);
            CHECK(attempts[run][i].proposedCoordinate.y==
                  attempts[0][i].proposedCoordinate.y);
            CHECK(attempts[run][i].proposedCoordinate.z==
                  attempts[0][i].proposedCoordinate.z);
        }
}

struct ParsedTetExport {
    std::uint64_t nodes=0,tets=0,triangles=0,lines=0;
    bool indicesValid=true,typesValid=true,classificationsValid=true;
};

ParsedTetExport parseVtkTet4(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("VTK parse open failed");
    ParsedTetExport result;
    std::string token;
    while (input>>token) {
        if (token=="POINTS") {
            input>>result.nodes>>token;
            double value=0.0;
            for (std::uint64_t i=0;i<result.nodes*3;++i) input>>value;
        } else if (token=="CELLS") {
            std::uint64_t cells=0,listLength=0;
            input>>cells>>listLength;
            result.tets=cells;
            std::uint64_t observedLength=0;
            for (std::uint64_t i=0;i<cells;++i) {
                std::uint64_t count=0;
                input>>count;
                observedLength+=1+count;
                result.typesValid=result.typesValid&&count==4;
                for (std::uint64_t j=0;j<count;++j) {
                    std::uint64_t node=0;
                    input>>node;
                    result.indicesValid=result.indicesValid&&node<result.nodes;
                }
            }
            result.indicesValid=result.indicesValid&&observedLength==listLength;
        } else if (token=="CELL_TYPES") {
            std::uint64_t cells=0;
            input>>cells;
            result.typesValid=result.typesValid&&cells==result.tets;
            for (std::uint64_t i=0;i<cells;++i) {
                int type=0; input>>type;
                result.typesValid=result.typesValid&&type==10;
            }
        }
    }
    result.indicesValid=result.indicesValid&&result.nodes>0&&result.tets>0;
    return result;
}

ParsedTetExport parseMshTet4(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("MSH parse open failed");
    ParsedTetExport result;
    std::string line;
    while (std::getline(input,line)) {
        if (line=="$Nodes") {
            std::getline(input,line);
            std::istringstream header(line);
            std::uint64_t blocks=0,minimum=0,maximum=0;
            header>>blocks>>result.nodes>>minimum>>maximum;
            std::uint64_t observed=0;
            for (std::uint64_t block=0;block<blocks;++block) {
                std::getline(input,line);
                std::istringstream row(line);
                int dim=0,tag=0,parametric=0;
                std::uint64_t count=0;
                row>>dim>>tag>>parametric>>count;
                result.classificationsValid=result.classificationsValid&&
                    dim>=0&&dim<=3&&tag>=0&&parametric==0;
                std::vector<std::uint64_t> ids(count);
                for (auto& id:ids) {
                    std::getline(input,line);
                    std::istringstream idRow(line); idRow>>id;
                    result.indicesValid=result.indicesValid&&
                        id>=1&&id<=result.nodes;
                }
                for (std::uint64_t i=0;i<count;++i) {
                    std::getline(input,line);
                    std::istringstream point(line);
                    double x=0.0,y=0.0,z=0.0;
                    point>>x>>y>>z;
                    result.indicesValid=result.indicesValid&&
                        std::isfinite(x)&&std::isfinite(y)&&std::isfinite(z);
                }
                observed+=count;
            }
            result.indicesValid=result.indicesValid&&observed==result.nodes;
        } else if (line=="$Elements") {
            std::getline(input,line);
            std::istringstream header(line);
            std::uint64_t blocks=0,total=0,minimum=0,maximum=0;
            header>>blocks>>total>>minimum>>maximum;
            std::uint64_t observed=0;
            for (std::uint64_t block=0;block<blocks;++block) {
                std::getline(input,line);
                std::istringstream blockRow(line);
                int dim=0,tag=0,type=0;
                std::uint64_t count=0;
                blockRow>>dim>>tag>>type>>count;
                const int nodeCount=type==1?2:(type==2?3:(type==4?4:0));
                result.typesValid=result.typesValid&&nodeCount!=0;
                result.classificationsValid=result.classificationsValid&&
                    dim==(type==1?1:(type==2?2:3))&&tag>=0;
                if (type==1) result.lines+=count;
                if (type==2) result.triangles+=count;
                if (type==4) result.tets+=count;
                for (std::uint64_t i=0;i<count;++i) {
                    std::getline(input,line);
                    std::istringstream element(line);
                    std::uint64_t id=0; element>>id;
                    for (int j=0;j<nodeCount;++j) {
                        std::uint64_t node=0; element>>node;
                        result.indicesValid=result.indicesValid&&
                            node>=1&&node<=result.nodes;
                    }
                }
                observed+=count;
            }
            result.indicesValid=result.indicesValid&&observed==total;
        }
    }
    result.indicesValid=result.indicesValid&&result.nodes>0&&result.tets>0;
    return result;
}

void auditedTet4Export() {
    TEST("T3-A audited Tet4 publication exports to VTK and MSH");
    auto kernel=geo::makeBuiltinKernel();
    auto model=kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(0.55);
    auto request=requestFor(provider);
    request.refinement.maxEdgeToTarget=1.8;
    mesh::Mesh published;
    const auto outcome=algo::meshModelRefinedTet4(
        *model,spec,request,published,provider);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult&&outcome.retainedResult->quality.inverted==0);
    CHECK(outcome.retainedResult&&
          outcome.retainedResult->metrics.committedInsertions==1);
    if (outcome.retainedResult) {
        const auto& quality=outcome.retainedResult->quality;
        CHECK(quality.exactPositive==quality.tetCount);
        CHECK(quality.nonfinite==0&&quality.inverted==0);
        CHECK(quality.sizeViolations==0);
        CHECK(quality.minimumSignedVolume>0.0);
        CHECK(quality.minimumRadiusRatio>0.0);
        CHECK(quality.p1RadiusRatio>=quality.minimumRadiusRatio);
        CHECK(quality.p5RadiusRatio>=quality.p1RadiusRatio);
        CHECK(quality.p99RadiusEdge>=quality.p95RadiusEdge);
        std::uint64_t edgeHistogram=0,ratioHistogram=0;
        for (const auto count:quality.radiusEdgeHistogram) edgeHistogram+=count;
        for (const auto count:quality.radiusRatioHistogram) ratioHistogram+=count;
        CHECK(edgeHistogram==quality.tetCount);
        CHECK(ratioHistogram==quality.tetCount);
    }
    const auto base="paacs_t3a_"+std::to_string(::getpid());
    const auto vtk=(std::filesystem::temp_directory_path()/(base+".vtk")).string();
    const auto msh=(std::filesystem::temp_directory_path()/(base+".msh")).string();
    io::writeVtk(published,vtk,true);
    io::writeMsh(published,msh);
    const auto vtkParsed=parseVtkTet4(vtk);
    const auto mshParsed=parseMshTet4(msh);
    CHECK(vtkParsed.nodes==static_cast<std::uint64_t>(published.nodeCount()));
    CHECK(vtkParsed.tets==outcome.retainedResult->quality.tetCount);
    CHECK(vtkParsed.indicesValid&&vtkParsed.typesValid);
    CHECK(mshParsed.nodes==static_cast<std::uint64_t>(published.nodeCount()));
    CHECK(mshParsed.tets==outcome.retainedResult->quality.tetCount);
    CHECK(mshParsed.triangles==outcome.retainedResult->immutableBoundary.size());
    CHECK(mshParsed.lines+mshParsed.triangles+mshParsed.tets==
          static_cast<std::uint64_t>(published.cellCount()));
    CHECK(mshParsed.indicesValid);
    CHECK(mshParsed.typesValid);
    CHECK(mshParsed.classificationsValid);
    const auto semantic=mesh::semanticHashEx(published);
    CHECK(semantic.canonical);
    CHECK(published.nodeCount()==27);
    CHECK(published.cellCount()==123);
    CHECK(outcome.retainedResult->quality.tetCount==51);
    CHECK(outcome.retainedResult->immutableBoundary.size()==48);
    CHECK(mshParsed.lines==24);
    CHECK(semantic.value==0x84f3e6e847108a8full);
#if defined(_MSC_VER)
    CHECK(mesh::meshHash(published)==0xb3753fee82bad073ull);
#endif
    std::printf("  t3a reference nodes=%lld cells=%lld tets=%llu tris=%llu "
                "strict=%016llx semantic=%016llx\n",
        static_cast<long long>(published.nodeCount()),
        static_cast<long long>(published.cellCount()),
        static_cast<unsigned long long>(outcome.retainedResult->quality.tetCount),
        static_cast<unsigned long long>(outcome.retainedResult->immutableBoundary.size()),
        static_cast<unsigned long long>(mesh::meshHash(published)),
        static_cast<unsigned long long>(semantic.value));
    std::filesystem::remove(vtk);
    std::filesystem::remove(msh);
}

} // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* t3aTests(std::size_t& count) {
    static const Entry tests[] = {
        {"t3a.predicates.power", checkedPowerPredicates},
        {"t3a.protected.certified_aabb", certifiedProtectedBounds},
        {"t3a.violation.identity_volume", violationIdentityAndVolumeFloor},
        {"t3a.transaction.staged", stagedTetRefinementTransaction},
        {"t3a.candidate.indexed_contracts", indexedLocationAndCandidateContracts},
        {"t3a.route.bounded", boundedPublicRoute},
        {"t3a.route.shapes_local_field", shapeAndLocalFieldCoverage},
        {"t3a.route.multi_insertion", multiInsertionPublicRoute},
        {"t3a.metrics.multi_insertion_scaling", multiInsertionScaling},
        {"t3a.metrics.public_fixture_scaling", publicFixtureScaling},
        {"t3a.route.sealed_cavity", sealedCavityBoundary},
        {"t3a.index.dynamic_scaling", dynamicLocationScaling},
        {"t3a.metrics.scaling", operationScaling},
        {"t3a.route.sizing_budgets", sizingTaxonomyAndBudgets},
        {"t3a.route.determinism", threadDeterminism},
        {"t3a.export.tet4", auditedTet4Export},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
} // namespace paacs_tests
