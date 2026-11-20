// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T3-A standalone allocation qualification. This target replaces global
// operator new and is intentionally excluded from sanitizer builds.
#include <algorithm>
#include <array>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <new>
#include <span>
#include <string_view>
#include <vector>

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#include <intrin.h>
#else
#include <unistd.h>
#endif

#include "paacs/algo/TetRefinementMesher.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "workspace_digest.hpp"

static long gAllocationOrdinal=0;
static long gFailOrdinal=0;
static long gPhase3Allocations=0;
static long gPublicationTailAllocations=0;
static bool gTrackAllocations=false;
static bool gFailureInjected=false;
static bool gProgressDiagnostics=false;
static std::uint8_t gFailurePhase=0;

namespace {
void progressDiagnostic(long ordinal) noexcept {
    char buffer[96];
    char* cursor=buffer;
    constexpr char prefix[]="probe_t3a begin ordinal=";
    std::memcpy(cursor,prefix,sizeof(prefix)-1);
    cursor+=sizeof(prefix)-1;
    const auto converted=std::to_chars(cursor,buffer+sizeof(buffer),ordinal);
    if (converted.ec==std::errc{}) cursor=converted.ptr;
    *cursor++='\n';
#if defined(_WIN32)
    DWORD written=0;
    (void)::WriteFile(::GetStdHandle(STD_ERROR_HANDLE),buffer,
                      static_cast<DWORD>(cursor-buffer),&written,nullptr);
#else
    (void)::write(STDERR_FILENO,buffer,
                  static_cast<std::size_t>(cursor-buffer));
#endif
}

void terminateDiagnostic() noexcept {
    char buffer[256];
    char* cursor=buffer;
    auto literal=[&](const char* value) noexcept {
        const std::size_t count=std::strlen(value);
        std::memcpy(cursor,value,count);
        cursor+=count;
    };
    auto number=[&](long value) noexcept {
        const auto converted=std::to_chars(cursor,buffer+sizeof(buffer),value);
        if (converted.ec==std::errc{}) cursor=converted.ptr;
    };
    literal("probe_t3a terminate ordinal="); number(gAllocationOrdinal);
    literal(" fail="); number(gFailOrdinal);
    literal(" phase="); number(gFailurePhase);
    literal(" phase3="); number(gPhase3Allocations);
    literal(" publication="); number(gPublicationTailAllocations);
    literal("\n");
#if defined(_WIN32)
    DWORD written=0;
    (void)::WriteFile(::GetStdHandle(STD_ERROR_HANDLE),buffer,
                      static_cast<DWORD>(cursor-buffer),&written,nullptr);
    __fastfail(7);
#else
    (void)::write(STDERR_FILENO,buffer,
                  static_cast<std::size_t>(cursor-buffer));
    std::_Exit(134);
#endif
}
} // namespace

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    if (paacs::topo::detail::phase3Active()) ++gPhase3Allocations;
    if (paacs::algo::boundary_conforming_detail::publicationTailActive())
        ++gPublicationTailAllocations;
    if (gTrackAllocations) {
        ++gAllocationOrdinal;
        if (gFailOrdinal>0&&gAllocationOrdinal==gFailOrdinal) {
            gFailureInjected=true;
            gFailurePhase=paacs::algo::boundary_conforming_detail::
                testingAllocationPhaseMarker();
            throw std::bad_alloc{};
        }
    }
    if (void* memory=std::malloc(bytes?bytes:1)) return memory;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory,std::size_t) noexcept { std::free(memory); }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;
int failures=0;
#define REQUIRE(condition)                                                   \
    do {                                                                     \
        if (!(condition)) {                                                  \
            std::printf("FAIL %d: %s\n",__LINE__,#condition);              \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

struct TransactionFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::array<topo::NodeHandle,4> nodes{};
    topo::CellHandle cell{};
};

TransactionFixture buildTransactionFixture() {
    TransactionFixture fixture;
    topo::NodePool nodes;
    fixture.nodes[0]=nodes.add({0,0,0},{3,9});
    fixture.nodes[1]=nodes.add({1,0,0},{3,9});
    fixture.nodes[2]=nodes.add({0,1,0},{3,9});
    fixture.nodes[3]=nodes.add({0,0,1},{3,9});
    topo::CellPool cells;
    fixture.cell=cells.add(topo::ElementType::Tet4,fixture.nodes,{3,9},nodes);
    topo::EdgePool edges;
    (void)edges.insertCellEdges(nodes,cells,fixture.cell);
    topo::FacePool faces;
    (void)faces.insertCellFaces(nodes,cells,edges,fixture.cell);
    const std::array<topo::NodeHandle,3> protectedNodes{
        fixture.nodes[0],fixture.nodes[2],fixture.nodes[1]};
    const auto face=faces.find(nodes,protectedNodes);
    if (!face) throw std::logic_error("probe protected face missing");
    faces.setClassification(*face,{2,77});
    fixture.workspace=std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    return fixture;
}

topo::TetRefinementCavitySpec transactionSpec(
    const TransactionFixture& fixture) {
    topo::TetRefinementCavitySpec spec;
    spec.stamp=topo::WorkspaceStamp::of(*fixture.workspace);
    spec.stagedNode={0,{0.2,0.2,0.2},{3,9}};
    spec.removed.push_back(fixture.cell);
    static constexpr std::array<std::array<std::size_t,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    for (const auto& face:faces) {
        topo::TetRefinementCellSpec cell;
        for (std::size_t i=0;i<3;++i)
            cell.vertices[i]=topo::TetRefinementNodeRef::existingNode(
                fixture.nodes[face[i]]);
        cell.vertices[3]=topo::TetRefinementNodeRef::stagedNode();
        cell.classification={3,9};
        spec.added.push_back(cell);
    }
    return spec;
}

void retryTransaction(TransactionFixture& fixture) {
    const auto spec=transactionSpec(fixture);
    auto tx=fixture.workspace->beginTetRefinementTransaction();
    (void)tx.stage(spec);
    const auto& predicted=tx.prepare();
    const auto result=tx.commit();
    REQUIRE(result.actualNode==predicted.stagedNode);
    REQUIRE(result.actualNode==result.predictedNode);
    for (const auto& row:result.cells) REQUIRE(row.actual==row.predicted);
    for (const auto& row:result.edges) REQUIRE(row.actual==row.predicted);
    for (const auto& row:result.faces) REQUIRE(row.actual==row.predicted);
}

void transactionSweeps() {
    {
        auto fixture=buildTransactionFixture();
        gPhase3Allocations=0;
        retryTransaction(fixture);
        REQUIRE(gPhase3Allocations==0);
    }

    long stageFaults=0;
    for (long ordinal=1;;++ordinal) {
        auto fixture=buildTransactionFixture();
        const auto before=paacs_tests::observableDigest(*fixture.workspace);
        bool faulted=false;
        const auto spec=transactionSpec(fixture);
        {
            try {
                auto tx=fixture.workspace->beginTetRefinementTransaction();
                gAllocationOrdinal=0; gFailOrdinal=ordinal;
                gTrackAllocations=true;
                (void)tx.stage(spec);
                gTrackAllocations=false; gFailOrdinal=0;
                tx.abandon();
            } catch (const std::bad_alloc&) {
                gTrackAllocations=false; gFailOrdinal=0; faulted=true;
            }
        }
        REQUIRE(paacs_tests::observableDigest(*fixture.workspace)==before);
        retryTransaction(fixture);
        if (!faulted) break;
        ++stageFaults;
    }

    long prepareFaults=0;
    for (long ordinal=1;;++ordinal) {
        auto fixture=buildTransactionFixture();
        const auto before=
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        bool faulted=false;
        const auto spec=transactionSpec(fixture);
        {
            try {
                auto tx=fixture.workspace->beginTetRefinementTransaction();
                (void)tx.stage(spec);
                gAllocationOrdinal=0; gFailOrdinal=ordinal;
                gTrackAllocations=true;
                (void)tx.prepare();
                gTrackAllocations=false; gFailOrdinal=0;
                tx.abandon();
            } catch (const std::bad_alloc&) {
                gTrackAllocations=false; gFailOrdinal=0; faulted=true;
            }
        }
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace)==before);
        retryTransaction(fixture);
        if (!faulted) break;
        ++prepareFaults;
    }

    long placeholderFaults=0;
    for (int ordinal=1;;++ordinal) {
        auto fixture=buildTransactionFixture();
        const auto before=
            paacs_tests::observableLogicalDigest(*fixture.workspace);
        bool faulted=false;
        const auto spec=transactionSpec(fixture);
        {
            try {
                auto tx=fixture.workspace->beginTetRefinementTransaction();
                (void)tx.stage(spec);
                topo::faultCountdown()=ordinal;
                (void)tx.prepare();
                topo::faultCountdown()=0;
                tx.abandon();
            } catch (const std::bad_alloc&) {
                topo::faultCountdown()=0; faulted=true;
            }
        }
        REQUIRE(paacs_tests::observableLogicalDigest(*fixture.workspace)==before);
        retryTransaction(fixture);
        if (!faulted) break;
        ++placeholderFaults;
    }
    std::printf("t3a transaction: stage=%ld prepare=%ld placeholders=%ld "
                "phase3=%ld\n",stageFaults,prepareFaults,placeholderFaults,
                gPhase3Allocations);
    REQUIRE(stageFaults>0);
    REQUIRE(prepareFaults>0);
    REQUIRE(placeholderFaults>0);
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
        result.providerStamp=current; result.facePersistentId=face;
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
    Vec3 lo{},hi{};
    const auto solids=model.solids();
    if (solids.size()!=1) throw std::logic_error("probe solid missing");
    solids.front()->bbox(lo,hi);
    const Vec3 centre=(lo+hi)*0.5;
    for (const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{}; face->paramBounds(uvLo,uvHi);
        const Vec2 uv=(uvLo+uvHi)*0.5;
        const int sense=face->normalAt(uv).dot(
            face->pointAt(uv)-centre)>0.0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]=
            {{0,sense}};
    }
}

algo::T3ARequest routeRequest(const Provider& provider,bool refine) {
    algo::T3ARequest request;
    request.conformingRequest.expectedSurfaceSetStamp={9801,1};
    request.conformingRequest.expectedOwnershipProviderStamp=provider.current;
    request.conformingRequest.envelopeMode=
        algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    request.refinement.maxEdgeToTarget=refine?0.5:2.0;
    request.refinement.maxRadiusEdge=100.0;
    request.refinement.maxInsertions=refine?2:8;
    request.refinement.maxWorkUnits=refine?100000:10000;
    request.refinement.maxOutputNodes=1000;
    request.refinement.maxOutputTets=10000;
    request.refinement.maxReportRows=45;
    request.refinement.maxReportWitnesses=10000;
    request.refinement.maxReportTopologyHandles=refine?100002:10008;
    return request;
}

algo::T3AOutcome runRoute(const geo::Model& model,Provider& provider,
                          bool refine,mesh::Mesh& destination) {
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(refine?0.55:100.0);
    return algo::meshModelRefinedTet4(
        model,spec,routeRequest(provider,refine),destination,provider);
}

long measureRoute(const geo::Model& model,Provider& provider,bool refine,
                  algo::T3AOutcome& outcome,mesh::Mesh& destination,
                  long failOrdinal) {
    sizing::MeshSpec spec;
    spec.field=std::make_shared<sizing::ConstantField>(refine?0.55:100.0);
    const auto request=routeRequest(provider,refine);
    gAllocationOrdinal=0; gFailOrdinal=failOrdinal;
    gFailureInjected=false; gTrackAllocations=true;
    outcome=algo::meshModelRefinedTet4(
        model,spec,request,destination,provider);
    gTrackAllocations=false; gFailOrdinal=0;
    return gAllocationOrdinal;
}

void normalRouteProof(const geo::Model& model,Provider& provider) {
    mesh::Mesh destination;
    gPublicationTailAllocations=0;
    auto outcome=runRoute(model,provider,false,destination);
    REQUIRE(outcome.published());
    REQUIRE(gPublicationTailAllocations==0);
}

struct CleanRetrySignature {
    algo::T3AStatus status{};
    algo::T3AReason reason{};
    algo::T3APhase phase{};
    bool destinationUnchanged=false,reportComplete=false,
         reportEnvelopePrepared=false,reportRowsComplete=false,
         hasTerminalViolation=false,boundaryPassed=false;
    std::uint64_t attempts=0,witnesses=0,blocked=0,unattempted=0,
        commits=0,topologyRows=0,committedInsertions=0,
        blockedRecorded=0,blockedInvalidated=0,blockedCurrent=0,
        work=0,globalScans=0,contextRebuilds=0,phase3Allocations=0,
        publicationTailAllocations=0,semanticTopology=0,revision=0,
        liveNodes=0,liveCells=0,prescribedExterior=0,matchedExterior=0,
        orphanExterior=0,unexpectedExterior=0;
    friend bool operator==(const CleanRetrySignature&,
                           const CleanRetrySignature&)=default;
};

CleanRetrySignature cleanRetrySignature(const algo::T3AOutcome& outcome) {
    CleanRetrySignature signature;
    signature.status=outcome.status;
    signature.reason=outcome.failure.reason;
    signature.phase=outcome.failure.phase;
    signature.destinationUnchanged=outcome.destinationUnchanged;
    signature.reportComplete=outcome.reportComplete;
    if (outcome.report) {
        signature.reportEnvelopePrepared=outcome.report->envelopePrepared;
        signature.reportRowsComplete=outcome.report->reportComplete;
        signature.hasTerminalViolation=outcome.report->hasTerminalViolation;
        signature.attempts=outcome.report->attemptsUsed;
        signature.witnesses=outcome.report->protectedWitnessesUsed;
        signature.blocked=outcome.report->blockedUsed;
        signature.unattempted=outcome.report->unattemptedUsed;
        signature.commits=outcome.report->committedPrefixUsed;
        signature.topologyRows=outcome.report->topologyHandlesUsed;
    }
    if (outcome.retainedResult) {
        const auto& metrics=outcome.retainedResult->metrics;
        signature.committedInsertions=metrics.committedInsertions;
        signature.blockedRecorded=metrics.blockedRowsRecorded;
        signature.blockedInvalidated=metrics.blockedRowsInvalidated;
        signature.blockedCurrent=metrics.blockedRowsCurrent;
        signature.work=metrics.workUnitsConsumed;
        signature.globalScans=metrics.perInsertionGlobalTetScans;
        signature.contextRebuilds=metrics.contextRebuilds;
        signature.phase3Allocations=metrics.phase3Allocations;
        signature.publicationTailAllocations=metrics.publicationTailAllocations;
        const auto& audit=outcome.retainedResult->boundaryAudit;
        signature.boundaryPassed=audit.passed;
        signature.prescribedExterior=audit.prescribedExterior;
        signature.matchedExterior=audit.matchedExterior;
        signature.orphanExterior=audit.orphanExterior;
        signature.unexpectedExterior=audit.unexpectedExterior;
        if (outcome.retainedResult->authoritative&&
            outcome.retainedResult->authoritative->workspace) {
            const auto& workspace=
                *outcome.retainedResult->authoritative->workspace;
            signature.semanticTopology=workspace.semanticTopologyHash();
            signature.revision=workspace.topologyRevision();
            signature.liveNodes=workspace.nodes().stats().live;
            signature.liveCells=workspace.cells().stats().live;
        }
    }
    return signature;
}

void routeRange(const geo::Model& model,Provider& provider,long first,long last) {
    mesh::Mesh warmDestination;
    const auto warm=runRoute(model,provider,true,warmDestination);
    const auto referenceSignature=cleanRetrySignature(warm);
    REQUIRE(warm.status==algo::T3AStatus::Incomplete);
    REQUIRE(warm.destinationUnchanged);
    REQUIRE(warm.report&&warm.report->reportComplete);
    REQUIRE(warm.retainedResult&&referenceSignature.boundaryPassed);
    REQUIRE(referenceSignature.orphanExterior==0);
    REQUIRE(referenceSignature.unexpectedExterior==0);
    REQUIRE(referenceSignature.globalScans==0);
    REQUIRE(referenceSignature.contextRebuilds==0);
    REQUIRE(referenceSignature.phase3Allocations==0);
    REQUIRE(referenceSignature.publicationTailAllocations==0);
    std::printf("t3a clean reference: status=%u reason=%u phase=%u "
                "commits=%llu blocked=%llu work=%llu semantic=%016llx\n",
        static_cast<unsigned>(referenceSignature.status),
        static_cast<unsigned>(referenceSignature.reason),
        static_cast<unsigned>(referenceSignature.phase),
        static_cast<unsigned long long>(referenceSignature.committedInsertions),
        static_cast<unsigned long long>(referenceSignature.blocked),
        static_cast<unsigned long long>(referenceSignature.work),
        static_cast<unsigned long long>(referenceSignature.semanticTopology));
    long faults=0;
    for (long ordinal=first;ordinal<=last;++ordinal) {
        if (gProgressDiagnostics) progressDiagnostic(ordinal);
        mesh::Mesh destination;
        destination.addNode({-9,-8,-7},{0,901});
        const auto before=mesh::meshHash(destination);
        algo::T3AOutcome outcome;
        (void)measureRoute(model,provider,true,outcome,destination,ordinal);
        REQUIRE(gFailureInjected);
        REQUIRE(outcome.status!=algo::T3AStatus::Published);
        REQUIRE(outcome.destinationUnchanged);
        REQUIRE(mesh::meshHash(destination)==before);
        if (outcome.report&&outcome.report->envelopePrepared) {
            REQUIRE(outcome.reportComplete);
            REQUIRE(outcome.report->reportComplete);
        }
        if (outcome.report&&outcome.report->attemptsUsed>0) {
            const auto& row=outcome.report->attemptStorage[
                static_cast<std::size_t>(outcome.report->attemptsUsed-1)];
            if (outcome.status==algo::T3AStatus::Failed&&
                outcome.failure.phase>=algo::T3APhase::CandidateGeneration&&
                outcome.failure.phase<=algo::T3APhase::CoordinatorPreparation)
                REQUIRE(row.disposition==algo::T3AAttemptDisposition::TerminalFailure||
                        row.disposition==algo::T3AAttemptDisposition::Committed);
        }
        mesh::Mesh retry;
        auto clean=runRoute(model,provider,true,retry);
        const auto cleanSignature=cleanRetrySignature(clean);
        REQUIRE(clean.status==algo::T3AStatus::Incomplete);
        REQUIRE(clean.destinationUnchanged);
        REQUIRE(clean.report&&clean.report->reportComplete);
        REQUIRE(clean.retainedResult&&cleanSignature.boundaryPassed);
        REQUIRE(cleanSignature.orphanExterior==0);
        REQUIRE(cleanSignature.unexpectedExterior==0);
        REQUIRE(cleanSignature.globalScans==0);
        REQUIRE(cleanSignature.contextRebuilds==0);
        REQUIRE(cleanSignature.phase3Allocations==0);
        REQUIRE(cleanSignature.publicationTailAllocations==0);
        REQUIRE(cleanSignature==referenceSignature);
        ++faults;
    }
    std::printf("t3a route range: %ld..%ld faults=%ld phase3=%ld tail=%ld\n",
                first,last,faults,gPhase3Allocations,
                gPublicationTailAllocations);
    REQUIRE(gPhase3Allocations==0);
    REQUIRE(gPublicationTailAllocations==0);
}

long parseLong(const char* text) {
    char* end=nullptr;
    const long value=std::strtol(text,&end,10);
    if (!end||*end!='\0'||value<1) throw std::invalid_argument("bad ordinal");
    return value;
}

} // namespace

int main(int argc,char** argv) {
    std::setvbuf(stdout,nullptr,_IONBF,0);
    std::set_terminate(&terminateDiagnostic);
    try {
        auto kernel=paacs::geo::makeBuiltinKernel();
        auto model=kernel->makeBox({0,0,0},{1,1,1});
        Provider provider; configureProvider(*model,provider);
        if (argc==1) {
            transactionSweeps();
            normalRouteProof(*model,provider);
        } else if (argc==2&&std::string_view(argv[1])=="--count") {
            mesh::Mesh warmDestination;
            (void)runRoute(*model,provider,true,warmDestination);
            paacs::mesh::Mesh destination;
            paacs::algo::T3AOutcome outcome;
            const long count=measureRoute(
                *model,provider,true,outcome,destination,0);
            std::printf("t3a route allocation count: %ld\n",count);
            REQUIRE(outcome.status==paacs::algo::T3AStatus::Incomplete);
        } else if (argc==3&&std::string_view(argv[1])=="--expect-route-count") {
            mesh::Mesh warmDestination;
            (void)runRoute(*model,provider,true,warmDestination);
            paacs::mesh::Mesh destination;
            paacs::algo::T3AOutcome outcome;
            const long count=measureRoute(
                *model,provider,true,outcome,destination,0);
            std::printf("t3a route allocation count: %ld\n",count);
            REQUIRE(count==parseLong(argv[2]));
            REQUIRE(outcome.status==paacs::algo::T3AStatus::Incomplete);
        } else if (argc==4&&
                  (std::string_view(argv[1])=="--range"||
                   std::string_view(argv[1])=="--range-progress")) {
            gProgressDiagnostics=std::string_view(argv[1])=="--range-progress";
            routeRange(*model,provider,parseLong(argv[2]),parseLong(argv[3]));
        } else {
            std::printf("usage: paacs_probe_t3a [--count|--expect-route-count N|"
                        "--range FIRST LAST|--range-progress FIRST LAST]\n");
            return 2;
        }
    } catch (const std::exception& error) {
        std::printf("UNEXPECTED: %s\n",error.what());
        return 3;
    } catch (...) {
        std::printf("UNEXPECTED NON-STANDARD EXCEPTION\n");
        return 3;
    }
    if (failures) {
        std::printf("T3A PROBE FAILED (%d)\n",failures);
        return 1;
    }
    std::printf("T3A PROBE PASSED\n");
    return 0;
}
