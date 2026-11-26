// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "../bnd/TetNumericalTopology.hpp"
#include "paacs/mesh/Quality.hpp"
#include <atomic>
#include <cmath>
#include <optional>
#include <queue>
#include <ranges>

namespace paacs::geo::detail {

struct OcctInitialTopologyQualityLimits {
    std::uint64_t maximumWork=1'000'000,maximumWorkingBytes=64ULL*1024*1024;
    std::size_t maximumEdits=64,maximumCandidates=128,maximumProposals=4096;
    // A construction trigger, not a publication quality requirement. The
    // deliberately small default targets numerical slivers rather than thin
    // cells whose anisotropy may represent the intended geometry.
    double searchBelowMeanRatio=1e-6;
    const std::atomic_bool* cancellation=nullptr;
};
struct OcctInitialTopologyQualityReport {
    bool attempted=false,completed=false,changed=false,cancelled=false,budgetLimited=false;
    bool summaryAvailable=false,candidateListTruncated=false;
    const char* diagnostic="not_attempted";
    std::uint64_t work=0,maximumModeledBytes=0,candidates=0,proposals=0,certifications=0;
    std::uint64_t stagedEdits=0,edits=0,flips23=0,flips32=0,flips44=0;
    std::uint64_t liveBefore=0,liveAfter=0,severeBefore=0,severeAfter=0;
    std::uint64_t protectedSourceFaces=0,protectedSourceEdges=0;
    double minimumMeanRatioBefore=1,minimumMeanRatioAfter=1;
    double minimumScaledJacobianBefore=1,minimumScaledJacobianAfter=1;
};

namespace occt_initial_topology_quality_detail {
using Complex=bnd::detail::BistellarComplex;
using Cell=bnd::TetFacetCell;
using Carrier=bnd::detail::TetNumericalTopology;
struct Limit {const char* reason;};
struct Finished {const char* reason;};
struct Cancelled {};
struct Invalid {};
struct Shape {double mean=0,jacobian=0;};

inline std::uint64_t add(std::uint64_t a,std::uint64_t b){
    if(b>UINT64_MAX-a)throw Limit{"memory_budget"};return a+b;
}
inline std::uint64_t multiply(std::uint64_t a,std::uint64_t b){
    if(b&&a>UINT64_MAX/b)throw Limit{"memory_budget"};return a*b;
}
inline bool guard(const Cell& cell,std::size_t original){
    return std::ranges::any_of(cell.vertices,[&](auto v){return v>=original&&v-original<8;});
}
template<class Charge>
Shape shape(const Complex& complex,const Cell& cell,Charge&& charge){
    charge(12);
    std::array<Vec3,4> p;
    for(unsigned i=0;i<4;++i){
        if(cell.vertices[i]>=complex.vertices.size())throw Invalid{};
        p[i]=complex.vertices[cell.vertices[i]].position;
        if(!std::isfinite(p[i].x)||!std::isfinite(p[i].y)||!std::isfinite(p[i].z))throw Invalid{};
    }
    // Sign is decided in the original coordinate system. Scaling is used only
    // for shape arithmetic; no coordinates, orientation decisions or CAD
    // constraints are snapped or perturbed.
    if(predicates::orient3dFiltered(p[0],p[1],p[2],p[3])!=predicates::Sign::Positive)throw Invalid{};
    std::array<Vec3,4> q{};double scale=0;
    for(unsigned i=1;i<4;++i){q[i]=p[i]-p[0];
        scale=std::max({scale,std::abs(q[i].x),std::abs(q[i].y),std::abs(q[i].z)});}
    if(!(scale>0)||!std::isfinite(scale))throw Invalid{};
    for(unsigned i=1;i<4;++i)q[i]=q[i]/scale;
    const auto det=predicates::orient3dDeterminantEstimate(q[0],q[1],q[2],q[3]);
    double sum=0;for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j)sum+=(q[i]-q[j]).norm2();
    const auto mean=det>0&&sum>0?12*std::pow(det/2.,2./3.)/sum:0.;
    const auto jacobian=mesh::tetScaledJacobian(q[0],q[1],q[2],q[3]);
    if(!std::isfinite(mean)||!std::isfinite(jacobian))throw Invalid{};
    return {std::max(0.,mean),std::max(0.,jacobian)};
}
} // namespace occt_initial_topology_quality_detail

// Optional finalization of a certified private numerical complex. Call only
// AFTER its last Delaunay point insertion: shape-driven bistellar flips do not
// establish or preserve the unweighted Delaunay property. Original coordinates,
// source facets/edges and body ownership are fixed. Existing recovery, material
// selection and final Native validity/CAD checks remain mandatory.
//
// This conservative prototype clones only when a severe candidate exists and
// admits both owners plus source/proposal scratch. Every refusal leaves the
// caller's carrier unchanged, including allocation or cancellation during an
// edit. A later optimization can stage individual patches without changing
// this contract. maximumWorkingBytes includes the retained input owner.
inline OcctInitialTopologyQualityReport finalizeOcctInitialTopologyQuality(
    bnd::detail::TetNumericalTopology& topology,const OcctInitialTopologyQualityLimits& limits={}) {
    using namespace occt_initial_topology_quality_detail;
    OcctInitialTopologyQualityReport report;
    if(!topology.available)return report;
    report.attempted=true;
    std::uint64_t workLimit=limits.maximumWork;
    const auto cancelled=[&]{if(limits.cancellation&&limits.cancellation->load(std::memory_order_relaxed))throw Cancelled{};};
    const auto charge=[&](std::uint64_t n=1){
        cancelled();if(report.work>workLimit||n>workLimit-report.work)throw Limit{"work_budget"};report.work+=n;
    };
    try {
        charge(0);
        if(!limits.maximumEdits||!limits.maximumCandidates||!limits.maximumProposals){
            report.completed=true;report.diagnostic="disabled";return report;
        }
        if(!(limits.searchBelowMeanRatio>0&&limits.searchBelowMeanRatio<1)||
           !std::isfinite(limits.searchBelowMeanRatio)||topology.complex.chargeWork||
           !topology.complex.walkConstraints||topology.sourceNodes>topology.complex.vertices.size()||
           topology.complex.vertices.size()-topology.sourceNodes<8||
           topology.complex.cells.size()>topology.complex.maximumStoredCells||
           topology.complex.cells.size()!=topology.complex.alive.size()||
           topology.complex.cells.size()!=topology.complex.neighbours.size())throw Invalid{};
        const auto retained=topology.retainedModeledBytes();
        const auto queueCount=std::min(limits.maximumCandidates,topology.complex.cells.size());
        const auto scanBytes=add(retained,add(65536,multiply(queueCount,64)));
        if(scanBytes>limits.maximumWorkingBytes)throw Limit{"memory_budget"};
        report.maximumModeledBytes=scanBytes;
        using Candidate=std::pair<double,std::size_t>;
        std::vector<Candidate> storage;storage.reserve(queueCount);
        std::priority_queue<Candidate,std::vector<Candidate>> worst(std::less<Candidate>{},std::move(storage));
        bool measured=false;
        for(std::size_t id=0;id<topology.complex.cells.size();++id){
            charge();if(!topology.complex.alive[id])continue;
            ++report.liveBefore;
            const auto& cell=topology.complex.cells[id];
            if(guard(cell,topology.sourceNodes))continue;
            const auto q=shape(topology.complex,cell,charge);measured=true;
            report.minimumMeanRatioBefore=std::min(report.minimumMeanRatioBefore,q.mean);
            report.minimumScaledJacobianBefore=std::min(report.minimumScaledJacobianBefore,q.jacobian);
            if(q.mean>=limits.searchBelowMeanRatio)continue;
            ++report.severeBefore;
            const Candidate candidate{q.mean,id};
            if(worst.size()<queueCount)worst.push(candidate);
            else if(candidate<worst.top()){worst.pop();worst.push(candidate);}
        }
        report.summaryAvailable=measured;
        report.liveAfter=report.liveBefore;report.severeAfter=report.severeBefore;
        report.minimumMeanRatioAfter=report.minimumMeanRatioBefore;
        report.minimumScaledJacobianAfter=report.minimumScaledJacobianBefore;
        if(report.liveBefore!=topology.liveCells)throw Invalid{};
        if(worst.empty()){report.completed=true;report.diagnostic="no_severe_candidates";return report;}
        report.candidateListTruncated=report.severeBefore>worst.size();
        if(report.candidateListTruncated){report.budgetLimited=true;report.diagnostic="candidate_budget";}
        // The source envelope covers temporary protection set nodes, bookkeeping
        // and proposals. Both copied owners retain their full admitted capacity.
        const auto scratch=add(1024ULL*1024,add(multiply(topology.sourceFaces.size(),2048),multiply(queueCount,128)));
        const auto peak=add(add(retained,retained),scratch);
        if(peak>limits.maximumWorkingBytes)throw Limit{"memory_budget"};
        report.maximumModeledBytes=peak;
        const auto retirement=add(multiply(topology.complex.cells.size(),11),
            add(topology.complex.vertices.size(),topology.sourceFaces.size()));
        if(retirement>(workLimit-report.work)/2)throw Limit{"work_budget"};
        // Reserve and eventually account destruction work even when a staged
        // edit is discarded. Cancellation does not erase completed work.
        workLimit-=retirement;
        struct Retirement {
            OcctInitialTopologyQualityReport& report;std::uint64_t work;
            ~Retirement(){report.work+=work;}
        } retire{report,retirement};
        charge(retirement);
        Complex staged=topology.complex;
        staged.chargeWork=charge;
        std::vector<Complex::Face> temporaryFaces;
        std::vector<Complex::Edge> temporaryEdges;
        temporaryFaces.reserve(topology.sourceFaces.size());temporaryEdges.reserve(3*topology.sourceFaces.size());
        for(const auto& [key,directed]:topology.sourceFaces){
            charge(4);(void)directed;
            if(staged.faces.contains(key)&&staged.protectedFaces.insert(key).second)temporaryFaces.push_back(key);
            for(unsigned j=0;j<3;++j){
                const auto edge=Complex::key(Complex::Edge{key[j],key[(j+1)%3]});
                if(staged.edges.contains(edge)&&staged.protectedEdges.insert(edge).second)temporaryEdges.push_back(edge);
            }
        }
        report.protectedSourceFaces=temporaryFaces.size();report.protectedSourceEdges=temporaryEdges.size();
        std::vector<Candidate> candidates;candidates.reserve(worst.size());
        while(!worst.empty()){candidates.push_back(worst.top());worst.pop();}
        std::ranges::sort(candidates);
        std::uint64_t live=report.liveBefore;
        std::array<std::uint64_t,3> stagedFlips{};
        const std::array thresholds{limits.searchBelowMeanRatio,.01,.1};
        struct Ranked {Complex::Proposal patch;double minimumMean,minimumJacobian;};
        for(const auto& [score,id]:candidates){
            charge();(void)score;
            if(report.stagedEdits>=limits.maximumEdits){report.budgetLimited=true;report.diagnostic="edit_budget";break;}
            if(!staged.alive[id])continue;
            ++report.candidates;
            std::optional<Ranked> best;
            const auto consider=[&](Complex::Proposal patch){
                if(report.proposals>=limits.maximumProposals)throw Limit{"proposal_budget"};
                ++report.proposals;
                if(patch.kind!=23&&patch.kind!=32&&patch.kind!=44)return;
                if(patch.added.size()>staged.maximumStoredCells-staged.cells.size())return;
                double oldMean=1,oldJacobian=1,newMean=1,newJacobian=1;
                std::array<unsigned,3> before{},after{};
                for(auto old:patch.removed){const auto q=shape(staged,staged.cells[old],charge);
                    oldMean=std::min(oldMean,q.mean);oldJacobian=std::min(oldJacobian,q.jacobian);
                    for(unsigned i=0;i<thresholds.size();++i)before[i]+=q.mean<thresholds[i];}
                for(const auto& cell:patch.added){const auto q=shape(staged,cell,charge);
                    newMean=std::min(newMean,q.mean);newJacobian=std::min(newJacobian,q.jacobian);
                    for(unsigned i=0;i<thresholds.size();++i)after[i]+=q.mean<thresholds[i];}
                if(!(newMean>oldMean*(1+1e-6))||newJacobian<oldJacobian||
                   std::ranges::any_of(std::views::iota(0U,3U),[&](unsigned i){return after[i]>before[i];}))return;
                if(best&&(newMean<best->minimumMean||(newMean==best->minimumMean&&
                    newJacobian<=best->minimumJacobian)))return;
                ++report.certifications;
                if(!staged.preservesSurvivors(patch))return;
                best=Ranked{std::move(patch),newMean,newJacobian};
            };
            const auto cell=staged.cells[id];
            for(const auto face:Complex::faceKeys(cell))if(auto p=staged.flip23(face,std::nullopt))consider(std::move(*p));
            for(const auto edge:Complex::edgeKeys(cell)){
                const auto owners=staged.edges.find(edge);
                if(owners==staged.edges.end()||owners->second.size()<3||owners->second.size()>4)continue;
                for(auto& p:staged.edgeFlips(edge,std::nullopt))consider(std::move(p));
            }
            if(!best)continue;
            staged.apply(best->patch);
            live=live-best->patch.removed.size()+best->patch.added.size();
            ++report.stagedEdits;
            ++stagedFlips[best->patch.kind==23?0:best->patch.kind==32?1:2];
        }
        if(!report.stagedEdits)throw Finished{"no_admissible_improvement"};
        double minimumMean=1,minimumJacobian=1;std::uint64_t severe=0;
        for(std::size_t id=0;id<staged.cells.size();++id){
            charge();if(!staged.alive[id]||guard(staged.cells[id],topology.sourceNodes))continue;
            const auto q=shape(staged,staged.cells[id],charge);
            minimumMean=std::min(minimumMean,q.mean);minimumJacobian=std::min(minimumJacobian,q.jacobian);
            severe+=q.mean<limits.searchBelowMeanRatio;
        }
        if(minimumMean<report.minimumMeanRatioBefore||minimumJacobian<report.minimumScaledJacobianBefore||
           severe>report.severeBefore)throw Finished{"no_global_improvement"};
        for(const auto face:temporaryFaces){charge();staged.protectedFaces.erase(face);}
        for(const auto edge:temporaryEdges){charge();staged.protectedEdges.erase(edge);}
        cancelled();
        staged.chargeWork={};
        static_assert(std::is_nothrow_move_assignable_v<Complex>);
        topology.complex=std::move(staged);topology.liveCells=live;
        report.changed=true;report.completed=true;report.edits=report.stagedEdits;
        report.flips23=stagedFlips[0];report.flips32=stagedFlips[1];report.flips44=stagedFlips[2];
        report.liveAfter=live;report.severeAfter=severe;
        report.minimumMeanRatioAfter=minimumMean;report.minimumScaledJacobianAfter=minimumJacobian;
        if(!report.budgetLimited)report.diagnostic=severe?"improved_with_remaining_candidates":"improved";
    }catch(const Finished& finished){report.completed=true;report.diagnostic=finished.reason;
    }catch(const Cancelled&){report.cancelled=true;report.diagnostic="cancelled";
    }catch(const Limit& limit){report.budgetLimited=true;report.diagnostic=limit.reason;
    }catch(const Invalid&){report.diagnostic="invalid_topology";
    }catch(const std::length_error&){report.budgetLimited=true;report.diagnostic="storage_budget";
    }catch(const std::bad_alloc&){report.diagnostic="allocation_failed";}
    return report;
}

} // namespace paacs::geo::detail
