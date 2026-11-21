// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "test_harness.hpp"

#include "paacs/algo/Tet10CurvingMesher.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

#if defined(_WIN32)
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

namespace {
using namespace paacs;

std::array<Vec3,10> straightTet() {
    const std::array<Vec3,4> c{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    return {{c[0],c[1],c[2],c[3],
             (c[0]+c[1])*.5,(c[1]+c[2])*.5,(c[2]+c[0])*.5,
             (c[0]+c[3])*.5,(c[1]+c[3])*.5,(c[2]+c[3])*.5}};
}

algo::Tet10CertificatePolicy policy() {
    algo::Tet10CertificatePolicy p;
    p.minimumJacobianRatio=0.01;
    p.maxDepth=4;
    p.maxVisitedSubcells=4096;
    p.maxCoefficientTests=100000;
    p.maxWorkUnits=1000000;
    return p;
}

class OwnershipProvider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{1701,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    mutable std::uint64_t stampCalls=0,ownershipCalls=0,envelopeCalls=0;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        ++stampCalls;
        return current;
    }
    algo::VolumeFaceOwnership ownership(algo::FacePersistentId face) const override {
        ++ownershipCalls;
        algo::VolumeFaceOwnership result;result.providerStamp=current;
        result.facePersistentId=face;
        const auto it=rows.find(face);if(it!=rows.end())result.owners=it->second;
        return result;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(std::span<const Vec3>) const override {
        ++envelopeCalls;
        algo::ProviderRecoveryEnvelope result;result.providerStamp=current;
        result.availability=algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;
        return result;
    }
};

class MidpointProvider final : public algo::Tet10PlacementProvider {
public:
    mutable algo::Tet10PlacementStamp current{4401,1};
    mutable std::uint64_t calls=0,stampCalls=0;
    mutable bool changeAfterPlace=false;
    mutable bool incidentFacetsCanonical=true;
    mutable bool endpointCoordinatesAligned=true;
    mutable std::vector<algo::Tet10SourceEdgeKey> requestedEdges;
    const std::map<algo::RefinementVertexKey,Vec3>* expectedCoordinates=nullptr;
    bool throwFailure=false,throwAllocation=false,nonfinite=false,
         foreignResult=false,staleResult=false,distortFirst=false,distortAll=false;
    int throwStampCall=0;
    bool throwStampAllocation=false;
    algo::Tet10PlacementStatus status=algo::Tet10PlacementStatus::Placed;
    algo::Tet10PlacementStamp stamp() const override{
        ++stampCalls;
        if(throwStampCall>0&&stampCalls==static_cast<std::uint64_t>(throwStampCall)) {
            if(throwStampAllocation)throw std::bad_alloc{};
            throw std::runtime_error("stamp failure");
        }
        return current;
    }
    algo::Tet10PlacementResult place(const algo::Tet10PlacementRequest& r) const override {
        ++calls;
        if(throwAllocation)throw std::bad_alloc{};
        if(throwFailure)throw std::runtime_error("placement failure");
        requestedEdges.push_back({r.endpoints});
        if(expectedCoordinates)for(unsigned i=0;i<2;++i) {
            const auto found=expectedCoordinates->find(r.endpoints[i]);
            endpointCoordinatesAligned=endpointCoordinatesAligned&&
                found!=expectedCoordinates->end()&&
                found->second.x==r.endpointCoordinates[i].x&&
                found->second.y==r.endpointCoordinates[i].y&&
                found->second.z==r.endpointCoordinates[i].z;
        }
        incidentFacetsCanonical=incidentFacetsCanonical&&
            std::is_sorted(r.incidentFacets.begin(),r.incidentFacets.end());
        const auto returned=foreignResult?algo::Tet10PlacementStamp{999,9}:
            (staleResult?algo::Tet10PlacementStamp{current.persistentId,current.revision+1}:current);
        const Vec3 coordinate=nonfinite?Vec3{std::numeric_limits<double>::quiet_NaN(),0,0}:
            ((distortFirst&&calls==1)||distortAll?
             Vec3{-10.0-double(calls),-10,-10}:
            Vec3{(r.endpointCoordinates[0].x+r.endpointCoordinates[1].x)*.5,
             (r.endpointCoordinates[0].y+r.endpointCoordinates[1].y)*.5,
             (r.endpointCoordinates[0].z+r.endpointCoordinates[1].z)*.5});
        if(changeAfterPlace)++current.revision;
        return {status,71,returned,coordinate};
    }
};

class CylinderProvider final : public algo::Tet10PlacementProvider {
public:
    algo::Tet10PlacementStamp current{4402,1};
    double curveFraction=.125;
    mutable std::uint64_t calls=0;
    mutable bool moved=false;
    algo::Tet10PlacementStamp stamp() const noexcept override{return current;}
    algo::Tet10PlacementResult place(const algo::Tet10PlacementRequest& r) const override {
        ++calls;
        Vec3 p{(r.endpointCoordinates[0].x+r.endpointCoordinates[1].x)*.5,
               (r.endpointCoordinates[0].y+r.endpointCoordinates[1].y)*.5,
               (r.endpointCoordinates[0].z+r.endpointCoordinates[1].z)*.5};
        const bool lateral=r.kind==algo::Tet10PlacementKind::Surface&&
                           (r.classification.tag==1||r.classification.tag==2);
        const bool circular=r.kind==algo::Tet10PlacementKind::Curve&&
                            r.endpointCoordinates[0].z==r.endpointCoordinates[1].z;
        if(lateral||circular) {
            const double length=std::hypot(p.x,p.y);
            if(length>0.0){
                const double target=length+curveFraction*(.5-length);
                if(target!=length)moved=true;
                p.x*=target/length;p.y*=target/length;
            }
        }
        return {algo::Tet10PlacementStatus::Placed,0,current,p};
    }
};

void configure(const geo::Model& model,OwnershipProvider& provider) {
    Vec3 lo{},hi{};model.solids().front()->bbox(lo,hi);const Vec3 center=(lo+hi)*.5;
    for(const auto* face:model.faces()) {
        Vec2 uvLo{},uvHi{};face->paramBounds(uvLo,uvHi);const Vec2 uv=(uvLo+uvHi)*.5;
        const int sense=face->normalAt(uv).dot(face->pointAt(uv)-center)>0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]={{0,sense}};
    }
}

algo::T4BRequest routeRequest(const OwnershipProvider& provider,
                              const MidpointProvider& placement) {
    algo::T4BRequest r;
    r.source.source.conformingRequest.expectedSurfaceSetStamp={17001,1};
    r.source.source.conformingRequest.expectedOwnershipProviderStamp=provider.current;
    r.source.source.conformingRequest.envelopeMode=algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    auto& refine=r.source.source.refinement;
    refine.maxEdgeToTarget=1e9;refine.maxRadiusEdge=1e9;refine.maxInsertions=0;
    refine.maxWorkUnits=200000;refine.maxOutputNodes=100000;refine.maxOutputTets=100000;
    refine.maxReportRows=5;refine.maxReportWitnesses=100000;
    refine.maxReportTopologyHandles=300000;
    r.source.quality.mode=algo::T4AQualityMode::ReportOnly;
    r.source.quality.minimumRadiusRatio=.2;r.source.quality.maxResultNodes=100000;
    r.source.quality.maxResultTets=100000;r.source.quality.maxReportRows=1000;
    r.source.quality.maxWorstTailRows=64;
    r.expectedPlacementStamp=placement.current;
    r.policy={.001,4,4096,100000,10000000,
              100000,100000,1000000,100000,100000,200000};
    return r;
}

algo::T4BRequest routeRequest(const OwnershipProvider& provider,
                              const CylinderProvider& placement) {
    MidpointProvider midpoint;midpoint.current=placement.current;
    return routeRequest(provider,midpoint);
}

void requireReportXor(const algo::T4BOutcome& result) {
    CHECK(result.report!=nullptr&&result.reportComplete);
    if(!result.report)return;
    std::set<algo::Tet10SourceEdgeKey> edges;
    for(const auto& row:result.report->placements)CHECK(edges.insert(row.edge).second);
    for(const auto& key:result.report->unattemptedEdges)CHECK(edges.insert(key).second);
    CHECK(std::is_sorted(result.report->unattemptedEdges.begin(),
                         result.report->unattemptedEdges.end()));
    std::vector<algo::Tet10SourceEdgeKey> ordered;
    for(const auto& row:result.report->placements)ordered.push_back(row.edge);
    ordered.insert(ordered.end(),result.report->unattemptedEdges.begin(),
                   result.report->unattemptedEdges.end());
    CHECK(std::is_sorted(ordered.begin(),ordered.end()));
    std::set<std::uint64_t> cells;
    for(const auto& row:result.report->certificates)
        CHECK(cells.insert(row.sourceCellOrdinal).second);
    for(const auto key:result.report->unattemptedCells)CHECK(cells.insert(key).second);
    CHECK(std::is_sorted(result.report->unattemptedCells.begin(),
                         result.report->unattemptedCells.end()));
    if(result.retainedResult&&result.retainedResult->privateMesh) {
        std::vector<std::uint64_t> expected;
        for(Index ci=0;ci<result.retainedResult->privateMesh->cellCount();++ci)
            if(result.retainedResult->privateMesh->cell(ci).type==mesh::CellType::Tet10)
                expected.push_back(static_cast<std::uint64_t>(ci));
        const std::vector<std::uint64_t> observed(cells.begin(),cells.end());
        CHECK(observed==expected);
    }
}

void straightCertificateGolden() {
    TEST("t4b.certificate.straight_golden");
    const auto result=algo::certifyTet10Jacobian(straightTet(),policy());
    CHECK(result.status==algo::Tet10CertificateStatus::Certified);
    CHECK(result.complete);
    CHECK(result.metrics.visitedSubcells==1);
    CHECK(result.metrics.coefficientTests==20);
    for(const auto& value:result.rootDeterminant) {
        CHECK(value.lo>0.0);
        CHECK(value.hi>=value.lo);
    }
}

void blossomConstantOracle() {
    TEST("t4b.certificate.blossom_children");
    std::array<algo::Tet10Interval,20> parent{};
    for(auto& value:parent)value={7.0,7.0};
    for(std::uint32_t child=0;child<8;++child) {
        std::array<algo::Tet10Interval,20> restricted{};
        CHECK(algo::restrictTet10BernsteinChildForTesting(
            parent,child,restricted));
        for(const auto& value:restricted) {
            CHECK(value.lo<=7.0);
            CHECK(value.hi>=7.0);
        }
    }
}

using Big=boost::multiprecision::cpp_int;
struct Rat {
    Big n=0,d=1;
};
Big bgcd(Big a,Big b){if(a<0)a=-a;if(b<0)b=-b;while(b!=0){Big r=a%b;a=b;b=r;}return a;}
Rat normalized(Big n,Big d){
    if(d<0){n=-n;d=-d;}const Big g=bgcd(n,d);if(g!=0){n/=g;d/=g;}return {n,d};
}
Rat plus(const Rat& a,const Rat& b){return normalized(a.n*b.d+b.n*a.d,a.d*b.d);}
Rat times(const Rat& a,const Rat& b){return normalized(a.n*b.n,a.d*b.d);}
Rat half(unsigned value){return {Big(value),Big(2)};}
Rat scaled(const Rat& a,long long n,long long d=1){return normalized(a.n*n,a.d*d);}

constexpr std::array<std::array<unsigned,4>,20> kAlpha{{
    {{3,0,0,0}},{{2,1,0,0}},{{2,0,1,0}},{{2,0,0,1}},
    {{1,2,0,0}},{{1,1,1,0}},{{1,1,0,1}},{{1,0,2,0}},
    {{1,0,1,1}},{{1,0,0,2}},{{0,3,0,0}},{{0,2,1,0}},
    {{0,2,0,1}},{{0,1,2,0}},{{0,1,1,1}},{{0,1,0,2}},
    {{0,0,3,0}},{{0,0,2,1}},{{0,0,1,2}},{{0,0,0,3}}}};
constexpr std::array<std::array<unsigned,4>,10> kBaryTwice{{
    {{2,0,0,0}},{{0,2,0,0}},{{0,0,2,0}},{{0,0,0,2}},
    {{1,1,0,0}},{{0,1,1,0}},{{1,0,1,0}},
    {{1,0,0,1}},{{0,1,0,1}},{{0,0,1,1}}}};
constexpr std::array<std::array<unsigned,4>,8> kChild{{
    {{0,4,6,7}},{{1,4,8,5}},{{2,6,5,9}},{{3,7,9,8}},
    {{4,9,6,7}},{{4,9,7,8}},{{4,9,8,5}},{{4,9,5,6}}}};

Rat oracleBlossom(const std::array<Rat,20>& parent,
                  const std::array<std::array<unsigned,4>,3>& args) {
    std::map<std::array<unsigned,4>,Rat> current;
    for(std::size_t i=0;i<20;++i)current.emplace(kAlpha[i],parent[i]);
    unsigned degree=3;
    for(const auto& argument:args) {
        std::map<std::array<unsigned,4>,Rat> next;
        for(unsigned a0=degree-1;;--a0) {
            for(unsigned a1=degree-1-a0;;--a1) {
                for(unsigned a2=degree-1-a0-a1;;--a2) {
                    const unsigned a3=degree-1-a0-a1-a2;
                    std::array<unsigned,4> beta{{a0,a1,a2,a3}};
                    Rat value{};
                    for(unsigned q=0;q<4;++q) {
                        auto key=beta;++key[q];
                        value=plus(value,times(half(argument[q]),current.at(key)));
                    }
                    next.emplace(beta,value);
                    if(a2==0)break;
                }
                if(a1==0)break;
            }
            if(a0==0)break;
        }
        current.swap(next);--degree;
    }
    return current.begin()->second;
}

void blossomExactRationalOracle() {
    TEST("t4b.certificate.blossom_exact_rational");
    std::array<algo::Tet10Interval,20> parent{};std::array<Rat,20> rational{};
    for(std::size_t i=0;i<20;++i) {
        const int value=static_cast<int>(i)*3-17;
        parent[i]={double(value),double(value)};rational[i]={Big(value),Big(1)};
    }
    for(unsigned child=0;child<8;++child) {
        std::array<algo::Tet10Interval,20> actual{};
        CHECK(algo::restrictTet10BernsteinChildForTesting(parent,child,actual));
        for(std::size_t ai=0;ai<20;++ai) {
            std::array<std::array<unsigned,4>,3> args{};unsigned at=0;
            for(unsigned v=0;v<4;++v)for(unsigned n=0;n<kAlpha[ai][v];++n)
                args[at++]=kBaryTwice[kChild[child][v]];
            const auto expected=oracleBlossom(rational,args);
            const double value=expected.n.convert_to<double>()/
                               expected.d.convert_to<double>();
            CHECK(actual[ai].lo<=value&&actual[ai].hi>=value);
        }
    }
}

struct RExp { unsigned r,s,t; };
constexpr std::array<RExp,20> kExp{{
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},{1,1,0},{1,0,1},
    {0,2,0},{0,1,1},{0,0,2},{3,0,0},{2,1,0},{2,0,1},{1,2,0},
    {1,1,1},{1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}}};
int rindex(unsigned r,unsigned s,unsigned t){
    for(std::size_t i=0;i<20;++i)if(kExp[i].r==r&&kExp[i].s==s&&kExp[i].t==t)return int(i);
    return -1;
}
struct RPoly { std::array<Rat,20> c{}; };
RPoly rpconst(long long n){RPoly p;p.c[0]={Big(n),Big(1)};return p;}
RPoly rpvar(unsigned axis){RPoly p;p.c[axis+1]={Big(1),Big(1)};return p;}
RPoly rpadd(const RPoly&a,const RPoly&b){RPoly p;for(unsigned i=0;i<20;++i)p.c[i]=plus(a.c[i],b.c[i]);return p;}
RPoly rpscale(const RPoly&a,long long n,long long d=1){RPoly p;for(unsigned i=0;i<20;++i)p.c[i]=scaled(a.c[i],n,d);return p;}
RPoly rpmul(const RPoly&a,const RPoly&b){
    RPoly p;
    for(unsigned i=0;i<20;++i)for(unsigned j=0;j<20;++j){
        const unsigned r=kExp[i].r+kExp[j].r,s=kExp[i].s+kExp[j].s,
                       t=kExp[i].t+kExp[j].t;
        if(r+s+t<=3)p.c[static_cast<std::size_t>(rindex(r,s,t))]=
            plus(p.c[static_cast<std::size_t>(rindex(r,s,t))],times(a.c[i],b.c[j]));
    }
    return p;
}
RPoly rpderiv(const RPoly&a,unsigned axis){
    RPoly p;
    for(unsigned i=0;i<20;++i){
        const unsigned n=axis==0?kExp[i].r:(axis==1?kExp[i].s:kExp[i].t);
        if (!n) continue;
        unsigned r=kExp[i].r,s=kExp[i].s,t=kExp[i].t;
        if(axis==0)--r;else if(axis==1)--s;else --t;
        const auto k=static_cast<std::size_t>(rindex(r,s,t));
        p.c[k]=plus(p.c[k],scaled(a.c[i],n));
    }
    return p;
}
unsigned smallChoose(unsigned n,unsigned k){
    if (k>n) return 0;
    if (k==0||k==n) return 1;
    if (n==2) return 2;
    return n==3?3:1;
}
unsigned smallFactorial(unsigned n){return n<2?1:(n==2?2:6);}
unsigned exactTetraDenominator(const RExp& e){
    return 6/(smallFactorial(e.r)*smallFactorial(e.s)*
              smallFactorial(e.t)*smallFactorial(3-e.r-e.s-e.t));
}

std::array<Rat,20> exactRootDeterminant() {
    const auto r=rpvar(0),s=rpvar(1),t=rpvar(2);
    const std::array<RPoly,4> l{{
        rpadd(rpconst(1),rpscale(rpadd(rpadd(r,s),t),-1)),r,s,t}};
    std::array<RPoly,10> shape{};
    for(unsigned i=0;i<4;++i)shape[i]=rpadd(rpscale(rpmul(l[i],l[i]),2),rpscale(l[i],-1));
    constexpr std::array<std::array<unsigned,2>,6> edge{{
        {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
    for(unsigned e=0;e<6;++e)shape[4+e]=rpscale(rpmul(l[edge[e][0]],l[edge[e][1]]),4);
    // Exact normalized positions for straightTet() with node 4 x=.25.
    // midpoint(min,max)=.5 and k=ilogb(.5)+1=0, hence q=p-.5.
    std::array<std::array<Rat,3>,10> q{};
    const auto nodes=straightTet();
    for(unsigned i=0;i<10;++i)for(unsigned a=0;a<3;++a){
        const double value=a==0?nodes[i].x:(a==1?nodes[i].y:nodes[i].z);
        const long long quarters=static_cast<long long>(std::llround(value*4.0-2.0));
        q[i][a]={Big(quarters),Big(4)};
    }
    q[4][0]={Big(-1),Big(4)};
    std::array<RPoly,3> x{};
    for(unsigned a=0;a<3;++a)for(unsigned n=0;n<10;++n){
        RPoly term=shape[n];for(auto& c:term.c)c=times(c,q[n][a]);x[a]=rpadd(x[a],term);
    }
    RPoly j[3][3];for(unsigned a=0;a<3;++a)for(unsigned d=0;d<3;++d)j[a][d]=rpderiv(x[a],d);
    const auto det=rpadd(rpadd(
        rpmul(j[0][0],rpadd(rpmul(j[1][1],j[2][2]),rpscale(rpmul(j[1][2],j[2][1]),-1))),
        rpscale(rpmul(j[0][1],rpadd(rpmul(j[1][0],j[2][2]),rpscale(rpmul(j[1][2],j[2][0]),-1))),-1)),
        rpmul(j[0][2],rpadd(rpmul(j[1][0],j[2][1]),rpscale(rpmul(j[1][1],j[2][0]),-1))));
    std::array<Rat,20> result{};
    for(unsigned ai=0;ai<20;++ai)for(unsigned bi=0;bi<20;++bi){
        const auto e=kExp[bi];
        if(e.r<=kAlpha[ai][1]&&e.s<=kAlpha[ai][2]&&e.t<=kAlpha[ai][3]){
            const unsigned num=smallChoose(kAlpha[ai][1],e.r)*
                smallChoose(kAlpha[ai][2],e.s)*smallChoose(kAlpha[ai][3],e.t);
            result[ai]=plus(result[ai],scaled(
                det.c[bi],num,exactTetraDenominator(e)));
        }
    }
    return result;
}

void initialCoefficientExactOracle() {
    TEST("t4b.certificate.initial_exact_rational");
    auto nodes=straightTet();nodes[4].x=.25;
    const auto result=algo::certifyTet10Jacobian(nodes,policy());
    const auto oracle=exactRootDeterminant();
    for(unsigned i=0;i<20;++i){
        const double value=oracle[i].n.convert_to<double>()/oracle[i].d.convert_to<double>();
        CHECK(result.rootDeterminant[i].lo<=value&&result.rootDeterminant[i].hi>=value);
    }
}

void taxonomyAndBudgets() {
    TEST("t4b.certificate.taxonomy");
    auto nodes=straightTet();
    auto p=policy();p.maxVisitedSubcells=0;
    auto result=algo::certifyTet10Jacobian(nodes,p);
    CHECK(result.status==algo::Tet10CertificateStatus::IncompleteSubcells);
    nodes[4].x=std::numeric_limits<double>::quiet_NaN();
    result=algo::certifyTet10Jacobian(nodes,policy());
    CHECK(result.status==algo::Tet10CertificateStatus::InvalidInput ||
          result.status==algo::Tet10CertificateStatus::RangeUnrepresentable);
    nodes=straightTet();std::swap(nodes[1],nodes[2]);
    result=algo::certifyTet10Jacobian(nodes,policy());
    CHECK(result.status==algo::Tet10CertificateStatus::InvalidInput);
}

void hardCertificateWorkBudget() {
    TEST("t4b.certificate.hard_work_budget");
    auto p=policy();p.maxWorkUnits=2048;
    const auto result=algo::certifyTet10Jacobian(straightTet(),p);
    CHECK(result.status==algo::Tet10CertificateStatus::IncompleteWork);
    CHECK(result.metrics.workUnits<=p.maxWorkUnits);
}

void curvedTaxonomy() {
    TEST("t4b.certificate.curved_taxonomy");
    auto below=straightTet();below[4].x=.3;
    auto threshold=policy();threshold.minimumJacobianRatio=.8;threshold.maxDepth=6;
    threshold.maxVisitedSubcells=100000;threshold.maxCoefficientTests=3000000;
    threshold.maxWorkUnits=10000000;
    const auto belowResult=algo::certifyTet10Jacobian(below,threshold);
    auto positivity=threshold;positivity.minimumJacobianRatio=0.0;
    CHECK(algo::certifyTet10Jacobian(below,positivity).status==
          algo::Tet10CertificateStatus::Certified);
    CHECK(belowResult.status==algo::Tet10CertificateStatus::PositiveBelowThreshold);
    auto inverted=straightTet();inverted[4].x=-1.0;
    CHECK(algo::certifyTet10Jacobian(inverted,threshold).status==
          algo::Tet10CertificateStatus::Nonpositive);
    auto unitThreshold=policy();unitThreshold.minimumJacobianRatio=1.0;
    const auto unitResult=algo::certifyTet10Jacobian(straightTet(),unitThreshold);
    CHECK(unitResult.status==algo::Tet10CertificateStatus::PositiveBelowThreshold||
          unitResult.status==algo::Tet10CertificateStatus::IncompleteDepth);
}

void extremeCertificateFrames() {
    TEST("t4b.certificate.extreme_frames");
    const auto base=straightTet();
    auto transformed=[&](double scale,double translation) {
        auto nodes=base;
        for(auto& p:nodes){p.x=translation+p.x*scale;
            p.y=translation+p.y*scale;p.z=translation+p.z*scale;}
        return nodes;
    };
    const auto translated=algo::certifyTet10Jacobian(
        transformed(1.0e90,1.0e100),policy());
    CHECK(translated.status==algo::Tet10CertificateStatus::Certified);
    const auto large=algo::certifyTet10Jacobian(
        transformed(1.0e300,0.0),policy());
    CHECK(large.status==algo::Tet10CertificateStatus::Certified);
    const double tiny=std::scalbn(1.0,-1060);
    const auto subnormal=algo::certifyTet10Jacobian(
        transformed(tiny,0.0),policy());
    CHECK(subnormal.status==algo::Tet10CertificateStatus::Certified);
    auto insufficient=policy();insufficient.maxWorkUnits=0;
    CHECK(algo::certifyTet10Jacobian(base,insufficient).status==
          algo::Tet10CertificateStatus::IncompleteWork);
}

void straightPublicRoute() {
    TEST("t4b.route.straight_box");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    mesh::Mesh destination;
    const auto outcome=algo::meshModelCurvedTet10(*box,spec,
        routeRequest(ownership,placement),destination,ownership,placement);
    CHECK(outcome.status==algo::T4BStatus::Published);
    CHECK(outcome.retainedResult&&outcome.retainedResult->certificationComplete);
    CHECK(destination.countCells(mesh::CellType::Tet10)>0);
    CHECK(destination.countCells(mesh::CellType::Tet4)==0);
    CHECK(outcome.retainedResult&&outcome.retainedResult->audit.boundaryOrphans==0);
    CHECK(outcome.retainedResult&&outcome.retainedResult->audit.boundaryUnexpected==0);
    CHECK(outcome.retainedResult&&outcome.retainedResult->metrics.perEdgeGlobalScans==0);
    CHECK(outcome.sourceReceipt!=nullptr&&outcome.sourceOutcome==nullptr);
    CHECK(outcome.reportComplete&&outcome.report!=nullptr);
    requireReportXor(outcome);
    if(outcome.retainedResult&&outcome.report) {
        const auto& result=*outcome.retainedResult;
        CHECK(result.edgeToMidside.size()==outcome.report->placements.size());
        CHECK(result.certificates.size()==outcome.report->certificates.size());
        CHECK(placement.calls==result.metrics.curveEdges+result.metrics.surfaceEdges);
        CHECK(placement.calls==result.metrics.providerCalls);
        CHECK(placement.stampCalls==result.metrics.providerStampCalls);
        CHECK(std::is_sorted(placement.requestedEdges.begin(),
                             placement.requestedEdges.end()));
        CHECK(std::adjacent_find(placement.requestedEdges.begin(),
                                 placement.requestedEdges.end())==
              placement.requestedEdges.end());
        CHECK(placement.incidentFacetsCanonical);
        std::set<Index> midsides;
        for(const auto& row:result.edgeToMidside) {
            CHECK(row.midsideIndex>=0&&row.midsideIndex<destination.nodeCount());
            CHECK(midsides.insert(row.midsideIndex).second);
            CHECK(destination.nodeClass(row.midsideIndex).dim==row.classification.dim);
            CHECK(destination.nodeClass(row.midsideIndex).tag==row.classification.tag);
        }
        CHECK(result.certificates.size()==
              static_cast<std::size_t>(destination.countCells(mesh::CellType::Tet10)));
        CHECK(result.sourceNodeToOutput.size()<static_cast<std::size_t>(destination.nodeCount()));
        CHECK(result.sourceCellToOutput.size()==static_cast<std::size_t>(destination.cellCount()));
        for(std::size_t i=0;i<result.sourceNodeToOutput.size();++i)
            CHECK(result.sourceNodeToOutput[i]==static_cast<Index>(i));
        for(std::size_t i=0;i<result.sourceCellToOutput.size();++i)
            CHECK(result.sourceCellToOutput[i]==static_cast<Index>(i));
        for(const auto& row:result.certificates)
            CHECK(row.conservativeMinimumLowerBound>0.0);
        CHECK(outcome.report->boundaryCompatibility==
              algo::T4BBoundaryCompatibility::ProviderTrusted);
        std::map<std::array<Index,2>,Index> midsideByCornerPair;
        auto record=[&](Index a,Index b,Index midside) {
            std::array<Index,2> key{{a,b}};std::sort(key.begin(),key.end());
            const auto [it,inserted]=midsideByCornerPair.emplace(key,midside);
            CHECK(inserted||it->second==midside);
        };
        for(Index ci=0;ci<destination.cellCount();++ci) {
            const auto cell=destination.cell(ci);const Index* n=destination.cellNodes(ci);
            if(cell.type==mesh::CellType::Line3)record(n[0],n[1],n[2]);
            else if(cell.type==mesh::CellType::Tri6) {
                record(n[0],n[1],n[3]);record(n[1],n[2],n[4]);
                record(n[2],n[0],n[5]);
            } else if(cell.type==mesh::CellType::Tet10) {
                static constexpr std::array<std::array<unsigned,2>,6> slots{{
                    {{0,1}},{{1,2}},{{2,0}},{{0,3}},{{1,3}},{{2,3}}}};
                for(unsigned e=0;e<6;++e)record(
                    n[slots[e][0]],n[slots[e][1]],n[4+e]);
            }
        }
        CHECK(midsideByCornerPair.size()==result.edgeToMidside.size());
        CHECK(result.strictHash==mesh::meshHash(destination));
    }
}

void sourceCellOrdinalAttribution() {
    TEST("t4b.route.source_cell_ordinals");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    mesh::Mesh destination;const auto outcome=algo::meshModelCurvedTet10(
        *box,spec,routeRequest(ownership,placement),destination,ownership,placement);
    CHECK(outcome.published()&&outcome.report);
    std::vector<std::uint64_t> expected,reported;
    for(Index ci=0;ci<destination.cellCount();++ci)
        if(destination.cell(ci).type==mesh::CellType::Tet10)
            expected.push_back(static_cast<std::uint64_t>(ci));
    if(outcome.report)for(const auto& row:outcome.report->certificates)
        reported.push_back(row.sourceCellOrdinal);
    CHECK(reported==expected);
}

void providerEndpointCoordinateAlignment() {
    TEST("t4b.route.provider_endpoint_coordinate_alignment");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider first;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    mesh::Mesh reference;const auto initial=algo::meshModelCurvedTet10(
        *box,spec,routeRequest(ownership,first),reference,ownership,first);
    CHECK(initial.published()&&initial.retainedResult&&
          initial.retainedResult->sourceArtifact);
    if(!initial.retainedResult||!initial.retainedResult->sourceArtifact)return;
    std::map<algo::RefinementVertexKey,Vec3> coordinates;
    const auto& artifact=*initial.retainedResult->sourceArtifact;
    for(const auto& [handle,key]:artifact.vertexKeyByHandle) {
        CHECK(handle.slot>=0);
        if(handle.slot<0)continue;
        const auto slot=static_cast<std::size_t>(handle.slot);
        CHECK(slot<artifact.publishedGlobalByAuthoritativeNodeSlot.size());
        if(slot>=artifact.publishedGlobalByAuthoritativeNodeSlot.size())continue;
        const Index output=artifact.publishedGlobalByAuthoritativeNodeSlot[slot];
        CHECK(output>=0&&output<reference.nodeCount());
        if(output>=0&&output<reference.nodeCount())coordinates.emplace(key,reference.node(output));
    }
    MidpointProvider checked;checked.expectedCoordinates=&coordinates;
    mesh::Mesh destination;const auto outcome=algo::meshModelCurvedTet10(
        *box,spec,routeRequest(ownership,checked),destination,ownership,checked);
    CHECK(outcome.published());
    CHECK(checked.calls>0);
    CHECK(checked.endpointCoordinatesAligned);
}

void curvedCylinderRoute() {
    TEST("t4b.route.curved_cylinder");
    auto kernel=geo::makeBuiltinKernel();auto cylinder=kernel->makeCylinder({0,0,0},.5,1.0);
    OwnershipProvider ownership;configure(*cylinder,ownership);CylinderProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.15);
    auto request=routeRequest(ownership,placement);
    request.policy.minimumCertifiedJacobianRatio=1e-12;
    request.policy.maxCertificateDepth=8;
    request.policy.maxCertificateCells=1000000;
    request.policy.maxCoefficientTests=20000000;
    request.policy.maxWorkUnits=100000000;
    mesh::Mesh destination;
    const auto result=algo::meshModelCurvedTet10(*cylinder,spec,request,
        destination,ownership,placement);
    CHECK(result.status==algo::T4BStatus::Published);
    CHECK(result.retainedResult&&result.retainedResult->certificationComplete);
    CHECK(result.retainedResult&&result.retainedResult->audit.boundaryOrphans==0);
    CHECK(result.retainedResult&&result.retainedResult->audit.boundaryUnexpected==0);
    bool curved=false;
    if(result.retainedResult)for(const auto& row:result.retainedResult->edgeToMidside) {
        if(row.classification.dim==2&&
           (row.classification.tag==1||row.classification.tag==2)) {
            const double radius=std::hypot(row.coordinate.x,row.coordinate.y);
            const auto& key=row.edge.vertices;
            (void)key;
            CHECK(radius>0.0&&radius<=.5);curved=true;
        }
    }
    CHECK(curved);
    CHECK(placement.moved);
}

void routeDeterminism() {
    TEST("t4b.route.determinism");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    std::uint64_t hash=0,semantic=0;std::vector<Index> midsides;
    for(const int threads:{1,4,0}) {
        OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
        sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
        spec.threads=threads;mesh::Mesh destination;
        const auto result=algo::meshModelCurvedTet10(*box,spec,
            routeRequest(ownership,placement),destination,ownership,placement);
        CHECK(result.status==algo::T4BStatus::Published);
        std::vector<Index> current;
        if(result.retainedResult)for(const auto& row:result.retainedResult->edgeToMidside)
            current.push_back(row.midsideIndex);
        const auto strict=mesh::meshHash(destination);
        const auto sem=mesh::semanticHashEx(destination);
        CHECK(sem.canonical);
        if(threads==1){hash=strict;semantic=sem.value;midsides=current;}
        else {CHECK(strict==hash);CHECK(sem.value==semantic);CHECK(current==midsides);}
    }
}

void publicRouteScaling() {
    TEST("t4b.route.public_scaling");
    // Each target-resolution row represents approximately 4x the preceding
    // three-dimensional problem size (h /= cbrt(4)).
    const std::array<double,4> sizes{{.500,.315,.198,.125}};
    std::array<std::uint64_t,4> cells{},work{};
    for(std::size_t row=0;row<sizes.size();++row) {
        auto kernel=geo::makeBuiltinKernel();
        auto box=kernel->makeBox({0,0,0},{1,1,1});
        OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
        sizing::MeshSpec spec;
        spec.field=std::make_shared<sizing::ConstantField>(sizes[row]);spec.threads=1;
        auto request=routeRequest(ownership,placement);
        request.source.source.conformingRequest.expectedSurfaceSetStamp={18001+row,1};
        request.policy.maxCertificateCells=1000000;
        request.policy.maxCoefficientTests=100000000;
        request.policy.maxWorkUnits=1000000000;
        request.policy.maxResultNodes=500000;
        request.policy.maxResultCells=500000;
        request.policy.maxResultConnectivity=5000000;
        request.policy.maxPlacementRows=500000;
        request.policy.maxCertificateRows=500000;
        request.policy.maxUnattemptedRows=1000000;
        mesh::Mesh output;const auto result=algo::meshModelCurvedTet10(
            *box,spec,request,output,ownership,placement);
        CHECK(result.published()&&result.retainedResult);
        if(!result.retainedResult)continue;
        const auto& metrics=result.retainedResult->metrics;
        cells[row]=static_cast<std::uint64_t>(output.countCells(mesh::CellType::Tet10));
        work[row]=metrics.workUnits;
        CHECK(metrics.perEdgeGlobalScans==0&&metrics.perCellGlobalScans==0&&
              metrics.contextRebuilds==0);
        CHECK(metrics.accountedHighWaterBytes>0&&metrics.auditRows>0);
        CHECK(result.retainedResult->audit.boundaryOrphans==0&&
              result.retainedResult->audit.boundaryUnexpected==0);
        std::printf("T4B scale h=%.3f cells=%llu work=%llu\n",sizes[row],
            static_cast<unsigned long long>(cells[row]),
            static_cast<unsigned long long>(work[row]));
    }
    for(std::size_t row=1;row<cells.size();++row) {
        CHECK(cells[row]>cells[row-1]);CHECK(work[row]>work[row-1]);
        CHECK(work[row]<12*work[row-1]);
    }
    CHECK(cells.back()>=4*cells.front());
}

void providerFailureContracts() {
    TEST("t4b.route.provider_failures");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    auto run=[&](MidpointProvider& placement) {
        mesh::Mesh destination;destination.addNode({-7,-11,-13},{0,77});
        const auto digest=mesh::meshHash(destination);
        auto result=algo::meshModelCurvedTet10(*box,spec,
            routeRequest(ownership,placement),destination,ownership,placement);
        CHECK(result.destinationUnchanged);
        CHECK(mesh::meshHash(destination)==digest);
        CHECK(result.sourceOutcome==nullptr);
        CHECK(result.sourceReceipt!=nullptr);
        CHECK(result.retainedResult==nullptr);
        CHECK(result.reportComplete);
        requireReportXor(result);
        CHECK(result.report&&result.report->placements.size()==1);
        if(!placement.throwFailure&&!placement.throwAllocation)
            CHECK(result.report&&result.report->placements.front().providerReason==71);
        return result;
    };
    for(const auto status:{algo::Tet10PlacementStatus::Unsupported,
                           algo::Tet10PlacementStatus::Rejected,
                           algo::Tet10PlacementStatus::Failed}) {
        MidpointProvider placement;placement.status=status;
        const auto result=run(placement);
        CHECK(result.status==(status==algo::Tet10PlacementStatus::Unsupported?
            algo::T4BStatus::Unsupported:(status==algo::Tet10PlacementStatus::Rejected?
            algo::T4BStatus::Rejected:algo::T4BStatus::Failed)));
    }
    {MidpointProvider p;p.nonfinite=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::NonfinitePlacement);}
    {MidpointProvider p;p.foreignResult=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderIdentityChanged);}
    {MidpointProvider p;p.staleResult=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::StalePlacementProvider);}
    {MidpointProvider p;p.changeAfterPlace=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::StalePlacementProvider);}
    {MidpointProvider p;p.throwFailure=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(p.stampCalls==3);}
    {MidpointProvider p;p.throwAllocation=true;const auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::AllocationFailed);
     CHECK(p.stampCalls==3);}
}

void providerStampCheckpoints() {
    TEST("t4b.route.provider_stamp_checkpoints");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    auto run=[&](MidpointProvider& provider) {
        mesh::Mesh destination;destination.addNode({59,61,67},{0,91});
        const auto before=mesh::meshHash(destination);
        auto result=algo::meshModelCurvedTet10(*box,spec,
            routeRequest(ownership,provider),destination,ownership,provider);
        CHECK(mesh::meshHash(destination)==before);CHECK(result.destinationUnchanged);
        return result;
    };
    {MidpointProvider p;p.throwStampCall=1;auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(r.sourceOutcome==nullptr&&r.sourceReceipt==nullptr);}
    {MidpointProvider p;p.throwStampCall=1;p.throwStampAllocation=true;auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::AllocationFailed);}
    {MidpointProvider p;p.throwStampCall=2;auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(r.sourceOutcome!=nullptr&&r.sourceReceipt==nullptr);}
    {MidpointProvider p;p.throwStampCall=3;auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(r.sourceOutcome==nullptr&&r.sourceReceipt!=nullptr);
     CHECK(p.calls==0);}
    {MidpointProvider p;p.throwStampCall=4;auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(p.calls==1);}
    MidpointProvider successful;mesh::Mesh published;
    const auto success=algo::meshModelCurvedTet10(*box,spec,
        routeRequest(ownership,successful),published,ownership,successful);
    CHECK(success.published());const auto total=successful.stampCalls;
    CHECK(total>4);
    {MidpointProvider p;p.throwStampCall=static_cast<int>(total-1);auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(r.sourceReceipt!=nullptr&&r.retainedResult==nullptr);}
    {MidpointProvider p;p.throwStampCall=static_cast<int>(total);auto r=run(p);
     CHECK(r.failure.reason==algo::T4BReason::PlacementProviderThrew);
     CHECK(r.sourceReceipt!=nullptr&&r.retainedResult!=nullptr);
     CHECK(r.retainedResult&&r.retainedResult->certificationComplete);
     CHECK(r.retainedResult&&!r.retainedResult->publicationReady);}

    auto initial=[&](bool foreign) {
        MidpointProvider p;auto request=routeRequest(ownership,p);
        if(foreign)++p.current.persistentId;else ++p.current.revision;
        mesh::Mesh destination;const auto before=mesh::meshHash(destination);
        auto r=algo::meshModelCurvedTet10(*box,spec,request,destination,ownership,p);
        CHECK(mesh::meshHash(destination)==before);
        CHECK(r.sourceOutcome==nullptr&&r.sourceReceipt==nullptr);
        return r;
    };
    CHECK(initial(true).failure.reason==algo::T4BReason::ForeignPlacementProvider);
    CHECK(initial(false).failure.reason==algo::T4BReason::StalePlacementProvider);
}

void certificationFailureRetainsCandidate() {
    TEST("t4b.route.certificate_retention");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    placement.distortAll=true;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    mesh::Mesh destination;destination.addNode({43,47,53},{0,89});
    const auto digest=mesh::meshHash(destination);
    const auto result=algo::meshModelCurvedTet10(*box,spec,
        routeRequest(ownership,placement),destination,ownership,placement);
    CHECK(result.status==algo::T4BStatus::Failed ||
          result.status==algo::T4BStatus::Unsupported ||
          result.status==algo::T4BStatus::Incomplete);
    CHECK(result.sourceReceipt!=nullptr&&result.sourceOutcome==nullptr);
    CHECK(result.retainedResult!=nullptr);
    CHECK(result.retainedResult&&!result.retainedResult->certificationComplete);
    CHECK(result.retainedResult&&!result.retainedResult->publicationReady);
    CHECK(result.reportComplete);
    requireReportXor(result);
    CHECK(mesh::meshHash(destination)==digest);
}

void hardCapacityBudgets() {
    TEST("t4b.route.hard_capacity_budgets");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    auto request=routeRequest(ownership,placement);request.policy.maxResultNodes=0;
    request.policy.maxResultCells=0;request.policy.maxResultConnectivity=0;
    mesh::Mesh destination;destination.addNode({31,37,41},{0,88});
    const auto digest=mesh::meshHash(destination);
    const auto result=algo::meshModelCurvedTet10(*box,spec,request,
        destination,ownership,placement);
    CHECK(result.status==algo::T4BStatus::Incomplete);
    CHECK(result.failure.reason==algo::T4BReason::ResultBudgetReached);
    CHECK(result.sourceOutcome!=nullptr);
    CHECK(result.sourceReceipt==nullptr&&result.retainedResult==nullptr);
    CHECK(mesh::meshHash(destination)==digest);
    CHECK(placement.calls==0);
}

void requestPreflightAndZeroWork() {
    TEST("t4b.route.request_preflight_zero_work");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    for(int discriminator=0;discriminator<3;++discriminator) {
        OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
        auto request=routeRequest(ownership,placement);
        if(discriminator==0)request.policy.minimumCertifiedJacobianRatio=
            std::numeric_limits<double>::quiet_NaN();
        else if(discriminator==1)request.policy.minimumCertifiedJacobianRatio=1.01;
        else request.policy.maxCertificateDepth=21;
        mesh::Mesh destination;destination.addNode({71,73,79},{0,97});
        const auto before=mesh::meshHash(destination);
        const auto result=algo::meshModelCurvedTet10(
            *box,spec,request,destination,ownership,placement);
        CHECK(result.status==algo::T4BStatus::Rejected&&
              result.failure.reason==algo::T4BReason::InvalidRequest);
        CHECK(placement.stampCalls==0&&placement.calls==0);
        CHECK(ownership.stampCalls==0&&ownership.ownershipCalls==0&&
              ownership.envelopeCalls==0);
        CHECK(mesh::meshHash(destination)==before);
    }
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    auto request=routeRequest(ownership,placement);request.policy.maxWorkUnits=0;
    mesh::Mesh destination;const auto result=algo::meshModelCurvedTet10(
        *box,spec,request,destination,ownership,placement);
    CHECK(result.status==algo::T4BStatus::Incomplete&&
          result.failure.reason==algo::T4BReason::JacobianCertificationBudgetReached);
    CHECK(result.sourceOutcome==nullptr&&result.sourceReceipt==nullptr);
    CHECK(result.reportComplete);requireReportXor(result);
    CHECK(placement.stampCalls==0&&placement.calls==0);

    OwnershipProvider boundaryOwnership;configure(*box,boundaryOwnership);
    MidpointProvider boundaryPlacement;
    auto boundaryRequest=routeRequest(boundaryOwnership,boundaryPlacement);
    boundaryRequest.policy.maxWorkUnits=1;
    mesh::Mesh boundaryDestination;
    const auto boundaryResult=algo::meshModelCurvedTet10(
        *box,spec,boundaryRequest,boundaryDestination,boundaryOwnership,
        boundaryPlacement);
    CHECK(boundaryResult.status==algo::T4BStatus::Incomplete&&
          boundaryResult.failure.reason==
              algo::T4BReason::JacobianCertificationBudgetReached);
    CHECK(boundaryResult.sourceOutcome!=nullptr&&
          boundaryResult.sourceReceipt==nullptr);
    CHECK(boundaryPlacement.stampCalls==1&&boundaryPlacement.calls==0);
}

void boundaryAuthorityTaxonomy() {
    TEST("t4b.source.boundary_authority_taxonomy");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    const auto t4Request=routeRequest(ownership,placement).source;
    mesh::Mesh source;auto t4=algo::meshModelQualityTet4(
        *box,spec,t4Request,source,ownership);
    CHECK(t4.published()&&t4.retainedResult&&t4.retainedResult->artifact);
    if(!t4.retainedResult||!t4.retainedResult->artifact)return;
    auto& artifact=*t4.retainedResult->artifact;
    const auto expected=t4Request.source.conformingRequest.expectedSurfaceSetStamp;
    auto validate=[&](const algo::ArtifactStamp& stamp=algo::ArtifactStamp{}) {
        return algo::t4b_detail::validateSourceAuthorityForTesting(
            source,artifact,stamp.persistentId==0?expected:stamp);
    };
    CHECK(validate()==algo::T4BReason::None);
    CHECK(validate({expected.persistentId+1,expected.revision})==
          algo::T4BReason::ForeignEdgeBinding);
    CHECK(validate({expected.persistentId,expected.revision+1})==
          algo::T4BReason::StaleEdgeBinding);
    auto& state=artifact.constraintState;
    const auto edgeSchedule=state.canonicalEdgeSchedule;
    state.canonicalEdgeSchedule.push_back(state.canonicalEdgeSchedule.front());
    std::sort(state.canonicalEdgeSchedule.begin(),state.canonicalEdgeSchedule.end());
    CHECK(validate()==algo::T4BReason::DuplicateEdgeBinding);
    state.canonicalEdgeSchedule=edgeSchedule;
    const auto edgeKey=state.canonicalEdgeSchedule.front();
    auto edgeBinding=state.edgeBindings.at(edgeKey);
    state.edgeBindings.erase(edgeKey);
    CHECK(validate()==algo::T4BReason::MissingEdgeBinding);
    state.edgeBindings.emplace(edgeKey,edgeBinding);
    auto& restored=state.edgeBindings.at(edgeKey);
    const auto originalClass=restored.classification;
    restored.classification.tag+=1000;
    CHECK(validate()==algo::T4BReason::ConflictingEdgeClassification);
    restored.classification=originalClass;
    const auto originalHandle=restored.edge;
    ++restored.edge.gen;
    CHECK(validate()==algo::T4BReason::StaleEdgeBinding);
    restored.edge=originalHandle;
    restored.required=false;
    CHECK(validate()==algo::T4BReason::ForeignEdgeBinding);
    restored.required=true;
    CHECK(validate()==algo::T4BReason::None);
}

struct ReloadCell {
    mesh::CellType type=mesh::CellType::Tet10;
    std::vector<Index> nodes;
    mesh::Classification cls{};
    friend bool operator==(const ReloadCell& a,const ReloadCell& b) noexcept {
        return a.type==b.type&&a.nodes==b.nodes&&a.cls.dim==b.cls.dim&&
               a.cls.tag==b.cls.tag;
    }
};
struct ReloadMesh {
    std::vector<Vec3> nodes;
    std::vector<mesh::Classification> nodeClasses;
    std::vector<ReloadCell> cells;
};

std::optional<mesh::CellType> vtkHighOrderType(int value) {
    if(value==21)return mesh::CellType::Line3;
    if(value==22)return mesh::CellType::Tri6;
    if(value==24)return mesh::CellType::Tet10;
    return std::nullopt;
}
std::optional<mesh::CellType> mshHighOrderType(int value) {
    if(value==8)return mesh::CellType::Line3;
    if(value==9)return mesh::CellType::Tri6;
    if(value==11)return mesh::CellType::Tet10;
    return std::nullopt;
}
int cellDimension(mesh::CellType type) noexcept {
    if(type==mesh::CellType::Line3)return 1;
    if(type==mesh::CellType::Tri6)return 2;
    return 3;
}

std::optional<ReloadMesh> reloadHighOrderVtk(const std::string& path) {
    std::ifstream input(path);if(!input)return std::nullopt;
    ReloadMesh out;std::string token,typeName;
    while(input>>token&&token!="POINTS"){}
    std::size_t nodeCount=0;if(token!="POINTS"||!(input>>nodeCount>>typeName))
        return std::nullopt;
    out.nodes.resize(nodeCount);out.nodeClasses.assign(nodeCount,{});
    for(auto& p:out.nodes)if(!(input>>p.x>>p.y>>p.z))return std::nullopt;
    std::size_t cellCount=0,listCount=0;
    if(!(input>>token>>cellCount>>listCount)||token!="CELLS")return std::nullopt;
    out.cells.resize(cellCount);std::size_t observed=0;
    for(auto& cell:out.cells){std::size_t n=0;if(!(input>>n))return std::nullopt;
        observed+=n+1;cell.nodes.resize(n);
        for(auto& node:cell.nodes)if(!(input>>node)||node<0||
            static_cast<std::size_t>(node)>=nodeCount)return std::nullopt;}
    if(observed!=listCount)return std::nullopt;
    std::size_t typeCount=0;
    if(!(input>>token>>typeCount)||token!="CELL_TYPES"||typeCount!=cellCount)
        return std::nullopt;
    for(auto& cell:out.cells){int code=0;if(!(input>>code))return std::nullopt;
        const auto type=vtkHighOrderType(code);if(!type||
            cell.nodes.size()!=static_cast<std::size_t>(mesh::nodesPerCell(*type)))
            return std::nullopt;
        cell.type=*type;cell.cls.dim=cellDimension(*type);}
    std::size_t dataCount=0;std::string name,integerType,components,lookup,table;
    if(!(input>>token>>dataCount)||token!="CELL_DATA"||dataCount!=cellCount||
       !(input>>token>>name>>integerType>>components)||token!="SCALARS"||
       name!="entity_tag"||!(input>>lookup>>table)||lookup!="LOOKUP_TABLE")
        return std::nullopt;
    for(auto& cell:out.cells)if(!(input>>cell.cls.tag))return std::nullopt;
    return out;
}

std::optional<ReloadMesh> reloadHighOrderMsh(const std::string& path) {
    std::ifstream input(path);if(!input)return std::nullopt;
    ReloadMesh out;std::string token;
    while(input>>token&&token!="$Nodes"){}
    std::size_t blocks=0,nodeCount=0;std::uint64_t minTag=0,maxTag=0;
    if(token!="$Nodes"||!(input>>blocks>>nodeCount>>minTag>>maxTag)||
       minTag!=1||maxTag!=nodeCount)return std::nullopt;
    out.nodes.resize(nodeCount);out.nodeClasses.resize(nodeCount);
    std::vector<unsigned char> seen(nodeCount,0);
    for(std::size_t block=0;block<blocks;++block){int dim=0,tag=0,parametric=0;
        std::size_t count=0;if(!(input>>dim>>tag>>parametric>>count)||parametric!=0)
            return std::nullopt;
        std::vector<std::uint64_t> ids(count);
        for(auto& id:ids)if(!(input>>id)||id==0||id>nodeCount)return std::nullopt;
        for(const auto id:ids){const auto index=static_cast<std::size_t>(id-1);
            if (seen[index]) return std::nullopt;
            seen[index]=1;
            out.nodeClasses[index]={dim,tag};auto& p=out.nodes[index];
            if(!(input>>p.x>>p.y>>p.z))return std::nullopt;}}
    if(!std::all_of(seen.begin(),seen.end(),[](unsigned char v){return v!=0;})||
       !(input>>token)||token!="$EndNodes")return std::nullopt;
    while(input>>token&&token!="$Elements"){}
    std::size_t elementBlocks=0,elementCount=0;
    std::uint64_t firstElement=0,lastElement=0;
    if(token!="$Elements"||!(input>>elementBlocks>>elementCount>>firstElement>>lastElement)||
       (elementCount&&(firstElement!=1||lastElement!=elementCount)))return std::nullopt;
    out.cells.reserve(elementCount);std::set<std::uint64_t> elementIds;
    for(std::size_t block=0;block<elementBlocks;++block){int dim=0,tag=0,code=0;
        std::size_t count=0;if(!(input>>dim>>tag>>code>>count))return std::nullopt;
        const auto type=mshHighOrderType(code);if(!type||dim!=cellDimension(*type))
            return std::nullopt;
        for(std::size_t row=0;row<count;++row){std::uint64_t element=0;
            if(!(input>>element)||!elementIds.insert(element).second)return std::nullopt;
            ReloadCell cell;cell.type=*type;cell.cls={dim,tag};
            cell.nodes.resize(static_cast<std::size_t>(mesh::nodesPerCell(*type)));
            for(auto& node:cell.nodes){std::uint64_t id=0;
                if(!(input>>id)||id==0||id>nodeCount)return std::nullopt;
                node=static_cast<Index>(id-1);}out.cells.push_back(std::move(cell));}}
    if(out.cells.size()!=elementCount||!(input>>token)||token!="$EndElements")
        return std::nullopt;
    return out;
}

double exportedValue(double value) {
    std::ostringstream text;text<<value;std::istringstream input(text.str());
    double result=0.0;input>>result;return result;
}
ReloadMesh expectedReload(const mesh::Mesh& source,bool volumeOnly,
                          bool classes,bool fullPrecision=false) {
    ReloadMesh out;out.nodes.reserve(static_cast<std::size_t>(source.nodeCount()));
    out.nodeClasses.reserve(static_cast<std::size_t>(source.nodeCount()));
    for(Index i=0;i<source.nodeCount();++i){const auto p=source.node(i);
        out.nodes.push_back(fullPrecision?p:Vec3{
            exportedValue(p.x),exportedValue(p.y),exportedValue(p.z)});
        out.nodeClasses.push_back(classes?source.nodeClass(i):mesh::Classification{});}
    for(Index i=0;i<source.cellCount();++i){const auto c=source.cell(i);
        if(volumeOnly&&c.type!=mesh::CellType::Tet10)continue;
        ReloadCell row;row.type=c.type;row.cls=c.onEntity;
        row.nodes.assign(source.cellNodes(i),source.cellNodes(i)+mesh::nodesPerCell(c.type));
        out.cells.push_back(std::move(row));}return out;
}
bool sameReload(const ReloadMesh& a,const ReloadMesh& b) noexcept {
    if(a.nodes.size()!=b.nodes.size()||a.nodeClasses.size()!=b.nodeClasses.size()||
       a.cells!=b.cells)return false;
    for(std::size_t i=0;i<a.nodes.size();++i)
        if(a.nodes[i].x!=b.nodes[i].x||a.nodes[i].y!=b.nodes[i].y||
           a.nodes[i].z!=b.nodes[i].z||
           a.nodeClasses[i].dim!=b.nodeClasses[i].dim||
           a.nodeClasses[i].tag!=b.nodeClasses[i].tag)return false;
    return true;
}
mesh::Mesh reloadAsMesh(const ReloadMesh& value) {
    mesh::Mesh out;std::size_t connectivity=0;
    for(const auto& cell:value.cells)connectivity+=cell.nodes.size();
    out.reserveAppend(value.nodes.size(),value.cells.size(),connectivity);
    for(std::size_t i=0;i<value.nodes.size();++i)
        out.addNode(value.nodes[i],value.nodeClasses[i]);
    for(const auto& cell:value.cells)out.addCell(cell.type,cell.nodes,cell.cls);
    return out;
}
std::set<std::array<Index,3>> tet10Boundary(const ReloadMesh& value) {
    static constexpr std::array<std::array<unsigned,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    std::map<std::array<Index,3>,unsigned> incidence;
    for(const auto& cell:value.cells)if(cell.type==mesh::CellType::Tet10)
        for(const auto& face:faces){std::array<Index,3> key{
            cell.nodes[face[0]],cell.nodes[face[1]],cell.nodes[face[2]]};
            std::sort(key.begin(),key.end());++incidence[key];}
    std::set<std::array<Index,3>> out;
    for (const auto& [key,count]:incidence)
        if (count==1) out.insert(key);
    return out;
}
std::set<std::array<Index,3>> tri6Boundary(const ReloadMesh& value) {
    std::set<std::array<Index,3>> out;
    for (const auto& cell:value.cells) if (cell.type==mesh::CellType::Tri6) {
        std::array<Index,3> key{cell.nodes[0],cell.nodes[1],cell.nodes[2]};
        std::sort(key.begin(),key.end());
        out.insert(key);
    }
    return out;
}

void auditedHighOrderExport() {
    TEST("t4b.export.reload_vtk_msh");
    auto kernel=geo::makeBuiltinKernel();auto box=kernel->makeBox({0,0,0},{1,1,1});
    OwnershipProvider ownership;configure(*box,ownership);MidpointProvider placement;
    sizing::MeshSpec spec;spec.field=std::make_shared<sizing::ConstantField>(.3);
    mesh::Mesh output;const auto result=algo::meshModelCurvedTet10(
        *box,spec,routeRequest(ownership,placement),output,ownership,placement);
    CHECK(result.published());
    const auto base=std::filesystem::temp_directory_path()/
        ("paacs_t4b_audited_"+std::to_string(::getpid()));
    const auto vtk=base.string()+".vtk",msh=base.string()+".msh";
    io::writeVtk(output,vtk,true);io::writeMsh(output,msh);
    const auto vtkRead=reloadHighOrderVtk(vtk),mshRead=reloadHighOrderMsh(msh);
    CHECK(vtkRead.has_value());CHECK(mshRead.has_value());
    const auto vtkExpected=expectedReload(output,true,false,true);
    const auto mshExpected=expectedReload(output,false,true);
    if(vtkRead){CHECK(sameReload(*vtkRead,vtkExpected));
        CHECK(tet10Boundary(*vtkRead)==tri6Boundary(mshExpected));
        CHECK(mesh::meshHash(reloadAsMesh(*vtkRead))==
              mesh::meshHash(reloadAsMesh(vtkExpected)));}
    if(mshRead){CHECK(sameReload(*mshRead,mshExpected));
        CHECK(tet10Boundary(*mshRead)==tri6Boundary(*mshRead));
        const auto actual=reloadAsMesh(*mshRead),expected=reloadAsMesh(mshExpected);
        CHECK(mesh::meshHash(actual)==mesh::meshHash(expected));
        const auto actualSemantic=mesh::semanticHashEx(actual);
        const auto expectedSemantic=mesh::semanticHashEx(expected);
        CHECK(actualSemantic.canonical&&expectedSemantic.canonical&&
              actualSemantic.value==expectedSemantic.value);}
    std::filesystem::remove(vtk);std::filesystem::remove(msh);
}

} // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* t4bTests(std::size_t& count) {
    static const Entry tests[]={{"t4b.certificate.straight_golden",straightCertificateGolden},
        {"t4b.certificate.blossom_children",blossomConstantOracle},
        {"t4b.certificate.blossom_exact_rational",blossomExactRationalOracle},
        {"t4b.certificate.initial_exact_rational",initialCoefficientExactOracle},
        {"t4b.certificate.taxonomy",taxonomyAndBudgets},
        {"t4b.certificate.hard_work_budget",hardCertificateWorkBudget},
        {"t4b.certificate.curved_taxonomy",curvedTaxonomy},
        {"t4b.certificate.extreme_frames",extremeCertificateFrames},
        {"t4b.route.straight_box",straightPublicRoute},
        {"t4b.route.source_cell_ordinals",sourceCellOrdinalAttribution},
        {"t4b.route.provider_endpoint_coordinate_alignment",providerEndpointCoordinateAlignment},
        {"t4b.route.curved_cylinder",curvedCylinderRoute},
        {"t4b.route.determinism",routeDeterminism},
        {"t4b.route.public_scaling",publicRouteScaling},
        {"t4b.route.provider_failures",providerFailureContracts},
        {"t4b.route.provider_stamp_checkpoints",providerStampCheckpoints},
        {"t4b.route.certificate_retention",certificationFailureRetainsCandidate},
        {"t4b.route.hard_capacity_budgets",hardCapacityBudgets},
        {"t4b.route.request_preflight_zero_work",requestPreflightAndZeroWork},
        {"t4b.source.boundary_authority_taxonomy",boundaryAuthorityTaxonomy},
        {"t4b.export.reload_vtk_msh",auditedHighOrderExport}};
    count=sizeof(tests)/sizeof(tests[0]);return tests;
}
} // namespace paacs_tests
