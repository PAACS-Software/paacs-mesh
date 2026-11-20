// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// T4-B route/publication allocation qualification. Excluded from sanitizer
// builds because it replaces global operator new.
#include "paacs/algo/Tet10CurvingMesher.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sizing/SizeField.hpp"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string_view>

static long gOrdinal=0,gFail=0,gTailAllocations=0;
static bool gTrack=false,gInjected=false,gPublicationOnly=false;

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif
void* operator new(std::size_t bytes) {
    if(paacs::algo::t4b_detail::publicationTailActive())++gTailAllocations;
    const bool selected=gPublicationOnly?
        paacs::algo::t4b_detail::publicationPreparationActive():
        paacs::algo::t4b_detail::routeActive();
    if(gTrack&&selected) {
        ++gOrdinal;
        if(gFail>0&&gOrdinal==gFail){gInjected=true;throw std::bad_alloc{};}
    }
    if(void* p=std::malloc(bytes?bytes:1))return p;
    throw std::bad_alloc{};
}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p) noexcept{std::free(p);}
void operator delete[](void* p) noexcept{std::free(p);}
void operator delete(void* p,std::size_t) noexcept{std::free(p);}
void operator delete[](void* p,std::size_t) noexcept{std::free(p);}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace {
using namespace paacs;
int failures=0;
#define REQUIRE(x) do{if(!(x)){std::printf("FAIL %d: %s\n",__LINE__,#x);++failures;}}while(0)

class Ownership final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{8701,1};
    std::map<algo::FacePersistentId,std::vector<bnd::VolumeRegionUse>> rows;
    algo::VolumeBoundaryProviderStamp stamp() const noexcept override{return current;}
    algo::VolumeFaceOwnership ownership(algo::FacePersistentId face) const override{
        algo::VolumeFaceOwnership out;out.providerStamp=current;out.facePersistentId=face;
        const auto i=rows.find(face);if(i!=rows.end())out.owners=i->second;return out;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(std::span<const Vec3>) const override{
        algo::ProviderRecoveryEnvelope out;out.providerStamp=current;
        out.availability=algo::VolumeProviderAvailability::RecoveryEnvelopeUnavailable;return out;
    }
};
class Placement final : public algo::Tet10PlacementProvider {
public:
    algo::Tet10PlacementStamp current{8801,1};
    algo::Tet10PlacementStamp stamp() const override{return current;}
    algo::Tet10PlacementResult place(const algo::Tet10PlacementRequest& r) const override{
        return {algo::Tet10PlacementStatus::Placed,0,current,
            {(r.endpointCoordinates[0].x+r.endpointCoordinates[1].x)*.5,
             (r.endpointCoordinates[0].y+r.endpointCoordinates[1].y)*.5,
             (r.endpointCoordinates[0].z+r.endpointCoordinates[1].z)*.5}};
    }
};
void configure(const geo::Model& model,Ownership& provider){
    Vec3 lo{},hi{};model.solids().front()->bbox(lo,hi);const Vec3 c=(lo+hi)*.5;
    for(const auto* face:model.faces()){
        Vec2 a{},b{};face->paramBounds(a,b);const Vec2 uv=(a+b)*.5;
        const int sense=face->normalAt(uv).dot(face->pointAt(uv)-c)>0?1:-1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1]={{0,sense}};
    }
}
algo::T4BRequest requestFor(const Ownership& o,const Placement& p){
    algo::T4BRequest r;
    r.source.source.conformingRequest.expectedSurfaceSetStamp={87001,1};
    r.source.source.conformingRequest.expectedOwnershipProviderStamp=o.current;
    r.source.source.conformingRequest.envelopeMode=algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    auto& x=r.source.source.refinement;x.maxEdgeToTarget=1e9;x.maxRadiusEdge=1e9;
    x.maxInsertions=0;x.maxWorkUnits=200000;x.maxOutputNodes=100000;x.maxOutputTets=100000;
    x.maxReportRows=5;x.maxReportWitnesses=100000;x.maxReportTopologyHandles=300000;
    r.source.quality.mode=algo::T4AQualityMode::ReportOnly;
    r.source.quality.minimumRadiusRatio=.2;r.source.quality.maxResultNodes=100000;
    r.source.quality.maxResultTets=100000;r.source.quality.maxReportRows=1000;
    r.source.quality.maxWorstTailRows=64;
    r.expectedPlacementStamp=p.current;
    r.policy={.001,4,4096,100000,10000000,100000,100000,1000000,
              100000,100000,200000};
    return r;
}
struct Fixture{
    std::unique_ptr<geo::Kernel> kernel;std::unique_ptr<geo::Model> model;
    Ownership ownership;Placement placement;sizing::MeshSpec spec;algo::T4BRequest request;
    Fixture(){kernel=geo::makeBuiltinKernel();model=kernel->makeBox({0,0,0},{1,1,1});
        configure(*model,ownership);spec.field=std::make_shared<sizing::ConstantField>(.3);
        spec.threads=1;request=requestFor(ownership,placement);}
};
algo::T4BOutcome run(Fixture& f,mesh::Mesh& destination,long fail,bool publication){
    gOrdinal=0;gFail=fail;gInjected=false;gPublicationOnly=publication;gTrack=true;
    auto out=algo::meshModelCurvedTet10(*f.model,f.spec,f.request,destination,
                                       f.ownership,f.placement);
    gTrack=false;gFail=0;return out;
}
struct Signature{std::uint64_t hash=0,edges=0,cells=0;
    friend bool operator==(const Signature&,const Signature&)=default;};
Signature signature(const algo::T4BOutcome& out){
    Signature s;if(out.retainedResult){s.hash=out.retainedResult->strictHash;
        s.edges=out.retainedResult->edgeToMidside.size();
        s.cells=out.retainedResult->certificates.size();}return s;
}
long measure(bool publication,Signature& reference){
    {Fixture warm;mesh::Mesh m;auto w=algo::meshModelCurvedTet10(*warm.model,warm.spec,
        warm.request,m,warm.ownership,warm.placement);REQUIRE(w.published());}
    Fixture f;mesh::Mesh m;auto out=run(f,m,0,publication);
    REQUIRE(out.published());REQUIRE(out.retainedResult!=nullptr);reference=signature(out);
    REQUIRE(gTailAllocations==0);return gOrdinal;
}
void sweep(bool publication,long first,long last,long expected){
    Signature reference;const long total=measure(publication,reference);
    REQUIRE(expected<=0||expected==total);if(last<=0||last>total)last=total;
    long faults=0;
    for(long ordinal=first;ordinal<=last;++ordinal){
        Fixture f;mesh::Mesh destination;destination.addNode({-71,-73,-79},{0,991});
        const auto before=mesh::meshHash(destination);auto out=run(f,destination,ordinal,publication);
        REQUIRE(gInjected);REQUIRE(!out.published());REQUIRE(out.destinationUnchanged);
        REQUIRE(mesh::meshHash(destination)==before);
        if(publication){REQUIRE(out.failure.reason==algo::T4BReason::AllocationFailed);
            REQUIRE(out.failure.phase==algo::T4BPhase::PublicationPreparation);
            REQUIRE(out.sourceReceipt!=nullptr&&out.retainedResult!=nullptr);
            REQUIRE(out.retainedResult&&out.retainedResult->certificationComplete);}
        Fixture retry;auto clean=run(retry,destination,0,publication);
        REQUIRE(clean.published());REQUIRE(signature(clean)==reference);++faults;
    }
    std::printf("probe_t4b %s range=%ld-%ld total=%ld faults=%ld tail=%ld\n",
        publication?"publication":"route",first,last,total,faults,gTailAllocations);
}
long number(const char* text){long out=0;const char* end=text+std::strlen(text);
    const auto p=std::from_chars(text,end,out);if(p.ec!=std::errc{}||p.ptr!=end||out<0)
        throw std::invalid_argument("ordinal");return out;}
bool publicationMode(std::string_view value){
    if(value=="publication")return true;
    if(value=="route")return false;
    throw std::invalid_argument("mode");
}
}

int main(int argc,char** argv){
    try{
        if(argc==3&&std::string_view(argv[1])=="--count"){
            const bool publication=publicationMode(argv[2]);
            Signature s;std::printf("t4b %s allocation count: %ld\n",argv[2],measure(publication,s));
        }else if(argc==4&&std::string_view(argv[1])=="--expect-count"){
            const bool publication=publicationMode(argv[2]);
            Signature s;const long observed=measure(publication,s);
            REQUIRE(observed==number(argv[3]));
            std::printf("t4b %s allocation count: %ld expected=%s\n",
                        argv[2],observed,argv[3]);
        }else if(argc==5&&std::string_view(argv[1])=="--range"){
            const bool publication=publicationMode(argv[2]);
            sweep(publication,number(argv[3]),number(argv[4]),0);
        }else if(argc==6&&std::string_view(argv[1])=="--expect-range"){
            const bool publication=publicationMode(argv[2]);
            sweep(publication,number(argv[3]),number(argv[4]),number(argv[5]));
        }else if(argc==1){Signature s;const long total=measure(false,s);sweep(false,1,std::min(16L,total),total);}
        else return 2;
    }catch(const std::exception& e){gTrack=false;std::printf("probe_t4b exception: %s\n",e.what());++failures;}
    REQUIRE(gTailAllocations==0);
    if(failures){std::printf("probe_t4b FAILED failures=%d\n",failures);return 1;}
    std::puts("probe_t4b PASSED");return 0;
}
