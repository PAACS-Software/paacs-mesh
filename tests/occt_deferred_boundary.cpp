// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctBoundaryInternal.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <atomic>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace paacs;
using namespace paacs::geo::detail;
namespace {
void check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(model),s);return s.str();}
}
int main(){try{
    auto box=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2,3,4).Shape());check(box.imported(),"native box fixture");
    const auto boxBefore=snapshot(*box.model);geo::OcctBoundaryOptions options;options.linearDeflection=.05;options.maximumChordError=.01;
    mesh::Mesh direct,candidate;const auto full=geo::buildOcctBoundary(*box.model,1,options,direct);
    const auto pending=buildOcctBoundaryCandidate(*box.model,1,options,candidate);
    check(full.published&&pending.published&&pending.closed&&pending.oriented&&pending.vertexManifold&&pending.sourceModelUnmodified,"native topology/node/source checks still mandatory while staging");
    check(mesh::meshHash(direct)==mesh::meshHash(candidate)&&full.nodeKeys==pending.nodeKeys,"candidate changes neither mesh nor native classifications/keys");
    const auto candidateHash=mesh::meshHash(candidate);
    const auto audit=auditOcctBoundaryFacets(*box.model,1,candidate,options);
    check(audit.completed&&audit.sourceModelUnmodified&&audit.reason==geo::OcctBoundaryReason::None,"deferred complete physical facet/source audit");
    check(audit.samples==4*candidate.cellCount()&&audit.maximumFacetSampleError==full.maximumFacetSampleError,"all four original CAD samples checked with original error threshold");
    auto capped=options;capped.maximumWork=1;const auto limited=auditOcctBoundaryFacets(*box.model,1,candidate,capped);
    check(!limited.completed&&limited.reason==geo::OcctBoundaryReason::BudgetExceeded&&limited.work==1,"shared-work refusal before incomplete facet acceptance");
    const auto tiny=auditOcctBoundaryFacets(*box.model,1,candidate,options,1);
    check(!tiny.completed&&tiny.reason==geo::OcctBoundaryReason::BudgetExceeded,"encoded source snapshot cap enforced");
    std::atomic_bool cancellation{true};capped=options;capped.cancellation=&cancellation;const auto cancelled=auditOcctBoundaryFacets(*box.model,1,candidate,capped);
    check(!cancelled.completed&&cancelled.reason==geo::OcctBoundaryReason::Cancelled&&cancelled.work==0,"cancelled audit cannot complete");
    check(mesh::meshHash(candidate)==candidateHash&&snapshot(*box.model)==boxBefore,"deferred audits are read-only on success and resource failure");
    auto sphere=OcctModelAccess::fromShape(BRepPrimAPI_MakeSphere(2).Shape());check(sphere.imported(),"native sphere fixture");
    const auto sphereBefore=snapshot(*sphere.model);options.linearDeflection=.1;options.maximumChordError=1e-6;
    mesh::Mesh sentinel;sentinel.addNode({-1,-2,-3},{3,99});const auto sentinelHash=mesh::meshHash(sentinel);
    const auto refused=geo::buildOcctBoundary(*sphere.model,1,options,sentinel);
    check(!refused.published&&refused.reason==geo::OcctBoundaryReason::GeometryError&&mesh::meshHash(sentinel)==sentinelHash,"public boundary still immediately refuses excessive spherical chord error");
    mesh::Mesh sphereCandidate;const auto spherePending=buildOcctBoundaryCandidate(*sphere.model,1,options,sphereCandidate);
    check(spherePending.published&&spherePending.closed&&spherePending.vertexManifold,"coarse spherical boundary may be privately staged before facet feasibility");
    const auto sphereCandidateHash=mesh::meshHash(sphereCandidate);
    const auto rejected=auditOcctBoundaryFacets(*sphere.model,1,sphereCandidate,options);
    check(!rejected.completed&&rejected.reason==geo::OcctBoundaryReason::GeometryError&&rejected.samples>0&&rejected.maximumFacetSampleError>options.maximumChordError,"postponed facet check still rejects actual measured CAD error");
    check(mesh::meshHash(sphereCandidate)==sphereCandidateHash&&snapshot(*sphere.model)==sphereBefore,"geometry refusal preserves source and staged payload");
    std::cout<<"deferred native candidate: identical payload/native topology; every accepted facet sample/source guard; physical chord/resource/cancel refusal passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
