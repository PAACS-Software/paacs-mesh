// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctInteriorSites.hpp"
#include "../src/geometry/OcctBoundaryInternal.hpp"
#include "../src/geometry/OcctFrontMesher.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include "../src/geometry/OcctSolidClassifier.hpp"
#include "../src/geometry/OcctSpatialSizing.hpp"
#include "../src/geometry/OcctStarMesher.hpp"
#include "../src/geometry/OcctStageTrace.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

using namespace paacs;
using namespace paacs::geo::detail;

namespace {
void check(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}

bool exactPoint(Vec3 a,Vec3 b) {
    return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x) &&
           std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y) &&
           std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);
}

void sameBoundary(const mesh::Mesh& before,const mesh::Mesh& after) {
    check(after.nodeCount()>=before.nodeCount() && after.cellCount()==before.cellCount(),
          "interior sites changed the original cell count or removed nodes");
    for (Index i=0;i<before.nodeCount();++i)
        check(exactPoint(before.node(i),after.node(i)) && before.nodeClass(i)==after.nodeClass(i),
              "interior sites changed original binary64 coordinates or CAD node identity");
    for (Index i=0;i<before.cellCount();++i) {
        const auto a=before.cell(i),b=after.cell(i);
        check(a.type==b.type && a.onEntity==b.onEntity &&
              std::equal(before.cellNodes(i),before.cellNodes(i)+mesh::nodesPerCell(a.type),
                         after.cellNodes(i)),
              "interior sites changed an original cell, winding or CAD face identity");
    }
}

void sameMesh(const mesh::Mesh& before,const mesh::Mesh& after) {
    check(before.nodeCount()==after.nodeCount(),"deterministic replay changed the node count");
    sameBoundary(before,after);
}

std::string snapshot(const geo::OcctModel& model) {
    std::ostringstream stream;
    BRepTools::Write(OcctModelAccess::shape(model),stream);
    return stream.str();
}

bool inControl(Vec3 p) {
    return p.x>-.3 && p.x<.3 && p.y>-.3 && p.y<.3 && p.z>.7 && p.z<1.3;
}

void exerciseNativePlacement() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    check(imported.imported(),"native cylinder fixture import failed");
    const auto& model=*imported.model;
    const auto originalCad=snapshot(model);
    int wall=-1;
    for (const auto* face:model.faces()) {
        const auto* descriptor=model.surfaceDescriptor(face->tag());
        if (descriptor && descriptor->kind==geo::SurfaceKind::Cylinder) wall=face->tag();
    }
    check(wall>0,"native cylinder fixture needs a cylindrical wall");
    geo::NativeCadMeshOptions options;
    options.geometry.targetSize=.7;
    options.minimumSize=.12;
    options.faceSizes.push_back({wall,.16,1.55});
    sizing::RegionSizeSpecification specification;
    specification.baseSize=.7;
    specification.minimumSize=.12;
    specification.controls.push_back({sizing::SizeRegionKind::Box,
        {-.3,-.3,.7},{.3,.3,1.3},1.,.15,1.55});
    auto compiled=sizing::compileRegionSizeField(specification);
    check(compiled.compiled,"native cylinder sizing fixture failed");
    const OcctSpatialSizing field(model,options,std::move(compiled.field));
    geo::OcctBoundaryOptions boundaryOptions;
    boundaryOptions.maximumChordError=.009;
    boundaryOptions.maximumNodes=6000;
    boundaryOptions.maximumTriangles=10000;
    mesh::Mesh boundary;
    OcctBoundaryAliasReceipt receipt;
    const auto built=buildSpacedOcctBoundary(model,1,boundaryOptions,.7,
        64ULL*1024*1024,boundary,receipt,true,false,0,&field);
    check(built.published && built.sourceModelUnmodified && built.closed && built.oriented &&
          certifyOcctBoundarySurface(boundary,OcctStarOptions{}).certified,
          "native fixture needs a certified oriented closed CAD boundary");
    const auto unchanged=boundary;
    for (const auto p:boundary.nodes())
        check(!inControl(p),"enclosed control unexpectedly intersects an original boundary vertex");
    constexpr mesh::Classification region{3,1};
    constexpr std::size_t cap=96;
    OcctSolidClassifier classifier(OcctModelAccess::solid(model,region.tag));
    std::size_t nativeQueries=0;
    const auto strictInterior=[&](Vec3 p) {
        ++nativeQueries;
        classifier.Perform(gp_Pnt(p.x,p.y,p.z),model.tolerance(3,region.tag));
        return classifier.State()==TopAbs_IN;
    };
    std::uint64_t work=0;
    const auto charge=[&](std::uint64_t n) {
        check(n<=2'000'000-work,"interior placement exceeded the unit work bound");
        work+=n;
    };
    OcctInteriorSiteReport firstReport;
    const auto sites=appendOcctInteriorSites(boundary,field,region,cap,
        strictInterior,charge,&firstReport);
    const auto added=std::size_t(sites.nodeCount()-boundary.nodeCount());
    check(added>0 && added<=cap && firstReport.boxSites>0 && firstReport.faceSites>0 &&
          added==firstReport.boxSites+firstReport.faceSites,
          "first placement must support both an enclosed control and the curved wall within the cap");
    check(nativeQueries>=added && nativeQueries<=4*cap+std::size_t(boundary.cellCount()),
          "native proposal admission escaped the bounded candidate set");
    sameBoundary(unchanged,sites);
    sameMesh(unchanged,boundary);
    std::vector<double> cachedSizes;
    cachedSizes.reserve(std::size_t(sites.nodeCount()));
    bool enclosedSite=false;
    for (Index i=0;i<sites.nodeCount();++i) {
        const auto p=sites.node(i);
        const auto classification=sites.nodeClass(i);
        cachedSizes.push_back(field.eval(p,i<boundary.nodeCount() && classification.dim==2?
            classification.tag:-1));
        if (i<boundary.nodeCount()) continue;
        check(classification==region && strictInterior(p),
              "appended site must retain the requested solid identity and be native strictly IN");
        enclosedSite|=inControl(p);
        for (Index j=0;j<i;++j) {
            const auto delta=p-sites.node(j);
            const auto separation=.65*std::min(cachedSizes[std::size_t(i)],cachedSizes[std::size_t(j)]);
            check(std::hypot(delta.x,delta.y,delta.z)>=separation*(1.-1e-12),
                  "interior spacing site violates separation from the fixed shell or earlier sites");
        }
    }
    check(enclosedSite,"a control with no original boundary vertex received no first-placement support");
    work=0;
    OcctInteriorSiteReport replayReport;
    const auto replay=appendOcctInteriorSites(boundary,field,region,cap,
        strictInterior,charge,&replayReport);
    sameMesh(sites,replay);
    check(firstReport.boxSites==replayReport.boxSites && firstReport.faceSites==replayReport.faceSites &&
          firstReport.spacingRejects==replayReport.spacingRejects &&
          firstReport.nativeRejects==replayReport.nativeRejects,
          "first-placement replay changed deterministic admission decisions");

    std::size_t zeroCalls=0;
    const auto noNative=[&](Vec3) {++zeroCalls;return true;};
    const auto noCharge=[&](std::uint64_t) {++zeroCalls;};
    const auto zero=appendOcctInteriorSites(boundary,field,region,0,noNative,noCharge);
    check(zeroCalls==0,"zero site cap must avoid field, work and native proposal queries");
    sameMesh(boundary,zero);
    work=0;
    const auto one=appendOcctInteriorSites(boundary,field,region,1,strictInterior,charge);
    check(one.nodeCount()==boundary.nodeCount()+1 && inControl(one.node(boundary.nodeCount())),
          "one-site cap must seed the enclosed control without overshooting the cap");
    sameBoundary(boundary,one);

    std::size_t refusedQueries=0;
    const auto refuse=[&](Vec3) {++refusedQueries;return false;};
    work=0;
    OcctInteriorSiteReport refusalReport;
    const auto refused=appendOcctInteriorSites(boundary,field,region,cap,refuse,charge,&refusalReport);
    check(refusedQueries>0 && refusalReport.nativeRejects==refusedQueries,
          "native refusal must guard every geometrically eligible site");
    sameMesh(boundary,refused);

    struct WorkLimit {};
    std::size_t beforeStopQueries=0;
    bool stopped=false;
    const auto admitThenStop=[&](Vec3 p) {++beforeStopQueries;return strictInterior(p);};
    const auto stopAfterAdmission=[&](std::uint64_t) {if (beforeStopQueries) throw WorkLimit{};};
    try {
        (void)appendOcctInteriorSites(boundary,field,region,cap,admitThenStop,stopAfterAdmission);
    } catch (const WorkLimit&) {stopped=true;}
    check(stopped && beforeStopQueries==1,
          "charged placement must propagate an exhausted budget after staging a site");
    sameMesh(unchanged,boundary);
    check(snapshot(model)==originalCad,"interior site construction modified original native CAD");
}

void exerciseRegularOffcentre() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2.,2.,3.).Shape());
    check(imported.imported(),"native offcentre admission fixture import failed");
    sizing::RegionSizeSpecification specification;
    specification.baseSize=1.;
    auto compiled=sizing::compileRegionSizeField(specification);
    check(compiled.compiled,"constant offcentre size fixture failed");
    const OcctSpatialSizing field(*imported.model,{},std::move(compiled.field));
    // An isolated oriented face tests local geometry independently of the
    // closed-surface native integration fixture above. Its inward apex is the
    // fourth vertex of a unit regular tetrahedron.
    mesh::Mesh face;
    face.addNode({.5,.5,1.5},{0,1});
    face.addNode({1.5,.5,1.5},{0,2});
    face.addNode({1.,.5+std::sqrt(3.)/2,1.5},{0,3});
    face.addCell(mesh::CellType::Tri3,{0,1,2},{2,1});
    const auto original=face;
    OcctSolidClassifier classifier(OcctModelAccess::solid(*imported.model,1));
    const auto strictInterior=[&](Vec3 p) {
        classifier.Perform(gp_Pnt(p.x,p.y,p.z),imported.model->tolerance(3,1));
        return classifier.State()==TopAbs_IN;
    };
    const auto charge=[](std::uint64_t) {};
    const auto sites=appendOcctInteriorSites(face,field,{3,1},1,strictInterior,charge);
    check(sites.nodeCount()==4,"regular face failed to propose an admitted inward apex");
    sameBoundary(face,sites);
    const auto apex=sites.node(3);
    check(std::abs(apex.x-1.)<1e-13 && std::abs(apex.y-(.5+std::sqrt(3.)/6))<1e-13 &&
          std::abs(apex.z-(1.5-std::sqrt(2./3.)))<1e-13,
          "regular face offcentre does not have the regular-tetrahedron inward height");
    for (Index i=0;i<3;++i)
        check(std::abs((apex-face.node(i)).norm()-1.)<1e-13,
              "regular offcentre must be one target edge length from every face vertex");
    mesh::Mesh reversed;
    for (Index i=0;i<face.nodeCount();++i) reversed.addNode(face.node(i),face.nodeClass(i));
    reversed.addCell(mesh::CellType::Tri3,{0,2,1},{2,1});
    const auto otherSide=appendOcctInteriorSites(reversed,field,{3,1},1,strictInterior,charge);
    check(otherSide.nodeCount()==4 &&
          std::abs(otherSide.node(3).z-(1.5+std::sqrt(2./3.)))<1e-13,
          "offcentre direction must follow the actual oriented face winding");
    sameMesh(original,face);
}

using CoordinateKey=std::array<std::uint64_t,3>;
using BoundarySignature=std::tuple<
    std::vector<std::tuple<CoordinateKey,int,int>>,
    std::vector<std::tuple<int,int,std::array<CoordinateKey,3>>>>;

BoundarySignature boundarySignature(const mesh::Mesh& value) {
    const auto coordinate=[&](Index id) {
        const auto p=value.node(id);
        return CoordinateKey{std::bit_cast<std::uint64_t>(p.x),
            std::bit_cast<std::uint64_t>(p.y),std::bit_cast<std::uint64_t>(p.z)};
    };
    BoundarySignature signature;
    auto& [nodes,triangles]=signature;
    std::set<Index> boundaryNodes;
    for (Index i=0;i<value.cellCount();++i) {
        const auto cell=value.cell(i);
        if (cell.type!=mesh::CellType::Tri3) continue;
        const auto* v=value.cellNodes(i);
        std::array<CoordinateKey,3> triangle{coordinate(v[0]),coordinate(v[1]),coordinate(v[2])};
        const auto first=std::min_element(triangle.begin(),triangle.end());
        std::rotate(triangle.begin(),first,triangle.end());
        triangles.emplace_back(cell.onEntity.dim,cell.onEntity.tag,triangle);
        boundaryNodes.insert(v,v+3);
    }
    for (const auto id:boundaryNodes) {
        const auto classification=value.nodeClass(id);
        nodes.emplace_back(coordinate(id),classification.dim,classification.tag);
    }
    std::ranges::sort(nodes);
    std::ranges::sort(triangles);
    return signature;
}

unsigned interiorSeedFallbacks=0,constraintSeedAttempts=0;
void observeFallback(const geo::NativeCadStageEvent& event) noexcept {
    if (!event.beginning || !event.stage) return;
    if (std::string_view(event.stage)=="interior_seed_fallback") ++interiorSeedFallbacks;
    if (std::string_view(event.stage)=="constraint-recovery") ++constraintSeedAttempts;
}

void exerciseProductionPipeline() {
    auto imported=OcctModelAccess::fromShape(BRepPrimAPI_MakeCylinder(1.,2.).Shape());
    check(imported.imported(),"production cylinder fixture import failed");
    const auto& model=*imported.model;
    const auto originalCad=snapshot(model);
    int wall=-1;
    for (const auto* face:model.faces()) {
        const auto* descriptor=model.surfaceDescriptor(face->tag());
        if (descriptor && descriptor->kind==geo::SurfaceKind::Cylinder) wall=face->tag();
    }
    check(wall>0,"production cylinder fixture needs its original curved wall tag");
    geo::NativeCadMeshOptions options;
    options.geometry.targetSize=.7;
    options.geometry.threads=1;
    options.geometry.maximumNodes=50'000;
    options.geometry.maximumCells=250'000;
    options.geometry.maximumOutputBytes=128ULL*1024*1024;
    options.geometry.maximumMeshingWork=100'000'000;
    options.geometry.maximumAuditTemporaryBytes=512ULL*1024*1024;
    options.geometry.minimumElementQuality=0.;
    options.geometry.warningElementQuality=.25;
    options.maximumValidationWork=100'000'000;
    options.maximumPrimitiveWorkingBytes=512ULL*1024*1024;
    options.maximumNativeValidationThreads=1;
    options.maximumChordError=.0086602540378443865;
    options.maximumRelativeVolumeError=.01;
    options.minimumSize=.12;
    options.faceSizes.push_back({wall,.16,1.55});
    options.sizeBoxes.push_back({{-.3,-.3,.7},{.3,.3,1.3},.15,1.55});
    options.useSpatialSizing=true;
    // This regression isolates the earlier quality-controlled initial sites.
    options.useAdvancingInteriorFront=false;
    options.allowCentroidsWithinChordTolerance=false;
    options.recoveryMode=geo::NativeCadRecoveryMode::Standard;
    mesh::Mesh ordinary;
    const auto baseline=geo::meshNativeCad(model,options,ordinary);
    if (!baseline.published)
        std::cerr<<"production cylinder baseline: "<<geo::nativeCadMeshReasonName(baseline.reason)<<'\n';
    check(baseline.published && baseline.audit.valid && baseline.quality.valid && baseline.bodies.size()==1,
          "production cylinder baseline must publish an audited valid Tet4 mesh");

    options.improveNativeTetQuality=true;
    mesh::Mesh improved;
    const auto built=geo::meshNativeCad(model,options,improved);
    if (!built.published)
        std::cerr<<"production cylinder first placement: "<<geo::nativeCadMeshReasonName(built.reason)<<'\n';
    check(built.published && built.reason==geo::NativeCadMeshReason::None &&
          built.audit.complete && built.audit.valid && built.quality.complete && built.quality.valid,
          "public native CAD pipeline must publish an audited valid first-placement mesh");
    check(built.bodies.size()==1 && built.bodies.front().solidTag==1 &&
          std::string(built.bodies.front().route)=="native_spatial_sizing" &&
          built.bodies.front().boundaryEmbeddingCertified,
          "production regression must exercise the original solid's spatial native volume route");
    const auto& body=built.bodies.front();
    const auto tetCount=std::uint64_t(improved.countCells(mesh::CellType::Tet4));
    check(tetCount>0 && built.quality.tet4==tetCount && built.quality.volumeCells==tetCount &&
          body.nativeCentroidsChecked==tetCount && body.nativeCentroidThreads==1 &&
          body.nativeCentroidSourceGuard && body.centroidOutsideDistanceLimit==0. &&
          body.maximumOutsideCentroidDistance==0. && body.centroidsWithinChordTolerance==0 &&
          body.nativeCentroidDistanceQueries==0,
          "every production tetrahedron must retain strict original-body centroid coverage");
    check(built.quality.medianQuality>.70 && built.quality.minimumQuality>.01,
          "production first placement regressed the cylinder's solver quality distribution");
    check(boundarySignature(ordinary)==boundarySignature(improved),
          "quality opt-in changed exact boundary coordinates, node tags, oriented Tri3 cells or face tags");

    using Edge=std::array<Index,2>;
    const auto edge=[](Index a,Index b) {return Edge{std::min(a,b),std::max(a,b)};};
    std::set<Index> boundaryNodes;
    std::set<Edge> protectedEdges,volumeEdges;
    std::uint64_t poor=0;
    for (Index i=0;i<improved.cellCount();++i) {
        const auto cell=improved.cell(i);
        const auto* v=improved.cellNodes(i);
        if (cell.type==mesh::CellType::Tri3) {
            check(cell.onEntity.dim==2 && model.face(cell.onEntity.tag),
                  "production boundary lost its original CAD face classification");
            boundaryNodes.insert(v,v+3);
            for (unsigned j=0;j<3;++j) protectedEdges.insert(edge(v[j],v[(j+1)%3]));
            continue;
        }
        check(cell.type==mesh::CellType::Tet4 && cell.onEntity==mesh::Classification{3,1},
              "production volume cell lost the original solid's Tet4 classification");
        const auto a=improved.node(v[0]),b=improved.node(v[1]),
                   c=improved.node(v[2]),d=improved.node(v[3]);
        double edgeSquared=0.;
        for (unsigned j=0;j<4;++j) for (unsigned k=j+1;k<4;++k) {
            edgeSquared+=(improved.node(v[j])-improved.node(v[k])).norm2();
            volumeEdges.insert(edge(v[j],v[k]));
        }
        const auto determinant=(b-a).cross(c-a).dot(d-a);
        check(determinant>0.,"production first placement retained a nonpositive tetrahedron");
        const auto meanRatio=12*std::pow(determinant/2,2./3.)/edgeSquared;
        poor+=meanRatio<.25;
    }
    check(double(poor)/double(tetCount)<.02,
          "production first placement left too many tetrahedra below mean ratio .25");
    for (const auto& e:volumeEdges) {
        if (protectedEdges.contains(e) || !boundaryNodes.contains(e[0]) ||
            !boundaryNodes.contains(e[1])) continue;
        check((improved.node(e[0])-improved.node(e[1])).norm()<=.7,
              "production first placement retained a long unconstrained boundary-to-boundary bridge");
    }
    check(snapshot(model)==originalCad,"public native volume pipeline modified original CAD");
    std::cout<<"production spatial cylinder: tets="<<tetCount<<" poor="<<poor
             <<" median="<<built.quality.medianQuality<<" minimum="<<built.quality.minimumQuality<<'\n';

    // Give the augmented point set less construction work than it needs,
    // while retaining enough for the measured boundary-only seed. The optional
    // first placement must then retry the established construction, and both
    // attempts must remain charged to the same enclosing work admission.
    const auto boundarySeedWork=baseline.bodies.front().constraintRecovery.work;
    check(boundarySeedWork>0,"fallback regression needs measured boundary-only constraint work");
    const auto earlyLimit=boundarySeedWork+boundarySeedWork/3;
    interiorSeedFallbacks=constraintSeedAttempts=0;
    const auto previousObserver=options.stageObserver;
    options.stageObserver=observeFallback;
    mesh::Mesh fallbackMesh;
    OcctFrontOutcome fallback;
    {
        OcctStageObserverScope observer(options.stageObserver);
        fallback=buildOcctFrontVolume(model,1,options,100'000'000,fallbackMesh,
            earlyLimit,OcctFrontSurface::SpatialSizing);
    }
    options.stageObserver=previousObserver;
    if (!fallback.published)
        std::cerr<<"production cylinder bounded seed retry: "<<fallback.diagnostic
                 <<" / "<<fallback.constraintRecovery.diagnostic<<'\n';
    check(fallback.published && fallback.reason==geo::NativeCadMeshReason::None &&
          interiorSeedFallbacks==1 && constraintSeedAttempts==2,
          "limited interior seed must actually retry and publish the boundary-only construction");
    check(fallback.constraintRecovery.work>earlyLimit && fallback.workUsed<=100'000'000 &&
          fallback.workUsed>=fallback.constraintRecovery.work,
          "fallback must charge both constraint attempts within the enclosing work bound");
    const auto fallbackTets=std::uint64_t(fallbackMesh.countCells(mesh::CellType::Tet4));
    const auto& validation=fallback.nativeValidation;
    check(fallbackTets>0 && validation.completed && validation.sourceModelUnmodified &&
          validation.threadsUsed==1 && validation.points==fallbackTets &&
          validation.nativeQueries+validation.strictAdmissionReuses+validation.trimmedRayAdmissions==fallbackTets &&
          validation.pointsWithinDistanceBand==0 && validation.nativeDistanceQueries==0 &&
          validation.maximumOutsideDistance==0. && validation.workUsed<=options.maximumValidationWork &&
          validation.maximumWorkingPayloadEstimate<=options.maximumPrimitiveWorkingBytes,
          "bounded seed fallback must retain complete strict native validation and resource admissions");
    check(mesh::auditLinearVolumeBoundary(fallbackMesh).valid &&
          algo::analyzeMixedMeshQuality(fallbackMesh).valid &&
          boundarySignature(fallbackMesh)==boundarySignature(ordinary),
          "bounded seed fallback must retain a valid mesh and the exact original classified boundary");
    check(snapshot(model)==originalCad,"bounded seed retry modified original native CAD");
    std::cout<<"bounded spatial seed retry: limit="<<earlyLimit
             <<" cumulative_constraint_work="<<fallback.constraintRecovery.work
             <<" total_work="<<fallback.workUsed<<'\n';
}
}

int main() {
    try {
        exerciseNativePlacement();
        exerciseRegularOffcentre();
        exerciseProductionPipeline();
        std::cout<<"native first-placement sites, fixed boundary, separation and budget regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    } catch (...) {
        std::cerr<<"unexpected native interior site regression failure\n";
        return 2;
    }
}
