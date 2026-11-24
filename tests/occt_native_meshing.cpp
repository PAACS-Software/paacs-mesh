// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <chrono>
#include <atomic>
#include <sstream>

using namespace paacs;
namespace {
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Case {const char* name;TopoDS_Shape shape;bool supported;};
std::string nativeSnapshot(const geo::OcctModel& model) {
    std::ostringstream stream;
    BRepTools::Write(geo::detail::OcctModelAccess::shape(model),stream);
    check(stream.good(),"native snapshot serialization");
    return stream.str();
}
void checkResourcePolicyPresets(const geo::OcctModel& model,const std::string& sourceSnapshot) {
    check(model.solids().size()==1,"resource preset fixture must contain one original body");
    const auto solidTag=model.solids().front()->tag();
    const std::array presets{geo::NativeCadMeshOptions{},geo::normalNativeCadMeshOptions()};
    for(std::size_t i=0;i<presets.size();++i) {
        auto request=presets[i];request.geometry.targetSize=.7;
        const bool normal=i==1;
        const auto expectedPolicy=normal?geo::CadResourcePolicy::Normal:geo::CadResourcePolicy::Bounded;
        const auto expectedWork=normal?geo::cadNormalWorkAllowance:std::uint64_t{50'000'000};
        check(request.construction.resourcePolicy==expectedPolicy&&
              request.maximumValidationWork==expectedWork&&request.geometry.maximumMeshingWork==expectedWork,
              "aggregate/factory resource policy or inherited work allowance changed");
        check(request.construction.effort==geo::CadConstructionEffort::Legacy&&
              request.maximumPrimitiveWorkingBytes==512ULL*1024*1024,
              "resource preset changed construction effort or inherited working memory");
        check(!request.useAdvancingInteriorFront.has_value()&&request.advancingInteriorFrontEnabled(1),
              "resource preset no longer selects automatic Tet4 construction");
        mesh::Mesh destination;destination.addNode({-10,-10,-10},{3,99});
        const auto sentinel=mesh::meshHash(destination);
        const auto accepted=geo::meshNativeCad(model,request,destination);
        check(accepted.published&&!accepted.destinationUnchanged&&accepted.audit.valid&&accepted.quality.valid&&
              accepted.bodies.size()==1&&accepted.bodies.front().solidTag==solidTag&&
              mesh::meshHash(destination)!=sentinel,
              "resource preset failed native publication or changed original body identity");
        check(accepted.validationWork>1&&accepted.validationWork<=request.maximumValidationWork&&
              accepted.validationWork<=request.geometry.maximumMeshingWork,
              "resource preset did not perform charged work within its effective allowance");
        check(accepted.bodies.front().initialFront.work>0&&
              std::string_view(accepted.bodies.front().initialFront.diagnostic)!="not_attempted",
              "resource preset bypassed automatic initial-point construction");
        std::uint64_t tetrahedra=0;
        for(Index cell=0;cell<destination.cellCount();++cell) {
            if(destination.cell(cell).type!=mesh::CellType::Tet4)continue;
            ++tetrahedra;
            check(destination.cell(cell).onEntity==mesh::Classification{3,solidTag},
                  "resource preset lost original tetrahedron body ownership");
        }
        check(tetrahedra>0&&nativeSnapshot(model)==sourceSnapshot,
              "resource preset produced no tetrahedra or changed the original native BRep");
        const auto retained=mesh::meshHash(destination);
        // Exercise each public work field separately; the other field keeps
        // the preset's allowance, so either ignored override is observable.
        for(unsigned field=0;field<2;++field) {
            auto limited=request;
            if(field==0)limited.geometry.maximumMeshingWork=1;
            else limited.maximumValidationWork=1;
            const auto refused=geo::meshNativeCad(model,limited,destination);
            check(!refused.published&&refused.destinationUnchanged&&
                  refused.reason==geo::NativeCadMeshReason::BudgetExceeded&&refused.validationWork<=1&&
                  mesh::meshHash(destination)==retained&&nativeSnapshot(model)==sourceSnapshot,
                  "resource preset bypassed explicit work or changed retained mesh/source on refusal");
        }
    }
}
}
int main(int argc,char** argv) {
    const bool keep=argc==3&&std::string_view(argv[1])=="--gallery";
    const auto dir=keep?std::filesystem::path(argv[2]):std::filesystem::temp_directory_path()/
        ("paacs-native-mesh-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(dir);
        const auto cylinder=BRepPrimAPI_MakeCylinder(1.,2.).Shape();
        const auto tube=BRepAlgoAPI_Cut(cylinder,BRepPrimAPI_MakeCylinder(.55,2.).Shape()).Shape();
        const auto hollowFrustum=BRepAlgoAPI_Cut(BRepPrimAPI_MakeCone(2.,1.8,2.).Shape(),
            BRepPrimAPI_MakeCone(1.9,1.7,2.).Shape()).Shape();
        gp_Trsf transform;transform.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,2,3)),.71);
        transform.SetTranslationPart(gp_Vec(4,-2,7));
        BRep_Builder builder;TopoDS_Compound assembly;builder.MakeCompound(assembly);
        builder.Add(assembly,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
        builder.Add(assembly,BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(3,0,0),gp_Dir(0,0,1)),.6,1.5).Shape());
        TopoDS_Compound touching;builder.MakeCompound(touching);
        builder.Add(touching,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
        builder.Add(touching,BRepPrimAPI_MakeBox(gp_Pnt(1,0,0),1.,1.,1.).Shape());
        std::vector<Case> cases={{"cylinder",cylinder,true},{"tube",tube,true},
            {"torus",BRepPrimAPI_MakeTorus(2.,.6).Shape(),true},
            {"elbow",BRepPrimAPI_MakeTorus(2.,.6,std::numbers::pi/2).Shape(),true},
            {"cone",BRepPrimAPI_MakeCone(1.,0.,2.).Shape(),true},
            {"frustum",BRepPrimAPI_MakeCone(1.,.5,2.).Shape(),true},
            {"rotated_cylinder",BRepBuilderAPI_Transform(cylinder,transform,true).Shape(),true},
            {"separated_bodies",assembly,true},{"touching_bodies",touching,false},
            {"sphere",BRepPrimAPI_MakeSphere(1.).Shape(),false},
            {"partial_cylinder",BRepPrimAPI_MakeCylinder(1.,2.,std::numbers::pi).Shape(),false},
            {"hollow_frustum",hollowFrustum,false}};
        int checks=0;std::ofstream manifest;
        if(keep){manifest.open(dir/"manifest.json");manifest<<"{\"title\":\"Native CAD meshes and selected-face refinement\",\"cases\":[";}
        bool first=true;
        for(const auto& fixture:cases) {
            const auto path=dir/(std::string(fixture.name)+".step");STEPControl_Writer writer;
            check(writer.Transfer(fixture.shape,STEPControl_AsIs)==IFSelect_RetDone,"STEP transfer");
            check(writer.Write(path.string().c_str())==IFSelect_RetDone,"STEP write");
            auto cad=geo::importOcctCad(path);check(cad.imported(),"native CAD import");
            const auto importedSnapshot=nativeSnapshot(*cad.model);
            check(nativeSnapshot(*cad.model)==importedSnapshot,"native snapshot serialization stable");
            if(std::string_view(fixture.name)=="cylinder") {
                checkResourcePolicyPresets(*cad.model,importedSnapshot);++checks;
            }
            geo::NativeCadMeshOptions options;options.geometry.targetSize=.3;
            // This compatibility matrix exercises analytic sweeps and established
            // star/partition routes. The preset tests above cover automatic Tet4 dispatch.
            options.useAdvancingInteriorFront=false;
            if(!fixture.supported)options.useNativeSurfaceFallback=false;
            mesh::Mesh mesh;mesh.addNode({-10,-10,-10},{3,99});
            const auto before=mesh::meshHash(mesh);
            auto result=geo::meshNativeCad(*cad.model,options,mesh);
            check(nativeSnapshot(*cad.model)==importedSnapshot,"recognition preserved original native BRep");++checks;
            std::cout<<fixture.name<<": "<<geo::nativeCadMeshReasonName(result.reason)<<" nodes="<<mesh.nodeCount()<<" cells="<<mesh.cellCount()<<'\n';
            if(!fixture.supported){
                check(!result.published&&result.destinationUnchanged&&mesh::meshHash(mesh)==before,"disabled fallback rollback");++checks;
                if(std::string_view(fixture.name)=="touching_bodies"){
                    auto contacts=options;contacts.geometry.targetSize=3.;
                    contacts.separateMatchedBoxFaceContacts=true;
                    result=geo::meshNativeCad(*cad.model,contacts,mesh);
                    check(result.published&&result.audit.valid&&result.quality.valid&&
                          result.separateBoxFaceContactPairs==1&&result.matchedBoxFaceTriangles==2&&
                          result.bodies.size()==2&&mesh.nodeCount()==16,"separate exact native box-face contacts");++checks;
                    for(const auto& body:result.bodies)check(body.boxContactVerticesEvaluated==8,"native box contact vertex identities evaluated");
                    const auto accepted=mesh::meshHash(mesh);
                    const auto repeat=geo::meshNativeCad(*cad.model,contacts,mesh);
                    check(repeat.published&&mesh::meshHash(mesh)==accepted,"box-face contact deterministic replay");++checks;
                    auto fine=contacts;fine.geometry.targetSize=.3;
                    check(!geo::meshNativeCad(*cad.model,fine,mesh).published&&mesh::meshHash(mesh)==accepted,"unsupported fine box-contact resolution rollback");++checks;
                    auto disabled=contacts;disabled.separateMatchedBoxFaceContacts=false;
                    check(geo::meshNativeCad(*cad.model,disabled,mesh).reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(mesh)==accepted,"box-face opt-in disabled rollback");++checks;
                    check(nativeSnapshot(*cad.model)==importedSnapshot,"box-face contact source preserved");++checks;
                }
                if(std::string_view(fixture.name)!="touching_bodies"){
                    options.useNativeSurfaceFallback=true;options.geometry.warningElementQuality=.25;
                    result=geo::meshNativeCad(*cad.model,options,mesh);
                    const bool partition=std::string_view(fixture.name)=="hollow_frustum";
                    check(result.published&&result.audit.valid&&result.quality.valid&&result.bodies.size()==1&&
                          std::string_view(result.bodies.front().route)==(partition?"native_partition_star":"native_surface_star"),"curved native surface fallback");++checks;
                    if(partition){
                        const auto accepted=mesh::meshHash(mesh);auto disabled=options;disabled.useNativePartitionFallback=false;
                        check(!geo::meshNativeCad(*cad.model,disabled,mesh).published&&mesh::meshHash(mesh)==accepted,"disabled partition preserves prior mesh");++checks;
                        auto limited=options;limited.maximumNativePartitionPieces=1;
                        const auto uncut=geo::meshNativeCad(*cad.model,limited,mesh);
                        if(!uncut.published||uncut.bodies.empty()||uncut.bodies.front().partitionPieces!=0){
                            std::cerr<<"uncut diagnostic: published="<<uncut.published<<" reason="<<geo::nativeCadMeshReasonName(uncut.reason)
                                     <<" work="<<uncut.validationWork;
                            if(!uncut.bodies.empty()){const auto& b=uncut.bodies.front();std::cerr<<" route="<<b.route<<" pieces="<<b.partitionPieces
                                <<" native_reason="<<b.nativeVolumeReason<<" boundary="<<b.boundaryEmbeddingCertified;}
                            std::cerr<<'\n';
                        }
                        check(uncut.published&&uncut.audit.valid&&uncut.quality.valid&&uncut.bodies.size()==1&&
                              std::string_view(uncut.bodies.front().route)=="native_plc_cavity"&&
                              uncut.bodies.front().partitionPieces==0&&uncut.bodies.front().boundaryEmbeddingCertified,
                              "partition limit permits certified uncut recovery");++checks;
                        const auto recoveredHash=mesh::meshHash(mesh);auto exhausted=limited;
                        exhausted.maximumValidationWork=1;
                        check(geo::meshNativeCad(*cad.model,exhausted,mesh).reason==geo::NativeCadMeshReason::BudgetExceeded&&
                              mesh::meshHash(mesh)==recoveredHash,"uncut recovery work budget rollback");++checks;
                        check(result.bodies.front().partitionInterfaces>0&&result.bodies.front().boundaryEmbeddingCertified,"partition joined interface certificate");++checks;
                        auto coarse=options;coarse.maximumChordError=.5;coarse.geometry.targetSize=1.;
                        // Keep the initial coarse volume approximation admitted
                        // so the independent native centroid gate drives repair.
                        coarse.maximumRelativeVolumeError=.1;
                        const auto recovered=geo::meshNativeCad(*cad.model,coarse,mesh);
                        check(recovered.published&&recovered.audit.valid&&recovered.quality.valid&&
                              recovered.validationWork<=coarse.maximumValidationWork,"bounded coarse curved boundary recovery");++checks;
                        const auto& attempts=recovered.bodies.front().partitionAttemptReports;
                        std::uint64_t consumed=0;for(const auto& attempt:attempts)consumed+=attempt.work;
                        check(!attempts.empty()&&consumed<=recovered.validationWork,
                              "all partition attempts remain charged to the common work budget");++checks;
                        result=geo::meshNativeCad(*cad.model,options,mesh);
                        check(result.published&&mesh::meshHash(mesh)==accepted,"recovery preserves normal deterministic mesh");++checks;
                        for(Index i=0;i<mesh.cellCount();++i){if(mesh.cell(i).type!=mesh::CellType::Tet4)continue;
                            const auto ids=mesh.cellNodes(i);Vec3 center{};for(unsigned j=0;j<4;++j)center=center+mesh.node(ids[j]);
                            check(cad.model->solid(mesh.cell(i).onEntity.tag)->isInside(center/4.),"hollow CAD retains empty bore");}
                    }
                    check(result.bodies.front().relativeVolumeError<=.01&&result.bodies.front().maximumBoundarySampleError<=result.bodies.front().surfaceChordErrorLimit*1.01&&
                          result.bodies.front().maximumBoundaryNodeError<1e-7,"fallback CAD geometry measurements");++checks;
                    const auto accepted=mesh::meshHash(mesh);auto strict=options;strict.geometry.minimumElementQuality=1;
                    check(geo::meshNativeCad(*cad.model,strict,mesh).reason==geo::NativeCadMeshReason::QualityRequirementNotMet&&mesh::meshHash(mesh)==accepted,"fallback strict quality rollback");++checks;
                    check(nativeSnapshot(*cad.model)==importedSnapshot,"fallback preserved native snapshot");++checks;
                    if(keep){const auto archive=dir/(std::string(fixture.name)+".paacsmesh");check(sdk::writeMeshArchive(mesh,archive.string()).completed,"fallback gallery archive");
                        if(!first)manifest<<',';first=false;manifest<<"{\"name\":\""<<fixture.name<<"\",\"mesh\":\""<<archive.filename().string()<<"\",\"source\":\""<<path.filename().string()<<"\",\"units\":\"mm\",\"section\":{\"axis\":\"y\",\"fraction\":0.5}}";}
                }
                continue;
            }
            check(result.published&&result.audit.valid&&result.quality.valid,"native CAD meshing");
            for(Index i=0;i<mesh.nodeCount();++i){const auto c=mesh.nodeClass(i);check(c.dim==0?cad.model->vertex(c.tag)!=nullptr:
                c.dim==1?cad.model->edge(c.tag)!=nullptr:c.dim==2?cad.model->face(c.tag)!=nullptr:c.dim==3&&cad.model->solid(c.tag)!=nullptr,"CAD node identity");}
            for(Index i=0;i<mesh.cellCount();++i){const auto c=mesh.cell(i).onEntity;check(c.dim==3?cad.model->solid(c.tag)!=nullptr:c.dim==2&&cad.model->face(c.tag)!=nullptr,"CAD cell identity");}
            for(const auto& b:result.bodies)check(b.relativeVolumeError<=.01&&b.maximumBoundaryNodeError<1e-7,"CAD volume/node error");
            ++checks;
            if(keep){const auto archive=dir/(std::string(fixture.name)+".paacsmesh");check(sdk::writeMeshArchive(mesh,archive.string()).completed,"gallery archive");
                if(!first)manifest<<',';first=false;manifest<<"{\"name\":\""<<fixture.name<<"\",\"mesh\":\""<<archive.filename().string()<<"\",\"source\":\""<<path.filename().string()<<"\",\"units\":\"mm\",\"section\":{\"axis\":\"y\",\"fraction\":0.5}}";}
            const auto publishedHash=mesh::meshHash(mesh);
            auto budget=options;budget.maximumValidationWork=0;const auto budgeted=geo::meshNativeCad(*cad.model,budget,mesh);
            check(!budgeted.published&&budgeted.reason==geo::NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(mesh)==publishedHash,"validation budget rollback");++checks;
            std::atomic<bool> cancel{true};auto cancelled=options;cancelled.geometry.cancellation=&cancel;
            check(geo::meshNativeCad(*cad.model,cancelled,mesh).reason==geo::NativeCadMeshReason::Cancelled&&mesh::meshHash(mesh)==publishedHash,"cancel rollback");++checks;
            if(std::string_view(fixture.name)=="separated_bodies") {
                // A control for a planar body must not become an unknown tag
                // inside a different body's view, or reject the cylinder.
                auto controlled=options;algo::GeometryBodyControl control;
                control.solidTag=result.bodies.front().solidTag;control.targetSize=.4;
                control.elements=algo::GeometryElements::Tet4;controlled.geometry.bodyControls={control};
                check(geo::meshNativeCad(*cad.model,controlled,mesh).published,"per-body control filtering");++checks;
                auto cap=options;
                for(const auto& b:result.bodies)if(std::string_view(b.route)=="cylinder")cap.faceSizes={{b.startFaceTag,.18,1.3}};
                cap.minimumSize=.15;
                check(geo::meshNativeCad(*cad.model,cap,mesh).published,"per-body face-control filtering");++checks;
            }
            if(std::string_view(fixture.name)=="cylinder"){
                check(result.bodies.size()==1&&std::string_view(result.bodies.front().route)=="cylinder"&&
                      cad.model->face(result.bodies.front().startFaceTag)!=nullptr&&
                      result.bodies.front().startFaceTag!=result.bodies.front().endFaceTag,
                      "profile-cap compatibility requires audited analytic cylinder face tags");++checks;
                options.faceSizes={{result.bodies.at(0).startFaceTag,.15,1.3}};options.minimumSize=.12;
                auto refined=geo::meshNativeCad(*cad.model,options,mesh);
                check(refined.published&&refined.bodies.at(0).primitive.profiledVolume,"imported cap sizing");
                const auto& p=refined.bodies.at(0).primitive;
                check(p.maximumSectionNodes>p.minimumSectionNodes&&p.maximumAxialGrowth<=1.3*(1+1e-10),"cap sizing section/growth measurements");++checks;
                if(keep){check(sdk::writeMeshArchive(mesh,(dir/"cylinder_refined_cap.paacsmesh").string()).completed,"refined archive");
                    manifest<<",{\"name\":\"Cylinder selected cap 0.15 / global 0.30\",\"mesh\":\"cylinder_refined_cap.paacsmesh\",\"source\":\"cylinder.step\",\"units\":\"mm\",\"section\":{\"axis\":\"y\",\"fraction\":0.5}}";}
                options.optimizeProfileInterior=true;options.geometry.warningElementQuality=.25;
                auto optimized=geo::meshNativeCad(*cad.model,options,mesh);
                check(optimized.published&&optimized.bodies.front().optimization.completed&&optimized.qualityWarning&&
                    optimized.quality.minimumQuality>refined.quality.minimumQuality,"bounded profile quality improvement");++checks;
                if(keep){check(sdk::writeMeshArchive(mesh,(dir/"cylinder_refined_cap_optimized.paacsmesh").string()).completed,"optimized archive");
                    manifest<<",{\"name\":\"Selected cap / bounded interior improvement\",\"mesh\":\"cylinder_refined_cap_optimized.paacsmesh\",\"source\":\"cylinder.step\",\"units\":\"mm\",\"section\":{\"axis\":\"y\",\"fraction\":0.5}}";}
                const auto optimizedHash=mesh::meshHash(mesh);auto strict=options;strict.geometry.minimumElementQuality=.25;
                check(geo::meshNativeCad(*cad.model,strict,mesh).reason==geo::NativeCadMeshReason::QualityRequirementNotMet&&
                    mesh::meshHash(mesh)==optimizedHash,"optimized strict quality rollback");++checks;
                auto tiny=options;tiny.maximumQualityNodeTrials=1;
                check(geo::meshNativeCad(*cad.model,tiny,mesh).reason==geo::NativeCadMeshReason::BudgetExceeded&&
                    mesh::meshHash(mesh)==optimizedHash,"profile optimization budget rollback");++checks;
                auto bad=options;bad.faceSizes={{cad.model->solids().front()->boundaryFaceTags().front(),.1,1.3}};
                // A curved-wall control is rejected instead of becoming an axial cap.
                for(auto face:cad.model->faces())if(face->kind()==geo::SurfaceKind::Cylinder)bad.faceSizes.front().faceTag=face->tag();
                const auto refinedHash=mesh::meshHash(mesh);auto refusal=geo::meshNativeCad(*cad.model,bad,mesh);
                check(refusal.reason==geo::NativeCadMeshReason::UnsupportedControl&&mesh::meshHash(mesh)==refinedHash,"wall control rollback");++checks;
            }
            check(nativeSnapshot(*cad.model)==importedSnapshot,"repeated meshing preserved original native BRep");++checks;
        }
        for(const auto [name,offset]:std::array<std::pair<const char*,Vec3>,4>{{
                {"point_contact",{1,1,1}},{"edge_contact",{1,1,0}},
                {"face_contact",{1,0,0}},{"overlap_contact",{.5,0,0}}}}) {
            TopoDS_Compound pair;builder.MakeCompound(pair);
            builder.Add(pair,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());
            builder.Add(pair,BRepPrimAPI_MakeBox(gp_Pnt(offset.x,offset.y,offset.z),1.,1.,1.).Shape());
            const auto path=dir/(std::string(name)+".step");STEPControl_Writer writer;
            check(writer.Transfer(pair,STEPControl_AsIs)==IFSelect_RetDone&&writer.Write(path.string().c_str())==IFSelect_RetDone,"contact STEP write");
            auto cad=geo::importOcctCad(path);check(cad.imported(),"contact CAD import");
            const auto snapshot=nativeSnapshot(*cad.model);mesh::Mesh mesh;mesh.addNode({-9,-9,-9},{3,99});const auto original=mesh::meshHash(mesh);
            geo::NativeCadMeshOptions options;options.geometry.targetSize=.5;
            check(geo::meshNativeCad(*cad.model,options,mesh).reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(mesh)==original,"default contact refusal");++checks;
            options.separateLowerDimensionalContacts=true;
            const auto result=geo::meshNativeCad(*cad.model,options,mesh);
            if(std::string_view(name)=="point_contact"||std::string_view(name)=="edge_contact") {
                check(result.published&&result.separateContactPairs==1&&result.audit.valid&&result.quality.valid&&
                      result.bodies.size()==2&&std::abs(result.quality.totalVolume-2)<1e-10,"separate contact mesh");
                const auto published=mesh::meshHash(mesh);options.maximumValidationWork=1;
                check(geo::meshNativeCad(*cad.model,options,mesh).reason==geo::NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(mesh)==published,"contact budget rollback");++checks;
            }else check(!result.published&&result.reason==geo::NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(mesh)==original,"face/overlap refusal despite separate-contact policy");
            check(nativeSnapshot(*cad.model)==snapshot,"contact tests preserve source BRep");++checks;
        }
        {
            // A native circle contact must reuse the accepted rim; independently
            // tessellating the hole and cap can introduce coplanar area overlap.
            TopoDS_Compound pair;builder.MakeCompound(pair);
            const auto outer=BRepPrimAPI_MakeCylinder(3.,2.).Shape();
            const auto recess=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,1.9),gp_Dir(0,0,1)),1.,.3).Shape();
            builder.Add(pair,BRepAlgoAPI_Cut(outer,recess).Shape());
            builder.Add(pair,BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,2),gp_Dir(0,0,1)),1.,1.).Shape());
            const auto path=dir/"circle-contact-cylinder.brep";
            check(BRepTools::Write(pair,path.string().c_str()),"contact cylinder fixture write");
            geo::NativeCadImportOptions importing;importing.format=geo::CadFileFormat::Brep;
            auto cad=geo::importOcctCad(path,importing);check(cad.imported(),"contact cylinder import");
            const auto snapshot=nativeSnapshot(*cad.model);mesh::Mesh mesh;
            geo::NativeCadMeshOptions options;options.geometry.targetSize=.6;
            options.geometry.maximumMeshingWork=options.maximumValidationWork=5'000'000;
            options.separateLowerDimensionalContacts=true;
            const auto result=geo::meshNativeCad(*cad.model,options,mesh);
            check(result.published&&result.audit.valid&&result.quality.valid&&result.bodies.size()==2&&
                  std::string_view(result.bodies.back().route)=="native_contact_cylinder"&&
                  result.bodies.back().contactRingNodes>=3&&result.contactTrianglePairs>0,
                  "integrated contact cylinder retains exact global checks");++checks;
            const auto accepted=mesh::meshHash(mesh);options.maximumPrimitiveWorkingBytes=1024;
            check(geo::meshNativeCad(*cad.model,options,mesh).reason==geo::NativeCadMeshReason::BudgetExceeded&&
                  mesh::meshHash(mesh)==accepted,"integrated contact working budget rollback");++checks;
            check(nativeSnapshot(*cad.model)==snapshot,"integrated contact preserves source");++checks;
        }
        if(keep)manifest<<"]}\n";
        std::cout<<"Native CAD meshing: PASS ("<<checks<<" checks)\n";
        if(!keep)std::filesystem::remove_all(dir);return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
