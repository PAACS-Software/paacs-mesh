// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pln.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Compound.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <sstream>
#include <iostream>
#include <set>
using namespace paacs;using namespace paacs::geo;using namespace paacs::geo::detail;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::string snapshot(const OcctModel& m){std::ostringstream s;BRepTools::Write(OcctModelAccess::shape(m),s);return s.str();}
int main(){try{
 BRep_Builder builder;TopoDS_Compound shape;builder.MakeCompound(shape);builder.Add(shape,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());builder.Add(shape,BRepPrimAPI_MakeBox(gp_Pnt(.5,0,0),1.,1.,1.).Shape());auto imported=OcctModelAccess::fromShape(shape);check(imported.imported()&&imported.model->solids().size()==2,"overlapping assembly import");const auto& model=*imported.model;const auto source=snapshot(model);
 NativeCadMeshOptions o;o.geometry.targetSize=.6;o.maximumValidationWork=10000000;o.geometry.maximumMeshingWork=10000000;o.maximumNativeValidationThreads=2;
 mesh::Mesh output;output.addNode({9,8,7});const auto untouched=mesh::meshHash(output);auto rejected=meshNativeCad(model,o,output);check(!rejected.published&&rejected.reason==NativeCadMeshReason::UnresolvedBodyInterface&&mesh::meshHash(output)==untouched,"default conforming overlap refusal");
 o.assemblyPolicy=NativeCadAssemblyPolicy::IndependentParts;auto result=meshNativeCad(model,o,output);if(!result.published)std::cerr<<nativeCadMeshReasonName(result.reason)<<' '<<result.nativeContactDiagnostic<<'\n';check(result.published&&result.assemblyPolicy==NativeCadAssemblyPolicy::IndependentParts&&result.validatedParts==2&&result.bodies.size()==2,"independent assembly publication");check(result.interferences.size()==1&&result.interferences[0].state==NativeCadInterferenceState::MaterialOverlap&&std::abs(result.interferences[0].nativeOverlapVolume-.5)<1e-10,"native interference report");check(result.audit.valid&&std::abs(result.quality.totalVolume-2.)<1e-10,"part sum is not union volume");check(result.nativeContactSourceGuard&&snapshot(model)==source,"assembly source preserved");
 std::set<Index> nodes[2];for(Index i=0;i<output.cellCount();++i){const auto& c=output.cell(i);if(c.type!=mesh::CellType::Tet4)continue;check(c.onEntity.tag==1||c.onEntity.tag==2,"original body tag");for(unsigned j=0;j<4;++j)nodes[c.onEntity.tag-1].insert(output.cellNodes(i)[j]);}for(auto i:nodes[0])check(!nodes[1].contains(i),"parts were welded");for(const auto& b:result.bodies)check(b.nativeCentroidSourceGuard&&b.nativeCentroidsChecked>0,"original body centroid certificate");
 auto adaptive=o;
 adaptive.useAdaptiveTetRecovery=true;
 adaptive.reuseNativeQueries=true;
 mesh::Mesh adaptiveMesh;
 const auto adaptiveResult=meshNativeCad(model,adaptive,adaptiveMesh);
 check(adaptiveResult.published&&adaptiveResult.validatedParts==2&&mesh::meshHash(adaptiveMesh)==mesh::meshHash(output),
       "adaptive assembly preserves independently meshed body identities");
 check(adaptiveResult.interferences.size()==1&&adaptiveResult.interferences[0].state==NativeCadInterferenceState::MaterialOverlap&&
       adaptiveResult.interferences[0].nativeOverlapVolume==result.interferences[0].nativeOverlapVolume&&snapshot(model)==source,
       "adaptive assembly preserves native interference reporting");
 const auto saved=mesh::meshHash(output);auto low=o;low.maximumValidationWork=1;auto budget=meshNativeCad(model,low,output);check(!budget.published&&budget.reason==NativeCadMeshReason::BudgetExceeded&&mesh::meshHash(output)==saved,"assembly budget transaction");
 std::atomic_bool cancellation{true};low=o;low.geometry.cancellation=&cancellation;auto cancelled=meshNativeCad(model,low,output);check(!cancelled.published&&cancelled.reason==NativeCadMeshReason::Cancelled&&mesh::meshHash(output)==saved,"assembly cancellation transaction");
 low=o;low.geometry.minimumElementQuality=1.;auto strict=meshNativeCad(model,low,output);check(!strict.published&&strict.reason==NativeCadMeshReason::QualityRequirementNotMet&&mesh::meshHash(output)==saved,"assembly strict quality gate");
 low=o;low.faceSizes.push_back({999999,.1,1.3});auto bad=meshNativeCad(model,low,output);check(!bad.published&&bad.reason==NativeCadMeshReason::InvalidRequest&&mesh::meshHash(output)==saved&&snapshot(model)==source,"unknown selection not silently discarded");
 gp_Trsf transform;transform.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),.4);transform.SetTranslationPart(gp_Vec(5,2,3));TopoDS_Compound located;builder.MakeCompound(located);builder.Add(located,BRepPrimAPI_MakeBox(1.,1.,1.).Shape());builder.Add(located,BRepPrimAPI_MakeBox(1.,1.,1.).Shape().Moved(TopLoc_Location(transform)));auto placed=OcctModelAccess::fromShape(located);check(placed.imported(),"located assembly import");mesh::Mesh placedMesh;const auto placementBytes=snapshot(*placed.model);auto placedResult=meshNativeCad(*placed.model,o,placedMesh);if(!placedResult.published)std::cerr<<nativeCadMeshReasonName(placedResult.reason)<<' '<<placedResult.nativeContactDiagnostic<<'\n';check(placedResult.published&&placedResult.validatedParts==2&&placedResult.interferences.size()==1&&placedResult.interferences[0].state==NativeCadInterferenceState::NoMaterialOverlap&&snapshot(*placed.model)==placementBytes,"located part geometry/history");
 TopoDS_Compound cylinders;builder.MakeCompound(cylinders);builder.Add(cylinders,BRepPrimAPI_MakeCylinder(1.,2.).Shape());builder.Add(cylinders,BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(5,0,0),gp_Dir(0,0,1)),1.,2.).Shape());auto sized=OcctModelAccess::fromShape(cylinders);check(sized.imported(),"sized assembly import");int selectedCap=-1;for(const auto* face:sized.model->faces())if(face->kind()==SurfaceKind::Plane){for(auto owner:sized.model->faceOwners(face->tag()))if(owner.solidTag==1)selectedCap=face->tag();if(selectedCap>0)break;}check(selectedCap>0,"original selected cap identity");auto sizing=o;sizing.minimumSize=.2;sizing.faceSizes.push_back({selectedCap,.25,1.3});mesh::Mesh sizedMesh;const auto sizedSource=snapshot(*sized.model);const auto sizedResult=meshNativeCad(*sized.model,sizing,sizedMesh);if(!sizedResult.published)std::cerr<<nativeCadMeshReasonName(sizedResult.reason)<<' '<<sizedResult.nativeContactDiagnostic<<'\n';check(sizedResult.published&&sizedResult.validatedParts==2&&snapshot(*sized.model)==sizedSource,"selected sizing must not reject unselected parts");unsigned cellsByPart[2]{};for(Index i=0;i<sizedMesh.cellCount();++i)if(sizedMesh.cell(i).type==mesh::CellType::Tet4){const int part=sizedMesh.cell(i).onEntity.tag;check(part==1||part==2,"sized part identity");++cellsByPart[part-1];}const auto& selectedProfile=sizedResult.bodies.at(0).primitive;check(cellsByPart[0]>0&&cellsByPart[1]>0&&selectedProfile.profiledVolume&&!sizedResult.bodies.at(1).primitive.profiledVolume&&selectedProfile.maximumSectionNodes>selectedProfile.minimumSectionNodes,"selected-face sizing applies to its original part only");const auto sizedHash=mesh::meshHash(sizedMesh);sizing.faceSizes.clear();const auto missingField=meshNativeCad(*sized.model,sizing,sizedMesh);check(!missingField.published&&missingField.reason==NativeCadMeshReason::UnsupportedControl&&mesh::meshHash(sizedMesh)==sizedHash,"minimum-size request without a face control remains explicit");
 TopoDS_Compound auxiliary;builder.MakeCompound(auxiliary);builder.Add(auxiliary,shape);builder.Add(auxiliary,BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(10,10,10),gp_Dir(0,0,1)),-1.,1.,-1.,1.).Face());NativeCadImportOptions auxiliaryOptions;auxiliaryOptions.retainAuxiliaryGeometry=true;auto withAuxiliary=OcctModelAccess::fromShape(auxiliary,auxiliaryOptions);check(withAuxiliary.imported(),"assembly with retained auxiliary face");int auxiliaryFace=-1;for(const auto* face:withAuxiliary.model->faces())if(withAuxiliary.model->faceOwners(face->tag()).empty())auxiliaryFace=face->tag();check(auxiliaryFace>0,"retained auxiliary face identity");auto auxiliaryControl=o;auxiliaryControl.faceSizes.push_back({auxiliaryFace,.25,1.3});const auto auxiliaryResult=meshNativeCad(*withAuxiliary.model,auxiliaryControl,sizedMesh);check(!auxiliaryResult.published&&auxiliaryResult.reason==NativeCadMeshReason::InvalidRequest&&mesh::meshHash(sizedMesh)==sizedHash,"auxiliary sizing selection must not be silently discarded");
 std::cout<<"Independent part identity, interference, sizing, native validation and transactional gates passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
