// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "../src/geometry/OcctContactCylinder.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax2.hxx>
#include <bit>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
using namespace paacs;
namespace {
void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
std::string snapshot(const geo::OcctModel& model){std::ostringstream s;BRepTools::Write(geo::detail::OcctModelAccess::shape(model),s);return s.str();}
bool same(Vec3 a,Vec3 b){return std::bit_cast<std::uint64_t>(a.x)==std::bit_cast<std::uint64_t>(b.x)&&
    std::bit_cast<std::uint64_t>(a.y)==std::bit_cast<std::uint64_t>(b.y)&&std::bit_cast<std::uint64_t>(a.z)==std::bit_cast<std::uint64_t>(b.z);}
std::uint64_t forbidden(const mesh::Mesh& a,const mesh::Mesh& b){std::uint64_t n=0;for(Index i=0;i<a.cellCount();++i){if(a.cell(i).onEntity.dim!=2)continue;const auto* u=a.cellNodes(i);
    const std::array<Vec3,3> x{a.node(u[0]),a.node(u[1]),a.node(u[2])};for(Index j=0;j<b.cellCount();++j){if(b.cell(j).onEntity.dim!=2)continue;const auto* v=b.cellNodes(j);
    const std::array<Vec3,3> y{b.node(v[0]),b.node(v[1]),b.node(v[2])};n+=bnd::tet_facet_detail::trianglesOpenIntersect(x,y)||bnd::tet_facet_detail::coplanarTrianglesHaveAreaOverlap(x,y);}}return n;}
mesh::Mesh select(const geo::OcctModel& model,const mesh::Mesh& input,int body){mesh::Mesh out;std::map<Index,Index> map;std::set<int> faces;for(auto use:model.faceUses(body))faces.insert(use.faceTag);
    for(Index i=0;i<input.cellCount();++i){const auto& cell=input.cell(i);if(!((cell.onEntity.dim==3&&cell.onEntity.tag==body)||(cell.onEntity.dim==2&&faces.contains(cell.onEntity.tag))))continue;
    std::array<Index,4> ids{};const int count=mesh::nodesPerCell(cell.type);check(count<=4,"nonlinear source");for(int j=0;j<count;++j){const auto old=input.cellNodes(i)[j];auto [it,inserted]=map.try_emplace(old,0);if(inserted)it->second=out.addNode(input.node(old),input.nodeClass(old));ids[j]=it->second;}out.addCell(cell.type,{ids.data(),static_cast<std::size_t>(count)},cell.onEntity);}return out;}
}
int main(int argc,char**argv){try{
    mesh::Mesh accepted;geo::NativeCadImportOutcome loaded;int target=2;
    if(argc==4){loaded=geo::importOcctCad(argv[1]);check(loaded.imported(),"real import");mesh::Mesh input;check(sdk::readMeshArchive(argv[2],input).completed,"real archive");accepted=select(*loaded.model,input,1);}
    else {const auto folder=argc==2?std::filesystem::path(argv[1]):std::filesystem::path("out/occt-contact-cylinder-tests");std::filesystem::create_directories(folder);
        const auto outer=BRepPrimAPI_MakeCylinder(3.,2.).Shape();const auto recess=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,1.9),gp_Dir(0,0,1)),1.,.3).Shape();
        BRepAlgoAPI_Cut cut(outer,recess);cut.SetNonDestructive(true);cut.Build();check(cut.IsDone()&&!cut.HasErrors(),"ring native cut");
        const auto cylinder=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0,0,2),gp_Dir(0,0,1)),1.,1.).Shape();
        TopoDS_Compound compound;BRep_Builder builder;builder.MakeCompound(compound);builder.Add(compound,cut.Shape());builder.Add(compound,cylinder);
        const auto path=folder/"ring-contact-cylinder.brep";check(BRepTools::Write(compound,path.string().c_str()),"native fixture write");geo::NativeCadImportOptions importing;importing.format=geo::CadFileFormat::Brep;loaded=geo::importOcctCad(path,importing);check(loaded.imported(),"native fixture import");
        geo::OcctBoundaryOptions boundary;boundary.linearDeflection=.005;boundary.maximumChordError=.02;const auto receipt=geo::buildOcctBoundary(*loaded.model,1,boundary,accepted);check(receipt.published,"prior actual ring boundary");}
    const auto before=snapshot(*loaded.model);const auto acceptedHash=mesh::meshHash(accepted);geo::detail::OcctContactCylinderOptions options;options.maximumChordError=.05;mesh::Mesh output;
    auto result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,options,output);
    std::cout<<geo::detail::occtContactCylinderReasonName(result.reason)<<" published="<<result.published<<" ring="<<result.ringNodes<<" nodes="<<output.nodeCount()<<" cells="<<output.cellCount()<<" volume_error="<<result.relativeVolumeError<<" sample="<<result.maximumBoundarySampleError<<" quality="<<result.quality.minimumQuality<<" work="<<result.workUsed<<'\n';
    check(result.published&&result.nativeCircleProved&&result.nativeCylinderProved&&result.sourceImmutabilityVerified,"generic contact route");
    check(result.audit.valid&&result.quality.valid&&result.relativeVolumeError<=.01&&result.maximumBoundarySampleError<=.05,"native/volume/audit requirements");
    check(forbidden(accepted,output)==0,"unchanged exact interbody predicates");check(mesh::meshHash(accepted)==acceptedHash&&snapshot(*loaded.model)==before,"source/accepted immutable");
    unsigned identical=0;for(Index i=0;i<output.nodeCount();++i){const auto cls=output.nodeClass(i);if(cls.dim!=1||cls.tag!=result.targetEdgeTag)continue;
        bool found=false;for(Index j=0;j<accepted.nodeCount();++j)if(same(accepted.node(j),output.node(i)))found=true;check(found,"actual rim coordinates retained bitwise");++identical;}
    check(identical+1==result.ringNodes,"classified native rim including seam vertex");
    const auto outputHash=mesh::meshHash(output);
    auto limited=options;limited.maximumWork=0;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::BudgetExceeded&&mesh::meshHash(output)==outputHash,"work rollback");
    limited=options;limited.maximumNodes=2;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::BudgetExceeded&&mesh::meshHash(output)==outputHash,"node rollback");
    limited=options;limited.maximumSnapshotBytes=8;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::BudgetExceeded&&mesh::meshHash(output)==outputHash,"snapshot rollback");
    limited=options;limited.maximumWorkingBytes=1024;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::BudgetExceeded&&mesh::meshHash(output)==outputHash,"working admission rollback");
    limited=options;std::atomic<bool> cancelled{true};limited.cancellation=&cancelled;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::Cancelled&&mesh::meshHash(output)==outputHash,"cancel rollback");
    limited=options;limited.maximumChordError=1e-14;result=geo::detail::buildOcctContactCylinder(*loaded.model,target,accepted,limited,output);check(!result.published&&result.reason==geo::detail::OcctContactCylinderReason::GeometryError&&mesh::meshHash(output)==outputHash,"measured geometry refusal rollback");
    check(mesh::meshHash(accepted)==acceptedHash&&snapshot(*loaded.model)==before,"all refusal source guards");
    if(argc==4)check(sdk::writeMeshArchive(output,argv[3]).completed,"generic archive");std::cout<<"PASS 13 contact/source/rollback checks\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
