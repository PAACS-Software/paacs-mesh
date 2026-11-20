// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/OcctMesher.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pln.hxx>
#include <cstdio>
#include <sstream>
#include <stdexcept>
using namespace paacs;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
std::string snapshot(const geo::OcctModel& m){std::ostringstream s;BRepTools::Write(geo::detail::OcctModelAccess::shape(m),s);return s.str();}
int main(){try{
    BRep_Builder builder;TopoDS_Compound document;builder.MakeCompound(document);
    builder.Add(document,BRepPrimAPI_MakeBox(1,1,1).Shape());
    builder.Add(document,BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0,0,3),gp_Dir(0,0,1)),-1,1,-1,1).Face());
    auto strict=geo::detail::OcctModelAccess::fromShape(document);
    require(!strict.imported()&&strict.reason==geo::CadImportReason::InvalidShape,"default rejects auxiliary geometry");
    geo::NativeCadImportOptions io;io.retainAuxiliaryGeometry=true;
    auto retained=geo::detail::OcctModelAccess::fromShape(document,io);
    require(retained.imported()&&retained.solids==1&&retained.faces==7&&retained.auxiliaryFaces==1&&
            retained.auxiliaryEdges==4&&retained.auxiliaryVertices==4&&retained.auxiliaryGeometryRetained,"retains full document topology");
    require(retained.model->faces().size()==7&&retained.model->faceOwners(7).empty(),"auxiliary remains queryable and unowned");
    const auto before=snapshot(*retained.model);
    geo::NativeCadMeshOptions options;options.geometry.targetSize=.5;
    mesh::Mesh mesh;const auto result=geo::meshNativeCad(*retained.model,options,mesh);
    require(result.published&&result.audit.valid&&result.bodies.size()==1&&result.bodies[0].solidTag==1,"meshes entire material body");
    require(std::abs(result.quality.totalVolume-1)<1e-12,"auxiliary does not become material");
    for(Index i=0;i<mesh.cellCount();++i){const auto c=mesh.cell(i);require(c.onEntity.dim!=2||c.onEntity.tag!=7,"auxiliary face excluded from volume boundary");}
    require(before==snapshot(*retained.model),"entire retained document source unchanged");
    io.maximumFaces=6;auto budget=geo::detail::OcctModelAccess::fromShape(document,io);
    require(!budget.imported()&&budget.reason==geo::CadImportReason::TopologyBudget,"full document topology counted for budgets");
    const char* path="paacs-auxiliary-fixture.brep";BRepTools::Write(document,path);
    io.maximumFaces=100;io.format=geo::CadFileFormat::Brep;auto imported=geo::importOcctCad(path,io);std::remove(path);
    require(imported.imported()&&imported.auxiliaryGeometryRetained&&imported.solids==1,"file import retains auxiliary geometry");
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
