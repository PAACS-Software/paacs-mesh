// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <Geom_BezierCurve.hxx>
#include <NCollection_Array1.hxx>
#include <TopoDS.hxx>
#include <stdexcept>

namespace paacs::test {
// The exact curved wire is simple. Joining only each curve's endpoints makes
// its two Bezier chords cross; one native midpoint per curve resolves the
// polygon topology. This is owned geometry, not a downloaded CAD derivative.
inline TopoDS_Shape narrowTrimPrism() {
    const auto bezier=[](gp_Pnt a,gp_Pnt b,gp_Pnt c){
        NCollection_Array1<gp_Pnt> poles(1,3);
        poles(1)=a;poles(2)=b;poles(3)=c;
        return occ::handle<Geom_BezierCurve>(new Geom_BezierCurve(poles));
    };
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(bezier(gp_Pnt(-1,0,0),gp_Pnt(0,2,0),gp_Pnt(1,0,0))).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1,0,0),gp_Pnt(0,-.1,0)).Edge());
    // Reverse this native curve's use to exercise decreasing parameters too.
    wire.Add(TopoDS::Edge(BRepBuilderAPI_MakeEdge(
        bezier(gp_Pnt(0,.1,0),gp_Pnt(.04,0,0),gp_Pnt(0,-.1,0))).Edge().Reversed()));
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0,.1,0),gp_Pnt(-1,0,0)).Edge());
    if(!wire.IsDone())throw std::runtime_error("narrow-trim prism wire construction failed");
    BRepBuilderAPI_MakeFace face(wire.Wire());
    if(!face.IsDone())throw std::runtime_error("narrow-trim prism face construction failed");
    return BRepPrimAPI_MakePrism(face.Face(),gp_Vec(0,0,1)).Shape();
}
} // namespace paacs::test
