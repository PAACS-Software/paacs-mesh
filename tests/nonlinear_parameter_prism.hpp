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
#include <stdexcept>

namespace paacs::test {
// Exact unit prism with C(t)=(.5*t+.5*t*t,0,0) on one boundary edge.
// The physical edge is straight but its positive speed is nonuniform. Native
// pcurve correspondence flags and tolerances are left to the ordinary builder.
inline TopoDS_Shape nonlinearParameterPrism() {
    NCollection_Array1<gp_Pnt> poles(1,3);
    poles(1)=gp_Pnt(0,0,0);poles(2)=gp_Pnt(.25,0,0);poles(3)=gp_Pnt(1,0,0);
    occ::handle<Geom_BezierCurve> curve=new Geom_BezierCurve(poles);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(curve).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1,0,0),gp_Pnt(1,1,0)).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1,1,0),gp_Pnt(0,1,0)).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0,1,0),gp_Pnt(0,0,0)).Edge());
    if(!wire.IsDone())throw std::runtime_error("nonlinear-parameter prism wire construction failed");
    BRepBuilderAPI_MakeFace face(wire.Wire());
    if(!face.IsDone())throw std::runtime_error("nonlinear-parameter prism face construction failed");
    return BRepPrimAPI_MakePrism(face.Face(),gp_Vec(0,0,1)).Shape();
}
} // namespace paacs::test
