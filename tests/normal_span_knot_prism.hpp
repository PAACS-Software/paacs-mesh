// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace paacs::test {
// An owned unit prism with one straight edge represented by a cubic B-spline.
// Greville-abscissa poles reproduce C(u)=(u,0,0) on an ordinary [0,1] domain.
// The mandatory C2 knot split is adjacent to an artificial sampler grid cut;
// it changes the representation, not the physical geometry or its tolerances.
inline occ::handle<Geom_BSplineCurve> normalSpanKnotCurve() {
    const double k=std::nextafter(.5,std::numeric_limits<double>::infinity());
    const double x[]{0.,k/3.,2.*k/3.,(2.*k+1.)/3.,(k+2.)/3.,1.};
    NCollection_Array1<gp_Pnt> poles(1,6);
    for(int i=1;i<=6;++i)poles(i)=gp_Pnt(x[i-1],0.,0.);
    NCollection_Array1<double> knots(1,3);knots(1)=0.;knots(2)=k;knots(3)=1.;
    NCollection_Array1<int> multiplicities(1,3);multiplicities(1)=multiplicities(3)=4;multiplicities(2)=2;
    return new Geom_BSplineCurve(poles,knots,multiplicities,3,false);
}
// Two mandatory continuity knots have a three-ULP interval, too short for a
// distinct five-point integration stencil. Neither may be discarded as an
// artificial integration seed. The Native boundary remains a valid unit prism.
inline occ::handle<Geom_BSplineCurve> adjacentMandatoryKnotCurve() {
    double k=.5;
    for(unsigned i=0;i<3;++i)k=std::nextafter(k,std::numeric_limits<double>::infinity());
    const double expanded[]{0.,0.,0.,0.,.5,.5,k,k,1.,1.,1.,1.};
    NCollection_Array1<gp_Pnt> poles(1,8);
    for(int i=1;i<=8;++i)poles(i)=gp_Pnt((expanded[i]+expanded[i+1]+expanded[i+2])/3.,0.,0.);
    NCollection_Array1<double> knots(1,4);knots(1)=0.;knots(2)=.5;knots(3)=k;knots(4)=1.;
    NCollection_Array1<int> multiplicities(1,4);multiplicities(1)=multiplicities(4)=4;multiplicities(2)=multiplicities(3)=2;
    return new Geom_BSplineCurve(poles,knots,multiplicities,3,false);
}
inline TopoDS_Shape normalSpanKnotPrism(bool adjacentMandatory=false) {
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(adjacentMandatory?adjacentMandatoryKnotCurve():normalSpanKnotCurve()).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1,0,0),gp_Pnt(1,1,0)).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1,1,0),gp_Pnt(0,1,0)).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0,1,0),gp_Pnt(0,0,0)).Edge());
    if(!wire.IsDone())throw std::runtime_error("normal-span knot prism wire construction failed");
    BRepBuilderAPI_MakeFace face(wire.Wire());
    if(!face.IsDone())throw std::runtime_error("normal-span knot prism face construction failed");
    return BRepPrimAPI_MakePrism(face.Face(),gp_Vec(0,0,1)).Shape();
}
} // namespace paacs::test
