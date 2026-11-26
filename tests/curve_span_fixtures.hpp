// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <array>

namespace paacs::test {
// Owned C0 parameterizations of the unit straight segment. The first has
// constant one-sided speeds 0.5 and 1.5; the second has second derivatives 1 and 2.
// Both have exact arclength 1 and an ordinary [0,1] parameter domain.
inline occ::handle<Geom_BSplineCurve> oneSidedSpeedCurve(bool varyingSecond=false){
    const auto x=varyingSecond?std::array{0.,1./24.,1./8.,.25,5./12.,2./3.,1.}:
        std::array{0.,1./12.,1./6.,.25,.5,.75,1.};
    NCollection_Array1<gp_Pnt> poles(1,7);
    for(int i=1;i<=7;++i)poles(i)=gp_Pnt(x[i-1],0.,0.);
    NCollection_Array1<double> knots(1,3);knots(1)=0.;knots(2)=.5;knots(3)=1.;
    NCollection_Array1<int> mult(1,3);mult(1)=mult(3)=4;mult(2)=3;
    return new Geom_BSplineCurve(poles,knots,mult,3,false);
}
// Ordinary smooth, non-collinear B-spline with multiple native knot intervals
// but one C2 span. Its native evaluation must not receive local overrides.
inline occ::handle<Geom_BSplineCurve> smoothSpanCurve(){
    NCollection_Array1<gp_Pnt> poles(1,7);
    const std::array x{0.,1./12.,.25,.5,.75,11./12.,1.};
    const std::array y{0.,.3,-.15,.4,-.2,.25,0.};
    for(int i=1;i<=7;++i)poles(i)=gp_Pnt(x[i-1],y[i-1],0.);
    NCollection_Array1<double> knots(1,5);for(int i=1;i<=5;++i)knots(i)=double(i-1)/4.;
    NCollection_Array1<int> mult(1,5);mult(1)=mult(5)=4;mult(2)=mult(3)=mult(4)=1;
    return new Geom_BSplineCurve(poles,knots,mult,3,false);
}
} // namespace paacs::test
