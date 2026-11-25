// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <BndLib_AddSurface.hxx>
#include <Bnd_Box.hxx>
#include <IntCurvesFace_Intersector.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Lin.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace paacs::geo::detail {

// A positive-only original-CAD witness. For a complete, transverse ray
// through a closed oriented shell, alternating exits and entries identify
// material at its origin. Both directions must establish that result.
// Ambiguous, on-boundary and outside queries use the ordinary classifier.
class OcctTrimmedRayWitness {
public:
    bool load(const TopoDS_Shape& solid,double tolerance){
        faces_.clear();
        if(solid.IsNull()||solid.ShapeType()!=TopAbs_SOLID||
           solid.Orientation()!=TopAbs_FORWARD||!std::isfinite(tolerance)||tolerance<0.)return false;
        try{
            for(TopoDS_Iterator child(solid);child.More();child.Next()){
                const auto& shell=child.Value();
                if(shell.ShapeType()!=TopAbs_SHELL||!BRep_Tool::IsClosed(shell)||
                   (shell.Orientation()!=TopAbs_FORWARD&&shell.Orientation()!=TopAbs_REVERSED))
                    return false;
            }
            std::vector<Face> faces;
            Bnd_Box body;
            double nativeTolerance=std::max(tolerance,Precision::Confusion());
            for(TopExp_Explorer explorer(solid,TopAbs_FACE);explorer.More();explorer.Next()){
                if(faces.size()>=128)return false;
                const auto face=TopoDS::Face(explorer.Current());
                if((face.Orientation()!=TopAbs_FORWARD&&face.Orientation()!=TopAbs_REVERSED)||
                   BRep_Tool::Surface(face).IsNull())return false;
                const auto faceTolerance=BRep_Tool::Tolerance(face);
                if(!std::isfinite(faceTolerance)||faceTolerance<0.)return false;
                nativeTolerance=std::max(nativeTolerance,faceTolerance);
                BRepAdaptor_Surface surface(face,true);
                // OCCT's Add uses analytic extrema for primitives and
                // pole hulls for Bezier and B-spline surfaces. Other types
                // may use sampled bounds, which cannot prove rejection.
                Bnd_Box bound;
                if(!surfaceBounds(surface,faceTolerance,bound))return false;
                const auto safe=box(bound);
                if(!safe)return false;
                body.Add(bound);
                faces.push_back({*safe,std::make_unique<IntCurvesFace_Intersector>(face,tolerance,true,true),surface.GetType()});
            }
            const auto bounds=box(body);
            if(!bounds)return false;
            const auto& lo=bounds->lo;
            const auto& hi=bounds->hi;
            const auto diagonal=std::hypot(hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]);
            const auto roundoff=std::max({std::abs(lo[0]),std::abs(lo[1]),std::abs(lo[2]),
                std::abs(hi[0]),std::abs(hi[1]),std::abs(hi[2]),diagonal})*
                (64.*std::numeric_limits<double>::epsilon());
            const auto margin=std::max(nativeTolerance,roundoff);
            if(!(diagonal>0.)||!std::isfinite(diagonal*3.)||
               !std::isfinite(margin)||margin>=diagonal*.25)return false;
            for(auto& face:faces){
                for(unsigned i=0;i<3;++i){
                    auto& left=face.bounds.lo[i];
                    auto& right=face.bounds.hi[i];
                    left=std::nextafter(left-margin,-std::numeric_limits<double>::infinity());
                    right=std::nextafter(right+margin,std::numeric_limits<double>::infinity());
                    if(!std::isfinite(left)||!std::isfinite(right))return false;
                }
            }
            bounds_=*bounds;
            faces_.swap(faces);
            margin_=margin;
            reach_=3.*diagonal;
            return true;
        }catch(const Standard_OutOfMemory&){throw std::bad_alloc{};}
         catch(const Standard_Failure&){return false;}
    }

    bool admits(const gp_Pnt& point,unsigned requiredPairs=2){
        if(faces_.empty()||!std::isfinite(point.X())||!std::isfinite(point.Y())||
           !std::isfinite(point.Z()))return false;
        const std::array coordinates{point.X(),point.Y(),point.Z()};
        for(unsigned i=0;i<3;++i)
            if(coordinates[i]<bounds_.lo[i]||coordinates[i]>bounds_.hi[i])return false;
        constexpr std::array<std::array<double,3>,3> directions{{
            {1.,.347,.113},{.113,1.,.347},{.347,.113,1.}}};
        try{
            unsigned pairs=0;
            for(const auto axis:directions){
                const auto forward=gp_Dir(axis[0],axis[1],axis[2]);
                if(!pairExits(point,forward))continue;
                if(++pairs==requiredPairs)return true;
            }
        }catch(const Standard_OutOfMemory&){throw std::bad_alloc{};}
         catch(const Standard_Failure&){return false;}
        return false;
    }

private:
    struct Box {std::array<double,3> lo,hi;};
    struct Face {Box bounds;std::unique_ptr<IntCurvesFace_Intersector> surface;GeomAbs_SurfaceType type;};
    struct Crossing {double distance;IntCurveSurface_TransitionOnCurve transition;};

    static bool surfaceBounds(const BRepAdaptor_Surface& surface,double tolerance,Bnd_Box& bounds){
        switch(surface.GetType()){
            case GeomAbs_Plane: case GeomAbs_Cylinder: case GeomAbs_Cone:
            case GeomAbs_Sphere: case GeomAbs_Torus:
            case GeomAbs_BezierSurface: case GeomAbs_BSplineSurface:
                BndLib_AddSurface::Add(surface,tolerance,bounds);
                return true;
            case GeomAbs_SurfaceOfExtrusion: {
                // Use the native location-aware basis curve and axis.
                // A finite extrusion is contained in the union bound
                // of the two translated basis-curve endpoint boxes.
                const auto curve=surface.BasisCurve();
                if(curve.IsNull()||curve->GetType()>GeomAbs_BSplineCurve||
                   surface.Trsf().ScaleFactor()!=1.)return false;
                const auto u0=surface.FirstUParameter(),u1=surface.LastUParameter();
                const auto v0=surface.FirstVParameter(),v1=surface.LastVParameter();
                if(!(u0<u1&&v0<v1)||
                   !std::isfinite(u0)||!std::isfinite(u1)||!std::isfinite(v0)||!std::isfinite(v1))return false;
                const auto axis=surface.Direction();
                Bnd_Box section;
                BndLib_Add3dCurve::Add(*curve,u0,u1,tolerance,section);
                const auto limits=box(section);
                if(!limits)return false;
                for(double parameter:{v0,v1}){
                    const auto x=axis.X()*parameter,y=axis.Y()*parameter,z=axis.Z()*parameter;
                    bounds.Update(limits->lo[0]+x,limits->lo[1]+y,limits->lo[2]+z,
                                  limits->hi[0]+x,limits->hi[1]+y,limits->hi[2]+z);
                }
                return true;
            }
            default:
                return false;
        }
    }

    static std::optional<Box> box(const Bnd_Box& source){
        if(source.IsVoid()||source.IsOpen()||source.IsWhole())return {};
        Box result;
        source.Get(result.lo[0],result.lo[1],result.lo[2],result.hi[0],result.hi[1],result.hi[2]);
        for(unsigned i=0;i<3;++i)
            if(!std::isfinite(result.lo[i])||!std::isfinite(result.hi[i])||result.lo[i]>result.hi[i])
                return {};
        return result;
    }

    bool intersects(const std::array<double,3>& point,const std::array<double,3>& direction,
                    const Box& bounds)const{
        double first=-margin_,last=reach_;
        for(unsigned i=0;i<3;++i){
            const auto a=(bounds.lo[i]-point[i])/direction[i];
            const auto b=(bounds.hi[i]-point[i])/direction[i];
            first=std::max(first,std::nextafter(std::min(a,b),-std::numeric_limits<double>::infinity()));
            last=std::min(last,std::nextafter(std::max(a,b),std::numeric_limits<double>::infinity()));
            if(first>last)return false;
        }
        return true;
    }

    bool validCrossings(std::vector<Crossing>& crossings)const{
        if(crossings.size()%2==0)return false;
        std::ranges::sort(crossings,std::less<>{},&Crossing::distance);
        double previous=0.;
        for(std::size_t i=0;i<crossings.size();++i){
            const auto& crossing=crossings[i];
            const auto expected=i%2==0?IntCurveSurface_Out:IntCurveSurface_In;
            if(crossing.transition!=expected||crossing.distance-previous<=margin_)return false;
            previous=crossing.distance;
        }
        return true;
    }

    bool pairExits(const gp_Pnt& point,const gp_Dir& direction){
        const gp_Lin line(point,direction);
        const std::array origin{point.X(),point.Y(),point.Z()};
        const std::array vector{direction.X(),direction.Y(),direction.Z()};
        const std::array negative{-direction.X(),-direction.Y(),-direction.Z()};
        std::array<std::size_t,128> negativeOnly{};std::size_t deferred=0;
        crossings_.clear();oppositeCrossings_.clear();
        const auto analytic=[&](std::size_t fi){
            switch(faces_[fi].type){
                case GeomAbs_Plane:case GeomAbs_Cylinder:case GeomAbs_Cone:
                case GeomAbs_Sphere:case GeomAbs_Torus:return true;
                default:return false;
            }
        };
        const auto directed=[&](std::size_t fi,bool positive){
            auto& surface=*faces_[fi].surface;
            const gp_Lin ray(point,positive?direction:direction.Reversed());
            surface.Perform(ray,-margin_,reach_);
            if(!surface.IsDone()||surface.IsParallel())return false;
            auto& hits=positive?crossings_:oppositeCrossings_;
            for(int i=1;i<=surface.NbPnt();++i){
                const auto distance=surface.WParameter(i);const auto transition=surface.Transition(i);
                if(!(distance>margin_&&distance<reach_)||surface.State(i)!=TopAbs_IN||
                   transition==IntCurveSurface_Tangent||hits.size()>=256)return false;
                hits.push_back({distance,transition});
            }
            return true;
        };
        const auto consume=[&](std::size_t fi,bool positive,bool negativeHalf){
            auto& surface=*faces_[fi].surface;
            surface.Perform(line,negativeHalf?-reach_:-margin_,positive?reach_:margin_);
            if(!surface.IsDone()||surface.IsParallel())return false;
            for(int i=1;i<=surface.NbPnt();++i){
                const auto distance=surface.WParameter(i),inverse=-distance;
                if(!std::isfinite(distance))return false;
                const auto transition=surface.Transition(i);
                const auto add=[&](double t,IntCurveSurface_TransitionOnCurve tr,std::vector<Crossing>& hits){
                    if(!(t>margin_&&t<reach_)||surface.State(i)!=TopAbs_IN||
                       tr==IntCurveSurface_Tangent||hits.size()>=256)return false;
                    hits.push_back({t,tr});return true;
                };
                if(positive&&distance>=-margin_&&distance<=reach_&&!add(distance,transition,crossings_))return false;
                if(negativeHalf&&inverse>=-margin_&&inverse<=reach_&&!add(inverse,(transition==IntCurveSurface_Out?IntCurveSurface_In:transition==IntCurveSurface_In?IntCurveSurface_Out:transition),oppositeCrossings_))return false;
            }
            return true;
        };
        for(std::size_t fi=0;fi<faces_.size();++fi){
            const bool positive=intersects(origin,vector,faces_[fi].bounds);
            const bool negativeHalf=intersects(origin,negative,faces_[fi].bounds);
            const bool share=positive&&negativeHalf&&analytic(fi);
            if(positive){const bool accepted=share?consume(fi,true,true):directed(fi,true);if(!accepted)return false;}
            if(negativeHalf&&!share)negativeOnly[deferred++]=fi;
        }
        // Preserve the old forward-half early refusal. Negative hits from
        // shared faces are retained, but negative-only faces are not queried
        // when the first half cannot provide a complete certificate.
        if(!validCrossings(crossings_))return false;
        for(std::size_t i=0;i<deferred;++i)
            if(!directed(negativeOnly[i],false))return false;
        return validCrossings(oppositeCrossings_);
    }

    std::vector<Face> faces_;
    std::vector<Crossing> crossings_,oppositeCrossings_;
    Box bounds_{};
    double margin_=0.,reach_=0.;
};

}
