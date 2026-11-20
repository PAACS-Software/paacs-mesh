// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private OCCT containment helper. Recover an OUT/UNKNOWN answer only when
// independent trimmed-face rays agree, with alternating entry/exit hits.
// IN/ON retain the ordinary classifier fast path. No tolerance is enlarged.
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <gp_Lin.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>

namespace paacs::geo::detail {
class OcctSolidClassifier {
public:
    struct Metrics {std::uint64_t queries=0,rayQueries=0,recoveredInside=0,uncertain=0;};
    OcctSolidClassifier()=default;
    explicit OcctSolidClassifier(const TopoDS_Shape& shape){Load(shape);}
    void Load(const TopoDS_Shape& shape){
        shape_=shape;primary_.Load(shape);state_=TopAbs_UNKNOWN;
        ray_.reset();bounds_.reset();rayTolerance_=-1;
    }
    void PerformInfinitePoint(double tolerance){primary_.PerformInfinitePoint(tolerance);state_=primary_.State();}
    void Perform(const gp_Pnt& p,double tolerance){
        ++metrics_.queries;state_=TopAbs_UNKNOWN;
        if(!std::isfinite(p.X())||!std::isfinite(p.Y())||!std::isfinite(p.Z())||!std::isfinite(tolerance)||tolerance<0||shape_.IsNull())return;
        primary_.Perform(p,tolerance);state_=primary_.State();
        if(state_==TopAbs_IN||state_==TopAbs_ON)return;
        try {
            if(!bounds_){Bnd_Box box;BRepBndLib::AddOptimal(shape_,box,false,false);if(box.IsVoid()||box.IsOpen())return;std::array<double,6> values;box.Get(values[0],values[1],values[2],values[3],values[4],values[5]);bounds_=values;}
            const auto b=*bounds_;
            if(p.X()<b[0]-tolerance||p.X()>b[3]+tolerance||p.Y()<b[1]-tolerance||p.Y()>b[4]+tolerance||p.Z()<b[2]-tolerance||p.Z()>b[5]+tolerance){state_=TopAbs_OUT;return;}
            const auto diagonal=std::hypot(b[3]-b[0],b[4]-b[1],b[5]-b[2]);if(!(diagonal>0&&std::isfinite(diagonal)))return;
            if(!ray_){ray_=std::make_unique<IntCurvesFace_ShapeIntersector>();}
            if(rayTolerance_!=tolerance){ray_->Load(shape_,tolerance);rayTolerance_=tolerance;}
            const auto proximity=std::max(tolerance,diagonal*1e-12);
            constexpr std::array<std::array<double,3>,8> directions{{
                {1,.347,.113},{-1,-.347,-.113},{.113,1,.347},{-.113,-1,-.347},
                {.347,.113,1},{-.347,-.113,-1},{.517,1,.193},{-.517,-1,-.193}}};
            std::optional<bool> consensus;unsigned clear=0,pairs=0;bool previousClear=false;
            for(unsigned direction=0;direction<directions.size();++direction){
                const auto d=directions[direction];++metrics_.rayQueries;
                ray_->Perform(gp_Lin(p,gp_Dir(d[0],d[1],d[2])),0.,diagonal*3.);
                bool ambiguous=!ray_->IsDone();unsigned crossings=0;double previous=-1;
                std::optional<IntCurveSurface_TransitionOnCurve> first;
                if(ray_->IsDone()){
                    ray_->SortResult();
                    for(int i=1;i<=ray_->NbPnt();++i){const auto t=ray_->WParameter(i);const auto transition=ray_->Transition(i);
                        if(!std::isfinite(t)||t<=proximity||ray_->State(i)!=TopAbs_IN||transition==IntCurveSurface_Tangent||(previous>=0&&t-previous<=proximity)){ambiguous=true;break;}
                        if(!first)first=transition;
                        const auto expected=crossings%2?(*first==IntCurveSurface_Out?IntCurveSurface_In:IntCurveSurface_Out):*first;
                        if(transition!=expected){ambiguous=true;break;}
                        previous=t;++crossings;
                    }
                }
                if(!ambiguous&&crossings&&((crossings%2&&*first!=IntCurveSurface_Out)||(!(crossings%2)&&*first!=IntCurveSurface_In)))ambiguous=true;
                if(!ambiguous){const bool inside=bool(crossings%2);
                    // The primary OUT answer can only be overturned by an
                    // unanimous inside consensus. A clear outside ray rules
                    // that out; later agreement, disagreement, or uncertainty
                    // would all retain OUT. Inside recovery still uses every
                    // direction and the same entry/exit and pair requirements.
                    if(!inside&&state_==TopAbs_OUT)return;
                    if(consensus&&*consensus!=inside){++metrics_.uncertain;return;}consensus=inside;++clear;if(direction%2&&previousClear)++pairs;}
                previousClear=!ambiguous;
            }
            if(clear<4||pairs<2||!consensus){++metrics_.uncertain;return;}
            if(*consensus){state_=TopAbs_IN;++metrics_.recoveredInside;}else state_=TopAbs_OUT;
        } catch(const Standard_OutOfMemory&){throw;}
          catch(const Standard_Failure&){++metrics_.uncertain;}
    }
    TopAbs_State State()const noexcept{return state_;}
    const Metrics& metrics()const noexcept{return metrics_;}
private:
    TopoDS_Shape shape_;
    BRepClass3d_SolidClassifier primary_;
    std::unique_ptr<IntCurvesFace_ShapeIntersector> ray_;
    std::optional<std::array<double,6>> bounds_;
    double rayTolerance_=-1;
    TopAbs_State state_=TopAbs_UNKNOWN;
    Metrics metrics_;
};
}
