// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/geometry/Kernel.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <Standard_Failure.hxx>
#include <Precision.hxx>
#include <Standard_OutOfMemory.hxx>
#include <algorithm>
#include <cmath>
#include <map>
#include <new>
#include <utility>

namespace paacs::geo::detail {

// One immutable original CAD face per chart. Cache only trim decisions.
// The mesher still checks physical chord distances using current geometry.
class OcctChartTrim {
public:
    OcctChartTrim(const GFace& original,const TopoDS_Face& face,double tolerance)
        : original_(original),face_(face) {
        const BRepAdaptor_Surface surface(face,false);
        tolerance=std::max(tolerance,Precision::Confusion());
        uvTolerance_=std::min(surface.UResolution(tolerance),surface.VResolution(tolerance));
    }

    bool contains(Vec2 parameter) {
        if(!std::isfinite(parameter.u)||!std::isfinite(parameter.v))return false;
        const auto key=std::pair{parameter.u,parameter.v};
        if(const auto found=decisions_.find(key);found!=decisions_.end())return found->second;
        const auto accepted=classify(parameter);
        constexpr std::size_t maximumDecisions=4096;
        if(decisions_.size()<maximumDecisions)decisions_.try_emplace(key,accepted);
        return accepted;
    }

private:
    bool classify(Vec2 parameter) const {
        if(!(uvTolerance_>0)||!std::isfinite(uvTolerance_))return false;
        try {
            // OCCT uses pcurve bounds to avoid intersections with irrelevant
            // trim edges. It retains the original wires and UV tolerance.
            const BRepClass_FaceClassifier classifier(face_,
                gp_Pnt2d(parameter.u,parameter.v),uvTolerance_,true);
            const auto state=classifier.State();
            if(state==TopAbs_IN||state==TopAbs_ON)return true;
        } catch(const Standard_OutOfMemory&) {
            throw std::bad_alloc{};
        } catch(const Standard_Failure&) {
        }
        // Ambiguous or refused accelerated queries retain the stock behavior.
        return original_.containsParam(parameter);
    }

    const GFace& original_;
    TopoDS_Face face_;
    double uvTolerance_;
    std::map<std::pair<double,double>,bool> decisions_;
};

} // namespace paacs::geo::detail
