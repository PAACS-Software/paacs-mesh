// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctModelInternal.hpp"
#include "OcctSpatialSizing.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace paacs::geo::detail {

using OcctCurveUses = std::map<int, std::vector<std::pair<int, std::uint64_t>>>;
using OcctVertexEdges = std::map<int, std::vector<int>>;

struct OcctCurveSpacing {
    unsigned pointsPerTurn = 32;
    double maximumChordError = 0;
};

template<class Tick>
std::vector<double> sampleOcctCurveMetric(const OcctModel& model, int tag, double sizeCap,
    const OcctCurveUses& uses, const OcctVertexEdges& incident, Tick tick,
    OcctCurveSpacing spacing = {},const OcctSpatialSizing* sizeField = nullptr) {
    constexpr std::size_t intervals = 256;
    const auto localCurvature = [&](int edge, double t) {
        tick(1);
        double curvature = 0;
        const auto include = [&](double value) {
            if (std::isfinite(value)) curvature = std::max(curvature, value);
        };
        include(model.edge(edge)->curvatureAt(t));
        for (const auto& [face, use] : uses.at(edge)) {
            tick(2);
            Vec2 uv;
            if (model.parameterOnUse(use, t, uv)) include(model.face(face)->curvatureAt(uv));
        }
        return curvature;
    };
    const auto vertexCurvature = [&](int vertex) {
        double curvature = 0;
        for (int edgeTag : incident.at(vertex)) {
            tick(1);
            const auto* edge = model.edge(edgeTag);
            double first, last;
            edge->parameterBounds(first, last);
            if (edge->startVertexTag() == vertex)
                curvature = std::max(curvature, localCurvature(edgeTag, first));
            if (edge->endVertexTag() == vertex)
                curvature = std::max(curvature, localCurvature(edgeTag, last));
        }
        return curvature;
    };
    const auto* edge = model.edge(tag);
    double first, last;
    edge->parameterBounds(first, last);
    BRepAdaptor_Curve curve(OcctModelAccess::edge(model, tag));
    std::array<double, intervals+1> t{}, speed{}, size{}, metric{};
    std::array<double, intervals> length{};
    for (std::size_t i = 0; i <= intervals; ++i) {
        tick(1);
        t[i] = i == intervals ? last : std::lerp(first, last, double(i)/intervals);
        gp_Pnt point;
        gp_Vec derivative;
        curve.D1(t[i], point, derivative);
        speed[i] = derivative.Magnitude();
        auto curvature = localCurvature(tag, t[i]);
        if (i == 0) curvature = std::max(curvature, vertexCurvature(edge->startVertexTag()));
        if (i == intervals) curvature = std::max(curvature, vertexCurvature(edge->endVertexTag()));
        size[i] = curvature > 0 ?
            std::min(sizeCap, 2*std::numbers::pi/(spacing.pointsPerTurn*curvature)) : sizeCap;
        if (curvature > 0 && spacing.maximumChordError > 0)
            size[i] = std::min(size[i], std::sqrt(8*spacing.maximumChordError/curvature));
        if (sizeField) {
            tick(sizeField->queryWork());
            size[i] = std::min(size[i], sizeField->eval({point.X(), point.Y(), point.Z()}));
        }
        if (!(std::isfinite(speed[i]) && std::isfinite(size[i]) && size[i] > 0)) return {};
        if (i) length[i-1] = (t[i]-t[i-1])*(speed[i]+speed[i-1])*.5;
    }
    const double slope = std::sqrt(1.8)-1;
    for (std::size_t i = 0; i < intervals; ++i)
        size[i+1] = std::min(size[i+1], size[i]+slope*length[i]);
    for (std::size_t i = intervals; i > 0; --i)
        size[i-1] = std::min(size[i-1], size[i]+slope*length[i-1]);
    for (std::size_t i = 1; i <= intervals; ++i)
        metric[i] = metric[i-1]+(t[i]-t[i-1])*(speed[i]/size[i]+speed[i-1]/size[i-1])*.5;
    if (!(std::isfinite(metric.back()) && metric.back() > 0 && metric.back() <= 256)) return {};
    const auto minimum = spacing.maximumChordError > 0 && edge->startVertexTag()==edge->endVertexTag()?3u:1u;
    const auto count = std::max(minimum, static_cast<unsigned>(std::ceil(metric.back()-.01)));
    std::vector<double> schedule{first};
    for (unsigned i = 1; i < count; ++i) {
        const auto target = metric.back()*i/count;
        const auto j = std::size_t(std::lower_bound(metric.begin(), metric.end(), target)-metric.begin());
        const auto parameter = std::lerp(t[j-1], t[j], (target-metric[j-1])/(metric[j]-metric[j-1]));
        if (!(parameter > schedule.back() && parameter < last)) return {};
        schedule.push_back(parameter);
    }
    schedule.push_back(last);
    return schedule;
}

struct OcctChartMetric {
    double u, shear, v;

    Vec2 physical(Vec2 uv) const { return {u*uv.u+shear*uv.v,v*uv.v}; }
    Vec2 parameter(Vec2 p) const {
        const double y=p.v/v;
        return {(p.u-shear*y)/u,y};
    }
};

inline std::optional<OcctChartMetric> occtChartMetric(double a,double b,double c) {
    if (!(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(c)&&a>0&&c>0)) return {};
    const auto u=std::sqrt(a),shear=b/u,v=std::sqrt(c-shear*shear);
    if (!(std::isfinite(v)&&v>0)) return {};
    return OcctChartMetric{u,shear,v};
}

inline std::optional<OcctChartMetric> occtChartMetricAt(const BRepAdaptor_Surface& surface,Vec2 uv) {
    gp_Pnt point;
    gp_Vec u,v;
    surface.D1(uv.u,uv.v,point,u,v);
    return occtChartMetric(u.SquareMagnitude(),u.Dot(v),v.SquareMagnitude());
}

template<class Tick>
bool metricOcctChart(const OcctModel& model, int face, std::span<const Vec2> samples,
                     Tick tick, std::vector<Vec2>& coordinates) {
    BRepAdaptor_Surface surface(OcctModelAccess::face(model, face));
    double a = 0, b = 0, c = 0;
    for (const auto uv : samples) {
        tick(1);
        gp_Pnt point;
        gp_Vec u, v;
        surface.D1(uv.u, uv.v, point, u, v);
        a += u.SquareMagnitude();
        b += u.Dot(v);
        c += v.SquareMagnitude();
    }
    if (samples.empty() || !(std::isfinite(a) && std::isfinite(b) && std::isfinite(c) &&
                             a > 0 && c > 0 && a*c > b*b)) return false;
    const auto scaleU=std::sqrt(a/samples.size());
    const auto shearV=b/(samples.size()*scaleU);
    const auto scaleV=std::sqrt(c/samples.size()-shearV*shearV);
    if (!(std::isfinite(scaleV)&&scaleV>0)) return false;
    const OcctChartMetric metric{scaleU,shearV,scaleV};
    coordinates.reserve(samples.size());
    for (const auto uv : samples) {
        const Vec2 transformed=metric.physical(uv);
        if (!(std::isfinite(transformed.u) && std::isfinite(transformed.v))) return false;
        coordinates.push_back(transformed);
    }
    return true;
}

template<class Triangles, class Tick>
std::optional<std::vector<int>> labelOcctChartRegions(const Triangles& triangles,
    const std::set<std::array<Index,2>>& constraints, Tick tick) {
    using Edge = std::array<Index,2>;
    std::map<Edge, std::vector<std::size_t>> owners;
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        tick(3);
        for (unsigned side = 0; side < 3; ++side) {
            const auto [lo,hi] = std::minmax(triangles[i].v[side], triangles[i].v[(side+1)%3]);
            owners[{lo,hi}].push_back(i);
        }
    }
    for (const auto& edge : constraints) {
        tick(1);
        if (!owners.contains(edge)) return {};
    }
    std::vector<std::vector<std::pair<std::size_t, int>>> neighbours(triangles.size());
    std::vector<int> parity(triangles.size(), -1);
    std::vector<std::size_t> pending;
    const auto assign = [&](std::size_t i, int label) {
        if (parity[i] < 0) {
            parity[i] = label;
            pending.push_back(i);
        }
        return parity[i] == label;
    };
    for (const auto& [edge, uses] : owners) {
        tick(1);
        const int crosses = constraints.contains(edge);
        if (uses.size() == 1) {
            if (!assign(uses.front(), crosses)) return {};
        } else if (uses.size() == 2) {
            neighbours[uses[0]].push_back({uses[1], crosses});
            neighbours[uses[1]].push_back({uses[0], crosses});
        } else return {};
    }
    for (std::size_t next = 0; next < pending.size(); ++next)
        for (const auto& [other, crosses] : neighbours[pending[next]]) {
            tick(1);
            if (!assign(other, parity[pending[next]] ^ crosses)) return {};
        }
    if (pending.size() != triangles.size()) return {};
    return parity;
}

} // namespace paacs::geo::detail
