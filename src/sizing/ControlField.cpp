// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/ControlField.hpp"
#include <algorithm>
#include <cmath>
#include <new>

namespace paacs::sizing {
namespace {
constexpr double coordinateLimit = 1e100;
constexpr double minimumGeometryLength = 1e-100;
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
constexpr double eps = std::numeric_limits<double>::epsilon();
double norm(Vec3 p) noexcept {
    return std::hypot(p.x, p.y, p.z);
}
double magnitude(Vec3 p) noexcept {
    return std::max({std::abs(p.x), std::abs(p.y), std::abs(p.z)});
}
bool coordinate(Vec3 p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) &&
           magnitude(p) <= coordinateLimit;
}
bool ordered(Vec3 a, Vec3 b) noexcept {
    return a.x <= b.x && a.y <= b.y && a.z <= b.z;
}
bool bounded(double x) noexcept {
    return std::isfinite(x) && std::abs(x) <= coordinateLimit;
}
unsigned arity(FieldOperation op) noexcept {
    switch (op) {
    case FieldOperation::Minimum:
    case FieldOperation::Maximum:
        return 2;
    case FieldOperation::ThresholdLinear:
    case FieldOperation::ThresholdSmooth:
    case FieldOperation::Scale:
    case FieldOperation::Clamp:
        return 1;
    default:
        return 0;
    }
}
double distanceBox(Vec3 p, Vec3 lo, Vec3 hi) noexcept {
    return norm({std::max({lo.x - p.x, p.x - hi.x, 0.0}), std::max({lo.y - p.y, p.y - hi.y, 0.0}),
                 std::max({lo.z - p.z, p.z - hi.z, 0.0})});
}
double distanceSegment(Vec3 p, Vec3 a, Vec3 b) noexcept {
    const Vec3 d = b - a;
    const double length = norm(d);
    const Vec3 unit = d * (1 / length);
    const double t = std::clamp((p - a).dot(unit) / length, 0.0, 1.0);
    return norm(p - (a + d * t));
}
double distanceTriangle(Vec3 p, const FieldNode &n) noexcept {
    const double length = std::max(norm(n.b - n.a), norm(n.c - n.a));
    const Vec3 u = (n.b - n.a) * (1 / length), v = (n.c - n.a) * (1 / length);
    const Vec3 normal = u.cross(v), w = (p - n.a) * (1 / length);
    const double denominator = normal.dot(normal);
    const double beta = w.cross(v).dot(normal) / denominator;
    const double gamma = u.cross(w).dot(normal) / denominator;
    if (beta >= 0 && gamma >= 0 && beta + gamma <= 1)
        return std::abs((p - n.a).dot(normal * (1 / norm(normal))));
    return std::min(
        {distanceSegment(p, n.a, n.b), distanceSegment(p, n.b, n.c), distanceSegment(p, n.c, n.a)});
}
double primitive(const FieldNode &n, Vec3 p) noexcept {
    switch (n.operation) {
    case FieldOperation::Constant:
        return n.parameters[0];
    case FieldOperation::DistancePoint:
        return norm(p - n.a);
    case FieldOperation::DistanceSegment:
        return distanceSegment(p, n.a, n.b);
    case FieldOperation::DistanceTriangle:
        return distanceTriangle(p, n);
    case FieldOperation::DistanceBox:
        return distanceBox(p, n.a, n.b);
    case FieldOperation::DistanceBall:
        return std::max(0.0, norm(p - n.a) - n.parameters[0]);
    default:
        return nan;
    }
}
double transform(const FieldNode &n, double a, double b = 0) noexcept {
    const auto &p = n.parameters;
    switch (n.operation) {
    case FieldOperation::Minimum:
        return std::min(a, b);
    case FieldOperation::Maximum:
        return std::max(a, b);
    case FieldOperation::Scale:
        return a * p[0];
    case FieldOperation::Clamp:
        return std::clamp(a, p[0], p[1]);
    case FieldOperation::ThresholdLinear:
    case FieldOperation::ThresholdSmooth: {
        if (a <= p[0])
            return p[2];
        if (a >= p[1])
            return p[3];
        double t = (a - p[0]) / (p[1] - p[0]);
        if (n.operation == FieldOperation::ThresholdSmooth)
            t = t * t * (3 - 2 * t);
        return std::lerp(p[2], p[3], t);
    }
    default:
        return nan;
    }
}
FieldRange outward(double lo, double hi, double error) noexcept {
    return {std::max(0.0, std::nextafter(lo - error, -INFINITY)),
            std::nextafter(hi + error, INFINITY)};
}
FieldRange primitiveRange(const FieldNode &n, Vec3 lo, Vec3 hi) noexcept {
    if (n.operation == FieldOperation::Constant)
        return {n.parameters[0], n.parameters[0]};
    double lower = 0, upper = 0;
    if (n.operation == FieldOperation::DistancePoint ||
        n.operation == FieldOperation::DistanceBall) {
        lower = distanceBox(n.a, lo, hi);
        if (n.operation == FieldOperation::DistanceBall)
            lower = std::max(0.0, lower - n.parameters[0]);
    } else if (n.operation == FieldOperation::DistanceBox) {
        lower = norm({std::max({n.a.x - hi.x, lo.x - n.b.x, 0.0}),
                      std::max({n.a.y - hi.y, lo.y - n.b.y, 0.0}),
                      std::max({n.a.z - hi.z, lo.z - n.b.z, 0.0})});
    } else {
        Vec3 a{std::min(n.a.x, n.b.x), std::min(n.a.y, n.b.y), std::min(n.a.z, n.b.z)};
        Vec3 b{std::max(n.a.x, n.b.x), std::max(n.a.y, n.b.y), std::max(n.a.z, n.b.z)};
        if (n.operation == FieldOperation::DistanceTriangle) {
            a = {std::min(a.x, n.c.x), std::min(a.y, n.c.y), std::min(a.z, n.c.z)};
            b = {std::max(b.x, n.c.x), std::max(b.y, n.c.y), std::max(b.z, n.c.z)};
        }
        const double boxGap =
            norm({std::max({a.x - hi.x, lo.x - b.x, 0.0}), std::max({a.y - hi.y, lo.y - b.y, 0.0}),
                  std::max({a.z - hi.z, lo.z - b.z, 0.0})});
        const Vec3 center = lo + (hi - lo) * .5;
        lower = std::max(boxGap, primitive(n, center) - norm((hi - lo) * .5));
    }
    // Distance to each supported primitive is convex: its maximum over a box
    // is bounded by its eight vertices. Lower bounds also cover the interior.
    for (unsigned mask = 0; mask < 8; ++mask)
        upper = std::max(upper, primitive(n, {mask & 1 ? hi.x : lo.x, mask & 2 ? hi.y : lo.y,
                                              mask & 4 ? hi.z : lo.z}));
    const double scale =
        std::max({magnitude(lo), magnitude(hi), magnitude(n.a), magnitude(n.b), magnitude(n.c),
                  n.operation == FieldOperation::DistanceBall ? n.parameters[0] : 0, upper});
    double conditioning = 1;
    if (n.operation == FieldOperation::DistanceTriangle) {
        const double length = std::max(norm(n.b - n.a), norm(n.c - n.a));
        const auto normal = ((n.b - n.a) * (1 / length)).cross((n.c - n.a) * (1 / length));
        conditioning = 1 / normal.dot(normal);
    }
    return outward(lower, upper, 128 * eps * scale * std::max(1.0, conditioning));
}
bool cancelled(const ControlFieldSpecification &s) noexcept {
    return s.cancellation && s.cancellation->load(std::memory_order_relaxed);
}
ControlFieldReason validateNode(const FieldNode &n) noexcept {
    if (static_cast<unsigned>(n.operation) > static_cast<unsigned>(FieldOperation::Clamp))
        return ControlFieldReason::InvalidNode;
    if (!coordinate(n.a) || !coordinate(n.b) || !coordinate(n.c))
        return ControlFieldReason::InvalidNode;
    for (double p : n.parameters)
        if (!bounded(p))
            return ControlFieldReason::InvalidNode;
    const auto &p = n.parameters;
    switch (n.operation) {
    case FieldOperation::Constant:
        if (p[0] < 0)
            return ControlFieldReason::InvalidNode;
        break;
    case FieldOperation::DistanceSegment:
        if (norm(n.b - n.a) < minimumGeometryLength)
            return ControlFieldReason::DegenerateGeometry;
        break;
    case FieldOperation::DistanceTriangle: {
        const double scale = std::max(norm(n.b - n.a), norm(n.c - n.a));
        if (scale < minimumGeometryLength ||
            norm(((n.b - n.a) * (1 / scale)).cross((n.c - n.a) * (1 / scale))) < 1e-10)
            return ControlFieldReason::DegenerateGeometry;
        break;
    }
    case FieldOperation::DistanceBox:
        if (!ordered(n.a, n.b))
            return ControlFieldReason::InvalidNode;
        break;
    case FieldOperation::DistanceBall:
        if (p[0] < 0)
            return ControlFieldReason::InvalidNode;
        break;
    case FieldOperation::ThresholdLinear:
    case FieldOperation::ThresholdSmooth:
        if (!(p[0] >= 0 && p[1] > p[0] && p[2] > 0 && p[3] > 0))
            return ControlFieldReason::InvalidNode;
        break;
    case FieldOperation::Scale:
        if (!(p[0] > 0))
            return ControlFieldReason::InvalidNode;
        break;
    case FieldOperation::Clamp:
        if (!(p[0] >= 0 && p[1] >= p[0]))
            return ControlFieldReason::InvalidNode;
        break;
    default:
        break;
    }
    return ControlFieldReason::None;
}
} // namespace
namespace field {
FieldNode constant(double v) noexcept {
    FieldNode n;
    n.parameters[0] = v;
    return n;
}
FieldNode distancePoint(Vec3 p) noexcept {
    FieldNode n;
    n.operation = FieldOperation::DistancePoint;
    n.a = p;
    return n;
}
FieldNode distanceSegment(Vec3 a, Vec3 b) noexcept {
    auto n = distancePoint(a);
    n.operation = FieldOperation::DistanceSegment;
    n.b = b;
    return n;
}
FieldNode distanceTriangle(Vec3 a, Vec3 b, Vec3 c) noexcept {
    auto n = distanceSegment(a, b);
    n.operation = FieldOperation::DistanceTriangle;
    n.c = c;
    return n;
}
FieldNode distanceBox(Vec3 a, Vec3 b) noexcept {
    auto n = distanceSegment(a, b);
    n.operation = FieldOperation::DistanceBox;
    return n;
}
FieldNode distanceBall(Vec3 c, double r) noexcept {
    auto n = distancePoint(c);
    n.operation = FieldOperation::DistanceBall;
    n.parameters[0] = r;
    return n;
}
FieldNode threshold(FieldId id, double d0, double d1, double h0, double h1, bool smooth) noexcept {
    FieldNode n;
    n.operation = smooth ? FieldOperation::ThresholdSmooth : FieldOperation::ThresholdLinear;
    n.first = id;
    n.parameters = {d0, d1, h0, h1};
    return n;
}
FieldNode minimum(FieldId a, FieldId b) noexcept {
    FieldNode n;
    n.operation = FieldOperation::Minimum;
    n.first = a;
    n.second = b;
    return n;
}
FieldNode maximum(FieldId a, FieldId b) noexcept {
    auto n = minimum(a, b);
    n.operation = FieldOperation::Maximum;
    return n;
}
FieldNode scale(FieldId a, double factor) noexcept {
    FieldNode n;
    n.operation = FieldOperation::Scale;
    n.first = a;
    n.parameters[0] = factor;
    return n;
}
FieldNode clamp(FieldId a, double lo, double hi) noexcept {
    FieldNode n;
    n.operation = FieldOperation::Clamp;
    n.first = a;
    n.parameters = {lo, hi, 0, 0};
    return n;
}
} // namespace field

double ControlField::eval(const Vec3 &point) const noexcept {
    if (!valid() || !coordinate(point))
        return nan;
    std::array<double, maximumControlFieldNodes> values;
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        const auto &n = nodes_[i];
        const auto count = arity(n.operation);
        values[i] = count ? transform(n, values[n.first], count == 2 ? values[n.second] : 0)
                          : primitive(n, point);
        if (!std::isfinite(values[i]))
            return nan;
    }
    return std::clamp(values[nodes_.size() - 1], minimum_, maximum_) * factor_;
}
FieldRange ControlField::rangeOnBox(const Vec3 &lo, const Vec3 &hi) const noexcept {
    if (!valid() || !coordinate(lo) || !coordinate(hi) || !ordered(lo, hi))
        return {nan, nan};
    std::array<FieldRange, maximumControlFieldNodes> values;
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        const auto &n = nodes_[i];
        const auto count = arity(n.operation);
        if (!count)
            values[i] = primitiveRange(n, lo, hi);
        else {
            const auto a = values[n.first];
            const auto b = count == 2 ? values[n.second] : FieldRange{};
            const double x = transform(n, a.lower, b.lower), y = transform(n, a.upper, b.upper);
            values[i] = x == y ? FieldRange{x, x}
                               : outward(std::min(x, y), std::max(x, y),
                                         16 * eps * std::max(std::abs(x), std::abs(y)));
        }
        if (!std::isfinite(values[i].lower) || !std::isfinite(values[i].upper))
            return {nan, nan};
    }
    const auto r = values[nodes_.size() - 1];
    const double loValue = std::clamp(r.lower, minimum_, maximum_) * factor_;
    const double hiValue = std::clamp(r.upper, minimum_, maximum_) * factor_;
    // Preserve exact envelopes, particularly constant fields at dyadic sizes.
    if (loValue == hiValue)
        return {loValue, hiValue};
    return {std::max(minimumSize(), std::nextafter(loValue, -INFINITY)),
            std::min(maximumSize(), std::nextafter(hiValue, INFINITY))};
}
double ControlField::minimumOnBox(const Vec3 &lo, const Vec3 &hi) const noexcept {
    return rangeOnBox(lo, hi).lower;
}
ControlFieldOutcome compileControlField(const ControlFieldSpecification &s) noexcept {
    ControlFieldOutcome out;
    out.inputNodes = s.nodes.size();
    auto fail = [&](ControlFieldReason reason) {
        out.reason = reason;
        return std::move(out);
    };
    if (!(s.minimumSize > 0 && s.maximumSize >= s.minimumSize && s.sizeFactor > 0) ||
        !bounded(s.minimumSize) || !bounded(s.maximumSize) || !bounded(s.sizeFactor) ||
        !(s.minimumSize * s.sizeFactor > 0) || !std::isfinite(s.maximumSize * s.sizeFactor))
        return fail(ControlFieldReason::InvalidEnvelope);
    if (cancelled(s))
        return fail(ControlFieldReason::Cancelled);
    if (s.nodes.size() > maximumControlFieldNodes || s.nodes.size() > s.maximumNodes)
        return fail(ControlFieldReason::BudgetExceeded);
    if (s.output >= s.nodes.size())
        return fail(ControlFieldReason::InvalidOutput);
    for (std::size_t i = 0; i < s.nodes.size(); ++i) {
        const auto &n = s.nodes[i];
        out.rejectedNode = static_cast<FieldId>(i);
        const auto reason = validateNode(n);
        if (reason != ControlFieldReason::None)
            return fail(reason);
        const auto count = arity(n.operation);
        if ((count && n.first >= i) || (count == 2 && n.second >= i))
            return fail(ControlFieldReason::InvalidDependency);
    }
    out.rejectedNode = invalidFieldId;
    std::array<bool, maximumControlFieldNodes> live{};
    live[s.output] = true;
    for (std::size_t i = s.nodes.size(); i-- > 0;)
        if (live[i]) {
            const auto count = arity(s.nodes[i].operation);
            if (count)
                live[s.nodes[i].first] = true;
            if (count == 2)
                live[s.nodes[i].second] = true;
            ++out.retainedNodes;
        }
    out.logicalBytes = out.retainedNodes * sizeof(FieldNode);
    if (out.logicalBytes > s.maximumBytes)
        return fail(ControlFieldReason::BudgetExceeded);
    try {
        auto &f = out.field;
        f.nodes_.reserve(static_cast<std::size_t>(out.retainedNodes));
        f.minimum_ = s.minimumSize;
        f.maximum_ = s.maximumSize;
        f.factor_ = s.sizeFactor;
        std::array<FieldId, maximumControlFieldNodes> map{};
        for (std::size_t i = 0; i < s.nodes.size(); ++i)
            if (live[i]) {
                if (cancelled(s))
                    return fail(ControlFieldReason::Cancelled);
                auto n = s.nodes[i];
                const auto count = arity(n.operation);
                if (count)
                    n.first = map[n.first];
                if (count == 2)
                    n.second = map[n.second];
                map[i] = static_cast<FieldId>(f.nodes_.size());
                f.nodes_.push_back(n);
                f.queryWork_ += (!count && n.operation != FieldOperation::Constant) ? 10 : 1;
            }
        f.valid_ = true;
        const auto range = f.rangeOnBox({-coordinateLimit, -coordinateLimit, -coordinateLimit},
                                        {coordinateLimit, coordinateLimit, coordinateLimit});
        if (!std::isfinite(range.lower) || !std::isfinite(range.upper)) {
            f.valid_ = false;
            return fail(ControlFieldReason::InvalidNode);
        }
        out.compiled = true;
        out.reason = ControlFieldReason::None;
        return out;
    } catch (const std::bad_alloc &) {
        return fail(ControlFieldReason::AllocationFailed);
    } catch (...) {
        return fail(ControlFieldReason::InternalError);
    }
}
const char *controlFieldReasonName(ControlFieldReason r) noexcept {
    switch (r) {
    case ControlFieldReason::None:
        return "none";
    case ControlFieldReason::InvalidEnvelope:
        return "invalid_envelope";
    case ControlFieldReason::InvalidNode:
        return "invalid_node";
    case ControlFieldReason::InvalidDependency:
        return "invalid_dependency";
    case ControlFieldReason::InvalidOutput:
        return "invalid_output";
    case ControlFieldReason::DegenerateGeometry:
        return "degenerate_geometry";
    case ControlFieldReason::BudgetExceeded:
        return "budget_exceeded";
    case ControlFieldReason::Cancelled:
        return "cancelled";
    case ControlFieldReason::AllocationFailed:
        return "allocation_failed";
    case ControlFieldReason::InternalError:
        return "internal_error";
    }
    return "unknown";
}
} // namespace paacs::sizing
