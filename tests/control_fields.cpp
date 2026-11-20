// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/ControlField.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>

using namespace paacs;
using namespace paacs::sizing;
namespace {
std::uint64_t checks = 0, refusals = 0;
void need(bool b, const char *message) {
    ++checks;
    if (!b)
        throw std::runtime_error(message);
}
void close(double a, long double b, const char *message) {
    need(std::abs(static_cast<long double>(a) - b) <= 3e-11L * std::max(1e-100L, std::abs(b)),
         message);
}
ControlField compile(std::vector<FieldNode> nodes, double lo = 1e-12, double hi = 1e90,
                     double factor = 1) {
    ControlFieldSpecification s;
    s.nodes = std::move(nodes);
    s.output = static_cast<FieldId>(s.nodes.size() - 1);
    s.minimumSize = lo;
    s.maximumSize = hi;
    s.sizeFactor = factor;
    auto r = compileControlField(s);
    need(r.compiled, "compile valid field");
    return std::move(r.field);
}
long double segment2(long double x, long double y, long double z, long double ax, long double ay,
                     long double bx, long double by) {
    const auto dx = bx - ax, dy = by - ay;
    const auto t = std::clamp(((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy), 0.L, 1.L);
    return std::sqrt((x - ax - t * dx) * (x - ax - t * dx) + (y - ay - t * dy) * (y - ay - t * dy) +
                     z * z);
}
void geometryOracles() {
    std::mt19937_64 rng(73519);
    std::uniform_real_distribution<double> unit(-1, 2);
    const Vec3 u{.6, .8, 0}, v{-.8, .6, 0}, normal{0, 0, 1};
    for (double scale : {1e-9, 1e-6, 1., 1e6, 1e9}) {
        const Vec3 origin{2 * scale, -3 * scale, 4 * scale};
        auto point = compile({field::distancePoint(origin)}, scale * 1e-12);
        auto segment = compile({field::distanceSegment(origin, origin + u * scale)}, scale * 1e-12);
        auto triangle =
            compile({field::distanceTriangle(origin, origin + u * scale, origin + v * scale)},
                    scale * 1e-12);
        auto ball = compile({field::distanceBall(origin, .2 * scale)}, scale * 1e-12);
        auto box = compile({field::distanceBox(origin, origin + Vec3{scale, 2 * scale, 3 * scale})},
                           scale * 1e-12);
        for (int i = 0; i < 300; ++i) {
            const long double x = unit(rng), y = unit(rng), z = unit(rng);
            const Vec3 p = origin + (u * static_cast<double>(x) + v * static_cast<double>(y) +
                                     normal * static_cast<double>(z)) *
                                        scale;
            const auto distance = std::sqrt(x * x + y * y + z * z);
            close(point.eval(p), scale * distance, "point distance oracle");
            close(segment.eval(p), scale * segment2(x, y, z, 0, 0, 1, 0),
                  "segment distance oracle");
            const auto tri =
                x >= 0 && y >= 0 && x + y <= 1
                    ? std::abs(z)
                    : std::min({segment2(x, y, z, 0, 0, 1, 0), segment2(x, y, z, 0, 0, 0, 1),
                                segment2(x, y, z, 1, 0, 0, 1)});
            close(triangle.eval(p), scale * tri, "triangle distance oracle");
            close(ball.eval(p), scale * std::max(1e-12L, distance - .2L), "ball distance oracle");
            const long double px = p.x - origin.x, py = p.y - origin.y, pz = p.z - origin.z;
            const long double dx = std::max({-px, px - scale, 0.L}),
                              dy = std::max({-py, py - 2 * scale, 0.L}),
                              dz = std::max({-pz, pz - 3 * scale, 0.L});
            close(box.eval(p),
                  std::max(static_cast<long double>(scale * 1e-12),
                           std::sqrt(dx * dx + dy * dy + dz * dz)),
                  "box distance oracle");
        }
    }
}
void composition() {
    auto linear =
        compile({field::distancePoint({0, 0, 0}), field::threshold(0, .1, .9, .05, .45)}, .01, 1);
    auto smooth = compile(
        {field::distancePoint({0, 0, 0}), field::threshold(0, .1, .9, .05, .45, true)}, .01, 1);
    auto decreasing = compile(
        {field::distancePoint({0, 0, 0}), field::threshold(0, .1, .9, .45, .05, true)}, .01, 1);
    for (int i = 0; i <= 1000; ++i) {
        const long double x = i / 1000.L, t = std::clamp((x - .1L) / .8L, 0.L, 1.L);
        close(linear.eval({double(x), 0, 0}), .05L + .4L * t, "linear threshold oracle");
        close(smooth.eval({double(x), 0, 0}), .05L + .4L * t * t * (3 - 2 * t),
              "smooth threshold oracle");
        close(decreasing.eval({double(x), 0, 0}), .45L - .4L * t * t * (3 - 2 * t),
              "decreasing threshold oracle");
    }
    std::vector<FieldNode> nodes = {field::distancePoint({0, 0, 0}),
                                    field::threshold(0, 0, 1, .05, .5),
                                    field::distancePoint({1, 0, 0}),
                                    field::threshold(2, 0, 1, .1, .6),
                                    field::minimum(1, 3),
                                    field::maximum(1, 3),
                                    field::scale(4, 2),
                                    field::clamp(6, .2, .7)};
    auto min = compile(std::vector<FieldNode>(nodes.begin(), nodes.begin() + 5), .01, 1);
    auto max = compile(std::vector<FieldNode>(nodes.begin(), nodes.begin() + 6), .01, 1);
    auto nested = compile(nodes, .15, .6, 1.5);
    need(nested.nodeCount() == 7, "unused maximum node pruned");
    for (int i = 0; i <= 100; ++i) {
        const long double x = i / 100.L, a = .05L + .45L * x, b = .1L + .5L * (1 - x);
        close(min.eval({double(x), 0, 0}), std::min(a, b), "minimum wins");
        close(max.eval({double(x), 0, 0}), std::max(a, b), "maximum wins");
        close(nested.eval({double(x), 0, 0}),
              1.5L * std::clamp(std::clamp(2 * std::min(a, b), .2L, .7L), .15L, .6L),
              "clamp before global factor");
    }
    auto constant = compile({field::constant(.25)}, .01, 1);
    need(constant.minimumOnBox({0, 0, 0}, {1, 1, 1}) == .25, "constant dyadic bound stays exact");
    auto moved = std::move(constant);
    need(moved.valid() && !constant.valid() && std::isnan(constant.eval({0, 0, 0})),
         "moved-from safe");
}
void intervalBounds() {
    std::vector<FieldNode> base = {field::distanceSegment({.1, .2, .1}, {.9, .8, .9}),
                                   field::threshold(0, .03, .4, .04, .5, true),
                                   field::distanceTriangle({0, 0, .7}, {1, 0, .7}, {.2, 1, .5}),
                                   field::threshold(2, .02, .3, .4, .06),
                                   field::minimum(1, 3),
                                   field::maximum(1, 3),
                                   field::scale(4, .7),
                                   field::clamp(6, .05, .3)};
    std::mt19937_64 rng(555012);
    std::uniform_real_distribution<double> random(0, 1);
    for (FieldId root : {0U, 1U, 2U, 3U, 4U, 5U, 7U}) {
        ControlFieldSpecification s;
        s.nodes = base;
        s.output = root;
        s.minimumSize = .001;
        s.maximumSize = 2;
        auto r = compileControlField(s);
        need(r.compiled, "interval program");
        for (int i = 0; i < 350; ++i) {
            Vec3 a{random(rng), random(rng), random(rng)}, b{random(rng), random(rng), random(rng)};
            Vec3 lo{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)},
                hi{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
            const auto bound = r.field.rangeOnBox(lo, hi);
            need(bound.lower > 0 && bound.upper >= bound.lower, "ordered range");
            for (int j = 0; j < 20; ++j) {
                const Vec3 p{std::lerp(lo.x, hi.x, random(rng)), std::lerp(lo.y, hi.y, random(rng)),
                             std::lerp(lo.z, hi.z, random(rng))};
                const double h = r.field.eval(p);
                need(h >= bound.lower && h <= bound.upper, "whole-box enclosure");
            }
        }
    }
    auto tiny = compile(
        {field::distanceBall({.137, .231, .419}, 1e-5), field::threshold(0, 0, .01, .03, .4)}, .01,
        1);
    close(tiny.minimumOnBox({0, 0, 0}, {1, 1, 1}), .03L, "tiny enclosed control found");
    for (int mask = 0; mask < 8; ++mask)
        need(tiny.eval({mask & 1 ? 1. : 0., mask & 2 ? 1. : 0., mask & 4 ? 1. : 0.}) == .4,
             "corner samples miss tiny region");
    auto f = compile(base, .001, 2);
    std::array<double, 2000> one{}, many{};
    for (std::size_t i = 0; i < one.size(); ++i)
        one[i] = f.eval({double(i % 37) / 37, double(i % 53) / 53, double(i % 71) / 71});
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < 4; ++t)
        workers.emplace_back([&, t] {
            for (std::size_t i = t; i < many.size(); i += 4)
                many[i] = f.eval({double(i % 37) / 37, double(i % 53) / 53, double(i % 71) / 71});
        });
    for (auto &thread : workers)
        thread.join();
    need(one == many, "thread-identical queries");
}
void invalid() {
    ControlFieldSpecification good;
    good.nodes = {field::constant(.2)};
    good.output = 0;
    auto reject = [&](ControlFieldSpecification s, ControlFieldReason reason) {
        auto r = compileControlField(s);
        need(!r.compiled && r.reason == reason, "explicit compile refusal");
        ++refusals;
    };
    for (double x : {0., -1., std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()}) {
        auto s = good;
        s.minimumSize = x;
        reject(s, ControlFieldReason::InvalidEnvelope);
        s = good;
        s.sizeFactor = x;
        reject(s, ControlFieldReason::InvalidEnvelope);
    }
    auto s = good;
    s.minimumSize = 2;
    reject(s, ControlFieldReason::InvalidEnvelope);
    s = good;
    s.output = 1;
    reject(s, ControlFieldReason::InvalidOutput);
    s = good;
    s.nodes.clear();
    reject(s, ControlFieldReason::InvalidOutput);
    s = good;
    s.nodes[0] = field::scale(0, 2);
    reject(s, ControlFieldReason::InvalidDependency);
    s = good;
    s.nodes.push_back(field::minimum(0, 2));
    reject(s, ControlFieldReason::InvalidDependency);
    s = good;
    s.nodes[0].operation = static_cast<FieldOperation>(255);
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes[0] = field::distanceSegment({0, 0, 0}, {0, 0, 0});
    reject(s, ControlFieldReason::DegenerateGeometry);
    s = good;
    s.nodes[0] = field::distanceTriangle({0, 0, 0}, {1, 0, 0}, {2, 0, 0});
    reject(s, ControlFieldReason::DegenerateGeometry);
    s = good;
    s.nodes[0] = field::distanceBox({1, 0, 0}, {0, 1, 1});
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes[0] = field::distanceBall({0, 0, 0}, -1);
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes[0] = field::distancePoint({1e200, 0, 0});
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes.push_back(field::threshold(0, 1, 1, .1, .2));
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes.push_back(field::threshold(0, 0, 1, -1, .2));
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes.push_back(field::scale(0, -1));
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.nodes.push_back(field::clamp(0, 1, .1));
    reject(s, ControlFieldReason::InvalidNode);
    s = good;
    s.maximumNodes = 0;
    reject(s, ControlFieldReason::BudgetExceeded);
    s = good;
    s.maximumBytes = sizeof(FieldNode) - 1;
    reject(s, ControlFieldReason::BudgetExceeded);
    s = good;
    s.maximumBytes = sizeof(FieldNode);
    need(compileControlField(s).compiled, "exact bytes accepted");
    s = good;
    s.nodes.resize(maximumControlFieldNodes + 1);
    reject(s, ControlFieldReason::BudgetExceeded);
    s = good;
    std::atomic<bool> stop = true;
    s.cancellation = &stop;
    reject(s, ControlFieldReason::Cancelled);
    auto f = compile(good.nodes);
    need(std::isnan(f.eval({INFINITY, 0, 0})), "nonfinite query rejected");
    need(std::isnan(f.minimumOnBox({1, 0, 0}, {0, 1, 1})), "reversed box rejected");
    need(std::isnan(f.minimumOnBox({0, 0, 0}, {1, 1, NAN})), "NaN box rejected");
}
} // namespace
int main() {
    try {
        geometryOracles();
        composition();
        intervalBounds();
        invalid();
        std::cout << "control fields: PASS (" << checks << " independent checks, " << refusals
                  << " compile refusals)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
