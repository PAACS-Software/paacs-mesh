// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/CappedSizeProfile.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
using namespace paacs::sizing;
namespace {
void check(bool ok, const char *why) {
    if (!ok)
        throw std::runtime_error(why);
}
// Independent minimum: min over the interval commutes with min over controls.
long double reference(const CappedSizeSpecification &s, double lo, double hi) {
    long double h = s.globalSize;
    for (auto c : s.controls) {
        const long double d = std::max(
            {0.L, static_cast<long double>(c.begin) - hi, static_cast<long double>(lo) - c.end});
        const long double g = c.growth;
        h = std::min(h, std::max<long double>(c.size, s.minimumSize) + 2 * ((g - 1) / (g + 1)) * d);
    }
    return h;
}
void verify(CappedSizeSpecification s) {
    auto compiled = compileCappedSizeProfile(s);
    check(compiled.compiled, "compile valid profile");
    const auto &p = compiled.profile;
    check(p.valid() && p.knots().front().position == 0 && p.knots().back().position == s.length,
          "domain");
    for (int i = 0; i <= 500; ++i) {
        const double x = s.length * (i / 500.);
        check(std::abs(p.eval(x) - reference(s, x, x)) <= 1e-12L * s.globalSize,
              "independent point field");
        const double b = std::min(s.length, x + s.length * .13);
        check(std::abs(p.minimumOnInterval(x, b) - reference(s, x, b)) <= 1e-12L * s.globalSize,
              "whole interval minimum");
    }
    const auto dist = planProfileDistribution(p);
    if (!dist.planned)
        std::cerr << "plan refused: " << profileReasonName(dist.reason) << '\n';
    check(dist.planned, "plan valid profile");
    check(dist.positions.front() == 0 && dist.positions.back() == s.length, "exact closure");
    double previous = 0;
    for (std::size_t i = 1; i < dist.positions.size(); ++i) {
        const double w = dist.positions[i] - dist.positions[i - 1];
        check(w > 0 && w <= reference(s, dist.positions[i - 1], dist.positions[i]) * (1 + 1e-10L),
              "independent interval cap");
        if (previous)
            check(std::max(w / previous, previous / w) <= p.maximumGrowth() * (1 + 1e-10),
                  "physical adjacent growth");
        previous = w;
    }
    auto reverse = s;
    std::reverse(reverse.controls.begin(), reverse.controls.end());
    const auto again = compileCappedSizeProfile(reverse);
    const auto repeated = planProfileDistribution(again.profile);
    check(repeated.planned && repeated.positions == dist.positions, "control ordering determinism");
    auto exact = s;
    exact.maximumBytes = compiled.workingBytes;
    check(compileCappedSizeProfile(exact).compiled, "exact compiler budget");
    exact.maximumBytes = compiled.workingBytes - 1;
    check(compileCappedSizeProfile(exact).reason == ProfileReason::BudgetExceeded,
          "one byte below compiler budget");
    ProfileDistributionLimits limits;
    limits.maximumBytes = dist.workingBytes;
    check(planProfileDistribution(p, limits).planned, "exact distribution budget");
    --limits.maximumBytes;
    check(planProfileDistribution(p, limits).reason == ProfileReason::BudgetExceeded,
          "one byte below distribution budget");
}
} // namespace
int main() {
    try {
        unsigned cases = 0, refusals = 0;
        for (double scale : {1e-6, 1., 1e6})
            for (int mode = 0; mode < 7; ++mode) {
                CappedSizeSpecification s;
                s.length = 6 * scale;
                s.globalSize = .5 * scale;
                if (mode == 1)
                    s.controls = {{0, 0, .1 * scale, 1.3}};
                if (mode == 2)
                    s.controls = {{s.length, s.length, .1 * scale, 1.3}};
                if (mode == 3)
                    s.controls = {{2 * scale, 3 * scale, .1 * scale, 1.3}};
                if (mode == 4)
                    s.controls = {{0, 0, .1 * scale, 1.2}, {s.length, s.length, .2 * scale, 1.5}};
                if (mode == 5)
                    s.controls = {{2 * scale, 4 * scale, .1 * scale, 1.25},
                                  {3 * scale, 5 * scale, .15 * scale, 1.1}};
                if (mode == 6)
                    s.controls = {{0, 0, .1 * scale, 1}};
                verify(s);
                ++cases;
            }
        std::mt19937_64 rng(0x706161637370726fULL);
        auto random = [&] { return std::generate_canonical<double, 53>(rng); };
        for (unsigned i = 0; i < 200; ++i) {
            CappedSizeSpecification s;
            s.length = std::pow(10., 12 * random() - 6);
            s.globalSize = s.length * (.1 + .2 * random());
            s.minimumSize = .005 * s.length;
            for (unsigned k = 0, n = 1 + unsigned(rng() % 12); k < n; ++k) {
                double a = random() * s.length, b = random() * s.length;
                if (a > b)
                    std::swap(a, b);
                if (k % 3 == 0)
                    b = a;
                s.controls.push_back({a, b, s.length * (.004 + .3 * random()), 1.01 + random()});
            }
            verify(s);
            ++cases;
        }
        CappedSizeSpecification base;
        base.length = 4;
        base.controls = {{0, 0, .05, 1.2}};
        const auto valid = compileCappedSizeProfile(base);
        const auto plan = planProfileDistribution(valid.profile);
        check(plan.planned && plan.plateauCellsRemoved > 0, "compact constant global plateau");
        const double last = plan.positions.back() - plan.positions[plan.positions.size() - 2];
        check(last > .9 * base.globalSize, "global plateau reaches useful size");
        ProfileDistributionLimits uniform;
        uniform.maximumAdjacentGrowth = 1;
        const auto allfine = planProfileDistribution(valid.profile, uniform);
        check(allfine.planned && allfine.maximumAdjacentGrowth < 1 + 1e-10,
              "explicit uniform refinement");
        for (int mode = 0; mode < 14; ++mode) {
            auto s = base;
            if (mode == 0)
                s.length = 0;
            if (mode == 1)
                s.globalSize = std::numeric_limits<double>::infinity();
            if (mode == 2)
                s.minimumSize = -1;
            if (mode == 3)
                s.minimumSize = 1;
            if (mode == 4)
                s.controls[0].begin = -.01;
            if (mode == 5)
                s.controls[0].end = 5;
            if (mode == 6)
                s.controls[0].size = 0;
            if (mode == 7)
                s.controls[0].growth = .99;
            if (mode == 8)
                s.controls[0].end = std::numeric_limits<double>::quiet_NaN();
            if (mode == 9)
                s.maximumControls = 0;
            if (mode == 10)
                s.maximumKnots = 1;
            if (mode == 11)
                s.maximumBytes = 1;
            if (mode == 12)
                s.maximumWork = 0;
            std::atomic<bool> stop{true};
            if (mode == 13)
                s.cancellation = &stop;
            const auto bad = compileCappedSizeProfile(s);
            check(!bad.compiled && !bad.profile.valid(),
                  "refusal cannot expose partially compiled field");
            ++refusals;
        }
        for (int mode = 0; mode < 5; ++mode) {
            ProfileDistributionLimits l;
            if (mode == 0)
                l.maximumCells = 1;
            if (mode == 1)
                l.maximumBytes = 1;
            if (mode == 2)
                l.maximumAdjacentGrowth = .9;
            if (mode == 3)
                l.maximumAdjacentGrowth = std::numeric_limits<double>::quiet_NaN();
            std::atomic<bool> stop{true};
            if (mode == 4)
                l.cancellation = &stop;
            const auto bad = planProfileDistribution(valid.profile, l);
            check(!bad.planned && bad.positions.empty(), "distribution refusal");
            ++refusals;
        }
        check(std::isnan(valid.profile.eval(-1)) &&
                  std::isnan(valid.profile.minimumOnInterval(2, 1)),
              "invalid query NaN");
        check(!planProfileDistribution(CappedSizeProfile{}).planned, "uncompiled profile refused");
        ++refusals;
        auto temporary = compileCappedSizeProfile(base);
        auto moved = std::move(temporary.profile);
        check(moved.valid() && !temporary.profile.valid() && std::isnan(temporary.profile.eval(0)), "moved-from profile remains safe");
        check(!planProfileDistribution(temporary.profile).planned, "moved-from planner refusal");
        std::cout << "capped profile: PASS (" << cases << " independent profiles, " << refusals
                  << " refusals)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
