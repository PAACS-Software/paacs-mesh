// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sizing/CappedSizeProfile.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <tuple>

namespace paacs::sizing {
namespace {
struct Stop {
    ProfileReason reason;
};
void need(bool ok, ProfileReason reason) {
    if (!ok)
        throw Stop{reason};
}
void cancel(const std::atomic<bool> *flag) {
    need(!flag || !flag->load(std::memory_order_relaxed), ProfileReason::Cancelled);
}
bool positive(double x) {
    return std::isfinite(x) && x > 0;
}
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
long double line(ProfileKnot a, ProfileKnot b, double x) {
    return static_cast<long double>(a.size) +
           (static_cast<long double>(b.size) - a.size) *
               ((static_cast<long double>(x) - a.position) /
                (static_cast<long double>(b.position) - a.position));
}
long double metric(ProfileKnot a, ProfileKnot b) {
    const long double d = static_cast<long double>(b.size) - a.size;
    const long double dx = static_cast<long double>(b.position) - a.position;
    return d == 0 ? dx / a.size : dx * std::log1p(d / a.size) / d;
}
double ratio(double a, double b) {
    return std::max(a / b, b / a);
}
} // namespace

double CappedSizeProfile::eval(double x) const noexcept {
    if (!valid() || !std::isfinite(x) || x < 0 || x > length_)
        return nan;
    auto it = std::lower_bound(knots_.begin(), knots_.end(), x,
                               [](auto k, double p) { return k.position < p; });
    if (it == knots_.begin() || it->position == x)
        return it->size;
    return static_cast<double>(line(*(it - 1), *it, x));
}
double CappedSizeProfile::minimumOnInterval(double a, double b) const noexcept {
    if (!valid() || !std::isfinite(a) || !std::isfinite(b) || a < 0 || b > length_ || b < a)
        return nan;
    double value = std::min(eval(a), eval(b));
    auto it = std::upper_bound(knots_.begin(), knots_.end(), a,
                               [](double p, auto k) { return p < k.position; });
    for (; it != knots_.end() && it->position < b; ++it)
        value = std::min(value, it->size);
    return value;
}
CappedSizeOutcome compileCappedSizeProfile(const CappedSizeSpecification &s) noexcept {
    CappedSizeOutcome out;
    try {
        cancel(s.cancellation);
        need(positive(s.length) && positive(s.globalSize) && std::isfinite(s.minimumSize) &&
                 s.minimumSize >= 0 && s.minimumSize <= s.globalSize,
             ProfileReason::InvalidRequest);
        need(s.controls.size() <= s.maximumControls && s.maximumKnots >= 2,
             ProfileReason::BudgetExceeded);
        for (const auto &c : s.controls) {
            cancel(s.cancellation);
            need(std::isfinite(c.begin) && std::isfinite(c.end) && c.begin >= 0 &&
                     c.begin <= c.end && c.end <= s.length && positive(c.size) &&
                     std::isfinite(c.growth) && c.growth >= 1,
                 ProfileReason::InvalidRequest);
        }
        auto budget = [&](long double bytes) {
            need(bytes <= s.maximumBytes && bytes <= SIZE_MAX, ProfileReason::BudgetExceeded);
            out.workingBytes = std::max(out.workingBytes, static_cast<std::uint64_t>(bytes));
        };
        budget(static_cast<long double>(s.controls.size()) * sizeof(ProfileSizeControl) +
               2 * sizeof(ProfileKnot));
        auto controls = s.controls;
        std::sort(controls.begin(), controls.end(), [](const auto &a, const auto &b) {
            return std::tie(a.begin, a.end, a.size, a.growth) <
                   std::tie(b.begin, b.end, b.size, b.growth);
        });
        CappedSizeProfile p;
        p.length_ = s.length;
        p.global_ = s.globalSize;
        p.knots_ = {{0, s.globalSize}, {s.length, s.globalSize}};
        for (const auto &c : controls) {
            cancel(s.cancellation);
            p.growth_ = std::max(p.growth_, c.growth);
            const long double slope = 2 * (1 - 2 / (static_cast<long double>(c.growth) + 1));
            const long double h = std::max(c.size, s.minimumSize);
            if (h >= s.globalSize)
                continue;
            p.gradient_ = std::max(p.gradient_, static_cast<double>(slope));
            const auto value = [&](double z) {
                const long double d = z < c.begin ? static_cast<long double>(c.begin) - z
                                      : z > c.end ? static_cast<long double>(z) - c.end
                                                  : 0;
                return h + slope * d;
            };
            std::vector<ProfileKnot> next;
            auto append = [&](double z, long double hnew) {
                if (!next.empty() && z == next.back().position) {
                    next.back().size = std::min(next.back().size, static_cast<double>(hnew));
                    return;
                }
                need(next.size() < s.maximumKnots, ProfileReason::BudgetExceeded);
                if (next.size() == next.capacity()) {
                    const auto capacity = std::min<std::uint64_t>(
                        s.maximumKnots, std::max<std::uint64_t>(2, 2 * next.capacity()));
                    budget(static_cast<long double>(controls.capacity()) *
                               sizeof(ProfileSizeControl) +
                           (static_cast<long double>(p.knots_.capacity()) + next.capacity() +
                            capacity) *
                               sizeof(ProfileKnot));
                    next.reserve(static_cast<std::size_t>(capacity));
                }
                need(positive(static_cast<double>(hnew)), ProfileReason::NumericallyUnresolved);
                next.push_back({z, static_cast<double>(hnew)});
            };
            for (std::size_t i = 1; i < p.knots_.size(); ++i) {
                cancel(s.cancellation);
                const auto a = p.knots_[i - 1], b = p.knots_[i];
                const std::array<double, 4> cuts{
                    a.position, std::clamp(c.begin, a.position, b.position),
                    std::clamp(c.end, a.position, b.position), b.position};
                for (unsigned j = 1; j < cuts.size(); ++j) {
                    const double lo = cuts[j - 1], hi = cuts[j];
                    if (lo == hi)
                        continue;
                    need(out.work < s.maximumWork, ProfileReason::BudgetExceeded);
                    ++out.work;
                    const auto oa = line(a, b, lo), ob = line(a, b, hi), ca = value(lo),
                               cb = value(hi);
                    append(lo, std::min(oa, ca));
                    const long double da = oa - ca, db = ob - cb;
                    if ((da < 0 && db > 0) || (da > 0 && db < 0)) {
                        const double z = static_cast<double>(static_cast<long double>(lo) +
                                                             (static_cast<long double>(hi) - lo) *
                                                                 da / (da - db));
                        if (z > lo && z < hi)
                            append(z, std::min(line(a, b, z), value(z)));
                    }
                    append(hi, std::min(ob, cb));
                }
            }
            p.knots_.swap(next);
        }
        p.minimum_ = s.globalSize;
        for (const auto k : p.knots_)
            p.minimum_ = std::min(p.minimum_, k.size);
        cancel(s.cancellation);
        p.valid_ = true;
        out.profile = std::move(p);
        out.compiled = true;
        out.reason = ProfileReason::None;
    } catch (const Stop &e) {
        out.reason = e.reason;
    } catch (const std::bad_alloc &) {
        out.reason = ProfileReason::AllocationFailed;
    } catch (...) {
        out.reason = ProfileReason::InternalError;
    }
    return out;
}

ProfileDistribution planProfileDistribution(const CappedSizeProfile &p,
                                            const ProfileDistributionLimits &limits) noexcept {
    ProfileDistribution out;
    try {
        cancel(limits.cancellation);
        const double growth =
            limits.maximumAdjacentGrowth == 0 ? p.maximumGrowth() : limits.maximumAdjacentGrowth;
        need(p.valid() && std::isfinite(growth) && growth >= 1, ProfileReason::InvalidRequest);
        const auto knots = p.knots();
        const long double slope = p.gradientBound();
        const bool uniform = growth == 1 || slope == 0;
        long double total = 0;
        if (uniform)
            total = static_cast<long double>(p.length()) / p.minimumSize();
        else
            for (std::size_t i = 1; i < knots.size(); ++i) {
                cancel(limits.cancellation);
                total += metric(knots[i - 1], knots[i]);
            }
        long double step = 1;
        if (!uniform)
            step = std::min(std::log1p(slope), std::log(static_cast<long double>(growth))) / slope;
        const long double requested = std::ceil(total / step);
        need(std::isfinite(requested) && requested >= 1 && requested <= limits.maximumCells &&
                 (requested + 1) * sizeof(double) <= limits.maximumBytes &&
                 requested + 1 <= SIZE_MAX / sizeof(double),
             ProfileReason::BudgetExceeded);
        const auto cells = static_cast<std::size_t>(requested);
        out.workingBytes = (cells + 1) * sizeof(double);
        out.positions.resize(cells + 1);
        std::size_t segment = 1;
        long double consumed = 0, span = metric(knots[0], knots[1]);
        for (std::size_t i = 1; i < cells; ++i) {
            if (i % 1024 == 0)
                cancel(limits.cancellation);
            if (uniform)
                out.positions[i] =
                    static_cast<double>(static_cast<long double>(p.length()) * i / cells);
            else {
                const long double target = total * i / cells;
                while (segment + 1 < knots.size() && consumed + span < target) {
                    consumed += span;
                    ++segment;
                    span = metric(knots[segment - 1], knots[segment]);
                }
                const auto a = knots[segment - 1], b = knots[segment];
                const long double s = (static_cast<long double>(b.size) - a.size) /
                                      (static_cast<long double>(b.position) - a.position);
                const long double dt = std::clamp(target - consumed, 0.L, span);
                const long double dx = s == 0 ? a.size * dt : a.size * std::expm1(s * dt) / s;
                out.positions[i] =
                    std::clamp(static_cast<double>(a.position + dx), a.position, b.position);
            }
        }
        out.positions.front() = 0;
        out.positions.back() = p.length();
        // Equalize fully constant runs, subject to their existing neighbours.
        // Read indices stay ahead of write indices: no additional allocation.
        std::size_t read = 0, write = 0;
        while (read < cells) {
            cancel(limits.cancellation);
            const double start = out.positions[read], cap = p.eval(start);
            std::size_t end = read;
            while (end < cells && p.eval(out.positions[end + 1]) == cap &&
                   p.minimumOnInterval(out.positions[end], out.positions[end + 1]) == cap)
                ++end;
            if (end > read + 1) {
                const double last = out.positions[end];
                const auto old = end - read;
                auto count = std::max<std::size_t>(
                    1, static_cast<std::size_t>(std::ceil((last - start) / cap)));
                count = std::min(count, old);
                for (; count < old; ++count) {
                    const double width = (last - start) / static_cast<double>(count);
                    if ((write == 0 || ratio(width, start - out.positions[write - 1]) <= growth) &&
                        (end == cells || ratio(width, out.positions[end + 1] - last) <= growth))
                        break;
                }
                out.positions[write] = start;
                for (std::size_t j = 1; j <= count; ++j)
                    out.positions[write + j] =
                        j == count ? last
                                   : static_cast<double>(static_cast<long double>(start) +
                                                         (static_cast<long double>(last) - start) *
                                                             j / count);
                write += count;
                read = end;
                out.plateauCellsRemoved += old - count;
            } else {
                out.positions[write++] = start;
                out.positions[write] = out.positions[++read];
            }
        }
        out.positions.resize(write + 1);
        out.minimumWidth = std::numeric_limits<double>::infinity();
        double previous = 0;
        for (std::size_t i = 1; i < out.positions.size(); ++i) {
            if (i % 1024 == 0)
                cancel(limits.cancellation);
            const double w = out.positions[i] - out.positions[i - 1];
            const double bound = p.minimumOnInterval(out.positions[i - 1], out.positions[i]);
            need(positive(w) && positive(bound), ProfileReason::NumericallyUnresolved);
            out.minimumWidth = std::min(out.minimumWidth, w);
            out.maximumWidth = std::max(out.maximumWidth, w);
            out.maximumWidthToTarget = std::max(out.maximumWidthToTarget, w / bound);
            if (previous > 0)
                out.maximumAdjacentGrowth = std::max(out.maximumAdjacentGrowth, ratio(w, previous));
            previous = w;
        }
        need(std::isfinite(out.maximumAdjacentGrowth) && out.maximumWidthToTarget <= 1 + 1e-10 &&
                 out.maximumAdjacentGrowth <= growth * (1 + 1e-10),
             ProfileReason::NumericallyUnresolved);
        cancel(limits.cancellation);
        out.planned = true;
        out.reason = ProfileReason::None;
    } catch (const Stop &e) {
        out.reason = e.reason;
    } catch (const std::bad_alloc &) {
        out.reason = ProfileReason::AllocationFailed;
    } catch (...) {
        out.reason = ProfileReason::InternalError;
    }
    if (!out.planned)
        out.positions.clear();
    return out;
}
const char *profileReasonName(ProfileReason reason) noexcept {
    switch (reason) {
    case ProfileReason::None:
        return "none";
    case ProfileReason::InvalidRequest:
        return "invalid_request";
    case ProfileReason::BudgetExceeded:
        return "budget_exceeded";
    case ProfileReason::Cancelled:
        return "cancelled";
    case ProfileReason::NumericallyUnresolved:
        return "numerically_unresolved";
    case ProfileReason::AllocationFailed:
        return "allocation_failed";
    case ProfileReason::InternalError:
        return "internal_error";
    }
    return "unknown";
}
} // namespace paacs::sizing
