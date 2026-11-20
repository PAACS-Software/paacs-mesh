// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sizing/ControlField.hpp"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <set>
#include <stdexcept>
namespace {
long remaining = -1;
}
void *operator new(std::size_t n) {
    if (remaining == 0)
        throw std::bad_alloc();
    if (remaining > 0)
        --remaining;
    if (auto p = std::malloc(n ? n : 1))
        return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t n) {
    return ::operator new(n);
}
// Keep standard-library temporary buffers on this same injection allocator.
// ASan's default nothrow interceptor otherwise bypasses our throwing overload.
void *operator new(std::size_t n, const std::nothrow_t &) noexcept {
    try {
        return ::operator new(n);
    } catch (...) {
        return nullptr;
    }
}
void *operator new[](std::size_t n, const std::nothrow_t &tag) noexcept {
    return ::operator new(n, tag);
}
void operator delete(void *p) noexcept {
    std::free(p);
}
void operator delete[](void *p) noexcept {
    std::free(p);
}
void operator delete(void *p, std::size_t) noexcept {
    std::free(p);
}
void operator delete[](void *p, std::size_t) noexcept {
    std::free(p);
}
void operator delete(void *p, const std::nothrow_t &) noexcept {
    std::free(p);
}
void operator delete[](void *p, const std::nothrow_t &) noexcept {
    std::free(p);
}
using namespace paacs;
void need(bool b) {
    if (!b)
        throw std::runtime_error("spatial allocation invariant");
}
std::vector<std::uint8_t> archive(const mesh::Mesh &m) {
    std::vector<std::uint8_t> b;
    need(sdk::serializeMeshArchive(m, b).completed);
    return b;
}
// Exhaust new compact paths. The inherited general Tet4 pipeline has many
// more allocations; use reproducible prefix/tail/stratified injection there.
template <class Call> void exercise(const char *name, Call call, bool exhaustive = true) {
    mesh::Mesh baseline;
    remaining = std::numeric_limits<long>::max();
    auto initial = call(baseline);
    const auto count = std::numeric_limits<long>::max() - remaining;
    remaining = -1;
    need(initial.published && count > 0);
    const auto expected = archive(baseline);
    std::set<long> points;
    if (exhaustive) {
        need(count < 20000);
        for (long i = 0; i <= count; ++i)
            points.insert(i);
    } else {
        for (long i = 0; i < 64 && i < count; ++i)
            points.insert(i);
        for (long i = std::max(0L, count - 128); i <= count; ++i)
            points.insert(i);
        for (long i = 0; i <= 256; ++i)
            points.insert(count * i / 256);
        for (long i = 1; i < count; i *= 2)
            points.insert(i);
    }
    long failures = 0, successes = 0;
    bool terminalSuccess = false;
    for (const long i : points) {
        mesh::Mesh held;
        held.addNode({9, 8, 7}, {0, 9});
        const auto before = archive(held);
        remaining = i;
        auto r = call(held);
        remaining = -1;
        if (r.published) {
            need(archive(held) == expected);
            ++successes;
            terminalSuccess = terminalSuccess || i == count;
        } else {
            need(r.destinationUnchanged && archive(held) == before);
            ++failures;
        }
    }
    need(failures > 0 && terminalSuccess);
    std::cout << name << ": PASS (" << count << " observed allocations; "
              << (exhaustive ? "exhaustive" : "stratified") << "; " << failures
              << " rejected injections, " << successes << " unchanged-output successes)\n"
              << std::flush;
}
int main() {
    try {
        sizing::ControlFieldSpecification s;
        s.minimumSize = .1;
        s.maximumSize = 2;
        s.nodes = {sizing::field::distancePoint({0, 0, 0}),
                   sizing::field::threshold(0, 0, 2, .3, .65)};
        s.output = 1;
        remaining = 0;
        auto failed = sizing::compileControlField(s);
        remaining = -1;
        need(!failed.compiled && !failed.field.valid() &&
             failed.reason == sizing::ControlFieldReason::AllocationFailed);
        remaining = 1;
        auto f = sizing::compileControlField(s);
        remaining = -1;
        need(f.compiled);
        std::cout << "field compile: PASS (1 allocation point)\n";
        bool okay = true;
        remaining = 0;
        for (int i = 0; i < 1000; ++i) {
            okay = okay && f.field.eval({i / 1000., .2, .3}) > 0 &&
                   f.field.minimumOnBox({0, 0, 0}, {.1, .2, .3}) > 0;
        }
        remaining = -1;
        need(okay);
        std::cout << "point and whole-box queries: zero allocations\n";
        algo::AdaptiveBoxRequest box;
        box.field = &f.field;
        box.threads = 1;
        exercise("adaptive graph field", [&](auto &m) { return algo::meshAdaptiveBox(box, m); });
        mesh::Mesh m;
        need(algo::meshAdaptiveBox(box, m).published);
        remaining = 0;
        auto audit = algo::auditMeshSize(m, f.field);
        remaining = -1;
        need(audit.complete);
        std::cout << "mesh size audit: zero allocations\n";
        s.nodes = {sizing::field::constant(.6)};
        s.output = 0;
        auto constant = sizing::compileControlField(s);
        need(constant.compiled);
        const std::vector<std::vector<Vec2>> loops{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
        auto model = geo::makePolyhedralModel(geo::extrudedPolygon(loops, {0, 0, 1}));
        need(model.built());
        algo::GeometryMeshingOptions o;
        o.targetSize = .6;
        o.spatialField = &constant.field;
        o.threads = 1;
        exercise(
            "geometry graph field",
            [&](auto &mesh) { return algo::meshGeometryModel(*model.model, o, mesh); }, false);
    } catch (const std::exception &e) {
        remaining = -1;
        std::cerr << e.what() << '\n';
        return 1;
    }
}
