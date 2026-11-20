// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
namespace {
long remaining = -1;
}
void *operator new(std::size_t n) {
    if (remaining == 0)
        throw std::bad_alloc();
    if (remaining > 0)
        --remaining;
    if (void *p = std::malloc(n ? n : 1))
        return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t n) {
    return ::operator new(n);
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
using namespace paacs;
void need(bool x) {
    if (!x)
        throw std::runtime_error("failure-injection invariant");
}
std::vector<std::uint8_t> bytes(const mesh::Mesh &m) {
    std::vector<std::uint8_t> b;
    need(sdk::serializeMeshArchive(m, b).completed);
    return b;
}
template <class Call> void exercise(const char *label, Call call) {
    bool succeeded = false;
    int failures = 0;
    for (int i = 0; i < 10000; ++i) {
        mesh::Mesh held;
        held.addNode({9, 8, 7}, {0, 9});
        const auto before = bytes(held);
        remaining = i;
        const auto r = call(held);
        remaining = -1;
        if (r.published) {
            succeeded = true;
            break;
        }
        need(r.destinationUnchanged && bytes(held) == before);
        ++failures;
    }
    need(succeeded && failures > 4);
    std::cout << label << ": PASS (" << failures << " allocation points)\n" << std::flush;
}
int main() {
    try {
        for (bool triangle : {false, true}) {
            std::vector<std::vector<Vec2>> loops =
                triangle ? std::vector<std::vector<Vec2>>{{{0, 0}, {1, 0}, {.2, 1}}}
                         : std::vector<std::vector<Vec2>>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
            auto model = geo::makePolyhedralModel(geo::extrudedPolygon(loops, {0, 0, 1}));
            need(model.built());
            sizing::CappedSizeSpecification s;
            s.length = 1;
            s.globalSize = .5;
            s.controls = {{0, 0, .2, 1.3}};
            auto f = sizing::compileCappedSizeProfile(s);
            need(f.compiled);
            algo::GeometryMeshingOptions o;
            o.targetSize = .5;
            o.threads = 1;
            algo::GeometryBodyControl c;
            c.solidTag = 1;
            c.sourceFaceTag = 1;
            c.sizeProfile = &f.profile;
            o.bodyControls = {c};
            for (bool tetra : {false, true}) {
                o.bodyControls[0].elements =
                    tetra ? algo::GeometryElements::Tet4 : algo::GeometryElements::Native;
                exercise(triangle ? (tetra ? "profiled prism to tet" : "profiled prism")
                                  : (tetra ? "profiled hex to tet" : "profiled hex"),
                         [&](auto &m) { return algo::meshGeometryModel(*model.model, o, m); });
            }
            o.bodyControls[0].elements = algo::GeometryElements::Native;
            mesh::Mesh native;
            need(algo::meshGeometryModel(*model.model, o, native).published);
            exercise(triangle ? "prism converter" : "hex converter",
                     [&](auto &m) { return algo::convertConvexCellsToTet4(native, {}, m); });
        }
    } catch (const std::exception &e) {
        remaining = -1;
        std::cerr << e.what() << '\n';
        return 1;
    }
}
