// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "geometry_checks.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/geometry/PolyhedralModel.hpp"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
using namespace paacs;
using geometry_test::require;
geo::PolyhedralInput voxels(const std::set<std::array<int, 3>> &cells) {
    geo::PolyhedralInput p;
    std::map<std::array<int, 3>, std::uint32_t> index;
    const std::array<std::array<int, 3>, 8> corners{
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    const std::array<std::array<int, 4>, 6> faces{
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}}};
    const std::array<std::array<int, 3>, 6> directions{
        {{0, 0, -1}, {0, 0, 1}, {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}}};
    for (const auto cell : cells)
        for (std::size_t f = 0; f < 6; ++f) {
            auto neighbor = cell;
            for (int j = 0; j < 3; ++j)
                neighbor[j] += directions[f][j];
            if (cells.contains(neighbor))
                continue;
            std::vector<std::uint32_t> ring;
            for (int c : faces[f]) {
                auto v = cell;
                for (int j = 0; j < 3; ++j)
                    v[j] += corners[c][j];
                auto [it, inserted] =
                    index.emplace(v, static_cast<std::uint32_t>(p.vertices.size()));
                if (inserted)
                    p.vertices.push_back({double(v[0]), double(v[1]), double(v[2])});
                ring.push_back(it->second);
            }
            p.faces.push_back({{ring}});
        }
    return p;
}
int main(int argc, char **argv) {
    try {
        const std::string name = argc > 1 ? argv[1] : "pocket_block";
        geo::PolyhedralInput input;
        double expected = 0;
        if (name == "tapered_hollow_duct") {
            const std::vector<std::vector<Vec2>> loops{{{0, 0}, {3, 0}, {3, 3}, {0, 3}},
                                                       {{1, 1}, {2, 1}, {2, 2}, {1, 2}}};
            input = geo::extrudedPolygon(loops, {0, 0, 3});
            for (std::size_t i = 8; i < 16; ++i) {
                input.vertices[i].x = .4 + .7 * input.vertices[i].x;
                input.vertices[i].y = .2 + .7 * input.vertices[i].y;
            }
            expected = 8 * 3 * (1 + .7 + .49) / 3;
        } else {
            std::set<std::array<int, 3>> cells;
            if (name == "pocket_block") {
                for (int k = 0; k < 2; ++k)
                    for (int j = 0; j < 3; ++j)
                        for (int i = 0; i < 3; ++i)
                            if (k == 0 || i != 1 || j != 1)
                                cells.insert({i, j, k});
            } else if (name == "branch_junction") {
                for (int k = 0; k < 3; ++k)
                    cells.insert({1, 1, k});
                for (int i = 0; i < 3; ++i)
                    cells.insert({i, 1, 2});
                for (int j = 0; j < 3; ++j)
                    cells.insert({1, j, 1});
            } else
                throw std::runtime_error("unknown fixture");
            expected = static_cast<double>(cells.size());
            input = voxels(cells);
        }
        auto model = geo::makePolyhedralModel(input);
        require(model.built(), std::string("model: ") + geo::polyhedralReasonName(model.reason));
        require(std::abs(model.volume - expected) < expected * 1e-9, "input volume");
        algo::GeometryMeshingOptions o;
        o.targetSize = .7;
        o.maximumMeshingWork = 500000;
        o.maximumNodes = 100000;
        o.maximumCells = 500000;
        o.threads = 1;
        mesh::Mesh m;
        m.addNode({99, 98, 97}, {0, 99});
        const auto before = geometry_test::digest(m);
        const auto start = std::chrono::steady_clock::now();
        const auto r = algo::meshGeometryModel(*model.model, o, m);
        if (argc > 2 && std::string(argv[2]) == "--require-published")
            require(r.published, "general tapered duct must publish");
        if (r.published) {
            require(r.audit.valid, "accepted unaudited mesh");
            require(std::abs(geometry_test::positiveVolume(m) - expected) < expected * 1e-8,
                    "accepted wrong volume");
            geometry_test::checkClassification(*model.model, m, 1e-7);
        } else
            require(r.destinationUnchanged && geometry_test::digest(m) == before,
                    "failed probe changed destination");
        std::cout << std::setprecision(17) << "{\"case\":\"" << name
                  << "\",\"published\":" << (r.published ? "true" : "false") << ",\"reason\":\""
                  << algo::geometryReasonName(r.reason) << "\",\"seconds\":"
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                  << ",\"nodes\":" << m.nodeCount() << ",\"cells\":" << m.cellCount()
                  << ",\"rollback_verified\":" << (!r.published ? "true" : "false");
        for (const auto &b : r.analysis.bodies)
            std::cout << ",\"strategy\":\"" << algo::geometryStrategyName(b.strategy)
                      << "\",\"engine_status\":" << b.engineStatus
                      << ",\"engine_reason\":" << b.engineReason
                      << ",\"engine_phase\":" << b.enginePhase
                      << ",\"boundary_status\":" << b.boundaryStatus
                      << ",\"boundary_reason\":" << b.boundaryReason
                      << ",\"boundary_phase\":" << b.boundaryPhase;
        std::cout << "}\n";
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
