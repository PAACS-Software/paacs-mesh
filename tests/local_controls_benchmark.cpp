// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/Delaunay2DScalable.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <charconv>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace paacs;
template <class T> T number(const char *s) {
    T value{};
    const std::string text = s;
    const auto r = std::from_chars(text.data(), text.data() + text.size(), value);
    if (r.ec != std::errc{} || r.ptr != text.data() + text.size())
        throw std::invalid_argument("invalid number");
    return value;
}
int main(int argc, char **argv) {
    try {
        if (argc != 6 && argc != 7)
            throw std::invalid_argument("usage: local-controls-benchmark mode base-size-or-grid "
                                        "local-size-or-first-height growth threads [archive]");
        const std::string mode = argv[1];
        const double base = number<double>(argv[2]), local = number<double>(argv[3]),
                     growth = number<double>(argv[4]);
        const int threads = number<int>(argv[5]);
        mesh::Mesh m;
        std::cout << std::setprecision(17);
        if (mode == "cube" || mode == "sphere" || mode == "face" || mode == "uniform-fine") {
            sizing::RegionSizeSpecification spec;
            spec.baseSize = mode == "uniform-fine" ? local : base;
            sizing::RegionSizeControl c;
            c.targetSize = local;
            c.growthRatio = growth;
            if (mode == "cube") {
                c.lower = {.375, .375, .375};
                c.upper = {.625, .625, .625};
            }
            if (mode == "sphere") {
                c.kind = sizing::SizeRegionKind::Ball;
                c.lower = {.5, .5, .5};
                c.radius = .16;
            }
            if (mode == "face")
                c = sizing::boxFaceSizeControl({0, 0, 0}, {1, 1, 1}, sizing::BoxFace::XMin, local,
                                               growth);
            if (mode != "uniform-fine")
                spec.controls.push_back(c);
            const auto field = sizing::compileRegionSizeField(spec);
            if (!field.compiled)
                throw std::invalid_argument("field refused");
            algo::AdaptiveBoxRequest q;
            q.field = &field.field;
            q.threads = threads;
            const auto begin = std::chrono::steady_clock::now();
            const auto r = algo::meshAdaptiveBox(q, m);
            const double wall =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            if (!r.published)
                throw std::runtime_error(algo::adaptiveBoxReasonName(r.reason));
            std::cout << "{\"case\":\"" << mode << "\",\"baseSize\":" << base
                      << ",\"localSize\":" << local << ",\"growth\":" << growth
                      << ",\"threads\":" << threads << ",\"nodes\":" << r.nodes
                      << ",\"volumeCells\":" << r.tetrahedra << ",\"leaves\":" << r.leaves
                      << ",\"balanceSplits\":" << r.balanceSplits
                      << ",\"minQuality\":" << r.minimumQuality
                      << ",\"actualGrowth\":" << r.audit.maximumAdjacentSizeRatio
                      << ",\"outputBytes\":" << r.outputBytes
                      << ",\"workingBytes\":" << r.workingBytes << ",\"wallSeconds\":" << wall
                      << ",\"planningSeconds\":" << r.planningSeconds
                      << ",\"generationSeconds\":" << r.generationSeconds
                      << ",\"validationSeconds\":" << r.validationSeconds
                      << ",\"auditSeconds\":" << r.auditSeconds << ",\"hash\":" << mesh::meshHash(m)
                      << "}\n";
        } else if (mode == "structured" || mode == "unstructured" || mode == "mixed") {
            if (!std::isfinite(base) || base < 1 || base > 500 || std::floor(base) != base)
                throw std::invalid_argument("grid must be 1..500");
            const auto n = static_cast<unsigned>(base);
            mesh::Mesh source;
            std::vector<Vec2> points;
            std::mt19937_64 rng(0x707269736d);
            std::uniform_real_distribution<double> jitter(-.28, .28);
            for (unsigned j = 0; j <= n; ++j)
                for (unsigned i = 0; i <= n; ++i) {
                    const bool move = mode == "unstructured" && i && j && i < n && j < n;
                    points.push_back(
                        {(i + (move ? jitter(rng) : 0)) / n, (j + (move ? jitter(rng) : 0)) / n});
                    source.addNode({points.back().u, points.back().v, 0}, {2, 1});
                }
            if (mode == "unstructured")
                for (const auto &t : algo::delaunay2DScalable(points))
                    source.addCell(mesh::CellType::Tri3, std::span<const Index>(t.v), {2, 1});
            else
                for (unsigned j = 0; j < n; ++j)
                    for (unsigned i = 0; i < n; ++i) {
                        const Index a = static_cast<Index>(j) * (n + 1) + i, b = a + 1,
                                    d = a + n + 1, c = d + 1;
                        if (mode == "mixed" && i < n / 2)
                            source.addCell(mesh::CellType::Quad4, {a, b, c, d}, {2, 1});
                        else {
                            source.addCell(mesh::CellType::Tri3, {a, b, c}, {2, 1});
                            source.addCell(mesh::CellType::Tri3, {a, c, d}, {2, 1});
                        }
                    }
            algo::LayeredPatchRequest q;
            q.threads = threads;
            q.schedule.firstLayerThickness = local;
            q.schedule.growthRatio = growth;
            q.schedule.maximumAllowedGrowth = std::max(1.3, growth);
            const auto begin = std::chrono::steady_clock::now();
            const auto r = algo::meshLayeredPatch(source, q, m);
            const double seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            if (!r.published)
                throw std::runtime_error("layer request refused: " + std::to_string(int(r.reason)));
            std::cout << "{\"case\":\"" << mode << "\",\"grid\":" << n
                      << ",\"firstHeight\":" << local << ",\"growth\":" << growth
                      << ",\"threads\":" << threads << ",\"nodes\":" << m.nodeCount()
                      << ",\"volumeCells\":" << r.sweep.metrics.wedges + r.sweep.metrics.hexes
                      << ",\"wedges\":" << r.sweep.metrics.wedges
                      << ",\"hexes\":" << r.sweep.metrics.hexes
                      << ",\"wallLayers\":" << r.schedule.wallLayers
                      << ",\"coreLayers\":" << r.schedule.coreLayers
                      << ",\"actualGrowth\":" << r.schedule.maximumAdjacentGrowth
                      << ",\"outputBytes\":" << r.outputBytes
                      << ",\"workingBytes\":" << r.plannedWorkingBytes
                      << ",\"sourcePairChecks\":" << r.sourcePairChecks
                      << ",\"wallSeconds\":" << seconds << ",\"hash\":" << mesh::meshHash(m)
                      << "}\n";
        } else
            throw std::invalid_argument("unknown mode");
        if (argc == 7 && !sdk::writeMeshArchive(m, argv[6]).completed)
            throw std::runtime_error("archive write failed");
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
