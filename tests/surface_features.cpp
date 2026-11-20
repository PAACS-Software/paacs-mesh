// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/SurfaceFeatures.hpp"
#include "primitive_checks.hpp"
#include <iostream>

using namespace primitive_test;
mesh::Mesh cube(bool triangles, double scale = 1) {
    mesh::Mesh m;
    for (auto p : {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{1, 1, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1},
                   Vec3{1, 0, 1}, Vec3{1, 1, 1}, Vec3{0, 1, 1}})
        m.addNode(p * scale);
    int tag = 0;
    for (auto q : {std::array<Index, 4>{0, 3, 2, 1},
                   {4, 5, 6, 7},
                   {0, 1, 5, 4},
                   {1, 2, 6, 5},
                   {2, 3, 7, 6},
                   {3, 0, 4, 7}}) {
        ++tag;
        if (triangles) {
            m.addCell(mesh::CellType::Tri3, {q[0], q[1], q[2]}, {2, tag});
            m.addCell(mesh::CellType::Tri3, {q[0], q[2], q[3]}, {2, tag});
        } else
            m.addCell(mesh::CellType::Quad4, q, {2, tag});
    }
    return m;
}
int main() {
    try {
        using R = algo::SurfaceFeatureReason;
        algo::SurfaceFeatureOptions o;
        unsigned successes = 0, refusals = 0;
        auto detect = [&](const mesh::Mesh &m, const auto &options) {
            const auto digest = hash(m);
            auto r = algo::detectSurfaceFeatures(m, options);
            require(r.completed && hash(m) == digest,
                    std::string("feature analysis: ") + algo::surfaceFeatureReasonName(r.reason));
            ++successes;
            return r;
        };
        auto reject = [&](const mesh::Mesh &m, const auto &options, R reason) {
            const auto digest = hash(m);
            auto r = algo::detectSurfaceFeatures(m, options);
            require(!r.completed && r.reason == reason && r.edges.empty() && hash(m) == digest,
                    std::string("feature refusal: ") + algo::surfaceFeatureReasonName(r.reason));
            ++refusals;
        };
        auto m = cube(true);
        auto r = detect(m, o);
        require(r.reason == R::Disabled && r.metrics.edgeUses == 0 && r.edges.empty(),
                "default off");
        o.enabled = true;
        for (bool triangles : {false, true})
            for (double scale : {1e-150, 1e-6, 1., 1e6, 1e150}) {
                auto model = cube(triangles, scale);
                r = detect(model, o);
                require(r.edges.size() == 12 && r.metrics.sharpEdges == 12 &&
                            r.metrics.tagBoundaryEdges == 12 && r.metrics.openBoundaryEdges == 0,
                        "cube exact features");
                std::pair<Index, Index> previous{-1, -1};
                for (const auto &edge : r.edges) {
                    require(std::abs(edge.length / scale - 1) < 1e-12 &&
                                std::abs(edge.angleDegrees - 90) < 1e-12,
                            "independent cube length/angle");
                    const std::pair key{edge.first, edge.second};
                    require(key > previous && edge.first < edge.second && edge.firstFace >= 0 &&
                                edge.secondFace >= 0,
                            "deterministic edge ordering");
                    previous = key;
                }
            }
        r = detect(m, o);
        auto exact = o;
        exact.maximumEdges = r.edges.size();
        exact.maximumOutputBytes = r.metrics.outputBytes;
        exact.maximumWorkingBytes = r.metrics.workingBytes;
        exact.maximumEdgeUses = r.metrics.edgeUses;
        require(detect(m, exact).edges.size() == 12, "exact budgets");
        auto bad = exact;
        --bad.maximumEdges;
        reject(m, bad, R::OutputBudget);
        bad = exact;
        --bad.maximumOutputBytes;
        reject(m, bad, R::OutputBudget);
        bad = exact;
        --bad.maximumWorkingBytes;
        reject(m, bad, R::WorkingBudget);
        bad = exact;
        --bad.maximumEdgeUses;
        reject(m, bad, R::WorkingBudget);
        bad = o;
        bad.minimumLength = 2;
        r = detect(m, bad);
        require(r.edges.empty() && r.metrics.filteredShortEdges == 12, "length filter");
        bad = o;
        bad.angleDegrees = 90;
        require(detect(m, bad).edges.empty(), "exact angle threshold excluded");
        bad.angleDegrees = 90 - 1e-9;
        require(detect(m, bad).edges.empty(), "angular tolerance band");
        bad.angleDegrees = 90 - 1e-6;
        require(detect(m, bad).edges.size() == 12, "outside tolerance band");
        bad.angleDegrees = 180;
        bad.includeTagBoundaries = true;
        r = detect(m, bad);
        require(r.edges.size() == 12 && r.metrics.sharpEdges == 0, "tag boundary option");
        std::atomic<bool> stop{true};
        bad = o;
        bad.cancellation = &stop;
        reject(m, bad, R::Cancelled);
        bad.enabled = false;
        bad.maximumWorkingBytes = 0;
        bad.angleDegrees = std::numeric_limits<double>::quiet_NaN();
        require(detect(m, bad).reason == R::Disabled, "disabled options ignored");
        for (double invalid : {-1., 181., std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
            bad = o;
            bad.angleDegrees = invalid;
            reject(m, bad, R::InvalidRequest);
            bad = o;
            bad.angularToleranceDegrees = invalid;
            reject(m, bad, R::InvalidRequest);
        }
        bad = o;
        bad.minimumLength = -1;
        reject(m, bad, R::InvalidRequest);
        bad.minimumLength = std::numeric_limits<double>::quiet_NaN();
        reject(m, bad, R::InvalidRequest);
        mesh::Mesh empty;
        reject(empty, o, R::EmptySurface);
        mesh::Mesh patch;
        for (auto p : {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{1, 1, 0}})
            patch.addNode(p);
        patch.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
        patch.addCell(mesh::CellType::Tri3, {1, 3, 2}, {2, 2});
        r = detect(patch, o);
        require(r.edges.size() == 4 && r.metrics.sharpEdges == 0 && r.metrics.tagBoundaryEdges == 1,
                "coplanar patch no diagonal feature");
        bad = o;
        bad.includeOpenBoundaries = false;
        require(detect(patch, bad).edges.empty(), "open-edge option");
        bad.includeTagBoundaries = true;
        require(detect(patch, bad).edges.size() == 1, "coplanar tag seam");
        patch.setNode(3, {1, 0, 1});
        r = detect(patch, o);
        require(r.metrics.sharpEdges == 1 && r.edges.size() == 5, "folded patch");
        auto broken = patch;
        broken.addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 3});
        bad = o;
        bad.minimumLength = 1e10;
        bad.angleDegrees = 180;
        reject(broken, bad, R::NonManifoldSurface); // topology cannot be filtered away
        mesh::Mesh winding;
        for (auto p : patch.nodes())
            winding.addNode(p);
        winding.addCell(mesh::CellType::Tri3, {0, 1, 2});
        winding.addCell(mesh::CellType::Tri3, {1, 2, 3});
        reject(winding, o, R::InconsistentWinding);
        mesh::Mesh collapsed;
        for (auto p : {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{2, 0, 0}})
            collapsed.addNode(p);
        collapsed.addCell(mesh::CellType::Tri3, {0, 1, 2});
        reject(collapsed, o, R::NumericallyUnresolved);
        collapsed.setNode(2, {0, 1, 0});
        collapsed.addCell(mesh::CellType::Tri3, {0, 0, 1});
        reject(collapsed, o, R::InvalidConnectivity);
        auto nonfinite = m;
        nonfinite.setNode(0, {std::numeric_limits<double>::infinity(), 0, 0});
        reject(nonfinite, o, R::NumericallyUnresolved);
        for (auto family : {algo::PrimitiveElements::HexDominant, algo::PrimitiveElements::Tet4,
                            algo::PrimitiveElements::Tet10})
            for (auto kind : {algo::PrimitiveKind::Cylinder, algo::PrimitiveKind::Torus,
                              algo::PrimitiveKind::Elbow}) {
                algo::PrimitiveShape s;
                s.kind = kind;
                algo::PrimitiveMeshingOptions p;
                p.elements = family;
                p.threads = 4;
                mesh::Mesh mesh;
                const auto meshed = algo::meshPrimitive(s, p, mesh);
                require(meshed.published, "feature fixture");
                r = detect(mesh, o);
                const auto expected =
                    kind == algo::PrimitiveKind::Torus ? 0 : 8 * meshed.segmentsPerQuarter;
                require(r.edges.size() == expected && r.metrics.sharpEdges == expected &&
                            r.metrics.openBoundaryEdges == 0,
                        "analytic rim features");
                require((r.metrics.highOrderFaces > 0) ==
                            (family == algo::PrimitiveElements::Tet10),
                        "corner high-order reporting");
                const auto again = detect(mesh, o);
                require(r.edges.size() == again.edges.size(), "repeat size");
                for (std::size_t i = 0; i < r.edges.size(); ++i)
                    require(r.edges[i].first == again.edges[i].first &&
                                r.edges[i].second == again.edges[i].second &&
                                r.edges[i].angleDegrees == again.edges[i].angleDegrees,
                            "repeat exact features");
            }
        std::cout << "{\"completed_checks\":" << successes << ",\"refusal_checks\":" << refusals
                  << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
