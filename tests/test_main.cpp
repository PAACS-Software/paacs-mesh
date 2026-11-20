// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// PAACS Mesh unit tests - minimal harness, no external framework.
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#if defined(_WIN32)
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include "paacs/accel/Backend.hpp"
#include "paacs/algo/Delaunay2D.hpp"
#include "paacs/algo/Delaunay3D.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/TetArena.hpp"
#include "paacs/topo/ElementTraits.hpp"
#include "paacs/topo/Handles.hpp"
#include "paacs/topo/Storage.hpp"
#include "paacs/topo/EdgeTopology.hpp"
#include "paacs/topo/FaceTopology.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/core/BuildInfo.hpp"
#include "paacs/algo/Meshers.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/core/PointGrid.hpp"
#include "paacs/algo/Optimize.hpp"
#include "paacs/algo/Smoothing.hpp"
#include "paacs/mesh/Quality.hpp"

#include "test_harness.hpp"

int gFailures = 0;   // shared harness counter (defined once, here)

// CHECK/TEST/canonU01 come from test_harness.hpp.

using namespace paacs;

// ---------------------------------------------------------------------------
static void testEdgeDiscretization() {
    TEST("edge discretization matches sizing field");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});

    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.25);

    mesh::Mesh m;
    algo::meshEdges(*model, spec, m);

    // 12 unit edges at h=0.25 -> 4 segments each = 48 Line2 cells,
    // 8 corner nodes + 12*3 interior edge nodes = 44 unique nodes.
    CHECK(m.countCells(mesh::CellType::Line2) == 48);
    CHECK(m.nodeCount() == 44);

    // Corner nodes classified on model vertices (dim 0).
    int dim0 = 0;
    for (Index i = 0; i < m.nodeCount(); ++i)
        if (m.nodeClass(i).dim == 0) ++dim0;
    CHECK(dim0 == 8);
}

// ---------------------------------------------------------------------------
static void testDelaunay2DEmptyCircumcircle() {
    TEST("2D Delaunay: empty circumcircle property");
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<Vec2> pts;
    for (int i = 0; i < 60; ++i) pts.push_back({U(rng), U(rng)});

    auto tris = algo::delaunay2D(pts);
    CHECK(!tris.empty());

    int violations = 0;
    for (const auto& t : tris)
        for (size_t k = 0; k < pts.size(); ++k) {
            if (k == size_t(t.v[0]) || k == size_t(t.v[1]) || k == size_t(t.v[2]))
                continue;
            if (algo::inCircumcircle(pts[size_t(t.v[0])], pts[size_t(t.v[1])],
                                     pts[size_t(t.v[2])], pts[k]))
                ++violations;
        }
    CHECK(violations == 0);

    // Convex hull of random points: E = 3n - 3 - h, F = 2n - 2 - h (Euler).
    // Just sanity-check the triangle count is in the plausible band.
    CHECK(tris.size() > pts.size());
    CHECK(tris.size() < 2 * pts.size());
}

// ---------------------------------------------------------------------------
static void testDelaunay3DValidity() {
    TEST("3D Delaunay: positive volumes, hull covered");
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<Vec3> pts;
    // Cube corners + random interior points.
    for (int i = 0; i < 8; ++i)
        pts.push_back({double(i & 1), double((i >> 1) & 1), double((i >> 2) & 1)});
    for (int i = 0; i < 40; ++i) pts.push_back({U(rng), U(rng), U(rng)});

    auto tets = algo::delaunay3D(pts);
    CHECK(!tets.empty());

    double vol = 0;
    int flat = 0;
    for (const auto& t : tets) {
        const double v = mesh::tetVolume(pts[size_t(t.v[0])], pts[size_t(t.v[1])],
                                         pts[size_t(t.v[2])], pts[size_t(t.v[3])]);
        if (std::abs(v) < 1e-12) ++flat;
        vol += std::abs(v);
    }
    CHECK(flat == 0);
    CHECK(std::abs(vol - 1.0) < 1e-6);  // Delaunay of points incl. cube = cube
}

// ---------------------------------------------------------------------------
static void testVolumeMeshBox() {
    TEST("full pipeline: box volume mesh");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});

    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.25);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);

    CHECK(m.countCells(mesh::CellType::Tet4) > 0);
    CHECK(std::abs(mesh::totalTetVolume(m) - 1.0) < 1e-6);

    auto q = mesh::tetQualityStats(m);
    CHECK(q.min > 0.05);   // no near-degenerate slivers
    CHECK(q.mean > 0.5);   // healthy average

    // Every tet classified on the solid.
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).type == mesh::CellType::Tet4)
            CHECK(m.cell(i).onEntity.dim == 3 && m.cell(i).onEntity.tag == 1);
}

// ---------------------------------------------------------------------------
static void testRefinementBall() {
    TEST("sizing: ball refinement produces locally smaller elements");
    auto field = std::make_shared<sizing::MinField>();
    field->add(std::make_shared<sizing::ConstantField>(0.3));
    field->add(std::make_shared<sizing::BallRefinement>(
        Vec3{0, 0, 0}, 0.25, 0.08, 0.3));
    CHECK(std::abs(field->eval({0, 0, 0}) - 0.08) < 1e-12);
    CHECK(std::abs(field->eval({5, 5, 5}) - 0.3) < 1e-12);

    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = field;

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);

    // Edges touching the refined corner get more segments than far edges.
    // Edge 1 runs (0,0,0)->(1,0,0) [refined end], edge 6 (1,0,1)->(1,1,1) far.
    // Count Line2 cells per edge tag.
    int segsEdge1 = 0, segsEdge7 = 0;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i);
        if (c.type != mesh::CellType::Line2) continue;
        if (c.onEntity.tag == 1) ++segsEdge1;
        if (c.onEntity.tag == 7) ++segsEdge7;
    }
    CHECK(segsEdge1 > segsEdge7);
    CHECK(std::abs(mesh::totalTetVolume(m) - 1.0) < 1e-6);
}

// ---------------------------------------------------------------------------
static void testQuadraticElevation() {
    TEST("order elevation: tet4 -> tet10 with shared midside dedup");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.5);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const Index nTet   = m.countCells(mesh::CellType::Tet4);
    const Index nNodes = m.nodeCount();

    // Count unique tet edges before elevation.
    std::set<std::pair<Index, Index>> edges;
    for (Index i = 0; i < m.cellCount(); ++i) {
        if (m.cell(i).type != mesh::CellType::Tet4) continue;
        const Index* v = m.cellNodes(i);
        static const int e[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
        for (const auto& ed : e)
            edges.insert(std::minmax(v[ed[0]], v[ed[1]]));
    }

    algo::elevateToQuadratic(*model, m);

    CHECK(m.countCells(mesh::CellType::Tet10) == nTet);
    CHECK(m.countCells(mesh::CellType::Tet4) == 0);
    // Midside dedup: new nodes <= unique edge count of the whole mesh
    // (tri/line midsides coincide with tet edge midsides on the shell).
    CHECK(m.nodeCount() > nNodes);
    CHECK(m.nodeCount() <= nNodes + static_cast<Index>(edges.size())
                          + 8);  // small slack for surface-only edges

    // Quadratic order flag coarsens the sizing (order-aware sizing).
    sizing::MeshSpec p2 = spec;
    p2.order = 2;
    CHECK(p2.h({0.5, 0.5, 0.5}) > spec.h({0.5, 0.5, 0.5}));
}

// ---------------------------------------------------------------------------
static void testVtkExport() {
    TEST("VTK export produces a valid file");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.5);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    // Portable temp path: /tmp does not exist on Windows (a cross-platform
    // test defect that fast-failed the MSVC io group). Use the platform's
    // temp directory and remove the artifact afterward.
    const std::string path =
        (std::filesystem::temp_directory_path() /
         ("paacs_test_out_" + std::to_string(::getpid()) + ".vtk")).string();
    io::writeVtk(m, path);

    {
        std::ifstream in(path);
        CHECK(in.good());
        std::string line;
        std::getline(in, line);
        CHECK(line.rfind("# vtk DataFile", 0) == 0);
        bool hasPoints = false, hasCells = false;
        while (std::getline(in, line)) {
            if (line.rfind("POINTS", 0) == 0) hasPoints = true;
            if (line.rfind("CELLS", 0) == 0)  hasCells = true;
        }
        CHECK(hasPoints && hasCells);
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
static void testQualityMetrics() {
    TEST("quality metrics: regular and degenerate tets");
    // Regular tetrahedron with unit edge.
    const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0.5, std::sqrt(3.0) / 2, 0},
               d{0.5, std::sqrt(3.0) / 6, std::sqrt(2.0 / 3.0)};

    CHECK(std::abs(mesh::tetQuality(a, b, c, d) - 1.0) < 1e-9);
    CHECK(std::abs(mesh::tetScaledJacobian(a, b, c, d) - 1.0) < 1e-9);
    CHECK(mesh::tetSkewness(a, b, c, d) < 1e-9);

    double dmin, dmax;
    mesh::tetDihedralAngles(a, b, c, d, dmin, dmax);
    CHECK(std::abs(dmin - 70.5288) < 0.01);   // arccos(1/3)
    CHECK(std::abs(dmax - 70.5288) < 0.01);

    // Near-degenerate sliver: 4th vertex almost in the base plane.
    const Vec3 dFlat{0.5, std::sqrt(3.0) / 6, 1e-4};
    CHECK(mesh::tetQuality(a, b, c, dFlat) < 0.01);
    CHECK(mesh::tetScaledJacobian(a, b, c, dFlat) < 0.01);
    CHECK(mesh::tetSkewness(a, b, c, dFlat) > 0.9);

    // Inverted tet: negative scaled Jacobian.
    CHECK(mesh::tetScaledJacobian(a, c, b, d) < 0);
}

// ---------------------------------------------------------------------------
static void testSmoothing() {
    TEST("smoothing: quality improves, volume and boundary preserved");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.2);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);

    // Perturb interior nodes to create bad elements smoothing must fix.
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> J(-0.03, 0.03);
    std::vector<Vec3> boundaryBefore;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        if (m.nodeClass(i).dim == 3)
            m.setNode(i, m.node(i) + Vec3{J(rng), J(rng), J(rng)});
        else
            boundaryBefore.push_back(m.node(i));
    }

    const auto before = mesh::analyzeTets(m);
    CHECK(before.nInverted == 0);   // perturbation kept mesh valid

    auto sr = algo::smoothInterior(m);
    const auto after = mesh::analyzeTets(m);

    CHECK(sr.acceptedMoves > 0);
    CHECK(after.radiusRatio.min >= before.radiusRatio.min - 1e-12);
    CHECK(after.radiusRatio.mean > before.radiusRatio.mean);
    CHECK(after.nInverted == 0);
    CHECK(std::abs(after.totalVolume - 1.0) < 1e-9);   // tiling preserved
    CHECK(after.minDihedralDeg > before.minDihedralDeg - 1e-9);

    // Boundary nodes must not have moved.
    size_t k = 0;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        if (m.nodeClass(i).dim == 3) continue;
        CHECK((m.node(i) - boundaryBefore[k]).norm() < 1e-15);
        ++k;
    }
}

// ---------------------------------------------------------------------------
static void testQualityReport() {
    TEST("quality report: consistent counts and histogram");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.3);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const auto r = mesh::analyzeTets(m);

    CHECK(r.nTets == m.countCells(mesh::CellType::Tet4));
    CHECK(r.nInverted == 0);   // mesher outputs positively oriented tets
    Index histSum = 0;
    for (Index h : r.radiusRatioHist) histSum += h;
    CHECK(histSum == r.nTets);
    CHECK(r.skewness.min >= 0 && r.skewness.max <= 1);
    CHECK(r.minDihedralDeg > 0 && r.maxDihedralDeg < 180);
    CHECK(!r.toString().empty());
}

// ---------------------------------------------------------------------------
static void testCylinderPipeline() {
    TEST("cylinder: curved geometry through full pipeline");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);

    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.15);
    spec.elementsPerCircle = 24;

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const auto r = mesh::analyzeTets(m);

    CHECK(r.nTets > 0);
    CHECK(r.nInverted == 0);
    // Chordal volume slightly under pi*R^2*H = 0.7854 (boundary sliver
    // peeling shaves ~0.5% more until constrained recovery lands).
    CHECK(r.totalVolume > 0.70 && r.totalVolume < 0.7854 + 1e-6);
    CHECK(r.radiusRatio.min > 5e-4);   // peeled: no flat hull slivers
    CHECK(r.radiusRatio.mean > 0.4);

    // Optimization (flips + smoothing) must lift the floor decisively.
    {
        algo::OptimizeOptions oo;
        oo.peelableFaceTags = algo::curvedFaceTags(*model);
        algo::optimizeTets(m, oo);
    }
    const auto ro = mesh::analyzeTets(m);
    CHECK(ro.radiusRatio.min > 0.05);
    CHECK(ro.radiusRatio.mean > 0.75);
    CHECK(ro.nInverted == 0);

    // Curvature adaptation: each bottom half-circle (edges 1,2) must carry
    // at least elementsPerCircle/2 segments.
    int segsArc1 = 0;
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).type == mesh::CellType::Line2 &&
            m.cell(i).onEntity.tag == 1)
            ++segsArc1;
    CHECK(segsArc1 >= 12);

    // All lateral-surface nodes lie exactly on radius R.
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto& c = m.nodeClass(i);
        if (c.dim == 2 && (c.tag == 1 || c.tag == 2)) {
            const Vec3& p = m.node(i);
            CHECK(std::abs(std::hypot(p.x, p.y) - 0.5) < 1e-9);
        }
    }
}

// ---------------------------------------------------------------------------
static void testCurvedQuadratic() {
    TEST("quadratic on cylinder: midside nodes snapped to true surface");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.2);
    spec.elementsPerCircle = 16;

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const Index nBefore = m.nodeCount();
    algo::elevateToQuadratic(*model, m);

    // Every new node classified on the lateral faces or circle arcs must be
    // on radius R, and at least one must differ from the linear midpoint
    // (i.e. the elements are genuinely curved).
    // Validity straightening may pull a few midsides slightly off the
    // surface (fidelity yields to a positive Jacobian); the vast majority
    // must sit exactly on radius R, the rest within half the local sagitta.
    int onSurface = 0, offSurface = 0, curved = 0;
    for (Index i = nBefore; i < m.nodeCount(); ++i) {
        const auto& c = m.nodeClass(i);
        const bool onCurve =
            (c.dim == 2 && (c.tag == 1 || c.tag == 2)) ||
            (c.dim == 1 && c.tag >= 1 && c.tag <= 4);
        if (!onCurve) continue;
        const double dr = std::abs(std::hypot(m.node(i).x, m.node(i).y) - 0.5);
        if (dr < 1e-9) ++onSurface;
        else { ++offSurface; CHECK(dr < 0.02); }
    }
    CHECK(onSurface > 0);
    CHECK(onSurface >= 3 * offSurface);   // strong majority exactly on CAD

    // Chord midpoints of boundary edges are strictly inside radius R, so a
    // snapped midside must have moved: verify via any tet10 midside on the
    // lateral surface.
    for (Index ci = 0; ci < m.cellCount() && !curved; ++ci) {
        if (m.cell(ci).type != mesh::CellType::Tet10) continue;
        const Index* v = m.cellNodes(ci);
        static const int mid[6][3] = {{4,0,1},{5,1,2},{6,2,0},
                                      {7,0,3},{8,1,3},{9,2,3}};
        for (const auto& e : mid) {
            const Vec3 linMid = (m.node(v[e[1]]) + m.node(v[e[2]])) * 0.5;
            if ((m.node(v[e[0]]) - linMid).norm() > 1e-6) { curved = 1; break; }
        }
    }
    CHECK(curved == 1);
}

// ---------------------------------------------------------------------------
static void testFlipOptimization() {
    TEST("optimization: flips lift minimum quality beyond smoothing");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.25);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const auto before = mesh::analyzeTets(m);

    auto res = algo::optimizeTets(m);
    const auto after = mesh::analyzeTets(m);

    CHECK(res.flips23 + res.flips32 > 0);            // topology changed
    CHECK(after.nInverted == 0);
    CHECK(std::abs(after.totalVolume - 1.0) < 1e-6); // tiling preserved
    CHECK(after.radiusRatio.min >= before.radiusRatio.min);
    CHECK(after.radiusRatio.mean >= before.radiusRatio.mean - 1e-9);
}

// ---------------------------------------------------------------------------
static void testMshExport() {
    TEST("MSH 4.1 export");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.5);

    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    // Portable temp path (see testVtkExport): /tmp is Unix-only.
    const std::string path =
        (std::filesystem::temp_directory_path() /
         ("paacs_test_out_" + std::to_string(::getpid()) + ".msh")).string();
    io::writeMsh(m, path);

    {
        std::ifstream in(path);
        CHECK(in.good());
        std::string line, all;
        bool fmt = false;
        while (std::getline(in, line)) {
            if (line == "4.1 0 8") fmt = true;
            all += line + "\n";
        }
        CHECK(fmt);
        for (const char* sec : {"$Entities", "$Nodes", "$Elements",
                                "$EndElements"})
            CHECK(all.find(sec) != std::string::npos);
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
static void testPointGrid() {
    TEST("point grid: radius queries match brute force");
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<Vec3> pts;
    core::PointGrid grid(0.13);
    for (int i = 0; i < 500; ++i) {
        Vec3 p{U(rng), U(rng), U(rng)};
        pts.push_back(p);
        grid.insert(p);
    }
    for (int i = 0; i < 200; ++i) {
        const Vec3 q{U(rng), U(rng), U(rng)};
        const double r = 0.02 + 0.3 * U(rng);
        bool brute = false;
        for (const auto& p : pts)
            brute = brute || (p - q).norm2() < r * r;
        CHECK(grid.anyWithin(q, r) == brute);
    }
}

// ---------------------------------------------------------------------------
static void testThreadDeterminism() {
    TEST("parallel: results bit-identical across thread counts");
    auto runOnce = [](int threads) {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.15);
        spec.threads = threads;
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        {
            algo::OptimizeOptions oo;
            oo.peelableFaceTags = algo::curvedFaceTags(*model);
            algo::optimizeTets(m, oo);
        }
        return std::make_tuple(m.nodeCount(),
                               m.countCells(mesh::CellType::Tet4),
                               mesh::analyzeTets(m).totalVolume,
                               mesh::analyzeTets(m).radiusRatio.min);
    };
    const auto serial = runOnce(1);
    const auto four   = runOnce(4);
    const auto autoN  = runOnce(0);
    CHECK(serial == four);   // exact tuple equality, doubles included
    CHECK(serial == autoN);
    core::setThreadCount(0);
}

// ---------------------------------------------------------------------------
static void testBakedField() {
    TEST("baked size field: trilinear grid tracks the analytic field");
    auto analytic = std::make_shared<sizing::BallRefinement>(
        Vec3{0.5, 0.5, 0.5}, 0.4, 0.05, 0.3);
    sizing::BakedGridField baked(*analytic, {0, 0, 0}, {1, 1, 1},
                                 64, 64, 64);
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(0.02, 0.98);
    double worst = 0;
    for (int i = 0; i < 500; ++i) {
        const Vec3 p{U(rng), U(rng), U(rng)};
        worst = std::max(worst,
                         std::abs(baked.eval(p) - analytic->eval(p)));
    }
    CHECK(worst < 0.01);   // grid resolution error only

    // Accelerator backend batch evaluation agrees with direct eval.
    auto backend = accel::makeCpuBackend();
    std::vector<Vec3> pts{{0.1, 0.2, 0.3}, {0.5, 0.5, 0.5}, {0.9, 0.9, 0.9}};
    std::vector<float> out;
    backend->evalField(baked, pts, out);
    CHECK(out.size() == 3);
    for (size_t k = 0; k < pts.size(); ++k)
        CHECK(std::abs(out[k] - baked.eval(pts[k])) < 1e-6);
}


// ---------------------------------------------------------------------------
static void testMidsideClassification() {
    TEST("elevation: straight-box midsides stay at linear midpoints");
    // Regression for the midside classification bug: an interior mesh edge
    // joining nodes on two different CAD entities must classify onto the
    // owning cell's entity, never get snapped onto one endpoint's entity.
    // On a box every entity is flat/straight, so EVERY midside must equal
    // its linear midpoint (within edge-projection tolerance ~1e-8).
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.2);
    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    algo::elevateToQuadratic(*model, m);

    static const int mid[6][3] = {{4,0,1},{5,1,2},{6,2,0},
                                  {7,0,3},{8,1,3},{9,2,3}};
    double worst = 0;
    for (Index ci = 0; ci < m.cellCount(); ++ci) {
        if (m.cell(ci).type != mesh::CellType::Tet10) continue;
        const Index* v = m.cellNodes(ci);
        for (const auto& e : mid)
            worst = std::max(worst,
                (m.node(v[e[0]]) -
                 (m.node(v[e[1]]) + m.node(v[e[2]])) * 0.5).norm());
    }
    CHECK(worst < 1e-6);
}

// ---------------------------------------------------------------------------
static void testTet10Jacobian() {
    TEST("tet10 Jacobian: straight valid, distorted detected, mesh clean");
    // Straight reference-style tet10: ratio must be ~1.
    Vec3 p[10] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},
                  {.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
    CHECK(std::abs(mesh::tet10MinJacobianRatio(p) - 1.0) < 1e-12);

    // Fold the map: pull the edge-01 midside toward corner 0. Corner tet
    // volume stays positive, but the edge parameterization reverses locally.
    // (A perpendicular bulge is shear-like and det-preserving, so it is
    // correctly NOT flagged - shape metrics, not validity, catch that.)
    p[4] = {0.1, 0.05, 0.0};
    CHECK(mesh::tet10MinJacobianRatio(p) <= 0.0);

    // Full curved cylinder pipeline: snapping + validity straightening must
    // leave zero invalid curved elements.
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.25);
    spec.order = 2;
    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    const auto er = algo::elevateToQuadratic(*model, m);
    CHECK(er.invalidAfter == 0);
    CHECK(mesh::countInvalidTet10(m, 0.0) == 0);
    CHECK(er.snapped > 0);          // elements are genuinely curved
}


// ---------------------------------------------------------------------------
static void testDedupLifecycle() {
    TEST("mesh dedup: fixed tolerance, index follows moved nodes");
    mesh::Mesh m;
    const Index a = m.addNodeUnique({0, 0, 0});
    CHECK(m.addNodeUnique({0, 0, 5e-10}) == a);        // dedups

    // Tolerance is latched: re-setting a different value must throw.
    bool threw = false;
    try { m.setDedupTolerance(1e-3); }
    catch (const std::logic_error&) { threw = true; }
    CHECK(threw);

    // Moving a node dirties the index; queries see the NEW position.
    m.setNode(a, {2, 2, 2});
    CHECK(m.addNodeUnique({2, 2, 2 + 5e-10}) == a);    // found at new spot
    const Index b = m.addNodeUnique({0, 0, 0});        // old spot now free
    CHECK(b != a);

    // addCell validates connectivity.
    threw = false;
    try { m.addCell(mesh::CellType::Tet4, {a, b, 99999, a}); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try { m.addCell(mesh::CellType::Tet4, {a, b}); }   // wrong count
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}


// ---------------------------------------------------------------------------
static void testQualityReportSafety() {
    TEST("quality report: empty, tet10-only and mixed meshes are safe");
    // Empty mesh: explicit state, no sentinel numbers (PM-QUAL-001).
    mesh::Mesh empty;
    const std::string es = mesh::analyzeTets(empty).toString();
    CHECK(es.find("no tetrahedra") != std::string::npos);
    CHECK(es.find("1e+300") == std::string::npos);

    // Tet10-only mesh: elevation replaces every Tet4, so the report must
    // describe the quadratic population, not claim zero tets with sentinels.
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.35);
    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    algo::elevateToQuadratic(*model, m);
    const auto r = mesh::analyzeTets(m);
    CHECK(r.nTets == 0);
    CHECK(r.nTet10 > 0);
    CHECK(r.nTet10Invalid == 0);
    CHECK(std::abs(r.tet10CornerVolume - 1.0) < 1e-6);   // straight box
    CHECK(r.tet10Jacobian.n == r.nTet10);
    CHECK(r.tet10Jacobian.min > 0.99);                   // straight elements
    const std::string s10 = r.toString();
    CHECK(s10.find("tet10") != std::string::npos);
    CHECK(s10.find("1e+300") == std::string::npos);

    // MIXED Tet4/Tet10: append a standalone straight Tet4 to the elevated
    // mesh; both populations must be reported, separately and sanely.
    {
        const Index a = m.addNode({5, 5, 5}), b = m.addNode({6, 5, 5}),
                    c = m.addNode({5, 6, 5}), d = m.addNode({5, 5, 6});
        m.addCell(mesh::CellType::Tet4, {a, b, c, d});
        const auto mr = mesh::analyzeTets(m);
        CHECK(mr.nTets == 1);
        CHECK(mr.nTet10 > 0);
        CHECK(mr.nInverted == 0 && mr.nTet10Invalid == 0);
        const std::string ms = mr.toString();
        CHECK(ms.find("tet4 1") != std::string::npos);
        CHECK(ms.find("tet10") != std::string::npos);
    }

    // NON-FINITE coordinates: a NaN node must be counted and excluded, not
    // silently poison the statistics or the formatted report.
    {
        mesh::Mesh bad;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const Index a = bad.addNode({0, 0, 0}), b = bad.addNode({1, 0, 0}),
                    c = bad.addNode({0, 1, 0}), d = bad.addNode({nan, 0, 1});
        bad.addCell(mesh::CellType::Tet4, {a, b, c, d});
        const Index e = bad.addNode({3, 0, 0}), f = bad.addNode({4, 0, 0}),
                    g = bad.addNode({3, 1, 0}), h = bad.addNode({3, 0, 1});
        bad.addCell(mesh::CellType::Tet4, {e, f, g, h});
        const auto br = mesh::analyzeTets(bad);
        CHECK(br.nNonFinite == 1);
        CHECK(br.nTets == 1);                    // only the finite tet
        CHECK(std::isfinite(br.totalVolume));
        CHECK(std::isfinite(br.radiusRatio.min));
        const std::string bs = br.toString();
        CHECK(bs.find("non-finite") != std::string::npos);
        CHECK(bs.find("nan") == std::string::npos);
    }
}


// ---------------------------------------------------------------------------
static void testScaleAndTranslationInvariance() {
    TEST("numerics: translated (1e9) and tiny (1e-10) boxes mesh correctly");
    // PM-NUM-001/002 regression. Before the normalized working frames these
    // produced radically different topology / wrong volume (1e9 offset) and
    // no usable mesh at all (1e-10 scale). Exact invariance needs exact
    // predicates (M3); the closure bar here is: correct volume, valid
    // elements, and tet counts in the same regime as the reference box.
    auto meshBox = [&](const Vec3& lo, double side, double h) {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeBox(lo, lo + Vec3{side, side, side});
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(h * side);
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        return mesh::analyzeTets(m);
    };

    const auto ref = meshBox({0, 0, 0}, 1.0, 0.3);
    CHECK(ref.nTets > 0 && ref.nInverted == 0);
    CHECK(std::abs(ref.totalVolume - 1.0) < 1e-9);

    // Large translation: coordinates near 1e9 (ulp ~1e-7 absolute).
    const auto far = meshBox({1e9, 1e9, 1e9}, 1.0, 0.3);
    CHECK(far.nInverted == 0);
    CHECK(std::abs(far.totalVolume - 1.0) < 1e-3);
    CHECK(far.nTets > ref.nTets / 2 && far.nTets < ref.nTets * 2);
    CHECK(far.radiusRatio.min > 0.01);

    // Tiny model: side 1e-10 (previously collapsed by absolute tolerances).
    const auto tiny = meshBox({0, 0, 0}, 1e-10, 0.3);
    CHECK(tiny.nTets > 0 && tiny.nInverted == 0);
    CHECK(std::abs(tiny.totalVolume - 1e-30) < 1e-33);
    CHECK(tiny.nTets > ref.nTets / 2 && tiny.nTets < ref.nTets * 2);

    // Large model: side 1e6.
    const auto big = meshBox({0, 0, 0}, 1e6, 0.3);
    CHECK(big.nTets > 0 && big.nInverted == 0);
    CHECK(std::abs(big.totalVolume - 1e18) < 1e12);
}


// ---------------------------------------------------------------------------
static void testThinDomains() {
    TEST("robustness: thin boxes 1x1x{1e-2,1e-3,1e-4} mesh end-to-end");
    // Reviewer finding: t <= 1e-3 silently produced ZERO tets. Two defects:
    // the finite-coordinate super-tet (fixed by the symbolic infinite
    // vertex) and quality-threshold peeling eating legitimately thin
    // domains (fixed by the curved-face peel gate). Quality is honestly
    // poor here - isotropic meshing of thin plates awaits anisotropy - but
    // the mesh must be complete, valid, and volume-exact.
    for (double t : {1e-2, 1e-3, 1e-4}) {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeBox({0, 0, 0}, {1, 1, t});
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.2);
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        const auto r = mesh::analyzeTets(m);
        CHECK(r.nTets > 0);
        CHECK(r.nInverted == 0);
        CHECK(std::abs(r.totalVolume - t) < 1e-9 * t);
    }
}

// ---------------------------------------------------------------------------
static void testFlattenedClouds() {
    TEST("robustness: flattened random clouds (z x 1e-6 and 1e-9) triangulate");
    // Reviewer stress: 59/80 failures at 1e-6, 80/80 at 1e-9 with the old
    // super-tet construction. The symbolic infinite vertex + exact
    // predicates must insert every point and produce a valid positive-
    // volume complex.
    for (double zf : {1e-6, 1e-9}) {
        for (unsigned seed : {11u, 12u, 13u}) {
            std::mt19937 rng(seed);
            std::uniform_real_distribution<double> U(0, 1);
            std::vector<Vec3> pts;
            for (int i = 0; i < 150; ++i)
                pts.push_back({U(rng), U(rng), U(rng) * zf});
            algo::Delaunay3DStats st;
            const auto tets = algo::delaunay3D(pts, &st);
            CHECK(st.pointsInserted + st.pointsDuplicate == 150);
            CHECK(st.pointsSkipped == 0);
            CHECK(st.auditPassed);
            CHECK(!tets.empty());
            double vol = 0;
            for (const auto& tt : tets)
                vol += mesh::tetVolume(pts[size_t(tt.v[0])], pts[size_t(tt.v[1])],
                                       pts[size_t(tt.v[2])], pts[size_t(tt.v[3])]);
            CHECK(vol > 0 && std::isfinite(vol));
        }
    }
}

// ---------------------------------------------------------------------------
static void testExactLattice() {
    TEST("robustness: exact 5^3 lattice tiles its hull (max cosphericality)");
    // Every insertion on a Cartesian grid hits exactly-cospherical and
    // exactly-coplanar configurations - the hardest degeneracy class. The
    // include-on-zero policy + star-shape trimming must tile the cube
    // exactly, insert every point, and keep reciprocal adjacency.
    std::vector<Vec3> pts;
    for (int k = 0; k < 5; ++k)
        for (int j = 0; j < 5; ++j)
            for (int i = 0; i < 5; ++i)
                pts.push_back({double(i), double(j), double(k)});
    algo::Delaunay3DStats st;
    algo::Delaunay3DOptions options;
    options.collectDetailedOperationCounts = true;
    const auto tets = algo::delaunay3D(pts, options, &st);
    CHECK(st.pointsInserted == 125);
    CHECK(st.pointsSkipped == 0);
    CHECK(st.auditPassed);
    double vol = 0;
    for (const auto& tt : tets)
        vol += std::abs(mesh::tetVolume(pts[size_t(tt.v[0])], pts[size_t(tt.v[1])],
                                        pts[size_t(tt.v[2])], pts[size_t(tt.v[3])]));
    CHECK(std::abs(vol - 64.0) < 1e-9);
    // The compact kernel uses instrumentation-free predicate entry points.
    // Request its local operation counters, rather than relying on unrelated
    // process-global calls left by earlier tests. Public predicate counters
    // are checked separately by predicates.conventions.
    CHECK(st.collectedDetailedOperationCounts);
    CHECK(st.conflictInsphereTests > 0);
    CHECK(st.bootstrapOrientationTests + st.walkOrientationTests +
          st.conflictOrientationTests + st.trimOrientationTests > 0);
}

// ---------------------------------------------------------------------------
static void testFullMeshHashDeterminism() {
    TEST("determinism: full mesh hash identical across thread counts");
    // The older determinism test compared four aggregates (counts, volume,
    // min quality), which cannot prove identical topology. Hash everything:
    // node coordinate bits, classifications, cell types, connectivity.
    // Uses the library reproducibility hash (M5): the audit and this test
    // must agree on what "identical mesh" means.
    auto hashMesh = [](const mesh::Mesh& m) { return mesh::meshHash(m); };
    std::uint64_t ref = 0;
    bool first = true;
    for (int threads : {1, 2, 4, 0}) {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.15);
        spec.threads = threads;
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        algo::OptimizeOptions oo;
        oo.peelableFaceTags = algo::curvedFaceTags(*model);
        algo::optimizeTets(m, oo);
        algo::elevateToQuadratic(*model, m);
        const std::uint64_t h = hashMesh(m);
        if (first) { ref = h; first = false; }
        CHECK(h == ref);
    }
}

// ---------------------------------------------------------------------------
static void testBoundaryAuditKnownGap() {
    TEST("boundary: cylinder facet audit via MeshAudit (PM-BND-001 pin)");
    // The M5 MeshAudit is now the single source of truth for boundary
    // conformity. This is a KNOWN open defect (PM-BND-001): without
    // constrained facet recovery a small number of prescribed triangles
    // are not volume-boundary faces (reviewer measured 8 of 1048 at
    // h=0.1). The pin prevents silent growth; M10 drives it to zero and
    // sets the budget to 0.
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.1);
    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    mesh::AuditOptions opt;
    opt.boundaryOrphanBudget = 16;       // PM-BND-001 pin, both sides:
    opt.boundaryUnexpectedBudget = 16;   // orphans AND unexpected facets
    opt.expectedComponents = 1;          // a single solid must stay one piece
    const auto a = mesh::auditMesh(m, opt);
    CHECK(a.boundaryPrescribed > 500);
    CHECK(a.boundaryMatched >= a.boundaryPrescribed - 16);   // gap pin
    CHECK(a.boundaryOrphans <= 16);          // must become == 0 at M10
    CHECK(a.boundaryUnexpected <= 16);
    CHECK(a.passed(diag::Severity::Warning));    // within declared budget
    mesh::AuditOptions zero;                     // and NOT with budget 0
    const auto z = mesh::auditMesh(m, zero);
    CHECK(z.boundaryOrphans == 0 ||
          !z.passed(diag::Severity::Warning));
}




// ---------------------------------------------------------------------------
static void testMeshAuditClean() {
    TEST("audit: clean cylinder pipeline passes with the declared gap budget");
    auto kernel = geo::makeBuiltinKernel();
    auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.15);
    mesh::Mesh m;
    algo::meshModel(*model, spec, m);
    mesh::AuditOptions opt;
    opt.boundaryOrphanBudget = 16;       // PM-BND-001 pin; 0 after M10
    opt.boundaryUnexpectedBudget = 16;   // separate policy, same pin
    opt.expectedComponents = 1;
    const auto a = mesh::auditMesh(m, opt);
    CHECK(a.passed(diag::Severity::Warning));
    CHECK(a.nInverted == 0 && a.nNonFinite == 0);
    CHECK(a.nonManifoldFaces == 0);
    CHECK(a.connectedComponents == 1);
    CHECK(a.unusedNodes == 0);
    CHECK(a.boundaryPrescribed > 100);
    CHECK(a.boundaryMatched + a.boundaryOrphans == a.boundaryPrescribed);
    CHECK(a.boundaryOrphans <= 16 && a.boundaryUnexpected <= 16);
    CHECK(a.meshHash != 0);
    // Region volume expectation: audit must AGREE with Quality and FAIL a
    // wrong expectation.
    const double vol = mesh::analyzeTets(m).totalVolume;
    double sum = 0;
    for (const auto& [tag, v] : a.regionVolumes) sum += v;
    CHECK(std::abs(sum - vol) < 1e-9 * vol);
    mesh::AuditOptions bad = opt;
    bad.expectedRegionVolumes[a.regionVolumes.begin()->first] = vol * 2;
    const auto abad = mesh::auditMesh(m, bad);
    CHECK(!abad.passed(diag::Severity::RecoverableFailure));
    CHECK(abad.report.count(diag::Severity::InvalidOutput) >= 1);
}

// ---------------------------------------------------------------------------
static void testMeshAuditDetects() {
    TEST("audit: detects NaN, inverted, non-manifold, unused, components");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // (a) NaN coordinates -> InvalidOutput
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({nan, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        const auto r = mesh::auditMesh(m);
        CHECK(!r.passed(diag::Severity::RecoverableFailure));
        CHECK(r.report.count(diag::Severity::InvalidOutput) >= 1);
    }
    // (b) inverted element -> InvalidOutput
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a2, a1, a3});   // negative vol
        const auto r = mesh::auditMesh(m);
        CHECK(r.nInverted == 1);
        CHECK(!r.passed());
    }
    // (c) non-manifold: three tets share one face -> InvalidOutput
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0});
        const Index d1 = m.addNode({0.3, 0.3, 1}), d2 = m.addNode({0.3, 0.3, -1}),
                    d3 = m.addNode({0.4, 0.4, 2});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, d1});
        m.addCell(mesh::CellType::Tet4, {a0, a2, a1, d2});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, d3});
        const auto r = mesh::auditMesh(m);
        CHECK(r.nonManifoldFaces >= 1);
        CHECK(!r.passed());
    }
    // (d) unused node -> Warning by default, failure when configured
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addNode({9, 9, 9});   // orphan
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        const auto r = mesh::auditMesh(m);
        CHECK(r.unusedNodes == 1);
        CHECK(r.passed(diag::Severity::Warning));
        mesh::AuditOptions strict;
        strict.unusedNodesAreErrors = true;
        CHECK(!mesh::auditMesh(m, strict).passed(diag::Severity::Warning));
    }
    // (e) two disjoint tets -> two components (Info only)
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        const Index b0 = m.addNode({5, 0, 0}), b1 = m.addNode({6, 0, 0}),
                    b2 = m.addNode({5, 1, 0}), b3 = m.addNode({5, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        m.addCell(mesh::CellType::Tet4, {b0, b1, b2, b3});
        const auto r = mesh::auditMesh(m);
        CHECK(r.connectedComponents == 2);
        CHECK(r.passed(diag::Severity::Warning));
    }
}

// ---------------------------------------------------------------------------
static void testMeshAuditTet10AndHardening() {
    TEST("audit: Tet10 corner topology + hardening battery");
    // (1) Tet10-only pipeline mesh: the review found it invisible to every
    // topology audit (1 Tet10, 0 components, empty volumes). Corner-based
    // topology must audit it identically to its Tet4 shadow, with Tri6
    // corner triangles as the prescribed boundary.
    {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeBox({0, 0, 0}, {1, 1, 1});
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.35);
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        mesh::AuditOptions o4;
        o4.expectedComponents = 1;
        const auto before = mesh::auditMesh(m, o4);
        CHECK(before.passed(diag::Severity::Warning));
        algo::elevateToQuadratic(*model, m);
        const auto a = mesh::auditMesh(m, o4);
        CHECK(a.nTet4 == 0 && a.nTet10 > 0 && a.nTri6 > 0);
        CHECK(a.connectedComponents == 1);
        CHECK(a.regionVolumes.size() == 1);
        double sum = 0;
        for (const auto& [tag, v] : a.regionVolumes) sum += v;
        CHECK(std::abs(sum - 1.0) < 1e-9);          // exact straight box
        CHECK(a.boundaryPrescribed > 0);
        CHECK(a.boundaryOrphans == 0 && a.boundaryUnexpected == 0);
        CHECK(a.nonManifoldFaces == 0);
        CHECK(a.passed(diag::Severity::Warning));
        // Region-volume expectations must also work on Tet10.
        mesh::AuditOptions bad = o4;
        bad.expectedRegionVolumes[a.regionVolumes.begin()->first] = 2.0;
        CHECK(!mesh::auditMesh(m, bad).passed(
            diag::Severity::RecoverableFailure));
    }
    // (2) Mixed Tet4/Tet10: both populations in one corner topology.
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3},
                  mesh::Classification{3, 7});
        // Tet10 stacked on the shared face (a1,a2,a3) with straight
        // midsides -> mixed complex is face-connected.
        const Index b3 = m.addNode({1, 1, 1});
        auto mid = [&](Index x, Index y) {
            return m.addNode((m.node(x) + m.node(y)) * 0.5);
        };
        m.addCell(mesh::CellType::Tet10,
                  {a1, a2, a3, b3, mid(a1, a2), mid(a2, a3), mid(a3, a1),
                   mid(a1, b3), mid(a2, b3), mid(a3, b3)},
                  mesh::Classification{3, 7});
        const auto a = mesh::auditMesh(m);
        CHECK(a.nTet4 == 1 && a.nTet10 == 1);
        CHECK(a.connectedComponents == 1);          // via the shared face
        CHECK(a.regionVolumes.count(7) == 1);
        CHECK(a.nInverted == 0 && a.nInvalidTet10 == 0);
    }
    // (3) NaN node in a CELL-LESS mesh must still be detected (previously
    // only EmptyMesh was reported).
    {
        mesh::Mesh m;
        m.addNode({std::numeric_limits<double>::quiet_NaN(), 0, 0});
        const auto a = mesh::auditMesh(m);
        CHECK(!a.passed(diag::Severity::RecoverableFailure));
        CHECK(a.report.count(diag::Severity::InvalidOutput) >= 1);
    }
    // (4) Duplicate prescribed facets are counted and diagnosed, not
    // silently collapsed.
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        m.addCell(mesh::CellType::Tri3, {a0, a1, a2});
        m.addCell(mesh::CellType::Tri3, {a1, a0, a2});   // same triangle
        const auto a = mesh::auditMesh(m);
        CHECK(a.boundaryDuplicates == 1);
        bool sawDup = false;
        for (const auto& d : a.report.items)
            sawDup = sawDup ||
                     d.code == diag::Code::DuplicateBoundaryFacet;
        CHECK(sawDup);
    }
    // (5) Expected component count: a "solid" split in two must FAIL when
    // one component is expected, and pass when two are declared.
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        const Index b0 = m.addNode({5, 0, 0}), b1 = m.addNode({6, 0, 0}),
                    b2 = m.addNode({5, 1, 0}), b3 = m.addNode({5, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        m.addCell(mesh::CellType::Tet4, {b0, b1, b2, b3});
        mesh::AuditOptions one;
        one.expectedComponents = 1;
        CHECK(!mesh::auditMesh(m, one).passed(diag::Severity::Warning));
        mesh::AuditOptions two;
        two.expectedComponents = 2;
        CHECK(mesh::auditMesh(m, two).passed(diag::Severity::Warning));
    }
    // (6) Per-cell attribution: an inverted tet's diagnostic names the cell.
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a2, a1, a3});
        const auto a = mesh::auditMesh(m);
        bool attributed = false;
        for (const auto& d : a.report.items)
            attributed = attributed ||
                         (d.code == diag::Code::InvertedElement &&
                          d.cellIndex == 0 && d.hasLocation);
        CHECK(attributed);
    }
    // (7b) Hex8/Wedge6/Quad4 cells are DIAGNOSED as audit-deferred, never
    // silently clean and never audited as truncated tets (review: a
    // one-Hex8 mesh previously passed with 0 components, empty volumes).
    {
        mesh::Mesh m;
        Index id[8];
        const Vec3 P[8] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                           {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (int i = 0; i < 8; ++i) id[i] = m.addNode(P[i], {3, 1});
        m.addCell(mesh::CellType::Hex8,
                  {id[0], id[1], id[2], id[3], id[4], id[5], id[6], id[7]},
                  {3, 1});
        const auto a = mesh::auditMesh(m);
        CHECK(a.deferredAuditCells == 1);
        CHECK(!a.report.items.empty());
        bool seen = false;
        for (const auto& dd : a.report.items)
            seen = seen || dd.code == diag::Code::AuditCoverageDeferred;
        CHECK(seen);
        CHECK(a.connectedComponents == 0);   // honestly: not computed
        // Requesting a component expectation on such a mesh fails loudly.
        mesh::AuditOptions one;
        one.expectedComponents = 1;
        CHECK(!mesh::auditMesh(m, one).passed(diag::Severity::Warning));
    }
    // (7c) Degenerate zero-volume Tet4 is InvalidOutput (Jacobian/volume
    // contract: non-positive is invalid, not merely negative).
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({1, 1, 0});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});   // coplanar
        const auto a = mesh::auditMesh(m);
        CHECK(!a.passed());
    }
    // (7) Malformed connectivity: currently unreachable through the public
    // API - Mesh::addCell is the enforcement point and must stay one. The
    // audit-side validation is defense-in-depth for future
    // deserialization/compaction paths.
    {
        mesh::Mesh m;
        m.addNode({0, 0, 0});
        bool threw = false;
        try { m.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}); }
        catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw);
    }
}

// ---------------------------------------------------------------------------
static void testMeshAuditBoundaryModes() {
    TEST("audit: boundary comparison modes (absent populations cannot pass)");
    auto tetOnly = [] {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        return m;
    };
    auto triOnly = [] {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0});
        m.addCell(mesh::CellType::Tri3, {a0, a1, a2});
        return m;
    };
    // (1) Tet-only: default mode is not engaged (no surface prescribed) and
    // passes; Required treats every hull face as unexpected and fails.
    {
        auto m = tetOnly();
        CHECK(mesh::auditMesh(m).passed(diag::Severity::Warning));
        mesh::AuditOptions req;
        req.boundaryMode = mesh::BoundaryComparisonMode::Required;
        const auto a = mesh::auditMesh(m, req);
        CHECK(a.boundaryUnexpected == 4);
        CHECK(!a.passed(diag::Severity::Warning));
    }
    // (2) Tri-only: the surface IS present, so the default mode engages,
    // and with no volume cells every prescribed facet is orphaned
    // (previously passed silently - review finding).
    {
        auto m = triOnly();
        const auto a = mesh::auditMesh(m);
        CHECK(a.boundaryOrphans == 1);
        CHECK(!a.passed(diag::Severity::Warning));
    }
    // (3) Empty mesh with boundary comparison Required must fail.
    {
        mesh::Mesh m;
        mesh::AuditOptions req;
        req.boundaryMode = mesh::BoundaryComparisonMode::Required;
        CHECK(!mesh::auditMesh(m, req).passed(diag::Severity::Warning));
    }
    // (4) Disabled: neither population raises boundary diagnostics.
    {
        mesh::AuditOptions off;
        off.boundaryMode = mesh::BoundaryComparisonMode::Disabled;
        auto mt = tetOnly();
        auto a1 = mesh::auditMesh(mt, off);
        CHECK(a1.boundaryUnexpected == 0 && a1.boundaryOrphans == 0);
        CHECK(a1.passed(diag::Severity::Warning));
        auto ms = triOnly();
        auto a2 = mesh::auditMesh(ms, off);
        CHECK(a2.boundaryOrphans == 0);
        CHECK(a2.passed(diag::Severity::Warning));
    }
    // (5) expectedComponents against a volume-less mesh: requesting one
    // component with zero found must fail.
    {
        auto m = triOnly();
        mesh::AuditOptions one;
        one.expectedComponents = 1;
        CHECK(!mesh::auditMesh(m, one).passed(diag::Severity::Warning));
    }
    // (6) duplicateFacetsAreErrors escalates duplicates for export gates.
    {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({1, 0, 0}),
                    a2 = m.addNode({0, 1, 0}), a3 = m.addNode({0, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        m.addCell(mesh::CellType::Tri3, {a0, a1, a2});
        m.addCell(mesh::CellType::Tri3, {a1, a0, a2});
        mesh::AuditOptions strict;
        strict.duplicateFacetsAreErrors = true;
        strict.boundaryOrphanBudget = 4;        // isolate the duplicate check
        strict.boundaryUnexpectedBudget = 4;
        CHECK(!mesh::auditMesh(m, strict).passed(diag::Severity::Warning));
    }
}

// ---------------------------------------------------------------------------
static void testPoolStragglerRace() {
    TEST("parallel: rapid consecutive pool jobs (straggler-capture regression)");
    // Regression for PM-PAR-001: a worker could capture the job pointer,
    // the job could complete and run() return (destroying the caller-owned
    // function), and the NEXT job's chunk-counter reset would hand the
    // straggler fresh indices to execute against the dangling previous
    // function - crash (std::bad_function_call) or silent chunk theft.
    //
    // REVIEW-CORRECTED: the first version of this test used 64 elements,
    // BELOW the pool's serial threshold, so it never touched the pool at
    // all. It now uses 4096 elements per job AND verifies that multiple
    // OS threads actually executed body chunks, so a future threshold
    // change cannot silently re-serialize it.
    core::setThreadCount(4);
    constexpr std::size_t kN = 4096;
    std::vector<int> data(kN);
    std::mutex idMx;
    std::set<std::thread::id> seenThreads;
    for (int round = 0; round < 1500; ++round) {
        const int stamp = round + 1;
        std::fill(data.begin(), data.end(), 0);
        core::parallelFor(0, kN, [&](std::size_t i) {
            data[i] += stamp;                  // body differs per job
            if ((i & 1023) == 0) {             // sample thread ids cheaply
                std::lock_guard<std::mutex> lk(idMx);
                seenThreads.insert(std::this_thread::get_id());
            }
        });
        for (const int v : data)
            if (v != stamp) { CHECK(v == stamp); return; }
    }
    // Multiple workers must have participated, or this test is not
    // exercising the pool (the original review finding).
    CHECK(seenThreads.size() >= 2);
    core::setThreadCount(0);
}


// ---------------------------------------------------------------------------
static void testFrozenBaselineHashes() {
    TEST("baseline: ALL NINE frozen cases machine-verified (M4-B.10 oracle)");
    // Frozen 2026-07-04. STRICT hashes are compared only on the reference
    // environment (gcc/linux/x86_64/optimized) per docs/BASELINE_FREEZE.md;
    // SEMANTIC v2 hashes (canonical physical policy) are compared on EVERY
    // supported environment - a semantic mismatch across compilers is a
    // finding, not noise. Review finding closed here: previously only 4 of
    // 9 baselines were automated and strict values ran unguarded off-
    // platform.
    const bool strictEnv = core::isStrictBaselineReferenceEnv();
    if (!strictEnv)
        std::puts("  (strict hashes skipped: not the reference environment)");
    auto verify = [&](mesh::Mesh& m, std::uint64_t strict,
                      std::uint64_t semantic) {
        if (strictEnv) CHECK(mesh::meshHash(m) == strict);
        // Review hardening: qualification FAILS on a non-canonical
        // relabeling - an ambiguous hash must never pass as a baseline
        // match (coincident-node cases could otherwise reintroduce
        // input-order dependence behind an apparently green comparison).
        const auto sh = mesh::semanticHashEx(m);
        CHECK(sh.canonical);
        CHECK(sh.value == semantic);
    };
    auto pipeline = [&](bool cylinder, bool elevate, double h, Vec3 lo,
                        Vec3 hi, int threads, std::uint64_t strict,
                        std::uint64_t semantic) {
        auto kernel = geo::makeBuiltinKernel();
        auto model = cylinder ? kernel->makeCylinder({0, 0, 0}, 0.5, 1.0)
                              : kernel->makeBox(lo, hi);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(h);
        spec.threads = threads;
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        if (elevate) algo::elevateToQuadratic(*model, m);
        verify(m, strict, semantic);
        return mesh::meshHash(m);
    };
    // 1-2: unit box; cylinder
    pipeline(false, false, 0.3, {0,0,0}, {1,1,1}, 0,
             0xfc6d2d242e745734ULL, 0xecc7fda58ef68700ULL);
    const auto cylStrict =
        pipeline(true, false, 0.15, {}, {}, 0,
                 0x3c58de9d246d2613ULL, 0x00c9d7259371b796ULL);
    // 3: thin box (PM-BND-001 22/22 boundary-gap case)
    pipeline(false, false, 0.2, {0,0,0}, {1,1,1e-3}, 0,
             0x407bcb430034ee78ULL, 0x88d4c774e346f80dULL);
    // 4-5: translated 1e9; scaled 1e6 (v2 must distinguish both from unit)
    pipeline(false, false, 0.3, {1e9,1e9,1e9}, {1e9+1,1e9+1,1e9+1}, 0,
             0x1d97492ad2ae2ae8ULL, 0xb2853015c6c3a2d7ULL);
    pipeline(false, false, 3e5, {0,0,0}, {1e6,1e6,1e6}, 0,
             0x68e41df16b5505d8ULL, 0x496eeb0c3c69bf2bULL);
    // 6: Tet10 elevated box
    pipeline(false, true, 0.35, {0,0,0}, {1,1,1}, 0,
             0xd809f64a9ba92494ULL, 0x001be401a71fad76ULL);
    // 7: threads=4 cylinder - strict hash must equal the default-thread
    // run IN THIS BUILD (environment-independent determinism statement).
    const auto cyl4Strict =
        pipeline(true, false, 0.15, {}, {}, 4,
                 0x3c58de9d246d2613ULL, 0x00c9d7259371b796ULL);
    CHECK(cyl4Strict == cylStrict);
    // 8-9: construction-level lattice + flattened cloud
    auto construct = [&](bool flat, std::uint64_t strict,
                         std::uint64_t semantic) {
        std::vector<Vec3> pts;
        if (flat) {
            std::mt19937 rng(778);
            for (int i = 0; i < 200; ++i)
                pts.push_back(
                    {canonU01(rng), canonU01(rng), canonU01(rng) * 1e-9});
        } else {
            for (int k = 0; k < 5; ++k)
                for (int jj = 0; jj < 5; ++jj)
                    for (int ii = 0; ii < 5; ++ii)
                        pts.push_back({double(ii), double(jj), double(k)});
        }
        auto tets = algo::delaunay3D(pts);
        mesh::Mesh m;
        m.setDedupTolerance(flat ? 1e-15 : 1e-12);
        for (auto& p : pts) m.addNode(p, {3, 1});
        for (auto& tt : tets)
            m.addCell(mesh::CellType::Tet4,
                      {tt.v[0], tt.v[1], tt.v[2], tt.v[3]}, {3, 1});
        verify(m, strict, semantic);
    };
    construct(false, 0xffc04ee273015e44ULL, 0x050061fc47f0350bULL);
    construct(true, 0x361514c645e05947ULL, 0x23f5fb22d60f7ef2ULL);
}


// ---------------------------------------------------------------------------
static void testSemanticHashCanonicality() {
    TEST("hash: semantic v3 canonical physical (relabel/translation/scale)");
    // Review repro (v1 era): identical physical mesh, different node
    // numbering -> v1 hashes DIFFERED. The canonical hash must be equal while
    // the strict hash (by contract) differs.
    auto build = [](const std::vector<int>& order) {
        mesh::Mesh m;
        const Vec3 P[5] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
        Index id[5];
        for (int i : order) id[i] = m.addNode(P[i], {3, 1});
        m.addCell(mesh::CellType::Tet4, {id[0], id[1], id[2], id[3]}, {3, 1});
        m.addCell(mesh::CellType::Tet4, {id[1], id[2], id[3], id[4]}, {3, 1});
        return m;
    };
    auto a = build({0, 1, 2, 3, 4});
    auto b = build({4, 2, 0, 3, 1});
    CHECK(mesh::semanticHash(a) == mesh::semanticHash(b));
    CHECK(mesh::meshHash(a) != mesh::meshHash(b));
    // Local connectivity rotation (same orientation class) is canonical.
    {
        mesh::Mesh c;
        Index id[4];
        const Vec3 P[4] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        for (int i = 0; i < 4; ++i) id[i] = c.addNode(P[i], {3, 1});
        c.addCell(mesh::CellType::Tet4, {id[1], id[2], id[0], id[3]}, {3, 1});
        mesh::Mesh d;
        for (int i = 0; i < 4; ++i) id[i] = d.addNode(P[i], {3, 1});
        d.addCell(mesh::CellType::Tet4, {id[0], id[1], id[2], id[3]}, {3, 1});
        CHECK(mesh::semanticHash(c) == mesh::semanticHash(d));
        // Orientation FLIP is semantic (inverted element != valid element).
        mesh::Mesh e;
        for (int i = 0; i < 4; ++i) id[i] = e.addNode(P[i], {3, 1});
        e.addCell(mesh::CellType::Tet4, {id[0], id[2], id[1], id[3]}, {3, 1});
        CHECK(mesh::semanticHash(e) != mesh::semanticHash(d));
    }
    // Review repro 2: translation near 1e9 must not be hidden (v1 clamp
    // overflow gave two DIFFERENT translated tets the SAME hash).
    auto tet = [](double off) {
        mesh::Mesh m;
        const Index a0 = m.addNode({off, 0, 0}), a1 = m.addNode({off + 1, 0, 0}),
                    a2 = m.addNode({off, 1, 0}), a3 = m.addNode({off, 0, 1});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        return m;
    };
    CHECK(mesh::semanticHash(tet(1e9)) != mesh::semanticHash(tet(1e9 + 1.0)));
    CHECK(mesh::semanticHash(tet(0)) != mesh::semanticHash(tet(1e9)));
    // Review repro 3: scale is semantic (v1 was scale-invariant: unit box
    // and 1e6 box collided).
    auto scaled = [](double s) {
        mesh::Mesh m;
        const Index a0 = m.addNode({0, 0, 0}), a1 = m.addNode({s, 0, 0}),
                    a2 = m.addNode({0, s, 0}), a3 = m.addNode({0, 0, s});
        m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        return m;
    };
    CHECK(mesh::semanticHash(scaled(1.0)) != mesh::semanticHash(scaled(1e6)));
    // FP-noise absorption: sub-2^-40 relative perturbation hashes equal
    // for bucket-interior values. HONEST LIMIT (documented in the
    // header): quantization absorbs noise only away from rounding
    // boundaries; a value straddling a bucket edge can flip on an
    // epsilon change (this exact test at 1e9 found such a boundary).
    // The engineered case below is provably interior: 1.0 has mantissa
    // exactly 0.5, so a 1e-14 relative nudge moves the 2^40-scaled
    // mantissa by ~5e-6 of a unit - nowhere near the 0.5 boundary.
    {
        auto unitTet = [](double f) {
            mesh::Mesh m;
            const Index a0 = m.addNode({0, 0, 0}),
                        a1 = m.addNode({1.0 * f, 0, 0}),
                        a2 = m.addNode({0, 1.0 * f, 0}),
                        a3 = m.addNode({0, 0, 1.0 * f});
            m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
            return m;
        };
        auto u0 = unitTet(1.0);
        auto u1 = unitTet(1.0 + 1e-14);
        CHECK(mesh::semanticHash(u0) == mesh::semanticHash(u1));
    }
    // v3 review repro 1: nodes 1e-13 apart share a quantization bucket;
    // the relabeling is then AMBIGUOUS and must be reported, not silently
    // input-order dependent.
    {
        auto near2 = [](bool swapOrder) {
            mesh::Mesh m;
            const Vec3 A{0.5, 0.5, 0.5}, B{0.5 + 1e-13, 0.5, 0.5};
            Index ia, ib;   // (were 'a','b': C4456/-Wshadow of outer locals)
            if (swapOrder) { ib = m.addNode(B); ia = m.addNode(A); }
            else           { ia = m.addNode(A); ib = m.addNode(B); }
            const Index c0 = m.addNode({1, 0, 0}), c1 = m.addNode({0, 1, 0});
            m.addCell(mesh::CellType::Tet4, {ia, c0, c1, ib});
            return m;
        };
        auto m1 = near2(false);
        CHECK(!mesh::semanticHashEx(m1).canonical);
        // And clearly separated nodes ARE canonical.
        auto m2 = mesh::Mesh();
        const Index a0 = m2.addNode({0,0,0}), a1 = m2.addNode({1,0,0}),
                    a2 = m2.addNode({0,1,0}), a3 = m2.addNode({0,0,1});
        m2.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
        CHECK(mesh::semanticHashEx(m2).canonical);
    }
    // v3 review repro 2: a 0.1 edge change at 1e12 must CHANGE the hash
    // (v2's absolute-relative encoding hid it).
    {
        auto farTet = [](double dx) {
            mesh::Mesh m;
            const double B = 1e12;
            const Index a0 = m.addNode({B, 0, 0}),
                        a1 = m.addNode({B + 1.0 + dx, 0, 0}),
                        a2 = m.addNode({B, 1, 0}), a3 = m.addNode({B, 0, 1});
            m.addCell(mesh::CellType::Tet4, {a0, a1, a2, a3});
            return m;
        };
        auto f0 = farTet(0.0);
        auto f1 = farTet(0.1);
        CHECK(mesh::semanticHash(f0) != mesh::semanticHash(f1));
        CHECK(mesh::semanticHashEx(f0).canonical);
    }
    // v3 review repro 3: Quad4 - a cyclic rotation is the SAME cell; a
    // twisted connectivity is a DIFFERENT cell (v2 had both backwards).
    {
        auto quad = [](std::array<int, 4> ord) {
            mesh::Mesh m;
            const Vec3 P[4] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0}};
            Index id[4];
            for (int i = 0; i < 4; ++i) id[i] = m.addNode(P[i], {2, 1});
            m.addCell(mesh::CellType::Quad4,
                      {id[ord[0]], id[ord[1]], id[ord[2]], id[ord[3]]},
                      {2, 1});
            return m;
        };
        auto valid = quad({0, 1, 2, 3});
        auto rotated = quad({1, 2, 3, 0});
        auto twisted = quad({0, 1, 3, 2});
        CHECK(mesh::semanticHash(valid) == mesh::semanticHash(rotated));
        CHECK(mesh::semanticHash(valid) != mesh::semanticHash(twisted));
    }
    // Hex8 rotation invariance through the generated 24-element group:
    // rotate the cube 90 degrees about z in CONNECTIVITY (same physical
    // element) -> same hash.
    {
        auto hex = [](std::array<int, 8> ord) {
            mesh::Mesh m;
            const Vec3 P[8] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                               {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
            Index id[8];
            for (int i = 0; i < 8; ++i) id[i] = m.addNode(P[i], {3, 1});
            std::vector<Index> conn;
            for (int i = 0; i < 8; ++i) conn.push_back(id[ord[i]]);
            m.addCell(mesh::CellType::Hex8, conn, {3, 1});
            return m;
        };
        auto h0 = hex({0, 1, 2, 3, 4, 5, 6, 7});
        auto h1 = hex({1, 2, 3, 0, 5, 6, 7, 4});
        CHECK(mesh::semanticHash(h0) == mesh::semanticHash(h1));
    }
    // Tet10 canonical form: midsides re-keyed under corner permutation.
    {
        auto t10 = [](const std::array<int, 4>& perm) {
            mesh::Mesh m;
            const Vec3 C[4] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
            Index cid[4];
            for (int i = 0; i < 4; ++i) cid[i] = m.addNode(C[i], {3, 1});
            // even permutations of corners with correctly re-derived
            // midsides describe the SAME physical element
            const int e[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
            std::array<Index, 10> conn;
            for (int i = 0; i < 4; ++i) conn[size_t(i)] = cid[perm[size_t(i)]];
            for (int q = 0; q < 6; ++q) {
                const Vec3 mid = (C[perm[size_t(e[q][0])]] +
                                  C[perm[size_t(e[q][1])]]) * 0.5;
                conn[size_t(4 + q)] = m.addNode(mid, {3, 1});
            }
            m.addCell(mesh::CellType::Tet10,
                      {conn[0], conn[1], conn[2], conn[3], conn[4], conn[5],
                       conn[6], conn[7], conn[8], conn[9]}, {3, 1});
            return m;
        };
        CHECK(mesh::semanticHash(t10({0, 1, 2, 3})) ==
              mesh::semanticHash(t10({1, 2, 0, 3})));   // even perm
    }
}

// ---------------------------------------------------------------------------
static void testPoolConcurrentRejection() {
    TEST("parallel: concurrent dispatch rejected WITHOUT poisoning caller");
    core::setThreadCount(4);
    constexpr std::size_t kN = 4096;
    std::atomic<bool> inside{false}, release{false};
    std::atomic<int> rejected{0};
    std::vector<int> data(kN, 0);
    std::thread holder([&] {
        core::parallelFor(0, kN, [&](std::size_t i) {
            if (i == 0) {
                inside.store(true);
                while (!release.load()) std::this_thread::yield();
            }
            data[i] = 1;
        });
    });
    while (!inside.load()) std::this_thread::yield();
    try {
        std::vector<int> other(kN, 0);
        core::parallelFor(0, kN, [&](std::size_t i) { other[i] = 1; });
    } catch (const std::logic_error&) { rejected.store(1); }
    release.store(true);
    holder.join();
    CHECK(rejected.load() == 1);
    std::size_t done = 0;
    for (const int v : data) done += static_cast<std::size_t>(v);
    CHECK(done == kN);                 // the holder's job was unharmed
    // Review finding: the rejection path previously leaked tlInPool=true,
    // permanently serializing THIS thread. Prove the guard was released:
    // the next dispatch from this thread must be genuinely parallel.
    {
        // Deterministic proof that the concurrent-dispatch rejection
        // restored THIS thread's re-entrancy guard (PM-PAR-001). A poisoned
        // guard leaves tlInPool == true on the caller, which we can now
        // observe directly. Distinguishes real release from serial
        // fallback, which the thread-count heuristic could not.
        CHECK(!core::inPoolContext());          // not in a pool here
        std::atomic<bool> guardInsideBody{false};
        std::atomic<bool> guardAfter{true};
        core::parallelFor(0, kN, [&](std::size_t i) {
            if (i == 0) guardInsideBody.store(core::inPoolContext());
        });
        guardAfter.store(core::inPoolContext());
        CHECK(guardInsideBody.load());          // guard set during dispatch
        CHECK(!guardAfter.load());               // and cleared afterward
    }
    core::setThreadCount(0);
}

// ---------------------------------------------------------------------------
static void testPoolNestedDispatch() {
    TEST("parallel: nested dispatch from a REAL worker degrades serial");
    // Independent of the rejection test (review: the poisoned guard made
    // the previous combined test run its outer loop serially, proving
    // nothing about workers). Here we additionally verify the outer
    // dispatch was genuinely parallel via participating thread ids.
    core::setThreadCount(4);
    constexpr std::size_t kOuter = 4096, kInner = 2048;   // both above threshold
    std::vector<long> row(kOuter, 0);
    std::mutex mx;
    std::set<std::thread::id> ids;
    core::parallelFor(0, kOuter, [&](std::size_t r) {
        if ((r & 255) == 0) {
            std::lock_guard<std::mutex> lk(mx);
            ids.insert(std::this_thread::get_id());
        }
        long acc = 0;
        core::parallelFor(0, kInner, [&](std::size_t i) {
            acc += static_cast<long>(i % 7);   // serial fallback: safe
        });
        row[r] = acc;
    });
    CHECK(ids.size() >= 2);            // outer WAS parallel: workers nested
    const long expect = [&] {
        long a = 0;
        for (std::size_t i = 0; i < kInner; ++i) a += long(i % 7);
        return a;
    }();
    std::size_t wrong = 0;
    for (const long v : row) wrong += v != expect ? 1 : 0;
    CHECK(wrong == 0);
    core::setThreadCount(0);
}


// ---------------------------------------------------------------------------
static void testElementTraits() {
    TEST("topo: M4-B.0 element traits golden reference (all 10 types)");
    using topo::ElementType;
    using topo::FaceShape;
    const auto V = [](const topo::RefCoord& r) {
        return Vec3{r.x, r.y, r.z};
    };
    int vtkSeen[256] = {};
    for (int ti = 0; ti < static_cast<int>(ElementType::Count_); ++ti) {
        const auto& tr = topo::traits(static_cast<ElementType>(ti));
        // Table sanity.
        CHECK(tr.type == static_cast<ElementType>(ti));
        CHECK(tr.cornerCount >= 2 && tr.cornerCount <= 8);
        CHECK(tr.nodeCount >= tr.cornerCount);
        CHECK(tr.dimension >= 1 && tr.dimension <= 3);
        CHECK(tr.vtkId != 0);
        CHECK(vtkSeen[tr.vtkId]++ == 0);           // export ids unique
        // Edges reference distinct corner-local indices.
        for (int e = 0; e < tr.edgeCount; ++e) {
            CHECK(tr.edges[e].a < tr.cornerCount);
            CHECK(tr.edges[e].b < tr.cornerCount);
            CHECK(tr.edges[e].a != tr.edges[e].b);
        }
        // Reference coordinates: one per node, corners distinct.
        for (int a = 0; a < tr.cornerCount; ++a)
            for (int b = a + 1; b < tr.cornerCount; ++b)
                CHECK((V(tr.refCoords[a]) - V(tr.refCoords[b])).norm() >
                      1e-12);
        // High-order convention: node cornerCount+q is the midpoint of
        // edges[q] in reference coordinates.
        for (int q = 0; q + tr.cornerCount < tr.nodeCount; ++q) {
            const Vec3 mid = (V(tr.refCoords[tr.edges[q].a]) +
                              V(tr.refCoords[tr.edges[q].b])) * 0.5;
            CHECK((V(tr.refCoords[tr.cornerCount + q]) - mid).norm() <
                  1e-12);
        }
        if (tr.dimension == 3) {
            CHECK(tr.faceCount > 0);
            // Cell centroid of corners.
            Vec3 cc{0, 0, 0};
            for (int a = 0; a < tr.cornerCount; ++a) cc = cc + V(tr.refCoords[a]);
            cc = cc * (1.0 / tr.cornerCount);
            int triF = 0, quadF = 0;
            for (int f = 0; f < tr.faceCount; ++f) {
                const auto& fc = tr.faces[f];
                CHECK((fc.shape == FaceShape::Tri && fc.n == 3) ||
                      (fc.shape == FaceShape::Quad && fc.n == 4));
                (fc.shape == FaceShape::Tri ? triF : quadF)++;
                Vec3 fcent{0, 0, 0};
                for (int k2 = 0; k2 < fc.n; ++k2) {
                    CHECK(fc.v[size_t(k2)] < tr.cornerCount);
                    fcent = fcent + V(tr.refCoords[fc.v[size_t(k2)]]);
                }
                fcent = fcent * (1.0 / fc.n);
                // Outward-orientation golden check: right-hand normal of
                // the first three corners points away from the centroid.
                const Vec3 n =
                    (V(tr.refCoords[fc.v[1]]) - V(tr.refCoords[fc.v[0]]))
                        .cross(V(tr.refCoords[fc.v[2]]) -
                               V(tr.refCoords[fc.v[0]]));
                CHECK(n.dot(fcent - cc) > 1e-12);
            }
            if (tr.type == ElementType::Tet4 || tr.type == ElementType::Tet10)
                CHECK(triF == 4 && quadF == 0);
            if (tr.type == ElementType::Hex8) CHECK(triF == 0 && quadF == 6);
            if (tr.type == ElementType::Wedge6) CHECK(triF == 2 && quadF == 3);
            if (tr.type == ElementType::Pyramid5) CHECK(triF == 4 && quadF == 1);
        }
    }
    // Positive-orientation convention agrees with the exact predicate.
    {
        const auto& t4 = topo::traits(ElementType::Tet4);
        CHECK(predicates::orient3d(V(t4.refCoords[0]), V(t4.refCoords[1]),
                                   V(t4.refCoords[2]), V(t4.refCoords[3])) ==
              predicates::Sign::Positive);
    }
    // Naming locks round-trip for every legacy CellType.
    {
        topo::ElementType e;
        CHECK(topo::elementTypeOf(mesh::CellType::Line2, e) &&
              e == ElementType::Segment2);
        CHECK(topo::elementTypeOf(mesh::CellType::Line3, e) &&
              e == ElementType::Segment3);
        CHECK(topo::elementTypeOf(mesh::CellType::Tri6, e) &&
              e == ElementType::Tri6);
        CHECK(topo::elementTypeOf(mesh::CellType::Tet10, e) &&
              e == ElementType::Tet10);
    }
    // Connectivity validation (std::span, table-driven over ALL TEN
    // types per the issue's own qualification criteria - review round 3):
    // correct arity accepted; short and long arity, out-of-range nodes,
    // and ANY repeated node (midsides included) rejected.
    {
        using ET = ElementType;
        for (int ti = 0; ti < static_cast<int>(ET::Count_); ++ti) {
            const auto et = static_cast<ET>(ti);
            const auto& tr = topo::traits(et);
            const Index N = Index(tr.nodeCount) + 4;   // node pool size
            std::vector<Index> good(tr.nodeCount);
            for (int q = 0; q < tr.nodeCount; ++q) good[size_t(q)] = q;
            CHECK(topo::validConnectivity(et, good, N));
            // Short and long arity.
            CHECK(!topo::validConnectivity(
                et, std::span<const Index>(good.data(), good.size() - 1),
                N));
            std::vector<Index> longer = good;
            longer.push_back(Index(tr.nodeCount));
            CHECK(!topo::validConnectivity(et, longer, N));
            // Out-of-range node (first slot).
            std::vector<Index> oob = good;
            oob[0] = N;
            CHECK(!topo::validConnectivity(et, oob, N));
            // Every possible duplicated pair, corners AND high-order.
            for (int a2 = 0; a2 < tr.nodeCount; ++a2)
                for (int b2 = a2 + 1; b2 < tr.nodeCount; ++b2) {
                    std::vector<Index> dup = good;
                    dup[size_t(b2)] = dup[size_t(a2)];
                    CHECK(!topo::validConnectivity(et, dup, N));
                }
        }
    }
    // Export mappings are the authoritative VTK/MSH ids.
    {
        using ET = ElementType;
        const std::pair<ET, std::pair<int, int>> ex[] = {
            {ET::Segment2, {3, 1}},  {ET::Segment3, {21, 8}},
            {ET::Tri3, {5, 2}},      {ET::Tri6, {22, 9}},
            {ET::Quad4, {9, 3}},     {ET::Tet4, {10, 4}},
            {ET::Tet10, {24, 11}},   {ET::Hex8, {12, 5}},
            {ET::Wedge6, {13, 6}},   {ET::Pyramid5, {14, 7}},
        };
        for (const auto& [et, ids] : ex) {
            CHECK(topo::traits(et).vtkId == ids.first);
            CHECK(topo::traits(et).mshId == ids.second);
        }
    }
    // Rotation-group sizes match theory (generated geometrically).
    {
        const std::pair<ElementType, std::size_t> sizes[] = {
            {ElementType::Segment2, 1}, {ElementType::Segment3, 1},
            {ElementType::Tri3, 3},     {ElementType::Tri6, 3},
            {ElementType::Quad4, 4},    {ElementType::Tet4, 12},
            {ElementType::Tet10, 12},   {ElementType::Hex8, 24},
            {ElementType::Wedge6, 6},   {ElementType::Pyramid5, 4},
        };
        for (const auto& [et, n] : sizes)
            CHECK(topo::rotations(et).size() == n);
    }
}

// ---------------------------------------------------------------------------
static void testWriterElementMappings() {
    TEST("io: writers emit traits ids + connectivity for ALL storable types");
    // Review round 3: the traits TABLE was tested but not what the
    // writers actually EMIT. Build one mesh containing every legacy-
    // storable type, write real files, parse them back, and verify the
    // element type ids and per-element connectivity ordering.
    mesh::Mesh m;
    std::vector<Index> id;
    // 12 well-spread nodes cover every element's needs.
    const Vec3 P[12] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
                        {2,0,0},{2,1,0},{2,0,1},{0.5,0.5,2}};
    for (const auto& p : P) id.push_back(m.addNode(p, {3, 1}));
    struct Want { mesh::CellType t; std::vector<Index> conn; };
    const std::vector<Want> cells = {
        {mesh::CellType::Line2,  {0, 1}},
        {mesh::CellType::Line3,  {0, 2, 1}},
        {mesh::CellType::Tri3,   {0, 1, 2}},
        {mesh::CellType::Tri6,   {0, 2, 4, 1, 3, 5}},
        {mesh::CellType::Quad4,  {0, 1, 2, 3}},
        {mesh::CellType::Tet4,   {0, 1, 2, 4}},
        {mesh::CellType::Tet10,  {0, 1, 2, 4, 3, 5, 6, 7, 8, 9}},
        {mesh::CellType::Prism6, {0, 1, 2, 4, 5, 6}},
        {mesh::CellType::Hex8,   {0, 1, 2, 3, 4, 5, 6, 7}},
    };
    for (const auto& w : cells)
        m.addCell(w.t, w.conn, {3, 1});
    // ---- VTK ---------------------------------------------------------------
    {
        const auto vtkPath =
            (std::filesystem::temp_directory_path() /
             ("pm_writer_" + std::to_string(::getpid()) + ".vtk"))
                .string();
        io::writeVtk(m, vtkPath, /*volumeOnly=*/false);
        std::ifstream in(vtkPath);
        CHECK(bool(in));
        std::string line;
        std::vector<std::vector<long>> conn;
        std::vector<int> types;
        while (std::getline(in, line)) {
            if (line.rfind("CELLS", 0) == 0) {
                for (std::size_t i = 0; i < cells.size(); ++i) {
                    std::getline(in, line);
                    std::istringstream ls(line);
                    long n2;
                    ls >> n2;
                    std::vector<long> c(static_cast<std::size_t>(n2), 0);
                    for (auto& v : c) ls >> v;
                    conn.push_back(c);
                }
            } else if (line.rfind("CELL_TYPES", 0) == 0) {
                for (std::size_t i = 0; i < cells.size(); ++i) {
                    std::getline(in, line);
                    types.push_back(std::stoi(line));
                }
            }
        }
        CHECK(conn.size() == cells.size() && types.size() == cells.size());
        for (std::size_t i = 0; i < cells.size(); ++i) {
            topo::ElementType et{};
            CHECK(topo::elementTypeOf(cells[i].t, et));
            CHECK(types[i] == topo::traits(et).vtkId);
            CHECK(conn[i].size() == cells[i].conn.size());
            for (std::size_t k = 0; k < conn[i].size(); ++k)
                CHECK(conn[i][k] == long(cells[i].conn[k]));
        }
    }
    // ---- MSH 4.1 ------------------------------------------------------------
    {
        const auto mshPath =
            (std::filesystem::temp_directory_path() /
             ("pm_writer_" + std::to_string(::getpid()) + ".msh"))
                .string();
        io::writeMsh(m, mshPath);
        std::ifstream in(mshPath);
        CHECK(bool(in));
        std::string line;
        std::map<int, int> typeCounts;          // MSH type ID -> elements
        std::map<int, std::vector<long>> firstConn;
        bool inElem = false;
        int blockType = -1, blockLeft = 0;
        while (std::getline(in, line)) {
            if (line == "$Elements") { inElem = true; std::getline(in, line); continue; }
            if (line == "$EndElements") inElem = false;
            if (!inElem) continue;
            std::istringstream ls(line);
            if (blockLeft == 0) {
                int dim, tag, ety, n2;
                if (ls >> dim >> tag >> ety >> n2) {
                    blockType = ety;
                    blockLeft = n2;
                }
            } else {
                long eid;
                ls >> eid;
                std::vector<long> c;
                long v;
                while (ls >> v) c.push_back(v - 1);   // MSH IDs are 1-based
                ++typeCounts[blockType];
                if (!firstConn.count(blockType)) firstConn[blockType] = c;
                --blockLeft;
            }
        }
        for (const auto& w : cells) {
            topo::ElementType et{};
            CHECK(topo::elementTypeOf(w.t, et));
            const int g = topo::traits(et).mshId;
            CHECK(typeCounts[g] == 1);
            const auto& c = firstConn[g];
            CHECK(c.size() == w.conn.size());
            for (std::size_t k = 0; k < c.size(); ++k)
                CHECK(c[k] == long(w.conn[k]));
        }
    }
}

// ---------------------------------------------------------------------------
static void testPoolExceptionPropagation() {
    TEST("parallel: body exceptions cancel, drain, rethrow on dispatcher");
    core::setThreadCount(4);
    constexpr std::size_t kN = 4096;
    // Throw from MANY indices so both the caller and workers hit throws
    // across rounds; the dispatcher must receive exactly one exception
    // per job and the process must not terminate.
    int caught = 0;
    std::vector<std::atomic<int>> ran(kN);
    for (int round = 0; round < 50; ++round) {
        for (auto& r : ran) r.store(0);
        try {
            core::parallelFor(0, kN, [&](std::size_t i) {
                if ((i % 97) == 13)
                    throw std::runtime_error("body failure");
                ran[i].fetch_add(1);
            });
        } catch (const std::runtime_error&) { ++caught; }
    }
    CHECK(caught == 50);
    // No index ever executed twice (cancellation consumes, never re-runs).
    // (ran was reset per round; the final round's data suffices.)
    for (const auto& r : ran) CHECK(r.load() <= 1);
    // The pool is healthy afterwards: a normal job completes fully and
    // in parallel.
    {
        std::mutex mx;
        std::set<std::thread::id> ids;
        std::vector<int> data(kN, 0);
        for (int round = 0; round < 300 && ids.size() < 2; ++round) {
            std::fill(data.begin(), data.end(), 0);
            core::parallelFor(0, kN, [&](std::size_t i) {
                data[i] = 1;
                if ((i & 255) == 0) {
                    std::lock_guard<std::mutex> lk(mx);
                    ids.insert(std::this_thread::get_id());
                }
            });
            std::size_t done = 0;
            for (const int v : data) done += size_t(v);
            if (done != kN) { CHECK(done == kN); break; }
        }
        CHECK(ids.size() >= 2);
    }
    core::setThreadCount(0);
}

// ---------------------------------------------------------------------------
template <class Tag>
static void handleTagBattery(const char* /*name*/) {
    topo::HandleAllocator<Tag> a;
    const auto h0 = a.alloc();
    CHECK(h0.isSpecified() && a.isCurrent(h0));
    const auto keep = h0;
    a.free(h0);
    CHECK(!a.isCurrent(h0) && !a.isCurrent(keep));
    CHECK(keep.isSpecified());
    const auto h1 = a.alloc();
    CHECK(h1.slot == h0.slot && h1.gen != h0.gen && a.isCurrent(h1));
    bool df = false;
    try { a.free(h0); } catch (const std::logic_error&) { df = true; }
    CHECK(df);
    bool sd = false;
    try { (void)a.deref(h0); } catch (const std::logic_error&) { sd = true; }
    CHECK(sd);
    CHECK(a.accountingConsistent());
}

static void testGenericHandles() {
    TEST("topo: M4-B.1 generic handles (state machine; all four tags)");
    using topo::HandleAllocator;
    using topo::NodeHandle;
    using topo::CellHandle;
    static_assert(!std::is_convertible_v<NodeHandle, CellHandle>);
    static_assert(!std::is_convertible_v<CellHandle, NodeHandle>);
    static_assert(!std::is_convertible_v<NodeHandle, std::int32_t>);
    static_assert(!std::is_convertible_v<std::int32_t, NodeHandle>);
    static_assert(
        !std::is_copy_constructible_v<HandleAllocator<topo::NodeTag>>);
    static_assert(
        std::is_move_constructible_v<HandleAllocator<topo::NodeTag>>);
    // Move ASSIGNMENT deleted (round-6 decision).
    static_assert(
        !std::is_move_assignable_v<HandleAllocator<topo::NodeTag>>);
    handleTagBattery<topo::NodeTag>("node");
    handleTagBattery<topo::EdgeTag>("edge");
    handleTagBattery<topo::FaceTag>("face");
    handleTagBattery<topo::CellTag>("cell");

    HandleAllocator<topo::NodeTag> a;
    NodeHandle def;
    CHECK(!def.isSpecified() && !a.isCurrent(def));
    const auto h0 = a.alloc();
    const auto h1 = a.alloc();
    CHECK(h0.isSpecified() && a.isCurrent(h0) && a.isCurrent(h1));
    CHECK(a.liveCount() == 2);
    const auto h0copy = h0;
    a.free(h0);
    CHECK(!a.isCurrent(h0) && !a.isCurrent(h0copy));
    const auto h2 = a.alloc();
    CHECK(h2.slot == h0.slot && h2.gen != h0.gen);
    CHECK(a.isCurrent(h2) && !a.isCurrent(h0));
    CHECK(std::hash<NodeHandle>{}(h0) == std::hash<NodeHandle>{}(h0copy));
    CHECK(std::hash<NodeHandle>{}(h0) != std::hash<NodeHandle>{}(h2));
    CHECK(h0 == h0copy && h0 != h2);
    CHECK((h0 < h2) != (h2 < h0));
    std::unordered_set<NodeHandle> set{h0, h2};
    CHECK(set.size() == 2 && set.count(h0copy) == 1);

    // Cross-pool identity + moved-from behavior.
    {
        HandleAllocator<topo::NodeTag> A, B;
        const auto ha = A.alloc();
        const auto hb = B.alloc();
        CHECK(A.isCurrent(ha) && B.isCurrent(hb));
        CHECK(!B.isCurrent(ha) && !A.isCurrent(hb) && ha != hb);
        bool foreign = false;
        try { (void)B.deref(ha); }
        catch (const std::logic_error&) { foreign = true; }
        CHECK(foreign);
        HandleAllocator<topo::NodeTag> C(std::move(A));
        CHECK(C.isCurrent(ha) && !A.isCurrent(ha));
        // Moved-from alloc() THROWS (review: it used to mint a
        // never-current handle and bump liveCount).
        bool movedAlloc = false;
        try { (void)A.alloc(); }
        catch (const std::logic_error&) { movedAlloc = true; }
        CHECK(movedAlloc);
        CHECK(A.liveCount() == 0 && A.accountingConsistent());
    }

    // Sequence-independent deterministic order under churn.
    {
        HandleAllocator<topo::CellTag> ca;
        std::vector<CellHandle> live;
        std::vector<std::uint64_t> seqs;
        std::mt19937 rng(42);
        for (int i = 0; i < 2000; ++i) {
            if (!live.empty() && (rng() & 3) == 0) {
                const std::size_t k = rng() % live.size();
                ca.free(live[size_t(k)]);
                live.erase(live.begin() + long(k));
            } else {
                const auto h = ca.alloc();
                seqs.push_back(ca.seqOf(h));
                live.push_back(h);
            }
        }
        for (std::size_t i = 1; i < seqs.size(); ++i)
            CHECK(seqs[i] > seqs[i - 1]);
        CHECK(ca.accountingConsistent());
        CHECK(ca.totalCreated() == seqs.size());
    }
}

template <class Tag>
static void forceSeqToMax(topo::HandleAllocator<Tag>& a) {
    a.testForceSequenceMax();
}

static void testHandleStateIntegrity() {
    TEST("topo: allocator state integrity under retirement + exhaustion");
    using topo::HandleAllocator;
    using topo::NodeHandle;
    static_assert(std::is_nothrow_move_constructible_v<
                  HandleAllocator<topo::EdgeTag>>);
    static_assert(std::is_nothrow_move_constructible_v<
                  HandleAllocator<topo::FaceTag>>);
    // Retirement accounting: materialized == live + reusableFree + retired.
    {
        HandleAllocator<topo::NodeTag> a;
        auto g = a.alloc();
        const auto slot = g.slot;
        a.testForceGeneration(g, std::numeric_limits<std::uint32_t>::max());
        NodeHandle forced{slot, std::numeric_limits<std::uint32_t>::max(),
                          a.poolId()};
        CHECK(a.isCurrent(forced));
        a.free(forced);
        CHECK(a.testSlotRetired(slot));
        CHECK(a.retiredSlots() == 1 && a.liveCount() == 0 &&
              a.reusableFreeSlots() == 0 && a.materializedSlots() == 1);
        CHECK(a.accountingConsistent());        // review invariant fix
        // A retired slot is never reissued.
        const auto again = a.alloc();
        CHECK(again.slot != slot);
        CHECK(a.accountingConsistent());
    }
    // Sequence exhaustion preserves state, WITH and WITHOUT a free slot.
    {
        // No free slot: a throw must not leave a dead materialized slot.
        HandleAllocator<topo::NodeTag> a;
        auto h = a.alloc();
        (void)h;
        const auto beforeMat = a.materializedSlots();
        const auto beforeLive = a.liveCount();
        forceSeqToMax(a);                       // test helper below
        bool threw = false;
        try { (void)a.alloc(); }
        catch (const std::length_error&) { threw = true; }
        CHECK(threw);
        CHECK(a.materializedSlots() == beforeMat);   // no leaked slot
        CHECK(a.liveCount() == beforeLive);
        CHECK(a.accountingConsistent());
    }
    {
        // With a free slot: a throw must not consume the free slot.
        HandleAllocator<topo::NodeTag> a;
        auto h1 = a.alloc();
        auto h2 = a.alloc();
        a.free(h1);
        CHECK(a.reusableFreeSlots() == 1);
        const auto beforeFree = a.reusableFreeSlots();
        (void)h2;
        forceSeqToMax(a);
        bool threw = false;
        try { (void)a.alloc(); }
        catch (const std::length_error&) { threw = true; }
        CHECK(threw);
        CHECK(a.reusableFreeSlots() == beforeFree);   // free slot intact
        CHECK(a.accountingConsistent());
    }
    // Moved-from allocator is an EMPTY accounting-consistent husk (review:
    // moved-from vectors are valid, not guaranteed empty; moveFrom now
    // clears them explicitly).
    {
        HandleAllocator<topo::NodeTag> src;
        for (int i = 0; i < 100; ++i) (void)src.alloc();
        std::vector<topo::NodeHandle> hs;
        HandleAllocator<topo::NodeTag> src2;
        for (int i = 0; i < 50; ++i) hs.push_back(src2.alloc());
        for (int i = 0; i < 20; ++i) src2.free(hs[size_t(i)]);
        CHECK(src2.materializedSlots() > 0);
        HandleAllocator<topo::NodeTag> dst(std::move(src2));
        CHECK(src2.materializedSlots() == 0);        // explicitly cleared
        CHECK(src2.liveCount() == 0 && src2.reusableFreeSlots() == 0 &&
              src2.retiredSlots() == 0);
        CHECK(src2.accountingConsistent());
        CHECK(dst.accountingConsistent());
    }
    // Reviewed tokenless adoption: direct move consumes no new pool id and
    // preserves exact churn/free-list/retirement/sequence state. The source
    // receives the private pool-zero empty state and cannot allocate.
    {
        using A = HandleAllocator<topo::EdgeTag>;
        A src;
        const auto sourcePool = src.poolId();
        const auto h0 = src.alloc();
        const auto h1 = src.alloc();
        const auto h2 = src.alloc();
        const auto h3 = src.alloc();
        src.free(h1);
        src.free(h3);                         // free-list order: h1, h3
        src.testForceGeneration(
            h0, std::numeric_limits<std::uint32_t>::max());
        const topo::EdgeHandle retired{
            h0.slot, std::numeric_limits<std::uint32_t>::max(), sourcePool};
        src.free(retired);
        const auto liveSeq = src.seqOf(h2);
        const auto totalBefore = src.totalCreated();
        CHECK(src.liveCount() == 1 && src.reusableFreeSlots() == 2);
        CHECK(src.retiredSlots() == 1 && src.materializedSlots() == 4);

        A dst(std::move(src));
        CHECK(dst.poolId() == sourcePool);
        CHECK(dst.isCurrent(h2) && dst.seqOf(h2) == liveSeq);
        CHECK(dst.testSlotRetired(h0.slot));
        CHECK(dst.liveCount() == 1 && dst.reusableFreeSlots() == 2);
        CHECK(dst.retiredSlots() == 1 && dst.materializedSlots() == 4);
        CHECK(dst.totalCreated() == totalBefore && dst.accountingConsistent());
        CHECK(src.poolId() == 0 && src.materializedSlots() == 0);
        CHECK(src.liveCount() == 0 && src.reusableFreeSlots() == 0 &&
              src.retiredSlots() == 0 && src.totalCreated() == 0);
        CHECK(src.accountingConsistent());
        bool unusable = false;
        try { (void)src.alloc(); }
        catch (const std::logic_error&) { unusable = true; }
        CHECK(unusable);

        const auto reuseFirst = dst.alloc();
        const auto reuseSecond = dst.alloc();
        CHECK(reuseFirst.slot == h3.slot);     // exact LIFO free-list order
        CHECK(reuseSecond.slot == h1.slot);
        CHECK(dst.seqOf(reuseFirst) == totalBefore);
        CHECK(dst.seqOf(reuseSecond) == totalBefore + 1);
        A nextOrdinary;
        CHECK(nextOrdinary.poolId() == sourcePool + 1);  // move minted none
    }
    {
        using A = HandleAllocator<topo::FaceTag>;
        A empty;
        const auto sourcePool = empty.poolId();
        A moved(std::move(empty));
        CHECK(moved.poolId() == sourcePool && moved.materializedSlots() == 0);
        CHECK(empty.poolId() == 0 && empty.accountingConsistent());
        A nextOrdinary;
        CHECK(nextOrdinary.poolId() == sourcePool + 1);
    }
    // Slot-space boundary via the pure constexpr predicate (review: the
    // INT32_MAX boundary was previously untestable without billions of
    // slots).
    {
        constexpr std::size_t kMax =
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
        static_assert(!HandleAllocator<topo::NodeTag>::wouldExhaustSlots(0));
        static_assert(
            !HandleAllocator<topo::NodeTag>::wouldExhaustSlots(kMax - 1));
        static_assert(
            HandleAllocator<topo::NodeTag>::wouldExhaustSlots(kMax));
        static_assert(
            HandleAllocator<topo::NodeTag>::wouldExhaustSlots(kMax + 1));
        CHECK(!HandleAllocator<topo::NodeTag>::wouldExhaustSlots(1));
        CHECK(HandleAllocator<topo::NodeTag>::wouldExhaustSlots(kMax));
    }
    // Pool-ID exhaustion via the pure transition helper (review: the
    // saturating path was claimed tested but was not). Prove:
    //   MAX-1 -> MAX ; MAX -> throws ; and it keeps throwing.
    {
        using A = HandleAllocator<topo::NodeTag>;
        constexpr std::uint64_t kMaxId =
            std::numeric_limits<std::uint64_t>::max();
        // Ordinary progression and the reserved-0 boundary.
        static_assert(A::poolIdStep(0) == 1);        // ids start at 1
        static_assert(A::poolIdStep(1) == 2);
        static_assert(A::poolIdStep(kMaxId - 1) == kMaxId);
        CHECK(A::poolIdStep(kMaxId - 1) == kMaxId);
        bool atMax = false;
        try { (void)A::poolIdStep(kMaxId); }
        catch (const std::length_error&) { atMax = true; }
        CHECK(atMax);
        // Saturating: a subsequent call from the same exhausted value keeps
        // throwing (never wraps to reissue id 1).
        bool stillMax = false;
        try { (void)A::poolIdStep(kMaxId); }
        catch (const std::length_error&) { stillMax = true; }
        CHECK(stillMax);
    }
}

// ---------------------------------------------------------------------------
static void testTetArena() {
    TEST("topo: arena slot reuse, generations, epochs, memory counters");
    topo::TetArena a;
    const auto s0 = a.alloc({0, 1, 2, 3});
    const auto s1 = a.alloc({1, 2, 3, 4});
    CHECK(a.alive(s0) && a.alive(s1));
    // Generation-protected handles: freeing invalidates.
    const topo::TetHandle h0 = a.handleOf(s0);
    CHECK(a.isCurrent(h0));
    a.free(s0);
    CHECK(!a.isCurrent(h0));
    bool threw = false;
    try { (void)a.deref(h0); } catch (const std::logic_error&) { threw = true; }
    CHECK(threw);
    // Double free is a structured failure.
    threw = false;
    try { a.free(s0); } catch (const std::logic_error&) { threw = true; }
    CHECK(threw);
    // Slot reuse: the freed slot comes back; the stale handle still fails.
    const auto s2 = a.alloc({5, 6, 7, 8});
    CHECK(s2 == s0);
    CHECK(!a.isCurrent(h0));                       // generation moved on
    CHECK(a.isCurrent(a.handleOf(s2)));
    CHECK(a.cell(s2).seq == 2);                    // creation sequence
    // Epoch marks: O(1) invalidation between rounds.
    a.beginEpoch();
    a.markVisited(s1);
    CHECK(a.visited(s1) && !a.visited(s2));
    a.beginEpoch();
    CHECK(!a.visited(s1));
    // Memory counters.
    const auto m = a.memory();
    CHECK(m.liveCells == 2);
    CHECK(m.peakLiveCells == 2);
    CHECK(m.totalCellsCreated == 3);
    CHECK(m.materializedSlots == m.liveCells + m.freeSlots);
    CHECK(m.bytes > 0);
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
static std::uint64_t canonTetHash(const std::vector<algo::Tet>& tets) {
    // Order-independent canonical hash: sorted node ids per tet, tets sorted.
    std::vector<std::array<Index, 4>> keys;
    for (const auto& tt : tets) {
        auto k = tt.v;
        std::sort(k.begin(), k.end());
        keys.push_back(k);
    }
    std::sort(keys.begin(), keys.end());
    std::uint64_t h = 1469598103934665603ULL;
    for (const auto& k : keys)
        for (Index v : k)
            for (int b = 0; b < 8; ++b) {
                h ^= (std::uint64_t(v) >> (8 * b)) & 0xFF;
                h *= 1099511628211ULL;
            }
    return h;
}

static void testDelaunayGoldenHashes() {
    TEST("delaunay: golden canonical output hashes (M3-A frozen behavior)");
    // These values were captured from the frozen M3-A implementation. The
    // M4 topology port (and any later refactor) must reproduce them
    // exactly: same input => same set of tetrahedra. If an intentional
    // algorithmic change alters them, the change must be justified in the
    // closure log and the goldens re-frozen in the same commit.
    {
        std::vector<Vec3> pts;
        for (int k = 0; k < 5; ++k)
            for (int j = 0; j < 5; ++j)
                for (int i = 0; i < 5; ++i)
                    pts.push_back({double(i), double(j), double(k)});
        const auto tets = algo::delaunay3D(pts);
        CHECK(tets.size() == 384);
        CHECK(canonTetHash(tets) == 0xea7550e52137ee03ULL);
    }
    {
        std::mt19937 rng(777);
        std::vector<Vec3> pts;
        for (int i = 0; i < 300; ++i)
            pts.push_back({canonU01(rng), canonU01(rng), canonU01(rng)});
        const auto tets = algo::delaunay3D(pts);
        CHECK(tets.size() == 1777);
        CHECK(canonTetHash(tets) == 0x4720fb7e5d0738ccULL);
    }
    {
        std::mt19937 rng(778);
        std::vector<Vec3> pts;
        for (int i = 0; i < 200; ++i)
            pts.push_back(
                {canonU01(rng), canonU01(rng), canonU01(rng) * 1e-9});
        const auto tets = algo::delaunay3D(pts);
        CHECK(tets.size() == 1109);
        CHECK(canonTetHash(tets) == 0xccf396a3b736f2daULL);
    }
}

// ---------------------------------------------------------------------------
static void testBoundaryPeelingSwitch() {
    TEST("volume: MeshSpec::boundaryPeeling=false prevents peeling");
    auto meshCyl = [&](bool peel) {
        auto kernel = geo::makeBuiltinKernel();
        auto model  = kernel->makeCylinder({0, 0, 0}, 0.5, 1.0);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.15);
        spec.boundaryPeeling = peel;
        mesh::Mesh m;
        algo::meshModel(*model, spec, m);
        return mesh::analyzeTets(m);
    };
    const auto on  = meshCyl(true);
    const auto off = meshCyl(false);
    // With peeling off, nothing may be removed: at least as many tets and
    // at least as much volume, and this cylinder resolution is known to
    // trigger the peeler, so the inequality must be strict.
    CHECK(off.nTets >= on.nTets);
    CHECK(off.totalVolume >= on.totalVolume - 1e-12);
    CHECK(off.nTets > on.nTets);
    CHECK(off.nInverted == 0);
}

// ---------------------------------------------------------------------------
static void testPredicateConventions() {
    TEST("predicates: sign conventions, exact zeros, fallback instrumentation");
    using namespace predicates;
    using predicates::Sign;
    // orient2d: CCW positive; exact collinear Zero.
    CHECK(orient2d({0,0},{1,0},{0,1}) == Sign::Positive);
    CHECK(orient2d({0,0},{1,0},{0,-1}) == Sign::Negative);
    CHECK(orient2d({0,0},{1,0},{2,0}) == Sign::Zero);
    // incircle: d strictly inside circle through CCW (a,b,c); exact
    // on-circle Zero.
    CHECK(incircle({1,0},{0,1},{-1,0},{0,0}) == Sign::Positive);
    CHECK(incircle({1,0},{0,1},{-1,0},{2,2}) == Sign::Negative);
    CHECK(incircle({1,0},{0,1},{-1,0},{0,-1}) == Sign::Zero);
    // orient3d: (b-a) x (c-a) . (d-a); exact coplanar Zero; resolves a
    // subnormal offset the naive determinant cannot.
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{0,0,1}) == Sign::Positive);
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{0,0,-1}) == Sign::Negative);
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{3,4,0}) == Sign::Zero);
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{0.3,0.3,1e-320}) == Sign::Positive);
    // insphere: e strictly inside sphere through positively oriented
    // (a,b,c,d). The corner (1,1,1) is exactly ON the circumsphere of the
    // unit corner tet (center (.5,.5,.5), r^2 = .75).
    const Vec3 p0{0,0,0}, p1{1,0,0}, p2{0,1,0}, p3{0,0,1};
    CHECK(orient3d(p0,p1,p2,p3) == Sign::Positive);
    CHECK(insphere(p0,p1,p2,p3,{0.5,0.5,0.5}) == Sign::Positive);
    CHECK(insphere(p0,p1,p2,p3,{9,9,9}) == Sign::Negative);
    CHECK(insphere(p0,p1,p2,p3,{1,1,1}) == Sign::Zero);
    CHECK(insphere(p0,p1,p2,p3,{1,1,1.0-1e-13}) == Sign::Positive);
    CHECK(insphere(p0,p1,p2,p3,{1,1,1.0+1e-13}) == Sign::Negative);
    // Fallback instrumentation: an exact-zero case MUST take the exact
    // path; a well-conditioned case MUST take the fast path.
    resetStats();
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{0.25,0.25,0}) == Sign::Zero);
    CHECK(stats().orient3dExact == 1);
    CHECK(orient3d({0,0,0},{1,0,0},{0,1,0},{0,0,2}) == Sign::Positive);
    CHECK(stats().orient3dFast == 1);
    resetStats();
    CHECK(insphere(p0,p1,p2,p3,{1,1,1}) == Sign::Zero);
    CHECK(stats().insphereExact == 1);
}

// ---------------------------------------------------------------------------
static void testPredicateDifferential() {
    TEST("predicates: differential vs exact int128 reference (dyadic grid)");
#if !defined(__SIZEOF_INT128__)
    // MSVC has no __int128; the direct convention/zero/fallback tests in
    // predicates.conventions still run there. An MSVC-compatible exact
    // reference (e.g. 64x64->128 intrinsics) is a portability follow-up.
    std::puts("  (skipped: no __int128 on this compiler)");
    return;
#else
    // Coordinates are integers in [-2^20, 2^20] scaled by 2^-10: exactly
    // representable doubles whose predicate determinants fit __int128
    // (orient3d ~ 2^66, insphere ~ 2^130 > int128!). To keep insphere in
    // range use a smaller grid for it: integers in [-2^14, 2^14], lifted
    // terms ~2^29, minors ~2^45, products ~2^75 - safe.
    __extension__ typedef __int128 I;   // GCC/Clang; MSVC needs _mul128 - test is POSIX-CI scoped
    auto sgn = [](I v) {
        return v > 0 ? predicates::Sign::Positive
             : v < 0 ? predicates::Sign::Negative : predicates::Sign::Zero;
    };
    std::mt19937_64 rng(2026);
    auto grid = [&](long lim) {
        const unsigned long long span = 2ULL * static_cast<unsigned long long>(lim) + 1ULL;
        return double(static_cast<long long>(rng() % span) - lim);
    };
    // orient3d differential (includes engineered degeneracies by drawing
    // from a coarse grid where coincidences are frequent).
    int zeros3 = 0;
    for (int it = 0; it < 4000; ++it) {
        const long lim = it % 2 ? 8 : (1L << 20);   // coarse & fine grids
        const double s = it % 2 ? 1.0 : 1.0 / 1024.0;
        double c[12];
        for (double& x : c) x = grid(lim) * s;
        Vec3 a{c[0],c[1],c[2]}, b{c[3],c[4],c[5]}, d{c[6],c[7],c[8]},
             e{c[9],c[10],c[11]};
        if (it % 5 == 0)   // engineered exact coplanarity: e = a + (b-a)+(d-a)
            e = {b.x + d.x - a.x, b.y + d.y - a.y, b.z + d.z - a.z};
        // Reference det[b-a; d-a; e-a] over exact integers (scale cancels
        // for the SIGN).
        auto ix = [&](double v) { return I(llround(v / s)); };
        const I r0x=ix(b.x)-ix(a.x), r0y=ix(b.y)-ix(a.y), r0z=ix(b.z)-ix(a.z);
        const I r1x=ix(d.x)-ix(a.x), r1y=ix(d.y)-ix(a.y), r1z=ix(d.z)-ix(a.z);
        const I r2x=ix(e.x)-ix(a.x), r2y=ix(e.y)-ix(a.y), r2z=ix(e.z)-ix(a.z);
        const I det = r0x*(r1y*r2z - r1z*r2y) - r0y*(r1x*r2z - r1z*r2x)
                    + r0z*(r1x*r2y - r1y*r2x);
        if (det == 0) ++zeros3;
        CHECK(predicates::orient3d(a, b, d, e) == sgn(det));
    }
    CHECK(zeros3 >= 4000 / 5);   // engineered coplanar cases must register as zeros
    // insphere differential on the small grid.
    int zerosS = 0, tested = 0;
    for (int it = 0; it < 2000; ++it) {
        const long lim = it % 2 ? 4 : (1L << 14);
        double c[15];
        for (double& x : c) x = grid(lim);
        Vec3 a{c[0],c[1],c[2]}, b{c[3],c[4],c[5]}, d{c[6],c[7],c[8]},
             e{c[9],c[10],c[11]}, q{c[12],c[13],c[14]};
        if (it % 7 == 0) {
            // Engineered exact cosphericality: unit-cube corner tets have
            // the antipodal corner exactly on their circumsphere; translate
            // by a random grid offset.
            const double ox = grid(lim / 2), oy = grid(lim / 2),
                         oz = grid(lim / 2);
            a = {ox, oy, oz}; b = {ox + 1, oy, oz};
            d = {ox, oy + 1, oz}; e = {ox, oy, oz + 1};
            q = {ox + 1, oy + 1, oz + 1};
        }
        auto ii = [&](double v) { return I(llround(v)); };
        // Orientation first (reference), skip degenerate tets.
        const I o = (ii(b.x)-ii(a.x))*((ii(d.y)-ii(a.y))*(ii(e.z)-ii(a.z))
                                     - (ii(d.z)-ii(a.z))*(ii(e.y)-ii(a.y)))
                  - (ii(b.y)-ii(a.y))*((ii(d.x)-ii(a.x))*(ii(e.z)-ii(a.z))
                                     - (ii(d.z)-ii(a.z))*(ii(e.x)-ii(a.x)))
                  + (ii(b.z)-ii(a.z))*((ii(d.x)-ii(a.x))*(ii(e.y)-ii(a.y))
                                     - (ii(d.y)-ii(a.y))*(ii(e.x)-ii(a.x)));
        if (o == 0) continue;
        // Reference lifted determinant, rows (v - q), for the tet ordered
        // positively.
        Vec3 A = a, B = b, C = d, D = e;
        if (o < 0) std::swap(C, D);
        auto row = [&](const Vec3& v, I r[4]) {
            r[0] = ii(v.x) - ii(q.x);
            r[1] = ii(v.y) - ii(q.y);
            r[2] = ii(v.z) - ii(q.z);
            r[3] = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        };
        I ra[4], rb[4], rc[4], rd[4];
        row(A, ra); row(B, rb); row(C, rc); row(D, rd);
        auto det3 = [](const I* u, const I* v, const I* w) {
            return u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0])
                 + u[2]*(v[0]*w[1]-v[1]*w[0]);
        };
        const I lifted = -ra[3]*det3(rb,rc,rd) + rb[3]*det3(ra,rc,rd)
                       - rc[3]*det3(ra,rb,rd) + rd[3]*det3(ra,rb,rc);
        // Locked convention: lifted determinant NEGATIVE when q inside.
        const predicates::Sign expect =
            lifted < 0 ? predicates::Sign::Positive
          : lifted > 0 ? predicates::Sign::Negative : predicates::Sign::Zero;
        if (lifted == 0) ++zerosS;
        CHECK(predicates::insphere(A, B, C, D, q) == expect);
        ++tested;
    }
    CHECK(tested > 1000);
    CHECK(zerosS >= 2000 / 7 - 5);   // engineered cospherical cases
#endif
}

// ---------------------------------------------------------------------------
static void testPredicateAffine() {
    TEST("predicates: exact affine orient2d + segment ordering (BND-2A)");
    using predicates::Sign;
    using predicates::AffinePoint2;
    auto S = [](double u, double v) { return AffinePoint2{{u, v}, {u, v}, 0, 0}; };
    auto threws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; } return false;
    };
    // ---- input validation: invalid tokens fail loudly ----------------------
    {
        // level 0 but numerator!=0, or lo!=hi -> invalid stored form.
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{0,0},1,0}, S(1,0), S(0,1)); }));
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{1,0},0,0}, S(1,0), S(0,1)); }));
        // level > 52.
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{1,0},1,53}, S(1,0), S(0,1)); }));
        // numerator == 0 or >= 2^level for a genuine affine level.
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{1,0},0,3}, S(1,0), S(0,1)); }));
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{1,0},8,3}, S(1,0), S(0,1)); }));
        // unreduced (even) numerator.
        CHECK(threws([&]{ predicates::orient2dAffine(
            AffinePoint2{{0,0},{1,0},2,3}, S(1,0), S(0,1)); }));
        // non-finite coordinate.
        const double inf = std::numeric_limits<double>::infinity();
        CHECK(threws([&]{ predicates::orient2dAffine(
            S(inf,0), S(1,0), S(0,1)); }));
        // degenerate (zero-length) ordering segment.
        CHECK(threws([&]{ predicates::orderOnSegmentAffine(
            {2,2}, {2,2}, S(2,2), S(2,2)); }));
        // non-finite ordering segment endpoint.
        CHECK(threws([&]{ predicates::orderOnSegmentAffine(
            {inf,0}, {1,0}, S(0,0), S(1,0)); }));
    }
    // ---- analytic level-52 (ODD numerators), all platforms -----------------
    {
        AffinePoint2 mid{{0,0},{4,4}, (1ULL<<52)-1, 52};   // t just below 1, on y=x
        CHECK(predicates::orient2dAffine(S(0,0), S(4,4), mid) == Sign::Zero);
        AffinePoint2 a{{0,0},{2,0}, 1, 52};                // ~origin on x-axis
        CHECK(predicates::orient2dAffine(a, S(10,0), S(0,5)) == Sign::Positive);
        CHECK(predicates::orient2dAffine(a, S(10,0), S(0,-5)) == Sign::Negative);
    }
    // ---- exact segment ordering: Stored/Affine, both directions ------------
    {
        Vec2 lo{0,0}, hi{8,0};
        auto Aff = [](std::uint64_t n, int L) {
            return AffinePoint2{{0,0},{8,0}, n, static_cast<std::uint8_t>(L)}; };
        CHECK(predicates::orderOnSegmentAffine(lo,hi,S(2,0),S(5,0)) == Sign::Negative);
        CHECK(predicates::orderOnSegmentAffine(lo,hi,S(5,0),S(2,0)) == Sign::Positive);
        CHECK(predicates::orderOnSegmentAffine(lo,hi,S(3,0),S(3,0)) == Sign::Zero);
        CHECK(predicates::orderOnSegmentAffine(lo,hi,Aff(1,1),S(3,0)) == Sign::Positive);
        CHECK(predicates::orderOnSegmentAffine(lo,hi,S(3,0),Aff(1,1)) == Sign::Negative);
        CHECK(predicates::orderOnSegmentAffine(lo,hi,Aff(1,2),Aff(3,2)) == Sign::Negative);
        CHECK(predicates::orderOnSegmentAffine(hi,lo,S(2,0),S(5,0)) == Sign::Positive); // reversed dir
        // level-52 odd-numerator tie-break.
        Vec2 l0{0,0}, h0{1,0};
        AffinePoint2 p{{0,0},{1,0}, (1ULL<<51)+1, 52};
        AffinePoint2 q{{0,0},{1,0}, (1ULL<<51)+3, 52};
        CHECK(predicates::orderOnSegmentAffine(l0,h0,p,q) == Sign::Negative);
        CHECK(predicates::orderOnSegmentAffine(l0,h0,q,p) == Sign::Positive);
    }
#if !defined(__SIZEOF_INT128__)
    std::puts("  (randomized oracle skipped: no __int128 on this compiler)");
    return;
#else
    __extension__ typedef __int128 I;
    auto isgn = [](I v) {
        return v > 0 ? Sign::Positive : v < 0 ? Sign::Negative : Sign::Zero; };
    std::mt19937_64 rng(4242);
    struct AP { long lo_u, lo_v, hi_u, hi_v; unsigned long num; int level; };
    auto coordNumOverL = [](long lo, long hi, unsigned long num, int level,
                            int Lmax) -> I {
        I pn = I(lo) * (I(1) << level) + I(num) * (I(hi) - lo);   // over 2^level
        return pn * (I(1) << (Lmax - level));                    // over 2^Lmax
    };
    auto toAP = [](const AP& a) {
        return AffinePoint2{{double(a.lo_u), double(a.lo_v)},
                            {double(a.hi_u), double(a.hi_v)},
                            static_cast<std::uint64_t>(a.num),
                            static_cast<std::uint8_t>(a.level)}; };
    // coords in [-2,2] keep even level-52 determinants inside __int128.
    auto rc = [&]() { return long(rng() % 5) - 2; };
    auto oddNum = [&](int L) -> unsigned long {
        // odd numerator in (0, 2^L); emphasize 1 and 2^L-1.
        const unsigned long hiN = (1UL << L) - 1;
        const unsigned r = unsigned(rng() % 8);
        if (r == 0) return 1UL;
        if (r == 1) return hiN;                 // odd (2^L even -> -1 odd)
        unsigned long n = 1UL + 2UL * (rng() % (((1UL << L) - 1 + 1) / 2));
        if (n >= (1UL << L)) n = hiN;
        return n | 1UL;
    };
    auto mk = [&]() -> AP {
        AP p; p.lo_u = rc(); p.lo_v = rc();
        if (rng() & 1) {
            int Lv = 1 + int(rng() % 52);        // FULL 1..52 range
            p.hi_u = rc(); p.hi_v = rc(); p.level = Lv; p.num = oddNum(Lv);
        } else { p.hi_u = p.lo_u; p.hi_v = p.lo_v; p.num = 0; p.level = 0; }
        return p; };
    auto oracle = [&](const AP& A, const AP& B, const AP& C) {
        const int L = std::max({A.level, B.level, C.level});
        const I ax = coordNumOverL(A.lo_u,A.hi_u,A.num,A.level,L),
                ay = coordNumOverL(A.lo_v,A.hi_v,A.num,A.level,L),
                bx = coordNumOverL(B.lo_u,B.hi_u,B.num,B.level,L),
                by = coordNumOverL(B.lo_v,B.hi_v,B.num,B.level,L),
                cx = coordNumOverL(C.lo_u,C.hi_u,C.num,C.level,L),
                cy = coordNumOverL(C.lo_v,C.hi_v,C.num,C.level,L);
        return isgn((bx-ax)*(cy-ay) - (by-ay)*(cx-ax)); };
    int naff[4] = {0,0,0,0}, zeros = 0;
    for (int it = 0; it < 200000; ++it) {
        AP A = mk(), B = mk(), C = mk();
        const Sign o = oracle(A, B, C);
        const AffinePoint2 pa = toAP(A), pb = toAP(B), pc = toAP(C);
        CHECK(predicates::orient2dAffine(pa, pb, pc) == o);
        // permutation/sign identities: one swap flips; cyclic preserves.
        const Sign flip = o == Sign::Positive ? Sign::Negative
                        : o == Sign::Negative ? Sign::Positive : Sign::Zero;
        CHECK(predicates::orient2dAffine(pb, pa, pc) == flip);
        CHECK(predicates::orient2dAffine(pa, pc, pb) == flip);
        CHECK(predicates::orient2dAffine(pb, pc, pa) == o);   // cyclic
        ++naff[(A.level>0)+(B.level>0)+(C.level>0)];
        if (o == Sign::Zero) ++zeros;
    }
    CHECK(naff[0] > 0 && naff[1] > 0 && naff[2] > 0 && naff[3] > 0);
    CHECK(zeros > 50);
    // stored-only agreement with ordinary orient2d.
    for (int it = 0; it < 20000; ++it) {
        Vec2 a{double(rc()),double(rc())}, b{double(rc()),double(rc())},
             c{double(rc()),double(rc())};
        CHECK(predicates::orient2dAffine(S(a.u,a.v),S(b.u,b.v),S(c.u,c.v))
              == predicates::orient2d(a,b,c));
    }
    // randomized exact ordering on a common integer segment, both directions,
    // vs an independent uint64 dyadic-alignment oracle (levels 1..52).
    int ord = 0;
    for (int it = 0; it < 40000; ++it) {
        const long sx = 1 + long(rng() % 4);         // segment (0,0)->(sx,0)
        Vec2 lo{0,0}, hi{double(sx),0};
        auto mkOn = [&](unsigned long& num, int& lvl) {
            if (rng() & 1) { lvl = 1 + int(rng() % 52); num = oddNum(lvl);
                return AffinePoint2{{0,0},{double(sx),0}, num,
                                    static_cast<std::uint8_t>(lvl)}; }
            lvl = 0; num = 0;                          // stored endpoint lo or hi
            return (rng() & 1) ? S(0,0) : S(double(sx),0);
        };
        unsigned long np=0,nq=0; int lp=0,lq=0;
        AffinePoint2 P = mkOn(np,lp), Q = mkOn(nq,lq);
        // oracle t as num/2^level; stored endpoints are t=0 (lo) or t=1 (hi).
        auto tAligned = [&](const AffinePoint2& X, unsigned long num, int lvl,
                            int L) -> unsigned long long {
            if (lvl == 0) return (X.lo.u == 0.0) ? 0ULL
                                 : (1ULL << L);       // hi -> t=1
            return (unsigned long long)num << (L - lvl);
        };
        const int L = std::max(std::max(lp,1), std::max(lq,1));
        const unsigned long long tp = tAligned(P,np,lp,L),
                                 tq = tAligned(Q,nq,lq,L);
        const Sign exp = tp < tq ? Sign::Negative
                       : tp > tq ? Sign::Positive : Sign::Zero;
        CHECK(predicates::orderOnSegmentAffine(lo,hi,P,Q) == exp);
        const Sign rexp = exp == Sign::Positive ? Sign::Negative
                        : exp == Sign::Negative ? Sign::Positive : Sign::Zero;
        CHECK(predicates::orderOnSegmentAffine(hi,lo,P,Q) == rexp);
        ++ord;
    }
    CHECK(ord == 40000);
#endif
}

// ---------------------------------------------------------------------------
static void testPredicateCompareOnAxisAffine() {
    TEST("predicates: exact compareOnAxisAffine (BND-2A frozen contact)");
    using predicates::Sign;
    using predicates::AffinePoint2;
    auto S = [](double u, double v) { return AffinePoint2{{u, v}, {u, v}, 0, 0}; };
    auto threws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; } return false;
    };
    // ---- input validation: malformed tokens + invalid axes fail loudly -----
    CHECK(threws([&]{ predicates::compareOnAxisAffine(S(0,0), S(1,0), 2); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(S(0,0), S(1,0), -1); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(       // level 0, num!=0
        AffinePoint2{{0,0},{0,0},1,0}, S(1,0), 0); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(       // lo!=hi at level 0
        AffinePoint2{{0,0},{1,0},0,0}, S(1,0), 0); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(       // level > 52
        AffinePoint2{{0,0},{1,0},1,53}, S(1,0), 0); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(       // num >= 2^level
        AffinePoint2{{0,0},{1,0},8,3}, S(1,0), 0); }));
    CHECK(threws([&]{ predicates::compareOnAxisAffine(       // even (unreduced)
        AffinePoint2{{0,0},{1,0},2,3}, S(1,0), 0); }));
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(threws([&]{ predicates::compareOnAxisAffine(S(inf,0), S(1,0), 0); }));
    // ---- both axes; stored/affine combinations -----------------------------
    CHECK(predicates::compareOnAxisAffine(S(1,0), S(3,0), 0) == Sign::Negative);
    CHECK(predicates::compareOnAxisAffine(S(3,0), S(1,0), 0) == Sign::Positive);
    CHECK(predicates::compareOnAxisAffine(S(0,2), S(0,7), 1) == Sign::Negative);
    CHECK(predicates::compareOnAxisAffine(S(5,5), S(5,5), 0) == Sign::Zero);  // eq
    CHECK(predicates::compareOnAxisAffine(S(5,5), S(5,5), 1) == Sign::Zero);
    const AffinePoint2 mid{{0,0},{4,0}, 1, 1};   // t=1/2 -> (2,0)
    CHECK(predicates::compareOnAxisAffine(S(0,0), mid, 0) == Sign::Negative);
    CHECK(predicates::compareOnAxisAffine(mid, S(4,0), 0) == Sign::Negative);  // 2 < 4
    CHECK(predicates::compareOnAxisAffine(mid, S(2,0), 0) == Sign::Zero);    // affine==stored
    // level-52 odd-numerator, both directions (antisymmetry).
    const AffinePoint2 p{{0,0},{1,0}, (1ULL<<51)+1, 52};
    const AffinePoint2 q{{0,0},{1,0}, (1ULL<<51)+3, 52};
    CHECK(predicates::compareOnAxisAffine(p, q, 0) == Sign::Negative);
    CHECK(predicates::compareOnAxisAffine(q, p, 0) == Sign::Positive);
    CHECK(predicates::compareOnAxisAffine(p, p, 0) == Sign::Zero);
#if defined(__SIZEOF_INT128__)
    // ---- differential vs exact int128 oracle, levels 1..52, both axes ------
    __extension__ typedef __int128 I;
    auto isgn = [](I v) {
        return v > 0 ? Sign::Positive : v < 0 ? Sign::Negative : Sign::Zero; };
    std::mt19937_64 rng(90210);
    auto rc = [&]() { return long(rng() % 5) - 2; };
    auto oddNum = [&](int L) -> unsigned long {
        const unsigned long hiN = (1UL << L) - 1;
        unsigned long n = 1UL + 2UL * (rng() % (((1UL << L)) / 2));
        if (n >= (1UL << L)) n = hiN;
        return n | 1UL;
    };
    struct AP { long lo_u, lo_v, hi_u, hi_v; unsigned long num; int level; };
    auto mk = [&]() -> AP {
        AP a; a.lo_u = rc(); a.lo_v = rc();
        if (rng() & 1) { int L = 1 + int(rng()%52); a.hi_u=rc(); a.hi_v=rc();
                         a.level=L; a.num=oddNum(L); }
        else { a.hi_u=a.lo_u; a.hi_v=a.lo_v; a.num=0; a.level=0; }
        return a; };
    auto toAP = [](const AP& a) {
        return AffinePoint2{{double(a.lo_u),double(a.lo_v)},
                            {double(a.hi_u),double(a.hi_v)},
                            static_cast<std::uint64_t>(a.num),
                            static_cast<std::uint8_t>(a.level)}; };
    // value on axis, scaled to a common level, as an exact int128.
    auto coord = [](long lo, long hi, unsigned long num, int level, int Lmax) -> I {
        I pn = I(lo) * (I(1) << level) + I(num) * (I(hi) - lo);
        return pn * (I(1) << (Lmax - level)); };
    int checks = 0, nz = 0;
    for (int it = 0; it < 120000; ++it) {
        AP a = mk(), b = mk();
        for (int axis = 0; axis < 2; ++axis) {
            const int L = std::max(a.level, b.level);
            const I av = axis==0 ? coord(a.lo_u,a.hi_u,a.num,a.level,L)
                                 : coord(a.lo_v,a.hi_v,a.num,a.level,L);
            const I bv = axis==0 ? coord(b.lo_u,b.hi_u,b.num,b.level,L)
                                 : coord(b.lo_v,b.hi_v,b.num,b.level,L);
            const Sign exp = isgn(av - bv);
            const AffinePoint2 pa = toAP(a), pb = toAP(b);
            CHECK(predicates::compareOnAxisAffine(pa, pb, axis) == exp);
            const Sign r = exp==Sign::Positive ? Sign::Negative
                         : exp==Sign::Negative ? Sign::Positive : Sign::Zero;
            CHECK(predicates::compareOnAxisAffine(pb, pa, axis) == r);  // antisym
            ++checks; if (exp == Sign::Zero) ++nz;
        }
    }
    CHECK(checks == 240000);
    CHECK(nz > 50);
#else
    std::puts("  (compareOnAxisAffine oracle skipped: no __int128)");
#endif
}

// ---------------------------------------------------------------------------
static void testDelaunayRejectsNonFinite() {
    TEST("delaunay: non-finite input is a structured failure, not a crash");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    auto throws = [](std::vector<Vec3> pts) {
        try { algo::delaunay3D(pts); }
        catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };
    CHECK(throws({{0,0,0},{1,0,0},{0,1,0},{0,0,1},{nan,0,0}}));
    CHECK(throws({{0,0,0},{1,0,0},{0,1,0},{0,0,1},{0,inf,0}}));
    CHECK(throws({{0,0,0},{1,0,0},{0,1,0},{0,0,-inf}}));
    // Review finding: validation previously ran AFTER the n<4 early return,
    // so short non-finite inputs were silently "skipped" instead of thrown.
    CHECK(throws({{nan,0,0}}));
    CHECK(throws({{0,0,0},{inf,0,0},{0,1,0}}));
}

// ---------------------------------------------------------------------------
// Topology tests live in tests/test_topology.cpp (split so the monolithic
// TU stays compilable in constrained review environments).
namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* topologyTests(std::size_t& count);
const Entry* adjacencyTests(std::size_t& count);
const Entry* transactionTests(std::size_t& count);
const Entry* bndTests(std::size_t& count);
const Entry* bndEmbeddingTests(std::size_t& count);
const Entry* bndSessionTests(std::size_t& count);
const Entry* bndCommitTests(std::size_t& count);
const Entry* bnd2Tests(std::size_t& count);
const Entry* faceAdapterTests(std::size_t& count);
const Entry* surfaceAdapterTests(std::size_t& count);
const Entry* delaunayMetricTests(std::size_t& count);
const Entry* delaunayScalableTests(std::size_t& count);
const Entry* productionSurfaceTests(std::size_t& count);
const Entry* constrainedV2Tests(std::size_t& count);
const Entry* tetHandoffTests(std::size_t& count);
const Entry* tetFacetTests(std::size_t& count);
const Entry* t2bTests(std::size_t& count);
const Entry* t3aTests(std::size_t& count);
const Entry* t4aTests(std::size_t& count);
const Entry* t4bTests(std::size_t& count);
const Entry* t5Tests(std::size_t& count);
const Entry* t2cCacheTests(std::size_t& count);
const Entry* volumeTypeTests(std::size_t& count);
const Entry* meshControlLayerTests(std::size_t& count);
} // namespace paacs_tests

int main(int argc, char** argv) {
    // Preserve the current fixture name if a fail-fast invariant terminates
    // the process while CTest is capturing redirected output.
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    // Registry: names group tests for CTest granularity. With no argument
    // every test runs; with an argument only tests whose registry name
    // contains it as a substring run (e.g. `paacs_mesh_tests delaunay`).
    struct Entry { const char* name; void (*fn)(); };
    static const Entry kTests[] = {
        {"edges.discretization",      testEdgeDiscretization},
        {"delaunay.2d",               testDelaunay2DEmptyCircumcircle},
        {"delaunay.3d",               testDelaunay3DValidity},
        {"pipeline.box",              testVolumeMeshBox},
        {"sizing.ball",               testRefinementBall},
        {"elevation.basic",           testQuadraticElevation},
        {"io.vtk",                    testVtkExport},
        {"quality.metrics",           testQualityMetrics},
        {"optimize.smoothing",        testSmoothing},
        {"quality.report",            testQualityReport},
        {"pipeline.cylinder",         testCylinderPipeline},
        {"elevation.curved",          testCurvedQuadratic},
        {"optimize.flips",            testFlipOptimization},
        {"io.msh",                    testMshExport},
        {"parallel.pointgrid",        testPointGrid},
        {"parallel.determinism",      testThreadDeterminism},
        {"parallel.straggler",        testPoolStragglerRace},
        {"parallel.reject",           testPoolConcurrentRejection},
        {"parallel.nested",           testPoolNestedDispatch},
        {"parallel.exceptions",       testPoolExceptionPropagation},
        {"sizing.baked",              testBakedField},
        {"elevation.midside",         testMidsideClassification},
        {"elevation.jacobian",        testTet10Jacobian},
        {"mesh.dedup",                testDedupLifecycle},
        {"quality.reportsafety",      testQualityReportSafety},
        {"numerics.invariance",       testScaleAndTranslationInvariance},
        {"robustness.thin",           testThinDomains},
        {"robustness.flattened",      testFlattenedClouds},
        {"robustness.lattice",        testExactLattice},
        {"determinism.meshhash",      testFullMeshHashDeterminism},
        {"boundary.audit",            testBoundaryAuditKnownGap},
        {"audit.clean",               testMeshAuditClean},
        {"audit.detects",             testMeshAuditDetects},
        {"audit.tet10hardening",      testMeshAuditTet10AndHardening},
        {"audit.boundarymodes",       testMeshAuditBoundaryModes},
        {"topo.traits",               testElementTraits},
        {"io.writers",                testWriterElementMappings},
        {"topo.handles",              testGenericHandles},
        {"topo.handles.integrity",    testHandleStateIntegrity},
        {"topo.arena",                testTetArena},
        {"delaunay.golden",           testDelaunayGoldenHashes},
        {"baseline.hashes",           testFrozenBaselineHashes},
        {"hash.semantic",             testSemanticHashCanonicality},
        {"volume.peelswitch",         testBoundaryPeelingSwitch},
        {"predicates.conventions",    testPredicateConventions},
        {"predicates.differential",   testPredicateDifferential},
        {"predicates.affine",         testPredicateAffine},
        {"predicates.axiscompare",    testPredicateCompareOnAxisAffine},
        {"delaunay.nonfinite",        testDelaunayRejectsNonFinite},
    };
    const std::string filter = argc > 1 ? argv[1] : "";
    int ran = 0;
    for (const auto& t : kTests) {
        if (!filter.empty() &&
            std::string(t.name).find(filter) == std::string::npos)
            continue;
        t.fn();
        ++ran;
    }
    // External test tables (tests/test_topology.cpp,
    // tests/test_adjacency.cpp), same filter rules.
    for (auto hook : {paacs_tests::topologyTests, paacs_tests::adjacencyTests,
                      paacs_tests::transactionTests,
                      paacs_tests::bndTests,
                      paacs_tests::bndEmbeddingTests,
                      paacs_tests::bndSessionTests,
                      paacs_tests::bndCommitTests,
                      paacs_tests::bnd2Tests,
                      paacs_tests::faceAdapterTests,
                      paacs_tests::surfaceAdapterTests,
                      paacs_tests::delaunayMetricTests,
                      paacs_tests::delaunayScalableTests,
                      paacs_tests::productionSurfaceTests,
                      paacs_tests::constrainedV2Tests,
                      paacs_tests::tetHandoffTests,
                      paacs_tests::tetFacetTests,
                       paacs_tests::t2bTests,
                       paacs_tests::t3aTests,
                       paacs_tests::t4aTests,
                       paacs_tests::t4bTests,
                       paacs_tests::t5Tests,
                       paacs_tests::t2cCacheTests,
                       paacs_tests::volumeTypeTests,
                       paacs_tests::meshControlLayerTests}) {
        std::size_t n = 0;
        const paacs_tests::Entry* ext = hook(n);
        for (std::size_t i = 0; i < n; ++i) {
            if (!filter.empty() &&
                std::string(ext[i].name).find(filter) == std::string::npos)
                continue;
            ext[i].fn();
            ++ran;
        }
    }
    if (!ran) {
        std::printf("no tests match filter '%s'\n", filter.c_str());
        return 2;
    }
    if (gFailures) {
        std::printf("\n%d FAILURE(S)\n", gFailures);
        return 1;
    }
    std::printf("\nALL TESTS PASSED (%d tests run)\n", ran);
    return 0;
}
