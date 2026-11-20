// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "complex_geometry_fixture.hpp"
#include "geometry_checks.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
using namespace paacs;
using geometry_test::digest;
using geometry_test::require;
namespace {
int meshes = 0, refusals = 0;
void boundaryEuler(const mesh::Mesh &m, std::size_t holes) {
    std::set<Index> vertices;
    std::map<std::pair<Index, Index>, int> edges;
    long faces = 0;
    for (Index i = 0; i < m.cellCount(); ++i)
        if (m.cell(i).onEntity.dim == 2) {
            ++faces;
            const auto *n = m.cellNodes(i);
            const int count = mesh::nodesPerCell(m.cell(i).type);
            for (int j = 0; j < count; ++j) {
                vertices.insert(n[j]);
                ++edges[std::minmax(n[j], n[(j + 1) % count])];
            }
        }
    for (const auto &[key, n] : edges) {
        (void)key;
        require(n == 2, "open/nonmanifold exterior edge");
    }
    require(static_cast<long>(vertices.size()) - static_cast<long>(edges.size()) + faces ==
                2 - 2 * static_cast<long>(holes),
            "lost tunnel / exterior Euler characteristic");
}
sizing::CappedSizeOutcome profile(double length, double scale, int mode) {
    sizing::CappedSizeSpecification s;
    s.length = length;
    s.globalSize = .45 * scale;
    if (mode == 2)
        s.controls = {{length * .38, length * .60, .12 * scale, 1.25}};
    else if (mode == 3)
        s.controls = {{0, 0, .12 * scale, 1.25}, {length, length, .16 * scale, 1.2}};
    else
        s.controls = {{0, 0, .12 * scale, 1.25}};
    auto c = sizing::compileCappedSizeProfile(s);
    require(c.compiled, "compile profile");
    return c;
}
void stations(const mesh::Mesh &m, const geo::Model &model, const algo::GeometryBodyReport &b,
              const sizing::CappedSizeProfile &p) {
    const auto *face = model.face(b.sourceFaceTag);
    Vec2 lo, hi;
    face->paramBounds(lo, hi);
    const auto origin = face->pointAt((lo + hi) * .5), normal = face->normalAt((lo + hi) * .5);
    const double denominator = normal.dot(b.sweep);
    std::vector<double> z;
    z.reserve(static_cast<std::size_t>(m.nodeCount()));
    for (Index i = 0; i < m.nodeCount(); ++i)
        z.push_back((m.node(i) - origin).dot(normal) / denominator * p.length());
    std::sort(z.begin(), z.end());
    std::vector<double> unique;
    for (double x : z)
        if (unique.empty() || std::abs(x - unique.back()) > p.length() * 1e-9)
            unique.push_back(x);
    require(unique.size() == static_cast<std::size_t>(b.layers) + 1, "wrong station count");
    require(std::abs(unique.front()) < p.length() * 1e-8 &&
                std::abs(unique.back() - p.length()) < p.length() * 1e-8,
            "wrong domain closure");
    unique.front() = 0;
    unique.back() = p.length();
    double previous = 0, ratio = 1;
    for (std::size_t i = 1; i < unique.size(); ++i) {
        const double w = unique[i] - unique[i - 1];
        double target = std::min(p.eval(unique[i - 1]), p.eval(unique[i]));
        for (const auto k : p.knots())
            if (k.position > unique[i - 1] && k.position < unique[i])
                target = std::min(target, k.size);
        require(w > 0 && w <= target * (1 + 1e-8),
                "stored station violates full-interval size cap");
        if (previous)
            ratio = std::max({ratio, w / previous, previous / w});
        previous = w;
    }
    require(ratio <= p.maximumGrowth() * (1 + 1e-8) && ratio <= 1.3 * (1 + 1e-8),
            "stored station violates growth");
}
void run(const complex_test::Shape &s, int mode, bool tets, double scale, bool rotated,
         const std::filesystem::path &gallery, bool compare) {
    const auto begin = std::chrono::steady_clock::now();
    auto built = geo::makePolyhedralModel(complex_test::input(s, scale, rotated));
    require(built.built(), std::string("geometry: ") + geo::polyhedralReasonName(built.reason));
    const double length = std::hypot(s.sweep.x, s.sweep.y, s.sweep.z) * scale;
    auto field = profile(length, scale, mode);
    algo::GeometryMeshingOptions o;
    o.targetSize = s.target * scale;
    o.threads = 1;
    o.warningElementQuality = .15;
    algo::GeometryBodyControl c;
    c.solidTag = 1;
    c.sourceFaceTag = mode == 1 ? 2 : 1;
    c.sizeProfile = &field.profile;
    c.elements = tets ? algo::GeometryElements::Tet4 : algo::GeometryElements::Native;
    o.bodyControls = {c};
    mesh::Mesh mesh;
    const auto r = algo::meshGeometryModel(*built.model, o, mesh);
    std::string failure = algo::geometryReasonName(r.reason);
    for (const auto &b : r.analysis.bodies)
        failure += " engine=" + std::to_string(b.engineStatus) + "/" +
                   std::to_string(b.engineReason) +
                   " conversion=" + algo::convexToTetReasonName(b.conversion.reason);
    require(r.published, s.name + " " + failure);
    require(r.audit.valid && r.qualityAssessed && r.quality.valid, "missing validation");
    const auto &b = r.analysis.bodies.front();
    require(b.usedSizeProfile, "profile ignored");
    require(!tets || b.conversion.published, "Tet4 conversion ignored");
    const double expected = complex_test::area(s) * s.sweep.z * scale * scale * scale;
    const double volume = geometry_test::positiveVolume(mesh);
    require(std::abs(volume - expected) <= 1e-8 * expected,
            "volume differs from polygon area times depth");
    geometry_test::checkClassification(*built.model, mesh, scale * 1e-7);
    boundaryEuler(mesh, s.loops.size() - 1);
    stations(mesh, *built.model, b, field.profile);
    require(mesh::referenceLinearVolumeBoundary(mesh).valid, "independent topology audit");
    for (Index i = 0; i < mesh.cellCount(); ++i)
        if (mesh.cell(i).onEntity.dim == 3 && tets)
            require(mesh.cell(i).type == mesh::CellType::Tet4, "wrong volume family");
    const auto hash = digest(mesh);
    if (compare) {
        o.threads = 4;
        mesh::Mesh other;
        require(algo::meshGeometryModel(*built.model, o, other).published && digest(other) == hash,
                "1/4-thread identity");
    }
    const std::string name = s.name + "-" +
                             (mode == 0   ? "start"
                              : mode == 1 ? "end"
                              : mode == 2 ? "band"
                                          : "both") +
                             (tets ? "-tet4" : "-native");
    if (!gallery.empty() && scale == 1 && !rotated) {
        std::filesystem::create_directories(gallery);
        require(sdk::writeMeshArchive(mesh, (gallery / (name + ".paacsmesh")).string()).completed,
                "archive export");
    }
    std::cout << std::setprecision(17) << "{\"case\":\"" << name << "\",\"scale\":" << scale
              << ",\"rotated\":" << (rotated ? "true" : "false")
              << ",\"nodes\":" << mesh.nodeCount() << ",\"volume_cells\":" << r.quality.volumeCells
              << ",\"layers\":" << b.layers << ",\"output_bytes\":" << r.outputBytes
              << ",\"conversion_working_bytes\":" << b.conversion.workingBytes
              << ",\"minimum_quality\":" << r.quality.minimumQuality
              << ",\"quality_warning\":" << (r.qualityWarning ? "true" : "false")
              << ",\"volume\":" << volume << ",\"generation_seconds\":" << r.generationSeconds
              << ",\"audit_seconds\":" << r.auditSeconds << ",\"wall_seconds\":"
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
              << ",\"hash\":\"" << hash << "\",\"passed\":true}\n"
              << std::flush;
    ++meshes;
}
void rejectedControls() {
    auto s = complex_test::shapes().back();
    auto built = geo::makePolyhedralModel(complex_test::input(s));
    auto f = profile(std::hypot(s.sweep.x, s.sweep.y, s.sweep.z), 1, 0);
    auto wrong = profile(99, 1, 0);
    sizing::CappedSizeProfile invalid;
    std::atomic<bool> cancel{true};
    algo::GeometryMeshingOptions base;
    algo::GeometryBodyControl control;
    control.solidTag = 1;
    control.sourceFaceTag = 1;
    control.sizeProfile = &f.profile;
    control.elements = algo::GeometryElements::Tet4;
    base.bodyControls = {control};
    const auto reject = [&](const auto &o, algo::GeometryReason reason) {
        mesh::Mesh held;
        held.addNode({9, 8, 7}, {0, 9});
        const auto before = digest(held);
        const auto r = algo::meshGeometryModel(*built.model, o, held);
        require(!r.published && r.destinationUnchanged && digest(held) == before &&
                    r.reason == reason,
                std::string("control refusal: ") + algo::geometryReasonName(r.reason));
        ++refusals;
    };
    using R = algo::GeometryReason;
    for (int i = 0; i < 18; ++i) {
        auto o = base;
        auto &c = o.bodyControls[0];
        R reason = R::InvalidRequest;
        switch (i) {
        case 0:
            c.sizeProfile = &invalid;
            break;
        case 1:
            c.sizeProfile = &wrong.profile;
            break;
        case 2:
            c.layers = 3;
            break;
        case 3:
            c.endToStartLayerWidthRatio = 2;
            break;
        case 4:
            c.elements = static_cast<algo::GeometryElements>(255);
            break;
        case 5:
            c.profileLimits.maximumAdjacentGrowth = NAN;
            break;
        case 6:
            c.profileLimits.maximumCells = 1;
            reason = R::ProfileRejected;
            break;
        case 7:
            c.profileLimits.maximumBytes = 1;
            reason = R::ProfileRejected;
            break;
        case 8:
            o.maximumProfileWorkingBytes = 1;
            reason = R::ProfileRejected;
            break;
        case 9:
            c.sourceFaceTag = 999;
            reason = R::NoMatchingStrategy;
            break;
        case 10:
            c.strategy = algo::GeometryStrategy::BoundaryConformingTet4;
            reason = R::NoMatchingStrategy;
            break;
        case 11:
            o.maximumConversionWorkingBytes = 1;
            reason = R::ConversionFailed;
            break;
        case 12:
            o.maximumOutputBytes = 1;
            reason = R::OutputBudget;
            break;
        case 13:
            o.cancellation = &cancel;
            reason = R::Cancelled;
            break;
        case 14:
            c.profileLimits.cancellation = &cancel;
            reason = R::Cancelled;
            break;
        case 15:
            o.minimumElementQuality = .999;
            reason = R::QualityRequirementNotMet;
            break;
        case 16:
            o.warningElementQuality = INFINITY;
            break;
        case 17:
            c.maximumAdjacentLayerRatio = NAN;
            break;
        }
        reject(o, reason);
    }
    mesh::Mesh warnings;
    auto warning = base;
    warning.warningElementQuality = .999;
    const auto w = algo::meshGeometryModel(*built.model, warning, warnings);
    require(w.published && w.qualityWarning && w.qualityAssessed, "warning must retain mesh");
}
void mappedFractions() {
    std::vector<double> f{0, .2, .45, .725, 1};
    algo::MappedHexBlockRequest r;
    r.geometry.cellsU = 2;
    r.geometry.cellsV = 3;
    r.geometry.cellsW = 4;
    r.growth.w.boundaryFractions = f;
    r.executionThreads = 1;
    algo::MappedHexWorkspace workspace;
    mesh::Mesh a, b;
    require(algo::generateMappedHexBlock(r, a).published(), "direct mapped fractions");
    require(algo::generateMappedHexBlock(r, workspace, b).published() && digest(a) == digest(b),
            "workspace mapped fractions");
    for (int i = 0; i < 6; ++i) {
        auto bad = f;
        switch (i) {
        case 0:
            bad[0] = .01;
            break;
        case 1:
            bad.back() = .99;
            break;
        case 2:
            bad[1] = 0;
            break;
        case 3:
            bad[2] = NAN;
            break;
        case 4:
            bad.pop_back();
            break;
        case 5:
            bad[1] = .01;
            break;
        }
        r.growth.w.boundaryFractions = bad;
        const auto before = digest(a);
        require(!algo::generateMappedHexBlock(r, a).published() && digest(a) == before,
                "bad mapped fractions published");
        ++refusals;
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path gallery;
        std::string selected;
        bool compare = true, small = false;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--gallery" && i + 1 < argc)
                gallery = argv[++i];
            else if (a == "--case" && i + 1 < argc)
                selected = argv[++i];
            else if (a == "--single")
                compare = false;
            else if (a == "--small")
                small = true;
            else
                throw std::runtime_error("invalid argument");
        }
        int found = 0;
        for (const auto &s : complex_test::shapes())
            if (selected.empty() || selected == s.name) {
                ++found;
                for (int mode = 0; mode < (small ? 1 : 4); ++mode)
                    for (bool tets : {false, true})
                        run(s, mode, tets, 1, false, gallery, compare);
                if (!small && (s.name == "multi_bore" || s.name == "mapped_block"))
                    for (double scale : {1e-6, 1e6})
                        for (bool tets : {false, true})
                            run(s, 2, tets, scale, true, {}, compare);
            }
        require(found > 0, "unknown case");
        rejectedControls();
        mappedFractions();
        std::cout << "complex geometry: PASS (" << meshes << " meshes, " << refusals
                  << " refusals; topology, volume, CAD tags, sizing, thread identity)\n";
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
