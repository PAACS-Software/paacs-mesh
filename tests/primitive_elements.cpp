// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Tet10FastCertificate.hpp"
#include "paacs/algo/Tet10JacobianCertificate.hpp"
#include "primitive_checks.hpp"
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <thread>
using namespace primitive_test;
namespace {
// Independent long-double FEM derivative evaluation at reference coordinates.
long double jacobian(const std::array<Vec3, 10> &p, long double r, long double s, long double t) {
    const long double l[4]{1 - r - s - t, r, s, t};
    const int dl[4][3]{{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const int edges[6][2]{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    long double d[10][3]{}, j[3][3]{};
    for (unsigned a = 0; a < 4; ++a)
        for (unsigned b = 0; b < 3; ++b)
            d[a][b] = (4 * l[a] - 1) * dl[a][b];
    for (unsigned e = 0; e < 6; ++e)
        for (unsigned b = 0; b < 3; ++b)
            d[e + 4][b] =
                4 * (l[edges[e][0]] * dl[edges[e][1]][b] + l[edges[e][1]] * dl[edges[e][0]][b]);
    for (unsigned a = 0; a < 10; ++a) {
        const long double x[3]{static_cast<long double>(p[a].x) - p[0].x,
                               static_cast<long double>(p[a].y) - p[0].y,
                               static_cast<long double>(p[a].z) - p[0].z};
        for (unsigned b = 0; b < 3; ++b)
            for (unsigned c = 0; c < 3; ++c)
                j[b][c] += x[b] * d[a][c];
    }
    return j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) -
           j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0]) +
           j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
}
void independentTopology(const mesh::Mesh &m, const algo::PrimitiveMeshingResult &r) {
    mesh::Mesh corners;
    for (Index i = 0; i < m.nodeCount(); ++i)
        corners.addNode(m.node(i), m.nodeClass(i));
    std::map<std::pair<Index, Index>, Index> midpoint;
    std::map<Index, std::pair<Index, Index>> reverse;
    std::set<Index> cornerIds;
    std::uint64_t tets = 0;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto type = m.cell(i).type;
        require(type == mesh::CellType::Tet4 || type == mesh::CellType::Tet10 ||
                    type == mesh::CellType::Tri3 || type == mesh::CellType::Tri6,
                "mixed family leak");
        const bool volume = type == mesh::CellType::Tet4 || type == mesh::CellType::Tet10;
        const bool high = type == mesh::CellType::Tet10 || type == mesh::CellType::Tri6;
        const unsigned n = volume ? 4 : 3;
        const auto *v = m.cellNodes(i);
        for (unsigned k = 0; k < n; ++k)
            cornerIds.insert(v[k]);
        corners.addCell(volume ? mesh::CellType::Tet4 : mesh::CellType::Tri3,
                        std::span<const Index>(v, n), m.cell(i).onEntity);
        if (volume)
            ++tets;
        if (high) {
            const auto &tr =
                topo::traits(volume ? topo::ElementType::Tet10 : topo::ElementType::Tri6);
            for (unsigned e = 0; e < tr.edgeCount; ++e) {
                const std::pair<Index, Index> key = std::minmax(v[tr.edges[e].a], v[tr.edges[e].b]);
                const Index id = v[n + e];
                const auto [it, inserted] = midpoint.emplace(key, id);
                require(inserted || it->second == id, "cracked quadratic edge");
                const auto [rev, added] = reverse.emplace(id, key);
                require(added || rev->second == key, "midside alias between edges");
                const auto a = m.node(key.first), b = m.node(key.second), p = m.node(id);
                require(p.x == std::midpoint(a.x, b.x) && p.y == std::midpoint(a.y, b.y) &&
                            p.z == std::midpoint(a.z, b.z),
                        "straight midside contract");
            }
        }
    }
    for (const auto &[id, key] : reverse) {
        (void)key;
        require(!cornerIds.contains(id), "midside aliases corner");
    }
    require(mesh::auditLinearVolumeBoundary(corners).valid, "independent emitted corner audit");
    require(tets == r.tet4s + r.tet10s, "tet counts");
    require(std::accumulate(r.quality.histogram.begin(), r.quality.histogram.end(),
                            std::uint64_t(0)) == tets,
            "quality histogram incomplete");
    if (r.tet10s)
        require(r.certifiedTet10s == tets && r.quadraticConnectivityVerified,
                "uncertified Tet10 publication");
}
} // namespace
int main() {
    try {
        using E = algo::PrimitiveElements;
        using K = algo::PrimitiveKind;
        using R = algo::PrimitiveReason;
        unsigned configurations = 0, refusals = 0;
        for (auto kind : {K::Cylinder, K::Cone, K::Frustum, K::Tube, K::Torus, K::Elbow})
            for (auto family : {E::Tet4, E::Tet10}) {
                algo::PrimitiveShape s;
                s.kind = kind;
                if (kind == K::Tube)
                    s.innerRadius = .8;
                algo::PrimitiveMeshingOptions o;
                o.elements = family;
                o.targetSize = .4;
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(s, o, m);
                require(r.published, std::string("family: ") + algo::primitiveKindName(kind) + " " +
                                         algo::primitiveReasonName(r.reason));
                independentTopology(m, r);
                ++configurations;
            }
        auto refuse = [&](algo::PrimitiveShape s, algo::PrimitiveMeshingOptions o, R reason) {
            mesh::Mesh m;
            m.addNode({1, 2, 3}, {0, 77});
            const auto before = hash(m);
            const auto r = algo::meshPrimitive(s, o, m);
            require(!r.published && r.destinationUnchanged && hash(m) == before,
                    "failure changed destination");
            require(r.reason == reason,
                    std::string("unexpected reason: ") + algo::primitiveReasonName(r.reason));
            ++refusals;
        };
        algo::PrimitiveShape s;
        algo::PrimitiveMeshingOptions o;
        o.targetSize = .4;
        for (auto e : {E::HexDominant, E::Tet4, E::Tet10}) {
            o.elements = e;
            o.warningQuality = .99;
            mesh::Mesh m;
            const auto warning = algo::meshPrimitive(s, o, m);
            require(warning.published && warning.qualityWarning &&
                        warning.quality.lowQualityCells > 32 && warning.quality.examplesTruncated,
                    "warning path");
            require(warning.quality.sharpFeatureLowQualityCells > 0, "sharp rim diagnostics");
            for (unsigned i = 0; i < warning.quality.exampleCount; ++i)
                require(warning.quality.examples[i].cell == static_cast<Index>(i),
                        "ordered quality examples");
            auto reject = o;
            reject.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
            refuse(s, reject, R::QualityRequirement);
            reject = o;
            reject.minimumQuality = .99;
            refuse(s, reject, R::QualityRequirement);
            o.warningQuality = .1;
        }
        for (auto e : {E::Tet4, E::Tet10})
            for (bool cone : {false, true}) {
                s.kind = cone ? K::Cone : K::Cylinder;
                o = {};
                o.elements = e;
                o.targetSize = .4;
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(s, o, m);
                require(r.published, "budget baseline");
                auto exact = o;
                exact.maximumNodes = r.nodes;
                exact.maximumCells = r.tet4s + r.tet10s + r.boundaryTriangles;
                exact.maximumOutputBytes = r.outputBytes;
                exact.maximumWorkingBytes = r.plannedWorkingBytes;
                // Conversion checks both its source and destination audits; final face
                // budget is larger. Temporary-byte budget remains the caller's bound.
                exact.auditLimits.maximumFaceUses = r.audit.faceUses;
                mesh::Mesh repeat;
                require(algo::meshPrimitive(s, exact, repeat).published && hash(repeat) == hash(m),
                        "exact budget success");
                auto bad = exact;
                --bad.maximumNodes;
                refuse(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumCells;
                refuse(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumOutputBytes;
                refuse(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumWorkingBytes;
                refuse(s, bad, R::WorkingBudget);
                std::atomic<bool> stop{true};
                bad = exact;
                bad.cancellation = &stop;
                refuse(s, bad, R::Cancelled);
            }
        s = {};
        o = {};
        o.elements = static_cast<E>(255);
        refuse(s, o, R::InvalidRequest);
        o = {};
        o.warningQuality = std::numeric_limits<double>::quiet_NaN();
        refuse(s, o, R::InvalidRequest);
        o = {};
        o.coneStrategy = algo::PrimitiveConeStrategy::Graded;
        refuse(s, o, R::InvalidRequest);
        o = {};
        o.qualityPolicy = static_cast<algo::PrimitiveQualityPolicy>(255);
        refuse(s, o, R::InvalidRequest);
        // A coarse cone may need no refinement at all. Preflight must still
        // run for the initial mesh (not only inside the refinement loop).
        for (auto family : {E::Tet4, E::Tet10}) {
            s = {};
            s.kind = K::Cone;
            o = {};
            o.elements = family;
            o.targetSize = 10;
            o.maximumChordError = .3;
            o.maximumRelativeVolumeError = .5;
            mesh::Mesh m;
            const auto coarse = algo::meshPrimitive(s, o, m);
            require(coarse.published && coarse.coneRefinementLevels == 0,
                    "unrefined cone baseline");
            auto bad = o;
            bad.maximumNodes = coarse.nodes - 1;
            refuse(s, bad, R::OutputBudget);
            bad = o;
            bad.maximumCells = coarse.tet4s + coarse.tet10s + coarse.boundaryTriangles - 1;
            refuse(s, bad, R::OutputBudget);
            bad = o;
            bad.maximumOutputBytes = coarse.outputBytes - 1;
            refuse(s, bad, R::OutputBudget);
            bad = o;
            bad.maximumWorkingBytes = coarse.plannedWorkingBytes - 1;
            refuse(s, bad, R::WorkingBudget);
        }
        // Exercise both new routes across scale, frame, taper and wall thickness.
        std::mt19937_64 shapes(0x26cad5);
        std::uniform_real_distribution<double> unit(0, 1);
        for (unsigned i = 0; i < 48; ++i) {
            s = {};
            s.kind = static_cast<K>(i % 6);
            const double scale = std::pow(10., -3 + 6 * unit(shapes));
            s.radius = scale;
            s.height = scale * (.4 + 3 * unit(shapes));
            s.endRadius = scale * (.05 + 1.5 * unit(shapes));
            s.majorRadius = scale * (1.1 + 2 * unit(shapes));
            if (s.kind == K::Tube || ((s.kind == K::Torus || s.kind == K::Elbow) && i % 2))
                s.innerRadius = scale * (.4 + .59 * unit(shapes));
            s.sweepAngle = .2 + 5.8 * unit(shapes);
            s.axis = {.2 + unit(shapes), .1 + unit(shapes), 1};
            s.origin = {10 * scale, -5 * scale, 2 * scale};
            o = {};
            o.elements = (i / 6) % 2 ? E::Tet10 : E::Tet4;
            o.targetSize = .6 * scale;
            o.maximumChordError = .01 * scale;
            mesh::Mesh m;
            const auto r = algo::meshPrimitive(s, o, m);
            require(r.published, std::string("random family ") + std::to_string(i) + ": " +
                                     algo::primitiveReasonName(r.reason));
            (void)check(s, m, r);
        }
        for (auto family : {E::Tet4, E::Tet10}) {
            s = {};
            o = {};
            o.elements = family;
            o.targetSize = .05;
            std::atomic<bool> stopped{false};
            o.cancellation = &stopped;
            mesh::Mesh m;
            m.addNode({9, 8, 7});
            const auto before = hash(m);
            std::jthread trigger([&] {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                stopped = true;
            });
            const auto result = algo::meshPrimitive(s, o, m);
            require(!result.published && result.reason == R::Cancelled && hash(m) == before,
                    "live tetrahedral cancellation");
        }
        s = {};
        // Stable cone-tip quality under refinement, and a fair same-family comparison.
        s.kind = K::Cone;
        o = {};
        o.elements = E::Tet4;
        double smallest = 1;
        std::uint64_t first = 0, last = 0;
        for (double h : {.4, .2, .1}) {
            o.targetSize = h;
            mesh::Mesh m;
            const auto r = algo::meshPrimitive(s, o, m);
            require(r.published, "graded cone refinement");
            require(r.minimumQuality > .2, "cone-tip quality regression");
            smallest = std::min(smallest, r.minimumQuality);
            if (!first)
                first = r.tet4s;
            last = r.tet4s;
        }
        o.targetSize = .1;
        o.coneStrategy = algo::PrimitiveConeStrategy::Swept;
        mesh::Mesh swept;
        const auto sr = algo::meshPrimitive(s, o, swept);
        require(sr.published && smallest > 5 * sr.minimumQuality,
                "graded tip did not improve quality");
        // Adversarial certificate oracle; sampling only tests the sufficient interval
        // certificate, it is never used by production to accept an element.
        std::mt19937_64 random(0x1bc526);
        std::uniform_real_distribution<double> move(-.22, .22);
        unsigned accepted = 0, rejected = 0;
        const auto &edges = topo::traits(topo::ElementType::Tet10);
        for (unsigned trial = 0; trial < 2000; ++trial) {
            std::array<Vec3, 10> p{};
            p[0] = {0, 0, 0};
            p[1] = {1, 0, 0};
            p[2] = {0, 1, 0};
            p[3] = {0, 0, 1};
            for (unsigned i = 0; i < 6; ++i) {
                p[i + 4] = (p[edges.edges[i].a] + p[edges.edges[i].b]) * .5;
                if (trial % 4)
                    p[i + 4] = p[i + 4] + Vec3{move(random), move(random), move(random)} *
                                              (trial % 4 == 1   ? .01
                                               : trial % 4 == 2 ? .1
                                                                : 1);
            }
            if (algo::positiveTet10JacobianFast(p)) {
                ++accepted;
                for (unsigned i = 0; i <= 6; ++i)
                    for (unsigned j = 0; j <= 6 - i; ++j)
                        for (unsigned k = 0; k <= 6 - i - j; ++k)
                            require(jacobian(p, i / 6.L, j / 6.L, k / 6.L) > 0,
                                    "interval certificate false acceptance");
                if (trial < 80)
                    require(algo::certifyTet10Jacobian(p, {}).status ==
                                algo::Tet10CertificateStatus::Certified,
                            "Bernstein cross-check");
            } else
                ++rejected;
        }
        require(accepted > 1000 && rejected > 400, "certificate adversaries ineffective");
        std::array<Vec3, 10> folded{};
        folded[1] = {1, 0, 0};
        folded[2] = {0, 1, 0};
        folded[3] = {0, 0, 1};
        for (unsigned e = 0; e < 6; ++e)
            folded[e + 4] = (folded[edges.edges[e].a] + folded[edges.edges[e].b]) * .5;
        for (double scale : {1e-90, 1., 1e90}) {
            auto scaled = folded;
            for (auto &p : scaled)
                p = p * scale;
            require(algo::positiveTet10JacobianFast(scaled), "certificate scaling");
        }
        folded[4] = {0, 0, 0};
        require(jacobian(folded, 0, 0, 0) < 0 && !algo::positiveTet10JacobianFast(folded),
                "positive-corner folded Tet10 accepted");
        folded[4].x = std::numeric_limits<double>::quiet_NaN();
        require(!algo::positiveTet10JacobianFast(folded), "nonfinite Tet10 accepted");
        std::cout << "{\"configurations\":" << configurations << ",\"rollback_cases\":" << refusals
                  << ",\"random_shapes\":48,\"certificate_trials\":2000,\"fast_accepted\":"
                  << accepted << ",\"unresolved\":" << rejected
                  << ",\"cone_minimum_quality\":" << smallest
                  << ",\"swept_cone_quality\":" << sr.minimumQuality
                  << ",\"coarse_cone_tets\":" << first << ",\"fine_cone_tets\":" << last << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
