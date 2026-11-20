// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/Hex8Jacobian.hpp"
#include "primitive_checks.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <thread>
using namespace primitive_test;
namespace {
// Independent long-double evaluation of FEM shape-function derivatives on
// [-1,1]^3. This oracle does not form Bernstein coefficients.
long double jacobian(const std::array<Vec3, 8> &p, long double u, long double v, long double w) {
    constexpr int signs[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                                 {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    long double j[3][3]{};
    for (int i = 0; i < 8; ++i) {
        const auto &s = signs[i];
        const long double d[3] = {s[0] * (1 + s[1] * v) * (1 + s[2] * w) / 8,
                                  s[1] * (1 + s[0] * u) * (1 + s[2] * w) / 8,
                                  s[2] * (1 + s[0] * u) * (1 + s[1] * v) / 8};
        const long double q[3] = {static_cast<long double>(p[i].x) - p[0].x,
                                  static_cast<long double>(p[i].y) - p[0].y,
                                  static_cast<long double>(p[i].z) - p[0].z};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                j[a][b] += q[a] * d[b];
    }
    return j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) -
           j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0]) +
           j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
}
} // namespace
int main() {
    try {
        using R = algo::PrimitiveReason;
        using K = algo::PrimitiveKind;
        algo::PrimitiveShape base;
        algo::PrimitiveMeshingOptions options;
        options.targetSize = .3;
        int refusals = 0;
        auto refuse = [&](algo::PrimitiveShape s, algo::PrimitiveMeshingOptions o, R expected) {
            mesh::Mesh m;
            m.addNode({9, 8, 7}, {0, 123});
            const auto before = hash(m);
            const auto r = algo::meshPrimitive(s, o, m);
            require(!r.published && r.destinationUnchanged && hash(m) == before,
                    "failure changed destination");
            require(r.reason == expected,
                    std::string("wrong refusal: ") + algo::primitiveReasonName(r.reason));
            ++refusals;
        };
        const double nan = std::numeric_limits<double>::quiet_NaN(),
                     inf = std::numeric_limits<double>::infinity();
        for (double v : {0., -1., nan, inf}) {
            auto s = base;
            s.radius = v;
            refuse(s, options, R::InvalidRequest);
            s = base;
            s.height = v;
            refuse(s, options, R::InvalidRequest);
            auto o = options;
            o.targetSize = v;
            refuse(base, o, R::InvalidRequest);
        }
        {
            auto s = base;
            s.kind = static_cast<K>(255);
            refuse(s, options, R::InvalidRequest);
            s = base;
            s.axis = {0, 0, 0};
            refuse(s, options, R::InvalidRequest);
            s = base;
            s.xDirection = {0, 0, 2};
            refuse(s, options, R::InvalidRequest);
            s = base;
            s.origin.x = nan;
            refuse(s, options, R::InvalidRequest);
        }
        for (double v : {0., 1., 2., nan}) {
            auto s = base;
            s.kind = K::Tube;
            s.innerRadius = v;
            refuse(s, options, R::InvalidRequest);
        }
        for (double v : {0., 1., nan}) {
            auto s = base;
            s.kind = K::Torus;
            s.majorRadius = v;
            refuse(s, options, R::InvalidRequest);
        }
        for (double v : {0., 2 * std::numbers::pi, nan}) {
            auto s = base;
            s.kind = K::Elbow;
            s.sweepAngle = v;
            refuse(s, options, R::InvalidRequest);
        }
        for (int v : {-1, 0, 1025}) {
            auto o = options;
            o.threads = v;
            refuse(base, o, R::InvalidRequest);
        }
        {
            auto o = options;
            o.maximumChordError = -1;
            refuse(base, o, R::InvalidRequest);
            o = options;
            o.maximumRelativeVolumeError = 0;
            refuse(base, o, R::InvalidRequest);
            o = options;
            o.minimumQuality = 2;
            refuse(base, o, R::InvalidRequest);
        }
        mesh::Mesh reference;
        const auto accepted = algo::meshPrimitive(base, options, reference);
        require(accepted.published, "reference failed");
        auto exact = options;
        exact.maximumNodes = accepted.nodes;
        exact.maximumCells = accepted.hexes + accepted.pyramids + accepted.boundaryQuads +
                             accepted.boundaryTriangles;
        exact.maximumOutputBytes = accepted.outputBytes;
        exact.maximumWorkingBytes = accepted.plannedWorkingBytes;
        exact.auditLimits.maximumFaceUses = accepted.audit.faceUses;
        exact.auditLimits.maximumTemporaryBytes = accepted.audit.temporaryBytes;
        {
            mesh::Mesh m;
            require(algo::meshPrimitive(base, exact, m).published && hash(m) == hash(reference),
                    "exact resource limits");
        }
        {
            auto o = exact;
            --o.maximumNodes;
            refuse(base, o, R::OutputBudget);
            o = exact;
            --o.maximumCells;
            refuse(base, o, R::OutputBudget);
            o = exact;
            --o.maximumOutputBytes;
            refuse(base, o, R::OutputBudget);
            o = exact;
            --o.maximumWorkingBytes;
            refuse(base, o, R::WorkingBudget);
            o = exact;
            --o.auditLimits.maximumFaceUses;
            refuse(base, o, R::WorkingBudget);
            o = exact;
            --o.auditLimits.maximumTemporaryBytes;
            refuse(base, o, R::WorkingBudget);
        }
        {
            auto o = options;
            o.targetSize = 1e-100;
            refuse(base, o, R::OutputBudget);
            o = options;
            o.minimumQuality = .999;
            refuse(base, o, R::QualityRequirement);
            o = options;
            o.maximumRelativeVolumeError = 1e-6;
            refuse(base, o, R::WorkingBudget);
        }
        {
            auto o = options;
            o.maximumRelativeVolumeError = 1e-4;
            mesh::Mesh m;
            const auto r = algo::meshPrimitive(base, o, m);
            require(r.published && r.relativeVolumeError <= o.maximumRelativeVolumeError &&
                        r.segmentsPerQuarter > accepted.segmentsPerQuarter,
                    "automatic volume-error refinement");
            check(base, m, r);
        }
        {
            auto s = base;
            s.origin = {1e16, 1e16, 1e16};
            refuse(s, options, R::GeometryTolerance);
        }
        {
            auto s = base;
            s.kind = K::Cone;
            s.origin = {0, 0, 1e16};
            s.height = .1;
            auto o = options;
            o.targetSize = 10;
            o.maximumChordError = 1000;
            refuse(s, o, R::InvalidCell); // The emitted apex rounds onto its base plane.
        }
        std::atomic<bool> stop{true};
        {
            auto o = options;
            o.cancellation = &stop;
            refuse(base, o, R::Cancelled);
            o = options;
            o.auditLimits.cancellation = &stop;
            refuse(base, o, R::Cancelled);
        }
        stop = false;
        {
            auto o = options;
            o.targetSize = .025;
            o.cancellation = &stop;
            std::jthread cancel([&] {
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
                stop = true;
            });
            refuse(base, o, R::Cancelled);
        }

        std::mt19937_64 random(0x504141435334ULL);
        auto uniform = [&](double a, double b) {
            return a + (b - a) * static_cast<double>(random() >> 11) / 9007199254740992.;
        };
        int randomPassed = 0;
        for (int i = 0; i < 128; ++i) {
            auto s = base;
            s.kind = static_cast<K>(i % 6);
            const double scale = std::pow(10., uniform(-3, 3));
            s.radius = scale;
            s.height = scale * uniform(.15, 6);
            s.endRadius = scale * uniform(.03, 2);
            s.majorRadius = scale * uniform(1.06, 4);
            s.innerRadius =
                s.kind == K::Tube
                    ? scale * uniform(.05, .98)
                    : (i % 2 && (s.kind == K::Torus || s.kind == K::Elbow) ? scale * .75 : 0);
            s.sweepAngle = uniform(.1, 6);
            s.axis = {uniform(-1, 1), uniform(-1, 1), uniform(.5, 1)};
            s.xDirection = {1, .3, -.2};
            s.origin = {scale * 3, scale * (-2), scale * 5};
            auto o = options;
            o.targetSize = scale * uniform(.2, .5);
            mesh::Mesh m;
            const auto r = algo::meshPrimitive(s, o, m);
            require(r.published, "random primitive " + std::to_string(i) + ": " +
                                     algo::primitiveReasonName(r.reason));
            check(s, m, r);
            ++randomPassed;
        }
        int convergence = 0;
        for (auto kind : {K::Cylinder, K::Cone, K::Frustum, K::Tube, K::Torus, K::Elbow}) {
            auto s = base;
            s.kind = kind;
            if (kind == K::Tube)
                s.innerRadius = .8;
            double previousVolume = 1, previousSurface = 1;
            for (double error : {.01, .0025, .000625}) {
                auto o = options;
                o.maximumChordError = error;
                o.maximumRelativeVolumeError = .05;
                o.targetSize = .3;
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(s, o, m);
                require(r.published, "convergence failed");
                const auto proof = check(s, m, r);
                require(r.relativeVolumeError < previousVolume * .6 &&
                            proof.maximumSurfaceError < previousSurface * .6,
                        "curved approximation did not converge");
                previousVolume = r.relativeVolumeError;
                previousSurface = proof.maximumSurfaceError;
                ++convergence;
                std::cout << "{\"convergence\":\"" << algo::primitiveKindName(kind)
                          << "\",\"chord_tolerance\":" << error
                          << ",\"volume_error\":" << r.relativeVolumeError
                          << ",\"surface_error\":" << proof.maximumSurfaceError
                          << ",\"cells\":" << r.hexes + r.pyramids << "}\n";
            }
        }
        const std::array<Vec3, 8> cube{{{0, 0, 0},
                                        {1, 0, 0},
                                        {1, 1, 0},
                                        {0, 1, 0},
                                        {0, 0, 1},
                                        {1, 0, 1},
                                        {1, 1, 1},
                                        {0, 1, 1}}};
        require(mesh::positiveHex8Jacobian(cube), "cube certificate");
        const std::array<Vec3, 8> folded{
            {{-.35087608027880646, .0261143510490367, .9089395824796755},
             {1.9835503585173442, -.7495904383613432, -1.0492740532258584},
             {1.7640165559082126, .7175817932028974, 1.0994783612164298},
             {-.6610601476918484, .36586759689771375, .9013842341629299},
             {-.3684522952035696, -.094571399410885, 1.6312070165615216},
             {.5099029164097465, .2867737378315247, 1.9578373139290217},
             {1.3977705922623176, .8106024191116803, 1.3725956320471633},
             {-.8408219437092093, 1.141783702146146, 1.0110664005949017}}};
        for (int w : {-1, 1})
            for (int v : {-1, 1})
                for (int u : {-1, 1})
                    require(jacobian(folded, u, v, w) > 0,
                            "folded fixture must fool corner-only checks");
        require(jacobian(folded, 0, 1, -1) < 0 && !mesh::positiveHex8Jacobian(folded),
                "fold between positive corners accepted");
        for (double scale : {1e-90, 1e90}) {
            auto p = cube;
            for (auto &v : p)
                v = v * scale;
            require(mesh::positiveHex8Jacobian(p), "scale invariant certificate");
        }
        {
            auto p = cube;
            for (auto &v : p)
                v.z = 0;
            require(!mesh::positiveHex8Jacobian(p), "collapsed hex accepted");
            p = cube;
            p[0].x = nan;
            require(!mesh::positiveHex8Jacobian(p), "nonfinite hex accepted");
        }
        int certified = 0, rejected = 0;
        for (int trial = 0; trial < 2000; ++trial) {
            auto p = cube;
            const double amplitude = trial % 2 ? .15 : 1.2;
            for (auto &v : p)
                v = v + Vec3{uniform(-amplitude, amplitude), uniform(-amplitude, amplitude),
                             uniform(-amplitude, amplitude)};
            const bool passed = mesh::positiveHex8Jacobian(p);
            if (passed) {
                ++certified;
                for (int k = 0; k <= 6; ++k)
                    for (int j = 0; j <= 6; ++j)
                        for (int i = 0; i <= 6; ++i)
                            require(jacobian(p, -1 + i / 3.L, -1 + j / 3.L, -1 + k / 3.L) > 0,
                                    "Bernstein certificate false positive");
            } else
                ++rejected;
        }
        require(certified >= 950 && rejected >= 500, "Jacobian adversaries were ineffective");
        std::cout << "primitive hardening PASS: " << refusals << " refusals, 128 seeded shapes, "
                  << convergence << " convergence meshes, 2000 Jacobian adversaries (" << certified
                  << " certified, " << rejected << " unresolved/rejected)\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
