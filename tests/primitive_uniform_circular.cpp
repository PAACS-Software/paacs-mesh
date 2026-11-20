// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "primitive_size_checks.hpp"
#include <iostream>
#include <random>

using namespace primitive_test;
int main() {
    try {
        using K = algo::PrimitiveKind;
        using E = algo::PrimitiveElements;
        using S = algo::PrimitiveCircularStrategy;
        using R = algo::PrimitiveReason;
        unsigned configurations = 0, refusals = 0;
        auto reject = [&](const auto &shape, const auto &options, R reason) {
            mesh::Mesh previous;
            previous.addNode({9, 8, 7}, {0, 9});
            const auto before = hash(previous);
            const auto r = algo::meshPrimitive(shape, options, previous);
            require(!r.published && r.destinationUnchanged && r.reason == reason &&
                        hash(previous) == before,
                    std::string("circular rollback: ") + algo::primitiveReasonName(r.reason));
            ++refusals;
        };
        std::mt19937_64 random(0x637572766564ULL);
        std::uniform_real_distribution<double> u(0, 1);
        for (auto family : {E::Tet4, E::Tet10})
            for (auto kind : {K::Torus, K::Elbow}) {
                for (unsigned i = 0; i < 8; ++i) {
                    algo::PrimitiveShape s;
                    s.kind = kind;
                    s.radius = std::pow(10., -6. + 12. * i / 7.);
                    s.majorRadius = (2 + 2 * u(random)) * s.radius;
                    s.sweepAngle = .4 + 5.7 * u(random);
                    s.axis = {1, 2, 3};
                    s.xDirection = {3, -1, .2};
                    s.origin = Vec3{17, -21, 31} * s.radius;
                    algo::PrimitiveMeshingOptions o;
                    o.elements = family;
                    o.targetSize = (.2 + .15 * u(random)) * s.radius;
                    o.maximumChordError = .01 * s.radius;
                    o.maximumRelativeVolumeError = .02;
                    mesh::Mesh m;
                    const auto r = algo::meshPrimitive(s, o, m);
                    require(r.published,
                            std::string("circular: ") + algo::primitiveReasonName(r.reason));
                    require(r.uniformCircular && !r.uniformCylinder && r.audit.valid,
                            "route/audit");
                    check(s, m, r);
                    const auto stats = sizes(m);
                    require(stats.minimumQuality > .4 && stats.cv < .32,
                            "balanced circular quality");
                    require(std::abs(stats.minimumQuality - r.minimumQuality) < 1e-8,
                            "independent quality");
                    require(r.circularRadialSpacing <= o.targetSize && r.circularInnerSpacing > 0 &&
                                r.circularOuterSpacing > r.circularInnerSpacing,
                            "reported spacings");
                    if (family == E::Tet10)
                        require(r.certifiedTet10s == r.tet10s && r.quadraticConnectivityVerified,
                                "quadratic certification");
                    o.circularStrategy = S::Uniform;
                    o.threads = 4;
                    mesh::Mesh parallel;
                    require(algo::meshPrimitive(s, o, parallel).published &&
                                hash(m) == hash(parallel),
                            "automatic/explicit and cross-thread hash");
                    ++configurations;
                }
                algo::PrimitiveShape s;
                s.kind = kind;
                algo::PrimitiveMeshingOptions o;
                o.elements = family;
                o.targetSize = 10;
                o.maximumChordError = .3;
                o.maximumRelativeVolumeError = .5;
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(s, o, m);
                require(r.published && r.radialLayers == 2, "coarsest circular section");
                check(s, m, r);
                auto exact = o;
                exact.maximumNodes = r.nodes;
                exact.maximumCells = m.cellCount();
                exact.maximumOutputBytes = r.outputBytes;
                exact.maximumWorkingBytes = r.plannedWorkingBytes;
                exact.auditLimits.maximumFaceUses = r.audit.faceUses;
                exact.auditLimits.maximumTemporaryBytes = r.audit.temporaryBytes;
                mesh::Mesh repeated;
                require(algo::meshPrimitive(s, exact, repeated).published &&
                            hash(m) == hash(repeated),
                        "exact budgets");
                auto bad = exact;
                --bad.maximumNodes;
                reject(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumCells;
                reject(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumOutputBytes;
                reject(s, bad, R::OutputBudget);
                bad = exact;
                --bad.maximumWorkingBytes;
                reject(s, bad, R::WorkingBudget);
                bad = exact;
                --bad.auditLimits.maximumFaceUses;
                reject(s, bad, R::WorkingBudget);
                bad = exact;
                --bad.auditLimits.maximumTemporaryBytes;
                reject(s, bad, R::WorkingBudget);
                bad = o;
                bad.targetSize = 1e-300;
                reject(s, bad, R::OutputBudget);
                bad = o;
                bad.maximumChordError = 1e-30;
                reject(s, bad, R::NumericallyUnresolved);
                std::atomic<bool> stop{true};
                bad = o;
                bad.cancellation = &stop;
                reject(s, bad, R::Cancelled);
                bad = o;
                bad.auditLimits.cancellation = &stop;
                reject(s, bad, R::Cancelled);
                o = {};
                o.elements = family;
                s.majorRadius = 1.0001;
                const auto tight = algo::meshPrimitive(s, o, m);
                require(tight.published && tight.qualityWarning &&
                            tight.quality.lowQualityCells > 0,
                        "tight-bend warning");
                check(s, m, tight);
                bad = o;
                bad.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
                reject(s, bad, R::QualityRequirement);
                bad = o;
                bad.minimumQuality = .5;
                reject(s, bad, R::QualityRequirement);
                ++configurations;
                if (kind == K::Elbow) {
                    s.majorRadius = 2;
                    s.sweepAngle = 1e-4;
                    const auto shortBend = algo::meshPrimitive(s, o, m);
                    require(shortBend.published && shortBend.qualityWarning, "short-bend warning");
                    check(s, m, shortBend);
                    ++configurations;
                }
            }
        for (auto kind : {K::Torus, K::Elbow}) {
            algo::PrimitiveShape s;
            s.kind = kind;
            algo::PrimitiveMeshingOptions o;
            o.elements = E::Tet4;
            mesh::Mesh balanced, old;
            const auto now = algo::meshPrimitive(s, o, balanced);
            o.circularStrategy = S::Swept;
            const auto before = algo::meshPrimitive(s, o, old);
            require(now.published && before.published && !before.uniformCircular, "legacy route");
            const auto a = sizes(balanced), b = sizes(old);
            require(a.cv < .8 * b.cv && a.minimumQuality > b.minimumQuality,
                    "measured quality improvement");
            require(balanced.countCells(mesh::CellType::Tet4) < old.countCells(mesh::CellType::Tet4),
                    "balanced default example should use fewer cells than legacy sweep");
            std::cout << "{\"shape\":\"" << algo::primitiveKindName(kind)
                      << "\",\"old_cv\":" << b.cv << ",\"new_cv\":" << a.cv
                      << ",\"old_quality\":" << b.minimumQuality
                      << ",\"new_quality\":" << a.minimumQuality << "}\n";
            o.circularStrategy = S::Uniform;
            s.innerRadius = .5;
            reject(s, o, R::InvalidRequest);
            o.circularStrategy = S::Automatic;
            require(algo::meshPrimitive(s, o, old).published, "hollow automatic compatibility");
            s.innerRadius = 0;
            o.circularStrategy = S::Uniform;
            o.elements = E::HexDominant;
            reject(s, o, R::InvalidRequest);
        }
        algo::PrimitiveShape s;
        algo::PrimitiveMeshingOptions o;
        o.elements = E::Tet4;
        o.circularStrategy = S::Uniform;
        reject(s, o, R::InvalidRequest);
        o.circularStrategy = static_cast<S>(255);
        reject(s, o, R::InvalidRequest);
        std::cout << "{\"configurations\":" << configurations << ",\"rollback_checks\":" << refusals
                  << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
