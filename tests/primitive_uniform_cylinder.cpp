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
        using E = algo::PrimitiveElements;
        using C = algo::PrimitiveCylinderStrategy;
        using R = algo::PrimitiveReason;
        unsigned configurations = 0, refusals = 0;
        auto reject = [&](const algo::PrimitiveShape &s, const algo::PrimitiveMeshingOptions &o,
                          R reason) {
            mesh::Mesh m;
            m.addNode({4, 5, 6}, {0, 77});
            const auto old = hash(m);
            const auto r = algo::meshPrimitive(s, o, m);
            require(!r.published && r.destinationUnchanged && hash(m) == old,
                    "transaction rollback");
            require(r.reason == reason,
                    std::string("wrong refusal: ") + algo::primitiveReasonName(r.reason));
            ++refusals;
        };
        std::mt19937_64 rng(0x26c71);
        std::uniform_real_distribution<double> unit(0, 1);
        for (auto family : {E::Tet4, E::Tet10}) {
            for (unsigned test = 0; test < 24; ++test) {
                algo::PrimitiveShape s;
                const double scale = std::pow(10., -6 + 12 * unit(rng));
                s.radius = scale;
                s.height = scale * (.7 + 3 * unit(rng));
                s.axis = {unit(rng) - .5, unit(rng) - .5, 1};
                s.xDirection = {1, .2, .1};
                s.origin = {17 * scale, -21 * scale, 31 * scale};
                algo::PrimitiveMeshingOptions o;
                o.elements = family;
                o.cylinderStrategy = C::Uniform;
                o.targetSize = scale * (.14 + .8 * unit(rng));
                if (test % 3 == 0) {
                    o.maximumChordError = .04 * scale;
                    o.maximumRelativeVolumeError = .04;
                }
                mesh::Mesh m;
                const auto r = algo::meshPrimitive(s, o, m);
                require(r.published,
                        std::string("uniform cylinder: ") + algo::primitiveReasonName(r.reason));
                require(r.uniformCylinder && r.audit.valid && r.hexes == 0 && r.pyramids == 0,
                        "route/family");
                check(s, m, r);
                const auto z = sizes(m);
                require(z.cv < .10 && z.maximumRatio < 1.311 && z.minimumQuality > .60,
                        "size/shape quality regression");
                require(std::abs(z.minimumQuality - r.minimumQuality) < 1e-8,
                        "independent quality");
                require(r.cylinderRadialSpacing <= o.targetSize &&
                            r.cylinderAxialSpacing <= r.cylinderRadialSpacing,
                        "nominal spacing");
                if (family == E::Tet10)
                    require(r.certifiedTet10s == r.tet10s && r.quadraticConnectivityVerified,
                            "Tet10 certification");
                auto parallel = o;
                parallel.threads = 4;
                mesh::Mesh repeated;
                require(algo::meshPrimitive(s, parallel, repeated).published &&
                            hash(m) == hash(repeated),
                        "cross-thread hash");
                ++configurations;
            }
            // The smallest permitted section must obey every budget too.
            algo::PrimitiveShape s;
            s.height = .7;
            algo::PrimitiveMeshingOptions o;
            o.elements = family;
            o.cylinderStrategy = C::Uniform;
            o.targetSize = 10;
            o.maximumChordError = .3;
            o.maximumRelativeVolumeError = .5;
            mesh::Mesh m;
            const auto r = algo::meshPrimitive(s, o, m);
            require(r.published && r.radialLayers == 2, "coarse baseline");
            auto exact = o;
            exact.maximumNodes = r.nodes;
            exact.maximumCells = m.cellCount();
            exact.maximumOutputBytes = r.outputBytes;
            exact.maximumWorkingBytes = r.plannedWorkingBytes;
            exact.auditLimits.maximumFaceUses = r.audit.faceUses;
            mesh::Mesh repeated;
            require(algo::meshPrimitive(s, exact, repeated).published && hash(m) == hash(repeated),
                    "exact resource limits");
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
            bad.auditLimits.maximumTemporaryBytes = 1;
            reject(s, bad, R::WorkingBudget);
            bad = o;
            bad.targetSize = 1e-300;
            reject(s, bad, R::OutputBudget);
            bad = o;
            bad.maximumChordError = 1e-30;
            reject(s, bad, R::NumericallyUnresolved);
            bad = o;
            std::atomic<bool> stop{true};
            bad.cancellation = &stop;
            reject(s, bad, R::Cancelled);
            bad = o;
            bad.auditLimits.cancellation = &stop;
            reject(s, bad, R::Cancelled);
            // A thin disk cannot be made isotropic just by changing its axial count.
            s.height = .001;
            o = {};
            o.elements = family;
            o.cylinderStrategy = C::Uniform;
            const auto thin = algo::meshPrimitive(s, o, m);
            require(thin.published && thin.qualityWarning && thin.quality.lowQualityCells > 0,
                    "thin-disk warning");
            bad = o;
            bad.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
            reject(s, bad, R::QualityRequirement);
            bad = o;
            bad.minimumQuality = .5;
            reject(s, bad, R::QualityRequirement);
            ++configurations;
        }
        algo::PrimitiveShape s;
        algo::PrimitiveMeshingOptions o;
        o.elements = E::Tet4;
        mesh::Mesh automatic, swept, uniform;
        const auto a = algo::meshPrimitive(s, o, automatic);
        o.cylinderStrategy = C::Uniform;
        const auto u = algo::meshPrimitive(s, o, uniform);
        require(a.published && u.published && hash(automatic) == hash(uniform),
                "automatic selection");
        o.cylinderStrategy = C::Swept;
        const auto old = algo::meshPrimitive(s, o, swept);
        require(old.published && !old.uniformCylinder, "legacy swept route");
        const auto before = sizes(swept), after = sizes(uniform);
        require(after.cv < before.cv * .25 && after.minimumQuality > before.minimumQuality,
                "comparison regression");
        o.cylinderStrategy = static_cast<C>(255);
        reject(s, o, R::InvalidRequest);
        o.cylinderStrategy = C::Uniform;
        o.elements = E::HexDominant;
        reject(s, o, R::InvalidRequest);
        o.elements = E::Tet4;
        s.kind = algo::PrimitiveKind::Cone;
        reject(s, o, R::InvalidRequest);
        s = {};
        o = {};
        o.elements = E::Tet4;
        o.maximumChordError = .01;
        o.maximumRelativeVolumeError = .1;
        o.targetSize = .06;
        const auto fine = algo::meshPrimitive(s, o, uniform);
        require(fine.published && !fine.cylinderGeometryLimited, "size-driven sizing");
        o.targetSize = .5;
        const auto coarse = algo::meshPrimitive(s, o, uniform);
        require(coarse.published && coarse.cylinderGeometryLimited, "geometry-driven sizing");
        require(fine.relativeVolumeError < coarse.relativeVolumeError, "geometry convergence");
        std::cout << "{\"passed\":true,\"configurations\":" << configurations
                  << ",\"rollback_cases\":" << refusals << ",\"swept_volume_cv\":" << before.cv
                  << ",\"uniform_volume_cv\":" << after.cv
                  << ",\"swept_quality\":" << before.minimumQuality
                  << ",\"uniform_quality\":" << after.minimumQuality << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
