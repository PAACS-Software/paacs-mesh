// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MixedMeshQuality.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/mesh/Quality.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace paacs;
namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
mesh::Mesh tetra(const std::array<Vec3, 4>& points, bool reversed = false) {
    mesh::Mesh out;
    for (auto p : points) out.addNode(p, {3, 1});
    const auto ids = reversed ? std::array<Index, 4>{0, 2, 1, 3}
                              : std::array<Index, 4>{0, 1, 2, 3};
    out.addCell(mesh::CellType::Tet4, ids, {3, 1});
    return out;
}

void positiveCancellationMagnitude() {
    // Owned decimal fixture: the ordinary determinant stays positive, so a
    // fallback limited to nonpositive volumes misses its magnitude error.
    const std::array<Vec3, 4> points{{{0, 0, 0}, {1.1, 2.2, 3.3},
                                   {4.4, 5.5, 6.6}, {5.5, 7.7, 9.899999999999997}}};
    // Fraction over the exact input binary64 coordinates gives determinant
    // 16105500875899657766089347794207 /
    // 1427247692705959881058285969449495136382746624. The mean ratio below
    // was independently evaluated from that rational and the exact squared
    // edge lengths using 90-digit decimal arithmetic.
    constexpr double expectedVolume = 1.8807178037150155555013302010948429e-15;
    constexpr double expectedMeanRatio = 8.634908748816873947012490638375e-12;
    for (int exponent : {-20, 0, 20}) {
        auto scaled = points;
        for (auto& p : scaled) {
            p.x = std::ldexp(p.x, exponent);
            p.y = std::ldexp(p.y, exponent);
            p.z = std::ldexp(p.z, exponent);
        }
        const double volume = std::ldexp(expectedVolume, 3 * exponent);
        const auto a = scaled[1] - scaled[0], b = scaled[2] - scaled[0],
                   c = scaled[3] - scaled[0];
        const double ordinaryVolume = a.cross(b).dot(c) / 6.;
        require(ordinaryVolume > 0 && std::abs(ordinaryVolume - volume) > .1 * volume,
                "owned fixture must expose a positive cancellation magnitude error");
        require(predicates::orient3d(scaled[0], scaled[1], scaled[2], scaled[3]) ==
                    predicates::Sign::Positive,
                "positive-cancellation fixture has exact positive orientation");
        const auto quality = algo::analyzeMixedMeshQuality(tetra(scaled));
        require(quality.complete && quality.valid && quality.invertedCells == 0,
                "positive-cancellation measurement preserves the exact validity gate");
        require(std::abs(quality.totalVolume - volume) <= volume * 1e-12 &&
                    std::abs(quality.minimumVolume - volume) <= volume * 1e-12,
                "positive ordinary volume must still agree with the exact reference");
        require(std::abs(quality.minimumQuality - expectedMeanRatio) <=
                    expectedMeanRatio * 1e-12 &&
                    std::abs(quality.meanQuality - expectedMeanRatio) <=
                    expectedMeanRatio * 1e-12,
                "mean ratio uses the independently measured positive volume");
        const auto reversed = algo::analyzeMixedMeshQuality(tetra(scaled, true));
        require(reversed.complete && !reversed.valid && reversed.invertedCells == 1,
                "magnitude correction must not accept an exactly inverted cell");
        std::array<Vec3, 4> flat{{{0, 0, 0}, {1, 2, 3}, {4, 5, 6}, {5, 7, 9}}};
        for (auto& p : flat) {
            p.x = std::ldexp(p.x, exponent);
            p.y = std::ldexp(p.y, exponent);
            p.z = std::ldexp(p.z, exponent);
        }
        const auto coplanar = algo::analyzeMixedMeshQuality(tetra(flat));
        require(coplanar.complete && !coplanar.valid && coplanar.invertedCells == 1,
                "magnitude correction must not accept an exactly flat cell");
    }
}
}
int main() {
    try {
        positiveCancellationMagnitude();
        // Reproduced from stored CAD recovery coordinates. Python Fraction over
        // the exact input doubles independently supplied these reference volumes.
        const std::array<std::array<Vec3, 4>, 2> cases{{
            {{{-10.798720329183794, -13.05341052813009, 30.615889574913982},
              {-7.2132074812276477, -15.32884617303554, 30.615889574913982},
              {-12.109151794544177, -14.637450059193773, 35.742526902246624},
              {-8.088534720143933, -17.189011242643666, 35.742526902246624}}},
            {{{-12.90098749253732, -22.345165804885198, 69.636806587820928},
              {-11.540819193361763, -19.989285203868594, 62.564589068209543},
              {-9.3430627509201631, -21.106141500706261, 62.564589068209543},
              {-10.444209693618932, -23.593651625090409, 69.636806587820928}}}
        }};
        constexpr std::array expected{1.68631142536503e-16, 2.955458494708311e-16};
        for (std::size_t i = 0; i < cases.size(); ++i) {
            const auto& p = cases[i];
            require(predicates::orient3d(p[0], p[1], p[2], p[3]) == predicates::Sign::Positive,
                    "reproduction has exact positive orientation");
            for (int exponent : {-20, 0, 20}) {
                auto scaled = p;
                for (auto& q : scaled) {
                    q.x = std::ldexp(q.x, exponent);
                    q.y = std::ldexp(q.y, exponent);
                    q.z = std::ldexp(q.z, exponent);
                }
                const auto q = algo::analyzeMixedMeshQuality(tetra(scaled));
                const double volume = std::ldexp(expected[i], 3 * exponent);
                require(q.complete && q.valid && q.invertedCells == 0,
                        "roundoff must not label an exact positive cell inverted");
                require(std::abs(q.totalVolume - volume) <= volume * 1e-12 &&
                            q.minimumQuality > 0 && std::isfinite(q.minimumQuality),
                        "positive volume agrees with the independent exact reference");
                const auto inverted = algo::analyzeMixedMeshQuality(tetra(scaled, true));
                require(inverted.complete && !inverted.valid && inverted.invertedCells == 1,
                        "true inversion remains invalid");
                const auto measured=mesh::analyzeTets(tetra(scaled));
                const auto audited=mesh::auditMesh(tetra(scaled));
                require(measured.nInverted==0 && audited.nInverted==0,
                        "archive and production audits agree on exact orientation");
                require(std::abs(mesh::tetVolume(scaled[0],scaled[1],scaled[2],scaled[3])-volume)
                            <=volume*1e-12,
                        "archive metrics measure the exact-input volume");
                require(mesh::analyzeTets(tetra(scaled,true)).nInverted==1,
                        "archive orientation still rejects real inversion");
            }
        }
        const auto flat = algo::analyzeMixedMeshQuality(tetra(
            {{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}}}));
        require(flat.complete && !flat.valid && flat.invertedCells == 1,
                "exactly coplanar cells remain invalid");
        const auto regular = algo::analyzeMixedMeshQuality(tetra(
            {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}}));
        require(regular.valid && regular.totalVolume == 1.0 / 6.0,
                "well-conditioned volume retains the existing fast result");
        const double nan = std::numeric_limits<double>::quiet_NaN();
        require(std::isnan(predicates::orient3dDeterminantEstimate(
                    {nan,0,0}, {1,0,0}, {0,1,0}, {0,0,1})),
                "nonfinite measurements remain unavailable");
        std::cout << "Exact positive cancellation cells accepted; true inversion and flat cells rejected.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
