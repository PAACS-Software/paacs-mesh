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
}
int main() {
    try {
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
