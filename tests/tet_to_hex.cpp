// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "geometry_checks.hpp"
#include "paacs/algo/TetToHexMesher.hpp"
#include "paacs/mesh/Hex8Jacobian.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <iostream>
using namespace paacs;
using geometry_test::require;
namespace {
mesh::Mesh tetra(double scale, int defect = -1) {
    mesh::Mesh m;
    for (Vec3 p : std::array<Vec3, 4>{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}})
        m.addNode(p * scale, {2, 1});
    if (defect == 4)
        m.addNode({3 * scale, 3 * scale, 3 * scale}, {3, 1});
    if (defect == 3)
        m.addCell(mesh::CellType::Tet4, {0, 2, 1, 3}, {3, 1});
    else
        m.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}, {3, 1});
    constexpr std::array<std::array<Index, 3>, 4> faces{
        {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}};
    for (unsigned i = 0; i < 4; ++i) {
        if (defect == 0 && i == 0)
            continue;
        auto f = faces[i];
        if ((defect == 1 && i == 0) || defect == 3)
            std::swap(f[1], f[2]);
        m.addCell(mesh::CellType::Tri3, std::span<const Index>(f),
                  {2, defect == 5 && i == 0 ? 2 : 1});
        if (defect == 2 && i == 0)
            m.addCell(mesh::CellType::Tri3, std::span<const Index>(f), {2, 1});
    }
    if (defect == 6) {
        for (Vec3 p : std::array<Vec3, 4>{{{3, 0, 0}, {4, 0, 0}, {3, 1, 0}, {3, 0, 1}}})
            m.addNode(p * scale, {2, 1});
        m.addCell(mesh::CellType::Tet4, {4, 5, 6, 7}, {3, 2});
        for (auto f : faces) {
            for (auto &id : f)
                id += 4;
            m.addCell(mesh::CellType::Tri3, std::span<const Index>(f), {2, 1});
        }
    }
    return m;
}
} // namespace
int main() {
    try {
        unsigned cases = 0, refusals = 0;
        for (double scale : {1e-6, 1., 1e6}) {
            const auto source = tetra(scale);
            const auto hash = mesh::meshHash(source);
            algo::TetToHexOptions o;
            mesh::Mesh m;
            const auto r = algo::convertTet4ToHex8(source, o, m);
            require(r.published && r.nodes == 15 && r.hexes == 4 && r.boundaryQuads == 12,
                    "canonical subdivision");
            require(mesh::meshHash(source) == hash, "immutable source");
            require(mesh::referenceLinearVolumeBoundary(m).valid,
                    "independent conforming boundary");
            for (Index i = 0; i < m.cellCount(); ++i)
                if (m.cell(i).onEntity.dim == 2) {
                    const auto *v = m.cellNodes(i);
                    const auto a = m.node(v[0]), b = m.node(v[1]), c = m.node(v[2]);
                    require((b - a).cross(c - a).dot(
                                a - Vec3{.25 * scale, .25 * scale, .25 * scale}) > 0,
                            "independent outward boundary winding");
                }
            const double volume = geometry_test::positiveVolume(m);
            require(std::abs(volume - scale * scale * scale / 6) < 1e-12 * volume,
                    "independent conserved volume");
            for (Index i = 0; i < 4; ++i) {
                std::array<Vec3, 8> p{};
                for (unsigned j = 0; j < 8; ++j)
                    p[j] = m.node(m.cellNodes(i)[j]);
                require(mesh::positiveHex8Jacobian(p), "full Jacobian");
            }
            auto parallel = o;
            parallel.threads = 4;
            mesh::Mesh other;
            require(algo::convertTet4ToHex8(source, parallel, other).published &&
                        mesh::meshHash(other) == mesh::meshHash(m),
                    "parallel determinism");
            auto inplace = source;
            require(algo::convertTet4ToHex8(inplace, o, inplace).published &&
                        mesh::meshHash(inplace) == mesh::meshHash(m),
                    "source/destination alias");
            auto exact = o;
            exact.maximumWorkingBytes = r.workingBytes;
            exact.maximumOutputBytes = r.outputBytes;
            require(algo::convertTet4ToHex8(source, exact, other).published,
                    "exact resource limits");
            ++cases;
            auto refuse = [&](const mesh::Mesh &src, const algo::TetToHexOptions &q) {
                mesh::Mesh old;
                old.addNode({9, 8, 7}, {0, 9});
                const auto before = mesh::meshHash(old);
                const auto result = algo::convertTet4ToHex8(src, q, old);
                require(!result.published && result.destinationUnchanged &&
                            mesh::meshHash(old) == before,
                        "conversion rollback case " + std::to_string(refusals) +
                            " scale=" + std::to_string(scale) +
                            " published=" + std::to_string(result.published) +
                            " reason=" + algo::primitiveReasonName(result.reason));
                ++refusals;
            };
            for (int defect = 0; defect < 7; ++defect)
                refuse(tetra(scale, defect), o);
            for (int which = 0; which < 10; ++which) {
                auto q = o;
                std::atomic<bool> stop{true};
                if (which == 0)
                    q.maximumWorkingBytes = r.workingBytes - 1;
                if (which == 1)
                    q.maximumOutputBytes = r.outputBytes - 1;
                if (which == 2)
                    q.maximumNodes = 14;
                if (which == 3)
                    q.maximumCells = 15;
                if (which == 4)
                    q.minimumQuality = .99;
                if (which == 5) {
                    q.warningQuality = .99;
                    q.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
                }
                if (which == 6)
                    q.cancellation = &stop;
                if (which == 7)
                    q.auditLimits.cancellation = &stop;
                if (which == 8)
                    q.auditLimits.maximumFaceUses = 1;
                if (which == 9)
                    q.warningQuality = std::numeric_limits<double>::quiet_NaN();
                refuse(source, q);
            }
            o.warningQuality = .99;
            const auto warned = algo::convertTet4ToHex8(source, o, other);
            require(warned.published && warned.quality.lowQualityCells == 4,
                    "bounded warning diagnostics");
        }
        std::cout << "Tet4 to Hex8: PASS (" << cases << " scaled meshes, " << refusals
                  << " refusal/rollback cases)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
