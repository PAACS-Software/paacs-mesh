// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "geometry_checks.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/algo/Tet10JacobianCertificate.hpp"
#include "paacs/algo/TetToHexMesher.hpp"
#include "paacs/mesh/Hex8Jacobian.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "profile_controls_fixture.hpp"
#include <filesystem>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
using namespace profile_test;
using geometry_test::require;
using E = algo::PrimitiveElements;
using S = algo::PrimitiveProfileStrategy;
using R = algo::PrimitiveReason;
namespace {
void inspect(const algo::PrimitiveShape &s, const mesh::Mesh &m,
             const algo::PrimitiveMeshingResult &r, const sizing::CappedSizeProfile &p) {
    mesh::Mesh corners;
    for (Index i = 0; i < m.nodeCount(); ++i)
        corners.addNode(m.node(i), m.nodeClass(i));
    std::map<std::array<Index, 2>, Index> mids;
    std::map<Index, std::array<Index, 2>> reverse;
    std::set<Index> boundaryVertices;
    std::set<std::array<Index, 2>> boundaryEdges;
    std::uint64_t boundaryFaces = 0, volumeCells = 0;
    const auto unit = [](Vec3 v) { return v / std::hypot(v.x, v.y, v.z); };
    const Vec3 z = unit(s.axis), hint = unit(s.xDirection), x = unit(hint - z * hint.dot(z)),
               y = z.cross(x);
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto c = m.cell(i);
        const auto *v = m.cellNodes(i);
        const bool high = c.type == mesh::CellType::Tet10 || c.type == mesh::CellType::Tri6;
        const bool volume = c.onEntity.dim == 3;
        const auto type = c.type == mesh::CellType::Tet10  ? mesh::CellType::Tet4
                          : c.type == mesh::CellType::Tri6 ? mesh::CellType::Tri3
                                                           : c.type;
        const int n = mesh::nodesPerCell(type);
        corners.addCell(type, std::span<const Index>(v, static_cast<std::size_t>(n)), c.onEntity);
        if (volume)
            ++volumeCells;
        else {
            ++boundaryFaces;
            for (int j = 0; j < n; ++j) {
                boundaryVertices.insert(v[j]);
                std::array<Index, 2> e{v[j], v[(j + 1) % n]};
                std::sort(e.begin(), e.end());
                boundaryEdges.insert(e);
                const Vec3 q = m.node(v[j]) - s.origin;
                const double a = q.dot(z), radius = std::hypot(q.dot(x), q.dot(y));
                const double eps = 2e-10 * std::max(s.height, s.radius);
                if (c.onEntity.tag == 2)
                    require(std::abs(a) <= eps, "start cap location");
                else if (c.onEntity.tag == 3)
                    require(std::abs(a - s.height) <= eps, "end cap location");
                else if (c.onEntity.tag == 1)
                    require(std::abs(radius - s.radius) <= r.chordErrorBound + eps,
                            "outer polygon geometry");
                else if (c.onEntity.tag == 4)
                    require(std::abs(radius - s.innerRadius) <= r.chordErrorBound + eps,
                            "inner polygon geometry");
                else
                    require(false, "unexpected semantic face");
            }
        }
        if (high) {
            constexpr std::array<std::array<unsigned, 2>, 6> edges{
                {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}};
            for (unsigned e = 0; e < (volume ? 6u : 3u); ++e) {
                std::array<Index, 2> key{v[edges[e][0]], v[edges[e][1]]};
                std::sort(key.begin(), key.end());
                const Index mid = v[n + e];
                const auto [it, added] = mids.emplace(key, mid);
                require(added || it->second == mid, "shared Tet10 midpoint");
                const auto [back, inserted] = reverse.emplace(mid, key);
                require(inserted || back->second == key, "no midpoint alias");
                const auto a = m.node(key[0]), b = m.node(key[1]), q = m.node(mid);
                require(q.x == std::midpoint(a.x, b.x) && q.y == std::midpoint(a.y, b.y) &&
                            q.z == std::midpoint(a.z, b.z),
                        "straight quadratic edge");
            }
            if (volume) {
                std::array<Vec3, 10> nodes{};
                for (unsigned j = 0; j < 10; ++j)
                    nodes[j] = m.node(v[j]);
                const auto certificate = algo::certifyTet10Jacobian(nodes, {});
                require(certificate.status == algo::Tet10CertificateStatus::Certified,
                        "full Tet10 reference-cell certificate");
            }
        }
        if (type == mesh::CellType::Hex8) {
            std::array<Vec3, 8> nodes{};
            for (unsigned j = 0; j < 8; ++j)
                nodes[j] = m.node(v[j]);
            require(mesh::positiveHex8Jacobian(nodes), "full Hex8 reference-cell certificate");
        }
    }
    const auto audit = mesh::referenceLinearVolumeBoundary(corners);
    require(audit.complete && audit.valid, "independent face incidence");
    const double volume = geometry_test::positiveVolume(corners);
    require(std::abs(volume - r.meshVolume) <= 1e-10 * r.meshVolume,
            "independent divergence volume");
    require(std::abs(volume - r.analyticVolume) / r.analyticVolume < .03, "analytic volume error");
    require(static_cast<long long>(boundaryVertices.size()) -
                    static_cast<long long>(boundaryEdges.size()) +
                    static_cast<long long>(boundaryFaces) ==
                (s.kind == algo::PrimitiveKind::Tube ? 0 : 2),
            "closed boundary Euler characteristic");
    require(std::accumulate(r.quality.histogram.begin(), r.quality.histogram.end(),
                            std::uint64_t(0)) == volumeCells,
            "complete quality histogram");
    const auto plan = sizing::planProfileDistribution(p);
    require(r.profiledSizing && r.sweepLayers + 1 == plan.positions.size(),
            "profile consumer counts");
    require(r.minimumAxialSpacing == plan.minimumWidth &&
                r.maximumAxialSpacing == plan.maximumWidth,
            "profile spacing diagnostics");
    if (r.tet10s)
        require(r.tet10s == r.certifiedTet10s && r.quadraticConnectivityVerified,
                "high order qualification");
    // Original station coordinates remain represented in all element families.
    for (double station : plan.positions) {
        bool found = false;
        for (Index i = 0; i < m.nodeCount() && !found; ++i)
            found = std::abs((m.node(i) - s.origin).dot(z) - station) <= 1e-10 * s.height;
        require(found, "planned axial station appears in output");
    }
}
mesh::Mesh patch(int mode) {
    mesh::Mesh m;
    for (Vec3 p :
         std::array<Vec3, 6>{{{0, 0, 0}, {.5, 0, 0}, {1, 0, 0}, {0, 1, 0}, {.5, 1, 0}, {1, 1, 0}}})
        m.addNode(p, {2, 1});
    for (int j = 0; j < 2; ++j) {
        if (mode == 0 || (mode == 2 && j == 0))
            m.addCell(mesh::CellType::Quad4, {j, j + 1, j + 4, j + 3}, {2, 1});
        else {
            m.addCell(mesh::CellType::Tri3, {j, j + 1, j + 4}, {2, 1});
            m.addCell(mesh::CellType::Tri3, {j, j + 4, j + 3}, {2, 1});
        }
    }
    return m;
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path gallery;
        if (argc == 3 && std::string(argv[1]) == "--gallery")
            gallery = argv[2];
        else
            require(argc == 1, "usage: test [--gallery directory]");
        if (!gallery.empty())
            std::filesystem::create_directories(gallery);
        unsigned cases = 0, refusals = 0;
        auto refuse = [&](algo::PrimitiveShape s, algo::PrimitiveMeshingOptions o, R why) {
            mesh::Mesh m;
            m.addNode({9, 8, 7}, {0, 9});
            const auto before = mesh::meshHash(m);
            const auto r = algo::meshPrimitive(s, o, m);
            require(!r.published && r.destinationUnchanged && mesh::meshHash(m) == before,
                    "strong rollback");
            require(r.reason == why,
                    std::string("wrong refusal: ") + algo::primitiveReasonName(r.reason));
            ++refusals;
        };
        for (auto kind : {algo::PrimitiveKind::Cylinder, algo::PrimitiveKind::Tube})
            for (auto strategy : {S::Sweep, S::Volume}) {
                if (kind == algo::PrimitiveKind::Tube && strategy == S::Volume)
                    continue;
                for (const std::string mode : {"start", "end", "both", "band", "overlap", "fine"})
                    for (auto family : {E::Tet4, E::Tet10, E::HexDominant}) {
                        auto s = shape();
                        s.kind = kind;
                        if (kind == algo::PrimitiveKind::Tube)
                            s.innerRadius = .3;
                        const auto spec = specification(mode);
                        const auto p = sizing::compileCappedSizeProfile(spec);
                        auto o = options(family, strategy);
                        o.sizeProfile = &p.profile;
                        mesh::Mesh m;
                        const auto r = algo::meshPrimitive(s, o, m);
                        const std::string name = std::string(algo::primitiveKindName(kind)) + "-" +
                                                 mode + "-" + algo::primitiveElementsName(family) +
                                                 (strategy == S::Volume ? "-volume" : "-sweep");
                        require(r.published, name + ": " + algo::primitiveReasonName(r.reason));
                        inspect(s, m, r, p.profile);
                        ++cases;
                        auto threaded = o;
                        threaded.threads = 4;
                        mesh::Mesh other;
                        const auto parallel = algo::meshPrimitive(s, threaded, other);
                        require(parallel.published && mesh::meshHash(other) == mesh::meshHash(m),
                                "thread-count determinism");
                        if (strategy == S::Volume) {
                            auto uncached = o;
                            uncached.reuseProfileSlabs = false;
                            mesh::Mesh reference;
                            const auto built = algo::meshPrimitive(s, uncached, reference);
                            require(built.published, "uncached comparison");
                            inspect(s, reference, built, p.profile);
                            require(r.profileReusedSlabs > 0 && built.profileReusedSlabs == 0,
                                    "slab reuse exercised");
                            require(std::abs(r.meshVolume - built.meshVolume) <
                                        1e-12 * r.meshVolume,
                                    "reuse preserves volume");
                            if (mode == "start" || mode == "end") {
                                std::uint64_t a = 0, b = 0;
                                for (Index i = 0; i < m.cellCount(); ++i) {
                                    a += m.cell(i).onEntity == mesh::Classification{2, 2};
                                    b += m.cell(i).onEntity == mesh::Classification{2, 3};
                                }
                                require(mode == "start" ? a > b : b > a,
                                        "selected cap actually finer than remote cap");
                            }
                        } else {
                            const auto plan = sizing::planProfileDistribution(p.profile);
                            auto direct = o;
                            direct.sizeProfile = nullptr;
                            direct.axialPositions = plan.positions;
                            direct.targetSize = p.profile.minimumSize();
                            direct.maximumAxialGrowth = p.profile.maximumGrowth();
                            mesh::Mesh preplanned;
                            require(algo::meshPrimitive(s, direct, preplanned).published &&
                                        mesh::meshHash(preplanned) == mesh::meshHash(m),
                                    "reusable plan identity");
                        }
                        if (!gallery.empty())
                            require(
                                sdk::writeMeshArchive(m, (gallery / (name + ".paacsmesh")).string())
                                    .completed,
                                "gallery archive");
                    }
            }
        for (double scale : {1e-6, 1., 1e6})
            for (auto strategy : {S::Sweep, S::Volume})
                for (auto family : {E::Tet4, E::Tet10, E::HexDominant}) {
                    auto s = shape(scale);
                    s.origin = {3 * scale, -2 * scale, 5 * scale};
                    s.axis = {1, 2, 3};
                    s.xDirection = {1, -2, 1};
                    const auto p = sizing::compileCappedSizeProfile(specification("start", scale));
                    auto o = options(family, strategy, scale);
                    o.sizeProfile = &p.profile;
                    mesh::Mesh m;
                    const auto r = algo::meshPrimitive(s, o, m);
                    require(r.published, "scaled rotated profile");
                    inspect(s, m, r, p.profile);
                    ++cases;
                }
        const auto p = sizing::compileCappedSizeProfile(specification("start"));
        const auto s = shape();
        for (auto strategy : {S::Sweep, S::Volume})
            for (auto family : {E::Tet4, E::Tet10, E::HexDominant}) {
                auto o = options(family, strategy);
                o.sizeProfile = &p.profile;
                mesh::Mesh m;
                const auto normal = algo::meshPrimitive(s, o, m);
                require(normal.published, "budget baseline");
                auto q = o;
                q.maximumWorkingBytes = normal.plannedWorkingBytes;
                mesh::Mesh exact;
                require(algo::meshPrimitive(s, q, exact).published,
                        "exact reported working budget");
                --q.maximumWorkingBytes;
                refuse(s, q, R::WorkingBudget);
                q = o;
                q.maximumOutputBytes = 1;
                refuse(s, q, R::OutputBudget);
                q = o;
                q.maximumNodes = 1;
                refuse(s, q, R::OutputBudget);
                q = o;
                q.profileLimits.maximumCells = 1;
                refuse(s, q, R::WorkingBudget);
                q = o;
                q.warningQuality = .99;
                const auto warning = algo::meshPrimitive(s, q, exact);
                require(warning.published && warning.qualityWarning &&
                            warning.quality.lowQualityCells > 0,
                        "quality warning permits valid mesh");
                q.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
                refuse(s, q, R::QualityRequirement);
                q = o;
                q.minimumQuality = .99;
                refuse(s, q, R::QualityRequirement);
                std::atomic<bool> stop{true};
                q = o;
                q.cancellation = &stop;
                refuse(s, q, R::Cancelled);
                q = o;
                q.profileLimits.cancellation = &stop;
                refuse(s, q, R::Cancelled);
            }
        auto o = options(E::Tet4, S::Volume);
        o.sizeProfile = &p.profile;
        auto tube = s;
        tube.kind = algo::PrimitiveKind::Tube;
        tube.innerRadius = .3;
        refuse(tube, o, R::InvalidRequest);
        auto q = o;
        q.sizeProfile = nullptr;
        refuse(s, q, R::InvalidRequest);
        q = o;
        q.maximumProfileSectionNodes = 1;
        refuse(s, q, R::WorkingBudget);
        q = o;
        q.maximumProfileSlabWork = 1;
        refuse(s, q, R::WorkingBudget);
        q = o;
        q.profileStrategy = static_cast<S>(99);
        refuse(s, q, R::InvalidRequest);
        sizing::CappedSizeProfile invalid;
        q = o;
        q.sizeProfile = &invalid;
        refuse(s, q, R::InvalidRequest);
        q = o;
        auto wrong = s;
        wrong.height = 5;
        refuse(wrong, q, R::InvalidRequest);
        q = options(E::Tet4, S::Sweep);
        std::vector<double> bad{0, 1, .5, 4.8};
        q.axialPositions = bad;
        refuse(s, q, R::InvalidRequest);
        q.axialPositions = {};
        q.sizeProfile = &p.profile;
        auto cone = s;
        cone.kind = algo::PrimitiveKind::Cone;
        refuse(cone, q, R::InvalidRequest);
        auto spec = specification("start");
        spec.controls = {algo::primitiveFaceSizeControl(s, 1, .1)};
        require(!sizing::compileCappedSizeProfile(spec).compiled,
                "curved wall must not silently become cap control");
        ++refusals;
        for (int mode = 0; mode < 3; ++mode) {
            const auto source = patch(mode);
            algo::LayeredPatchRequest request;
            request.sizeProfile = &p.profile;
            mesh::Mesh m;
            const auto r = algo::meshLayeredPatch(source, request, m);
            require(r.published && r.usedSizeProfile, "shared prism/hex profile");
            require(mesh::referenceLinearVolumeBoundary(m).valid, "independent layered incidence");
            require(std::abs(geometry_test::positiveVolume(m) - p.profile.length()) < 1e-10,
                    "independent layered volume");
            const auto plan = sizing::planProfileDistribution(p.profile);
            require(r.schedule.offsets == plan.positions, "one planner reused by layers");
            if (!gallery.empty())
                require(
                    sdk::writeMeshArchive(
                        m, (gallery / ("patch-" + std::to_string(mode) + ".paacsmesh")).string())
                        .completed,
                    "patch gallery");
            ++cases;
        }
        std::cout << "profiled primitives: PASS (" << cases << " meshes, " << refusals
                  << " refusals; independent geometry/topology, cached/uncached, 1/4 threads)\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
