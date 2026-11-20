// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Optimize.hpp"
#include "paacs/mesh/Quality.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <vector>

namespace paacs::algo {

using mesh::CellType;
using mesh::Classification;
using mesh::Mesh;

namespace {

struct WTet {
    std::array<Index, 4> v;
    Classification cls;
    bool alive = true;
};

double qual(const Mesh& m, const std::array<Index, 4>& t) {
    return mesh::tetQuality(m.node(t[0]), m.node(t[1]),
                            m.node(t[2]), m.node(t[3]));
}
double vol(const Mesh& m, const std::array<Index, 4>& t) {
    return mesh::tetVolume(m.node(t[0]), m.node(t[1]),
                           m.node(t[2]), m.node(t[3]));
}
// Orientation fix: returns tet with positive volume, or nullopt-ish flag.
bool fixOrient(const Mesh& m, std::array<Index, 4>& t, double epsV) {
    double V = vol(m, t);
    if (V < 0) { std::swap(t[2], t[3]); V = -V; }
    return V > epsV;
}

// Rebuild the mesh with the surviving tets (nodes and lower-dim cells kept).
void rebuild(Mesh& m, const std::vector<WTet>& tets) {
    Mesh out;
    for (Index i = 0; i < m.nodeCount(); ++i)
        out.addNode(m.node(i), m.nodeClass(i));
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i);
        if (c.type == CellType::Tet4) continue;
        const Index* v = m.cellNodes(i);
        out.addCell(c.type,
                    std::vector<Index>(v, v + mesh::nodesPerCell(c.type)),
                    c.onEntity);
    }
    for (const auto& t : tets)
        if (t.alive)
            out.addCell(CellType::Tet4, {t.v[0], t.v[1], t.v[2], t.v[3]},
                        t.cls);
    m = std::move(out);
}

std::vector<WTet> extract(const Mesh& m) {
    std::vector<WTet> tets;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i);
        if (c.type != CellType::Tet4) continue;
        const Index* v = m.cellNodes(i);
        tets.push_back({{v[0], v[1], v[2], v[3]}, c.onEntity, true});
    }
    return tets;
}

// Targeted sliver vertex perturbation ("exudation-lite"). Flips change
// connectivity and peeling removes all-boundary hull slivers; what remains
// are slivers owning at least one interior node in a geometry where the
// smart-Laplacian centroid move doesn't help (the centroid direction is not
// the volume-opening direction for a flat tet). Here we search a small
// deterministic set of directions/steps for that node and accept whichever
// move maximizes the cavity's WORST quality. Sequential and index-ordered,
// so results are reproducible on any thread count.
Index perturbSliverVertices(Mesh& m, double qBelow, double minGain) {
    std::vector<std::vector<Index>> nodeTets(
        static_cast<size_t>(m.nodeCount()));
    std::vector<Index> tetCells;
    for (Index ci = 0; ci < m.cellCount(); ++ci) {
        if (m.cell(ci).type != CellType::Tet4) continue;
        tetCells.push_back(ci);
        const Index* v = m.cellNodes(ci);
        for (int k = 0; k < 4; ++k)
            nodeTets[static_cast<size_t>(v[k])].push_back(ci);
    }

    // Worst quality over the cavity with `node` moved to `pos`; -1 if any
    // tet inverts.
    auto cavityMin = [&](const std::vector<Index>& cav, Index node,
                         const Vec3& pos) {
        double q = 1e300;
        for (Index ci : cav) {
            const Index* v = m.cellNodes(ci);
            Vec3 p[4];
            for (int k = 0; k < 4; ++k)
                p[k] = (v[k] == node) ? pos : m.node(v[k]);
            if (mesh::tetVolume(p[0], p[1], p[2], p[3]) <= 0) return -1.0;
            q = std::min(q, mesh::tetQuality(p[0], p[1], p[2], p[3]));
        }
        return q;
    };

    static const double s3 = 0.57735026918962576;
    static const Vec3 kDirs[10] = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        {s3, s3, s3}, {-s3, s3, -s3}, {s3, -s3, s3}, {-s3, -s3, -s3}};

    Index moved = 0;
    for (Index ci : tetCells) {
        const Index* v = m.cellNodes(ci);
        const Vec3 &a = m.node(v[0]), &b = m.node(v[1]), &c = m.node(v[2]),
                   &d = m.node(v[3]);
        if (mesh::tetQuality(a, b, c, d) >= qBelow) continue;

        const double hLoc =
            ((b - a).norm() + (c - a).norm() + (d - a).norm() +
             (c - b).norm() + (d - b).norm() + (d - c).norm()) / 6.0;

        for (int k = 0; k < 4; ++k) {
            const Index n = v[k];
            if (m.nodeClass(n).dim != 3) continue;   // interior nodes only
            const auto& cav = nodeTets[static_cast<size_t>(n)];
            const Vec3 cur = m.node(n);
            const double q0 = cavityMin(cav, n, cur);

            Vec3 best = cur;
            double bestQ = q0;
            for (const Vec3& dir : kDirs)
                for (double step : {0.25, 0.5}) {
                    const Vec3 cand = cur + dir * (step * hLoc);
                    const double qc = cavityMin(cav, n, cand);
                    if (qc > bestQ) { bestQ = qc; best = cand; }
                }
            if (bestQ > q0 + minGain) {
                m.setNode(n, best);
                ++moved;
            }
        }
    }
    return moved;
}

// Final-cleanup peel: remove sub-threshold tets exposing a hull face.
// Unlike the VolumeMesher peel this also allows interior nodes (shield
// points behind size-transition caps whose skinny base triangle no vertex
// move can fix); the exposed notch has near-zero volume by construction.
Index peelHullSlivers(std::vector<WTet>& tets, const Mesh& m,
                      double qBelow, const std::set<int>& peelableTags) {
    if (peelableTags.empty()) return 0;   // peeling disabled (planar model)
    Index peeled = 0;
    for (int round = 0; round < 5; ++round) {
        std::map<std::array<Index, 3>, int> faceCount;
        static const int ft[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
        for (const auto& t : tets) {
            if (!t.alive) continue;
            for (const auto& f : ft) {
                std::array<Index, 3> key{t.v[size_t(f[0])], t.v[size_t(f[1])],
                                         t.v[size_t(f[2])]};
                std::sort(key.begin(), key.end());
                ++faceCount[key];
            }
        }
        Index n = 0;
        for (auto& t : tets) {
            if (!t.alive) continue;
            if (qual(m, t.v) >= qBelow) continue;
            bool touchesCurved = false;
            for (Index vv : t.v) {
                const auto cls = m.nodeClass(vv);
                touchesCurved = touchesCurved ||
                    (cls.dim == 2 && peelableTags.count(cls.tag));
            }
            if (!touchesCurved) continue;
            bool hull = false;
            for (const auto& f : ft) {
                std::array<Index, 3> key{t.v[size_t(f[0])], t.v[size_t(f[1])],
                                         t.v[size_t(f[2])]};
                std::sort(key.begin(), key.end());
                hull = hull || faceCount[key] == 1;
            }
            if (!hull) continue;
            t.alive = false;
            ++n;
        }
        peeled += n;
        if (!n) break;
    }
    return peeled;
}

} // namespace

OptimizeResult optimizeTets(Mesh& m, const OptimizeOptions& opt) {
    OptimizeResult res;
    {
        const auto r = mesh::analyzeTets(m);
        res.minQualityBefore = r.radiusRatio.n ? r.radiusRatio.min : 0;
        res.meanQualityBefore = r.radiusRatio.mean;
    }
    const double epsV = 1e-14;

    for (int pass = 0; pass < opt.maxPasses; ++pass) {
        auto tets = extract(m);
        std::vector<bool> touched(tets.size(), false);
        Index applied = 0;

        // ------------------------------------------------- 3-2 flips ------
        // Edge shared by exactly three tets whose "ring" is three vertices:
        // replace by two tets sharing the ring face.
        {
            std::map<std::pair<Index, Index>, std::vector<size_t>> edgeMap;
            static const int etab[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
            for (size_t ti = 0; ti < tets.size(); ++ti)
                for (const auto& e : etab)
                    edgeMap[std::minmax(tets[ti].v[size_t(e[0])],
                                        tets[ti].v[size_t(e[1])])].push_back(ti);

            for (const auto& [edge, incident] : edgeMap) {
                if (incident.size() != 3) continue;
                bool skip = false;
                for (size_t ti : incident)
                    if (!tets[ti].alive || touched[ti]) skip = true;
                if (skip) continue;

                const Index a = edge.first, b = edge.second;
                // Ring vertices: each incident tet contributes its two
                // non-edge vertices; a proper ring has exactly 3 distinct.
                std::map<Index, int> ringCount;
                for (size_t ti : incident)
                    for (Index vv : tets[ti].v)
                        if (vv != a && vv != b) ++ringCount[vv];
                if (ringCount.size() != 3) continue;
                bool proper = true;
                std::array<Index, 3> ring{};
                int k = 0;
                for (const auto& [vv, cnt] : ringCount) {
                    if (cnt != 2) proper = false;
                    ring[size_t(k++)] = vv;
                }
                if (!proper) continue;

                std::array<Index, 4> t1{ring[0], ring[1], ring[2], a};
                std::array<Index, 4> t2{ring[0], ring[1], ring[2], b};
                if (!fixOrient(m, t1, epsV) || !fixOrient(m, t2, epsV))
                    continue;

                // Volume conservation rejects geometrically invalid flips.
                double oldV = 0, oldQ = 1e300;
                for (size_t ti : incident) {
                    oldV += std::abs(vol(m, tets[ti].v));
                    oldQ = std::min(oldQ, qual(m, tets[ti].v));
                }
                const double newV = vol(m, t1) + vol(m, t2);
                if (std::abs(newV - oldV) > 1e-9 * (1 + oldV)) continue;

                const double newQ = std::min(qual(m, t1), qual(m, t2));
                if (newQ < oldQ + opt.minGain) continue;

                for (size_t ti : incident) tets[ti].alive = false;
                const Classification cls = tets[incident[0]].cls;
                tets.push_back({t1, cls, true});
                tets.push_back({t2, cls, true});
                touched.resize(tets.size(), true);
                ++res.flips32;
                ++applied;
            }
        }

        // ------------------------------------------------- 2-3 flips ------
        // Interior face shared by two tets: replace by three tets around the
        // segment joining the two apexes (valid iff it pierces the face,
        // which the positive-volume + conservation checks enforce).
        {
            std::map<std::array<Index, 3>, std::vector<size_t>> faceMap;
            static const int ftab[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
            for (size_t ti = 0; ti < tets.size(); ++ti) {
                if (!tets[ti].alive) continue;
                for (const auto& fc : ftab) {
                    std::array<Index, 3> key{tets[ti].v[size_t(fc[0])],
                                             tets[ti].v[size_t(fc[1])],
                                             tets[ti].v[size_t(fc[2])]};
                    std::sort(key.begin(), key.end());
                    faceMap[key].push_back(ti);
                }
            }

            for (const auto& [face, incident] : faceMap) {
                if (incident.size() != 2) continue;
                const size_t ti = incident[0], tj = incident[1];
                if (!tets[ti].alive || !tets[tj].alive ||
                    touched[ti] || touched[tj]) continue;

                auto apex = [&](size_t t) {
                    for (Index vv : tets[t].v)
                        if (vv != face[0] && vv != face[1] && vv != face[2])
                            return vv;
                    return kInvalidIndex;
                };
                const Index a = apex(ti), b = apex(tj);
                if (a == kInvalidIndex || b == kInvalidIndex || a == b) continue;

                std::array<std::array<Index, 4>, 3> nt{{
                    {face[0], face[1], a, b},
                    {face[1], face[2], a, b},
                    {face[2], face[0], a, b}}};
                bool ok = true;
                for (auto& t : nt) ok = ok && fixOrient(m, t, epsV);
                if (!ok) continue;

                const double oldV = std::abs(vol(m, tets[ti].v)) +
                                    std::abs(vol(m, tets[tj].v));
                const double oldQ = std::min(qual(m, tets[ti].v),
                                             qual(m, tets[tj].v));
                double newV = 0, newQ = 1e300;
                for (const auto& t : nt) {
                    newV += vol(m, t);
                    newQ = std::min(newQ, qual(m, t));
                }
                if (std::abs(newV - oldV) > 1e-9 * (1 + oldV)) continue;
                if (newQ < oldQ + opt.minGain) continue;

                tets[ti].alive = tets[tj].alive = false;
                const Classification cls = tets[ti].cls;
                for (const auto& t : nt) tets.push_back({t, cls, true});
                touched.resize(tets.size(), true);
                ++res.flips23;
                ++applied;
            }
        }

        rebuild(m, tets);
        if (opt.perturbBelow > 0)
            res.perturbations +=
                perturbSliverVertices(m, opt.perturbBelow, opt.minGain);
        if (opt.smoothBetweenPasses) smoothInterior(m, opt.smooth);
        ++res.passes;
        if (applied == 0 && pass > 0) break;   // converged
    }

    // Final cleanup: peel irreparable hull caps, then one smoothing sweep
    // to relax the freshly exposed neighborhood.
    if (opt.peelBelow > 0) {
        auto tets = extract(m);
        res.peeledSlivers = peelHullSlivers(tets, m, opt.peelBelow,
                                            opt.peelableFaceTags);
        if (res.peeledSlivers) {
            rebuild(m, tets);
            if (opt.smoothBetweenPasses) smoothInterior(m, opt.smooth);
        }
    }

    {
        const auto r = mesh::analyzeTets(m);
        res.minQualityAfter = r.radiusRatio.n ? r.radiusRatio.min : 0;
        res.meanQualityAfter = r.radiusRatio.mean;
    }
    return res;
}

} // namespace paacs::algo
