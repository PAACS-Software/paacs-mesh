// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/mesh/Quality.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/ElementTraits.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace paacs::mesh {
namespace {

constexpr int kTetFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

// Corner connectivity via ElementTraits (M4-B.0): membership and corner
// counts come from the authoritative tables, not per-type branches, so a
// newly enabled dimension-3 or dimension-2 type is audited automatically.
// (Public mesh::CellType currently carries only simplex-family types, so
// the fixed-4/-3 output arrays are exact; generic corner counts arrive
// with generic storage at M4-B.2.)
// CRITICAL SCOPE RULE (second review round): only the SIMPLEX families
// participate in the sorted-triple face machinery below. Truncating a
// Hex8 to its first four corners would audit it as a WRONG TETRAHEDRON -
// worse than the earlier invisibility. Non-simplex volume/surface cells
// are counted as deferredAuditCells and diagnosed explicitly; their full
// audit arrives with shape-aware faces (M4-B.4) and the generic audit
// (M4-B.8).
bool volumeCorners(const Mesh& m, Index c, std::array<Index, 4>& out) {
    topo::ElementType et;
    if (!topo::elementTypeOf(m.cell(c).type, et)) return false;
    const auto& tr = topo::traits(et);
    if (tr.dimension != 3 || tr.cornerCount != 4) return false;
    const Index* v = m.cellNodes(c);
    for (int k = 0; k < 4; ++k) out[size_t(k)] = v[k];
    return true;
}
bool surfaceCorners(const Mesh& m, Index c, std::array<Index, 3>& out) {
    topo::ElementType et;
    if (!topo::elementTypeOf(m.cell(c).type, et)) return false;
    const auto& tr = topo::traits(et);
    if (tr.dimension != 2 || tr.cornerCount != 3) return false;
    const Index* v = m.cellNodes(c);
    for (int k = 0; k < 3; ++k) out[size_t(k)] = v[k];
    return true;
}
bool auditCovered(const Mesh& m, Index c) {
    topo::ElementType et;
    if (!topo::elementTypeOf(m.cell(c).type, et)) return false;
    const auto& tr = topo::traits(et);
    if (tr.dimension == 3) return tr.cornerCount == 4;
    if (tr.dimension == 2) return tr.cornerCount == 3;
    return true;                                 // lines: nothing to defer
}

// Cap on per-instance diagnostics of one kind: after this many, a summary
// diagnostic carries the remainder (keeps reports readable on badly broken
// meshes without hiding the total).
constexpr std::size_t kPerKindDetailCap = 8;

struct UnionFind {
    std::vector<Index> p;
    explicit UnionFind(std::size_t n) : p(n) {
        std::iota(p.begin(), p.end(), Index{0});
    }
    Index find(Index x) {
        while (p[size_t(x)] != x) {
            p[size_t(x)] = p[size_t(p[size_t(x)])];
            x = p[size_t(x)];
        }
        return x;
    }
    void unite(Index a, Index b) {
        a = find(a);
        b = find(b);
        if (a != b) p[size_t(b)] = a;
    }
};

} // namespace

namespace {
// Structured hash streams (review requirement): every hash begins with a
// domain marker + policy version, then explicit counts and section tags,
// and every cell record carries its connectivity length. Structurally
// different serializations can therefore never share an input stream,
// and policy changes are explicit in the stream itself rather than only
// in re-frozen constants.
constexpr std::uint64_t kSectionNodes = 0x4e4f444553ULL;   // "NODES"
constexpr std::uint64_t kSectionCells = 0x43454c4c53ULL;   // "CELLS"
void mixString(std::uint64_t& h, const char* s2) {
    for (const char* p = s2; *p; ++p) {
        h ^= static_cast<unsigned char>(*p);
        h *= 1099511628211ULL;
    }
}
} // namespace

std::uint64_t meshHash(const Mesh& m) {
    std::uint64_t h = 1469598103934665603ULL;
    auto mix = [&h](std::uint64_t v) {
        for (int b = 0; b < 8; ++b) {
            h ^= (v >> (8 * b)) & 0xFF;
            h *= 1099511628211ULL;
        }
    };
    mixString(h, "PAACS_STRICT_HASH_V2");   // V2 == structured stream
    mix(kSectionNodes);
    mix(static_cast<std::uint64_t>(m.nodeCount()));
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const Vec3& p = m.node(i);
        std::uint64_t bx, by, bz;
        std::memcpy(&bx, &p.x, 8);
        std::memcpy(&by, &p.y, 8);
        std::memcpy(&bz, &p.z, 8);
        mix(bx);
        mix(by);
        mix(bz);
        mix(static_cast<std::uint64_t>(m.nodeClass(i).dim));
        mix(static_cast<std::uint64_t>(m.nodeClass(i).tag));
    }
    mix(kSectionCells);
    mix(static_cast<std::uint64_t>(m.cellCount()));
    for (Index c = 0; c < m.cellCount(); ++c) {
        const auto& cell = m.cell(c);
        mix(static_cast<std::uint64_t>(cell.type));
        const int nv = nodesPerCell(cell.type);
        mix(static_cast<std::uint64_t>(nv));    // per-cell length
        const Index* v = m.cellNodes(c);
        for (int k = 0; k < nv; ++k)
            mix(static_cast<std::uint64_t>(v[size_t(k)]));
        mix(static_cast<std::uint64_t>(cell.onEntity.dim));
        mix(static_cast<std::uint64_t>(cell.onEntity.tag));
    }
    return h;
}

namespace {
struct Enc { std::int64_t s, e, mant; };
Enc encodeRel(double x) {
    // (sign, exponent, 40-bit mantissa): relative precision 2^-40 at any
    // magnitude, overflow-free. Used for the model ORIGIN and SCALE.
    if (!std::isfinite(x)) return {2, 0, 0};
    if (x == 0.0) return {0, 0, 0};
    int e = 0;
    const double mag = std::frexp(std::abs(x), &e);
    std::int64_t mi = std::llround(mag * 1099511627776.0);   // 2^40
    if (mi >= 1099511627776LL) { mi >>= 1; ++e; }
    return {x < 0 ? -1 : 1, e, mi};
}
double decodeRel(const Enc& en) {
    if (en.s == 0 || en.s == 2) return 0.0;
    return double(en.s) * std::ldexp(double(en.mant), int(en.e) - 40);
}
} // namespace

SemanticHash semanticHashEx(const Mesh& m) {
    // v3 CANONICAL PHYSICAL MESH HASH (second review redesign).
    //   * ORIGIN + SCALE HASHED SEPARATELY, LOCALS RELATIVE: the model
    //     origin (bbox center) and scale (bbox diagonal) are encoded with
    //     relative precision and hashed explicitly; node coordinates are
    //     encoded as (p - origin)/scale on an ABSOLUTE 2^-40 grid, so
    //     local geometric detail is preserved regardless of distance from
    //     the global origin (v2 hid a 0.1 edge change at 1e12 - review
    //     repro), while translation and scale remain semantic via the
    //     origin/scale terms. Origin/scale are decode-requantized before
    //     use so all platforms subtract identical values.
    //   * AMBIGUITY IS DETECTED, NOT HIDDEN: if two distinct nodes share
    //     an encoded key (same quantized coords + classification), the
    //     relabeling is not canonical; the result carries canonical=false
    //     instead of silently depending on input order (v2's "exact
    //     coincidence only" claim was wrong: any same-bucket pair broke
    //     invariance - review repro at 1e-13 separation).
    //   * AUTOMORPHISM-CANONICAL CELLS: each cell takes the
    //     lexicographically smallest relabeling over its type's
    //     orientation-preserving rotation group (topo::rotations,
    //     generated from reference geometry), with high-order nodes
    //     re-slotted through the edge mapping. This is generic: the v2
    //     sorted-corners+parity scheme collided a valid Quad4 with a
    //     twisted one and was rotation-sensitive (review repro).
    SemanticHash out;
    std::uint64_t h = 1469598103934665603ULL;
    auto mix = [&h](std::uint64_t v) {
        for (int b = 0; b < 8; ++b) {
            h ^= (v >> (8 * b)) & 0xFF;
            h *= 1099511628211ULL;
        }
    };
    mixString(h, "PAACS_SEMANTIC_HASH_V3");
    // ---- origin and scale ---------------------------------------------------
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    bool anyFinite = false;
    for (Index i2 = 0; i2 < m.nodeCount(); ++i2) {
        const Vec3& p = m.node(i2);
        if (!std::isfinite(p.x) || !std::isfinite(p.y) ||
            !std::isfinite(p.z))
            continue;
        anyFinite = true;
        lo.x = std::min(lo.x, p.x); hi.x = std::max(hi.x, p.x);
        lo.y = std::min(lo.y, p.y); hi.y = std::max(hi.y, p.y);
        lo.z = std::min(lo.z, p.z); hi.z = std::max(hi.z, p.z);
    }
    if (!anyFinite) { lo = {0,0,0}; hi = {0,0,0}; }
    const Enc ox = encodeRel((lo.x + hi.x) * 0.5);
    const Enc oy = encodeRel((lo.y + hi.y) * 0.5);
    const Enc oz = encodeRel((lo.z + hi.z) * 0.5);
    const Enc sc = encodeRel((hi - lo).norm());
    for (const Enc* e : {&ox, &oy, &oz, &sc}) {
        mix(std::uint64_t(e->s));
        mix(std::uint64_t(e->e));
        mix(std::uint64_t(e->mant));
    }
    const Vec3 origin{decodeRel(ox), decodeRel(oy), decodeRel(oz)};
    const double scale = std::max(decodeRel(sc), 1e-300);
    auto quantLocal = [&](double x, double o) -> std::int64_t {
        if (!std::isfinite(x)) return std::int64_t(1) << 62;
        const double l = (x - o) / scale;               // ~[-1, 1]
        return std::llround(l * 1099511627776.0);       // 2^-40 grid
    };
    // ---- canonical node relabeling with ambiguity detection ------------------
    struct NKey {
        std::array<std::int64_t, 3> q;
        int dim, tag;
        Index orig;
        bool sameKey(const NKey& o) const {
            return q == o.q && dim == o.dim && tag == o.tag;
        }
        bool operator<(const NKey& o) const {
            if (q != o.q) return q < o.q;
            if (dim != o.dim) return dim < o.dim;
            if (tag != o.tag) return tag < o.tag;
            return orig < o.orig;
        }
    };
    std::vector<NKey> keys(static_cast<std::size_t>(m.nodeCount()));
    for (Index i2 = 0; i2 < m.nodeCount(); ++i2) {
        const Vec3& p = m.node(i2);
        keys[size_t(i2)] = {{quantLocal(p.x, origin.x),
                             quantLocal(p.y, origin.y),
                             quantLocal(p.z, origin.z)},
                            m.nodeClass(i2).dim, m.nodeClass(i2).tag, i2};
    }
    std::sort(keys.begin(), keys.end());
    out.canonical = true;
    for (std::size_t r = 1; r < keys.size(); ++r)
        if (keys[r].sameKey(keys[r - 1])) out.canonical = false;
    std::vector<Index> label(static_cast<std::size_t>(m.nodeCount()));
    for (std::size_t r = 0; r < keys.size(); ++r)
        label[size_t(keys[r].orig)] = static_cast<Index>(r);
    mix(kSectionNodes);
    mix(static_cast<std::uint64_t>(m.nodeCount()));
    for (const NKey& k : keys) {
        for (const std::int64_t v : k.q) mix(std::uint64_t(v));
        mix(std::uint64_t(k.dim));
        mix(std::uint64_t(k.tag));
    }
    // ---- automorphism-canonical cells ----------------------------------------
    struct CKey {
        int type;
        std::vector<Index> t;
        int cdim, ctag;
        bool operator<(const CKey& o) const {
            if (type != o.type) return type < o.type;
            if (t != o.t) return t < o.t;
            if (cdim != o.cdim) return cdim < o.cdim;
            return ctag < o.ctag;
        }
    };
    std::vector<CKey> cells;
    cells.reserve(static_cast<std::size_t>(m.cellCount()));
    for (Index c = 0; c < m.cellCount(); ++c) {
        const auto& cell = m.cell(c);
        const Index* v = m.cellNodes(c);
        topo::ElementType et;
        const bool mapped = topo::elementTypeOf(cell.type, et);
        CKey k;
        k.type = static_cast<int>(cell.type);
        k.cdim = cell.onEntity.dim;
        k.ctag = cell.onEntity.tag;
        if (!mapped) {                       // future-proofing only
            for (int q = 0; q < nodesPerCell(cell.type); ++q)
                k.t.push_back(label[size_t(v[q])]);
            out.canonical = false;
        } else {
            const auto& tr = topo::traits(et);
            std::vector<Index> relabeled(size_t(tr.nodeCount));
            for (int q = 0; q < tr.nodeCount; ++q)
                relabeled[size_t(q)] = label[size_t(v[q])];
            std::vector<Index> best;
            std::vector<Index> cand(size_t(tr.nodeCount));
            for (const auto& rot : topo::rotations(et)) {
                for (int q = 0; q < tr.cornerCount; ++q)
                    cand[size_t(q)] = relabeled[size_t(rot[size_t(q)])];
                // High-order re-slotting: new slot q holds the node of
                // the edge whose endpoints are rot(edge q).
                for (int q = 0; q + tr.cornerCount < tr.nodeCount; ++q) {
                    const int a = rot[tr.edges[q].a], b = rot[tr.edges[q].b];
                    int src = -1;
                    for (int e2 = 0; e2 < tr.edgeCount; ++e2) {
                        const int ea = tr.edges[e2].a, eb = tr.edges[e2].b;
                        if ((ea == a && eb == b) || (ea == b && eb == a)) {
                            src = e2;
                            break;
                        }
                    }
                    cand[size_t(tr.cornerCount + q)] =
                        relabeled[size_t(tr.cornerCount + src)];
                }
                if (best.empty() || cand < best) best = cand;
            }
            k.t = std::move(best);
        }
        cells.push_back(std::move(k));
    }
    std::sort(cells.begin(), cells.end());
    mix(kSectionCells);
    mix(static_cast<std::uint64_t>(cells.size()));
    for (const CKey& k : cells) {
        mix(std::uint64_t(k.type));
        mix(std::uint64_t(k.t.size()));         // per-cell tuple length
        for (const Index v : k.t) mix(std::uint64_t(v));
        mix(std::uint64_t(k.cdim));
        mix(std::uint64_t(k.ctag));
    }
    out.value = h;
    return out;
}

std::uint64_t semanticHash(const Mesh& m) {
    // CONTRACT (review hardening): the value-only convenience MUST NOT
    // silently hand out a non-canonical (input-order-dependent) hash.
    // Callers that can handle ambiguity use semanticHashEx().
    const SemanticHash r = semanticHashEx(m);
    if (!r.canonical)
        throw std::logic_error(
            "mesh::semanticHash: relabeling is ambiguous (coincident "
            "quantized node keys); use semanticHashEx() and handle "
            "canonical=false explicitly");
    return r.value;
}

AuditResult auditMesh(const Mesh& m, const AuditOptions& opt) {
    AuditResult r;
    r.nNodes = m.nodeCount();
    r.meshHash = meshHash(m);

    auto add = [&](diag::Code code, diag::Severity sev, std::string msg,
                   std::string hint = {}, Index cell = kInvalidIndex,
                   Index node = kInvalidIndex,
                   Classification ent = Classification{}) {
        diag::Diagnostic d;
        d.code = code;
        d.severity = sev;
        d.stage = "audit";
        d.message = std::move(msg);
        d.suggestion = std::move(hint);
        d.cellIndex = cell;
        d.nodeIndex = node;
        d.entity = ent;
        if (node != kInvalidIndex && node < m.nodeCount()) {
            d.location = m.node(node);
            d.hasLocation = std::isfinite(d.location.x) &&
                            std::isfinite(d.location.y) &&
                            std::isfinite(d.location.z);
        }
        r.report.add(std::move(d));
    };

    // ---- connectivity validation (before ANY traversal) --------------------
    // Mesh::addCell rejects out-of-range indices today, so this cannot fire
    // via the public construction API; it is defense-in-depth for future
    // deserialization, compaction and transactional paths, so the audit
    // diagnoses malformed input instead of crashing while auditing it.
    bool connectivityOk = true;
    {
        std::size_t emitted = 0;
        for (Index c = 0; c < m.cellCount(); ++c) {
            const Index* v = m.cellNodes(c);
            const int nv = nodesPerCell(m.cell(c).type);
            bool ok = nv > 0;
            for (int k = 0; ok && k < nv; ++k)
                ok = v[size_t(k)] >= 0 && v[size_t(k)] < m.nodeCount();
            if (ok) continue;
            ++r.malformedCells;
            connectivityOk = false;
            if (emitted++ < kPerKindDetailCap)
                add(diag::Code::MalformedConnectivity,
                    diag::Severity::InvalidOutput,
                    "cell references an out-of-range node index",
                    "reject this mesh; corruption upstream", c);
        }
        if (r.malformedCells > static_cast<Index>(kPerKindDetailCap))
            add(diag::Code::MalformedConnectivity,
                diag::Severity::InvalidOutput,
                std::to_string(r.malformedCells) + " malformed cell(s) total");
    }

    // ---- non-finite node coordinates (runs even for cell-less meshes:
    //      review finding - a NaN node in an empty mesh went undetected) ----
    {
        std::size_t emitted = 0, total = 0;
        for (Index i = 0; i < m.nodeCount(); ++i) {
            const Vec3& p = m.node(i);
            if (std::isfinite(p.x) && std::isfinite(p.y) &&
                std::isfinite(p.z))
                continue;
            ++total;
            if (emitted++ < kPerKindDetailCap)
                add(diag::Code::NonFiniteCoordinates,
                    diag::Severity::InvalidOutput,
                    "node has non-finite coordinates",
                    "reject upstream input; see delaunay3D input gate",
                    kInvalidIndex, i, m.nodeClass(i));
        }
        if (total > kPerKindDetailCap)
            add(diag::Code::NonFiniteCoordinates,
                diag::Severity::InvalidOutput,
                std::to_string(total - kPerKindDetailCap) +
                    " further node(s) with non-finite coordinates");
    }

    if (m.cellCount() == 0) {
        add(diag::Code::EmptyMesh, diag::Severity::Warning,
            "mesh contains no cells");
        if (opt.boundaryMode == BoundaryComparisonMode::Required)
            add(diag::Code::BoundaryOrphanFacet,
                diag::Severity::RecoverableFailure,
                "boundary comparison required but the mesh has no cells");
        return r;
    }
    if (!connectivityOk) return r;   // traversals below assume valid indices

    // ---- audit-coverage accounting (loud, never silent) ----------------------
    for (Index c = 0; c < m.cellCount(); ++c)
        if (!auditCovered(m, c)) ++r.deferredAuditCells;
    if (r.deferredAuditCells > 0)
        add(diag::Code::AuditCoverageDeferred, diag::Severity::Warning,
            std::to_string(r.deferredAuditCells) +
                " cell(s) of types the topology audit does not yet cover "
                "(Hex8/Wedge6/Quad4): boundary, manifoldness, components "
                "and region volumes EXCLUDE them",
            "coverage lands with M4-B.4 shape-aware faces + M4-B.8 "
            "generic audit; do not treat this audit as complete for "
            "these cells");

    // ---- element validity via Quality ---------------------------------------
    {
        const QualityReport q = analyzeTets(m);
        r.nTet4 = q.nTets;
        r.nTet10 = q.nTet10;
        r.nInverted = q.nInverted;
        r.nInvalidTet10 = q.nTet10Invalid;
        r.nNonFinite = q.nNonFinite;
        // Per-cell attribution (review: aggregates alone waste the
        // vocabulary's attribution fields). The pass runs only when the
        // aggregate counts demand it.
        if (q.nInverted > 0 || q.nTet10Invalid > 0) {
            std::size_t invEmitted = 0, t10Emitted = 0;
            for (Index c = 0; c < m.cellCount(); ++c) {
                const auto type = m.cell(c).type;
                const Index* v = m.cellNodes(c);
                if (type == CellType::Tet4 &&
                    invEmitted < kPerKindDetailCap) {
                    const double vol =
                        tetVolume(m.node(v[0]), m.node(v[1]), m.node(v[2]),
                                  m.node(v[3]));
                    if (std::isfinite(vol) && predicates::orient3d(
                            m.node(v[0]),m.node(v[1]),m.node(v[2]),m.node(v[3]))
                            != predicates::Sign::Positive) {
                        diag::Diagnostic d;
                        d.code = diag::Code::InvertedElement;
                        d.severity = diag::Severity::InvalidOutput;
                        d.stage = "audit";
                        d.message = "inverted Tet4 (volume " +
                                    std::to_string(vol) + ")";
                        d.cellIndex = c;
                        d.entity = m.cell(c).onEntity;
                        d.location = (m.node(v[0]) + m.node(v[1]) +
                                      m.node(v[2]) + m.node(v[3])) * 0.25;
                        d.hasLocation = true;
                        r.report.add(std::move(d));
                        ++invEmitted;
                    }
                } else if (type == CellType::Tet10 &&
                           t10Emitted < kPerKindDetailCap) {
                    Vec3 p[10];
                    bool fin = true;
                    for (int k = 0; k < 10; ++k) {
                        p[k] = m.node(v[k]);
                        fin = fin && std::isfinite(p[k].x) &&
                              std::isfinite(p[k].y) && std::isfinite(p[k].z);
                    }
                    if (fin && tet10MinJacobianRatio(p) <= 0.0) {
                        diag::Diagnostic d;
                        d.code = diag::Code::InvalidTet10;
                        d.severity = diag::Severity::InvalidOutput;
                        d.stage = "audit";
                        d.message =
                            "Tet10 with non-positive sampled Jacobian";
                        d.cellIndex = c;
                        d.entity = m.cell(c).onEntity;
                        d.location =
                            (p[0] + p[1] + p[2] + p[3]) * 0.25;
                        d.hasLocation = true;
                        r.report.add(std::move(d));
                        ++t10Emitted;
                    }
                }
            }
        }
        if (q.nInverted > 0)
            add(diag::Code::InvertedElement, diag::Severity::InvalidOutput,
                std::to_string(q.nInverted) + " inverted Tet4 element(s) "
                "in total",
                "optimization gates should have prevented this; do not "
                "export");
        if (q.nTet10Invalid > 0)
            add(diag::Code::InvalidTet10, diag::Severity::InvalidOutput,
                std::to_string(q.nTet10Invalid) +
                    " Tet10 element(s) with non-positive sampled Jacobian "
                    "in total",
                "midside straightening failed; do not export");
        if (q.nNonFinite > 0)
            add(diag::Code::NonFiniteCoordinates,
                diag::Severity::InvalidOutput,
                std::to_string(q.nNonFinite) +
                    " element(s) reference non-finite coordinates");
    }

    // ---- faces (CORNER topology): manifoldness + boundary + components -----
    {
        std::map<std::array<Index, 3>, int> faceUse;
        std::vector<Index> volCells;
        for (Index c = 0; c < m.cellCount(); ++c) {
            std::array<Index, 4> tv;
            if (!volumeCorners(m, c, tv)) continue;
            volCells.push_back(c);
            for (const auto& f : kTetFace) {
                std::array<Index, 3> key{tv[size_t(f[0])], tv[size_t(f[1])],
                                         tv[size_t(f[2])]};
                std::sort(key.begin(), key.end());
                ++faceUse[key];
            }
        }

        std::size_t nmEmitted = 0;
        for (const auto& [key, uses] : faceUse) {
            if (uses <= 2) continue;
            ++r.nonManifoldFaces;
            if (nmEmitted++ < kPerKindDetailCap)
                add(diag::Code::NonManifoldFace,
                    diag::Severity::InvalidOutput,
                    "volume-cell corner face shared by " +
                        std::to_string(uses) + " cells",
                    "construction bug; complex is not a valid manifold");
        }
        if (r.nonManifoldFaces > static_cast<Index>(kPerKindDetailCap))
            add(diag::Code::NonManifoldFace, diag::Severity::InvalidOutput,
                std::to_string(r.nonManifoldFaces) +
                    " non-manifold faces in total");

        // Prescribed boundary: Tri3 + Tri6 corner triangles. Duplicates are
        // counted and diagnosed (a map would silently swallow them -
        // review finding).
        std::map<std::array<Index, 3>, int> prescribed;
        for (Index c = 0; c < m.cellCount(); ++c) {
            std::array<Index, 3> key;
            if (!surfaceCorners(m, c, key)) continue;
            if (m.cell(c).type == CellType::Tri3) ++r.nTri3;
            else ++r.nTri6;
            std::sort(key.begin(), key.end());
            ++prescribed[key];
        }
        r.boundaryPrescribed = static_cast<Index>(prescribed.size());
        for (const auto& [key, uses] : prescribed)
            r.boundaryDuplicates += static_cast<Index>(uses - 1);
        if (r.boundaryDuplicates > 0)
            add(diag::Code::DuplicateBoundaryFacet,
                opt.duplicateFacetsAreErrors
                    ? diag::Severity::RecoverableFailure
                    : diag::Severity::Warning,
                std::to_string(r.boundaryDuplicates) +
                    " prescribed facet(s) duplicate another facet's corner "
                    "triangle (boundaryPrescribed counts unique triangles)",
                "deduplicate the surface mesh upstream");

        const bool engaged =
            opt.boundaryMode == BoundaryComparisonMode::Required ||
            (opt.boundaryMode == BoundaryComparisonMode::WhenSurfacePresent &&
             !prescribed.empty());
        if (engaged && prescribed.empty() && !volCells.empty()) {
            // Required mode, no prescribed boundary: every hull face is
            // unexpected by definition.
            for (const auto& [key, uses] : faceUse)
                if (uses == 1) ++r.boundaryUnexpected;
            add(diag::Code::BoundaryUnexpectedFacet,
                static_cast<std::uint64_t>(r.boundaryUnexpected) <=
                        opt.boundaryUnexpectedBudget
                    ? diag::Severity::Warning
                    : diag::Severity::RecoverableFailure,
                "boundary comparison required but no prescribed facets "
                "exist; all " + std::to_string(r.boundaryUnexpected) +
                    " hull face(s) are unexpected",
                "mesh the surfaces first, or use "
                "BoundaryComparisonMode::Disabled");
        }
        if (engaged && !prescribed.empty() && volCells.empty()) {
            // Prescribed boundary with no volume complex: everything is
            // orphaned (previously passed silently - review finding).
            r.boundaryOrphans = r.boundaryPrescribed;
            add(diag::Code::BoundaryOrphanFacet,
                static_cast<std::uint64_t>(r.boundaryOrphans) <=
                        opt.boundaryOrphanBudget
                    ? diag::Severity::Warning
                    : diag::Severity::RecoverableFailure,
                "all " + std::to_string(r.boundaryOrphans) +
                    " prescribed facet(s) are orphaned: no volume cells "
                    "exist",
                "volume meshing missing or failed upstream");
        }
        if (engaged && !prescribed.empty() && !volCells.empty()) {
            for (const auto& [key, uses] : prescribed) {
                const auto it = faceUse.find(key);
                if (it != faceUse.end() && it->second == 1)
                    ++r.boundaryMatched;
                else
                    ++r.boundaryOrphans;
            }
            for (const auto& [key, uses] : faceUse) {
                if (uses != 1) continue;
                if (!prescribed.count(key)) ++r.boundaryUnexpected;
            }
            if (r.boundaryOrphans > 0) {
                const bool within =
                    static_cast<std::uint64_t>(r.boundaryOrphans) <=
                    opt.boundaryOrphanBudget;
                add(diag::Code::BoundaryOrphanFacet,
                    within ? diag::Severity::Warning
                           : diag::Severity::RecoverableFailure,
                    std::to_string(r.boundaryOrphans) + " of " +
                        std::to_string(r.boundaryPrescribed) +
                        " prescribed boundary facet(s) are not faces of "
                        "the volume complex" +
                        (within ? " (within the declared PM-BND-001 budget)"
                                : ""),
                    "constrained facet recovery (M10) eliminates these");
            }
            if (r.boundaryUnexpected > 0)
                add(diag::Code::BoundaryUnexpectedFacet,
                    static_cast<std::uint64_t>(r.boundaryUnexpected) <=
                            opt.boundaryUnexpectedBudget
                        ? diag::Severity::Warning
                        : diag::Severity::RecoverableFailure,
                    std::to_string(r.boundaryUnexpected) +
                        " volume boundary face(s) match no prescribed facet",
                    "constrained facet recovery (M10) eliminates these");
        }

        // Face-connected components over corner topology.
        if (volCells.empty() && opt.expectedComponents > 0)
            add(diag::Code::ComponentCountMismatch,
                diag::Severity::RecoverableFailure,
                "expected " + std::to_string(opt.expectedComponents) +
                    " volume component(s) but the mesh has no volume cells");
        if (!volCells.empty()) {
            std::map<std::array<Index, 3>, Index> firstOwner;
            UnionFind uf(volCells.size());
            for (std::size_t ti = 0; ti < volCells.size(); ++ti) {
                std::array<Index, 4> tv;
                volumeCorners(m, volCells[ti], tv);
                for (const auto& f : kTetFace) {
                    std::array<Index, 3> key{tv[size_t(f[0])],
                                             tv[size_t(f[1])],
                                             tv[size_t(f[2])]};
                    std::sort(key.begin(), key.end());
                    const auto it = firstOwner.find(key);
                    if (it == firstOwner.end())
                        firstOwner[key] = static_cast<Index>(ti);
                    else
                        uf.unite(it->second, static_cast<Index>(ti));
                }
            }
            std::map<Index, int> roots;
            for (std::size_t ti = 0; ti < volCells.size(); ++ti)
                roots[uf.find(static_cast<Index>(ti))] = 1;
            r.connectedComponents = static_cast<int>(roots.size());
            if (opt.expectedComponents > 0 &&
                r.connectedComponents != opt.expectedComponents)
                add(diag::Code::ComponentCountMismatch,
                    diag::Severity::RecoverableFailure,
                    "tet complex has " +
                        std::to_string(r.connectedComponents) +
                        " face-connected component(s); expected " +
                        std::to_string(opt.expectedComponents),
                    "a solid may have split or merged; investigate");
            else if (r.connectedComponents > 1)
                add(diag::Code::MultipleComponents, diag::Severity::Info,
                    "tet complex has " +
                        std::to_string(r.connectedComponents) +
                        " face-connected components",
                    "expected for multi-solid models; set "
                    "expectedComponents to enforce");
        }
    }

    // ---- unused nodes / point accounting -----------------------------------
    {
        std::vector<char> used(static_cast<std::size_t>(m.nodeCount()), 0);
        for (Index c = 0; c < m.cellCount(); ++c) {
            const Index* v = m.cellNodes(c);
            for (int k = 0; k < nodesPerCell(m.cell(c).type); ++k)
                used[size_t(v[size_t(k)])] = 1;
        }
        for (const char u : used) r.unusedNodes += u ? 0 : 1;
        if (r.unusedNodes > 0)
            add(diag::Code::UnusedNodes,
                opt.unusedNodesAreErrors
                    ? diag::Severity::RecoverableFailure
                    : diag::Severity::Warning,
                std::to_string(r.unusedNodes) +
                    " node(s) referenced by no cell",
                "compact the mesh or investigate upstream point loss");
    }

    // ---- classification sanity ----------------------------------------------
    {
        std::size_t bad = 0;
        for (Index i = 0; i < m.nodeCount() && bad < kPerKindDetailCap;
             ++i) {
            const auto c = m.nodeClass(i);
            if (c.dim < -1 || c.dim > 3) {
                ++bad;
                add(diag::Code::ClassificationInconsistent,
                    diag::Severity::Warning,
                    "node classification dim out of range (" +
                        std::to_string(c.dim) + ")",
                    {}, kInvalidIndex, i, c);
            }
        }
        for (Index c = 0; c < m.cellCount() && bad < kPerKindDetailCap;
             ++c) {
            const auto& cell = m.cell(c);
            const bool isVol = cell.type == CellType::Tet4 ||
                               cell.type == CellType::Tet10;
            const bool isSurf = cell.type == CellType::Tri3 ||
                                cell.type == CellType::Tri6;
            const bool isLine = cell.type == CellType::Line2 ||
                                cell.type == CellType::Line3;
            const int d = cell.onEntity.dim;
            const bool tetOk = !isVol || d == 3 || d == -1;
            const bool triOk = !isSurf || d == 2 || d == -1;
            const bool lineOk = !isLine || d == -1 || d == 0 || d == 1;
            if (!tetOk || !triOk || !lineOk) {
                ++bad;
                add(diag::Code::ClassificationInconsistent,
                    diag::Severity::Warning,
                    "cell classification dim inconsistent with cell type",
                    {}, c, kInvalidIndex, cell.onEntity);
            }
        }
    }

    // ---- region volumes -------------------------------------------------------
    {
        for (Index c = 0; c < m.cellCount(); ++c) {
            std::array<Index, 4> tv;
            if (!volumeCorners(m, c, tv)) continue;   // Tet4 AND Tet10
            const double vol = std::abs(
                tetVolume(m.node(tv[0]), m.node(tv[1]), m.node(tv[2]),
                          m.node(tv[3])));
            if (std::isfinite(vol))
                r.regionVolumes[m.cell(c).onEntity.tag] += vol;
        }
        for (const auto& [tag, expect] : opt.expectedRegionVolumes) {
            const auto it = r.regionVolumes.find(tag);
            const double got = it == r.regionVolumes.end() ? 0.0
                                                           : it->second;
            if (std::abs(got - expect) >
                opt.regionVolumeRelTol * std::max(std::abs(expect), 1e-300))
                add(diag::Code::RegionVolumeMismatch,
                    diag::Severity::InvalidOutput,
                    "region tag " + std::to_string(tag) + " volume " +
                        std::to_string(got) + " differs from expected " +
                        std::to_string(expect),
                    "domain not represented; investigate filtering/peeling",
                    kInvalidIndex, kInvalidIndex, Classification{3, tag});
        }
    }

    return r;
}

} // namespace paacs::mesh
