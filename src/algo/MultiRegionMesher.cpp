// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MultiRegionMesher.hpp"

#include "paacs/predicates/Predicates.hpp"

#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace paacs::algo {

static_assert(std::numeric_limits<double>::is_iec559,
              "T5 exact volume authority requires IEC-559 binary64");

T5Outcome::T5Outcome() = default;
T5Outcome::~T5Outcome() = default;
T5Outcome::T5Outcome(T5Outcome &&) noexcept = default;
T5Outcome &T5Outcome::operator=(T5Outcome &&) noexcept = default;

namespace {

bool sameClass(mesh::Classification a, mesh::Classification b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}
template <class T> bool checkedAdd(T a, T b, T &out) noexcept {
    if (b > (std::numeric_limits<T>::max)() - a)
        return false;
    out = a + b;
    return true;
}
bool checkedMul(std::size_t a, std::size_t b, std::size_t &out) noexcept {
    if (a != 0 && b > (std::numeric_limits<std::size_t>::max)() / a)
        return false;
    out = a * b;
    return true;
}
bool addBytes(std::size_t &total, std::size_t count, std::size_t width) noexcept {
    std::size_t bytes = 0;
    return checkedMul(count, width, bytes) && checkedAdd(total, bytes, total);
}

void fail(T5Outcome &out, T5Status status, T5Reason reason, T5Phase phase) noexcept {
    out.status = status;
    out.failure.reason = reason;
    out.failure.phase = phase;
    out.reportComplete = out.report && out.report->reportComplete;
}

T5Status sourceStatus(BoundaryConformingV3Status s) noexcept {
    switch (s) {
    case BoundaryConformingV3Status::Published:
        return T5Status::Published;
    case BoundaryConformingV3Status::Rejected:
        return T5Status::Rejected;
    case BoundaryConformingV3Status::Unsupported:
        return T5Status::Unsupported;
    case BoundaryConformingV3Status::Refused:
        return T5Status::Refused;
    case BoundaryConformingV3Status::Failed:
        return T5Status::Failed;
    }
    return T5Status::Failed;
}
T5Reason sourceReason(BoundaryConformingV3Status s) noexcept {
    switch (s) {
    case BoundaryConformingV3Status::Rejected:
        return T5Reason::SourceRejected;
    case BoundaryConformingV3Status::Unsupported:
        return T5Reason::SourceUnsupported;
    case BoundaryConformingV3Status::Refused:
        return T5Reason::SourceRefused;
    case BoundaryConformingV3Status::Failed:
        return T5Reason::SourceFailed;
    case BoundaryConformingV3Status::Published:
        return T5Reason::None;
    }
    return T5Reason::SourceFailed;
}

template <class H>
int permutationParity(const std::array<H, 3> &a, const std::array<H, 3> &canonical) noexcept {
    std::array<int, 3> p{-1, -1, -1};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            if (a[i] == canonical[j]) {
                if (p[i] != -1)
                    return 0;
                p[i] = static_cast<int>(j);
            }
    if (p[0] < 0 || p[1] < 0 || p[2] < 0 || p[0] == p[1] || p[0] == p[2] || p[1] == p[2])
        return 0;
    int inversions = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            inversions += p[i] > p[j];
    return (inversions & 1) ? -1 : 1;
}

struct NodePlan {
    Vec3 coordinate{};
    mesh::Classification classification{};
    Index source = kInvalidIndex;
};
struct SurfacePlan {
    bnd::VolumeFacetKey key{};
    std::array<Index, 3> nodes{};
    mesh::Classification classification{};
    Index source = kInvalidIndex;
};
struct TetPlan {
    topo::CellHandle cell{};
    std::array<Index, 4> nodes{};
    mesh::Classification classification{};
    Index source = kInvalidIndex;
    RegionPersistentId region = 0;
};

struct Expansion {
    std::vector<double> terms;
};
void twoDiff(double a, double b, double &x, double &y) noexcept {
    x = a - b;
    const double bv = a - x;
    const double av = x + bv;
    const double br = bv - b;
    const double ar = a - av;
    y = ar + br;
}
void split(double a, double &hi, double &lo) noexcept {
    constexpr double splitter = 134217729.0;
    const double c = splitter * a;
    const double abig = c - a;
    hi = c - abig;
    lo = a - hi;
}
void twoProduct(double a, double b, double &x, double &y) noexcept {
    x = a * b;
    double ahi = 0, alo = 0, bhi = 0, blo = 0;
    split(a, ahi, alo);
    split(b, bhi, blo);
    y = alo * blo - (((x - ahi * bhi) - alo * bhi) - ahi * blo);
}
Expansion grow(const Expansion &e, double b) {
    Expansion out;
    out.terms.reserve(e.terms.size() + 1);
    double q = b;
    for (double enow : e.terms) {
        const double x = q + enow;
        const double bv = x - q;
        const double hh = (q - (x - bv)) + (enow - bv);
        if (hh != 0)
            out.terms.push_back(hh);
        q = x;
    }
    if (q != 0 || out.terms.empty())
        out.terms.push_back(q);
    return out;
}
Expansion sumExpansion(const Expansion &a, const Expansion &b) {
    Expansion out = a;
    for (double v : b.terms)
        out = grow(out, v);
    return out;
}
Expansion negate(Expansion e) {
    for (double &v : e.terms)
        v = -v;
    return e;
}
Expansion scaleExpansion(const Expansion &e, double b) {
    Expansion out;
    for (double v : e.terms) {
        double x = 0, y = 0;
        twoProduct(v, b, x, y);
        if (y != 0)
            out = grow(out, y);
        if (x != 0)
            out = grow(out, x);
    }
    if (out.terms.empty())
        out.terms.push_back(0);
    return out;
}
Expansion productExpansion(const Expansion &a, const Expansion &b) {
    Expansion out;
    out.terms.push_back(0);
    for (double v : b.terms)
        out = sumExpansion(out, scaleExpansion(a, v));
    return out;
}
Expansion difference(double a, double b) {
    double x = 0, y = 0;
    twoDiff(a, b, x, y);
    Expansion e;
    if (y != 0)
        e.terms.push_back(y);
    if (x != 0 || e.terms.empty())
        e.terms.push_back(x);
    return e;
}
Expansion crossTerm(const Expansion &a, const Expansion &b, const Expansion &c,
                    const Expansion &d) {
    return sumExpansion(productExpansion(a, b), negate(productExpansion(c, d)));
}

bool exactSixVolume(const std::array<Vec3, 4> &p, Expansion &result, std::uint64_t &work) {
    static_assert(kT5SixVolumeExpansionTermsPerTet == 192);
    double magnitude = 0;
    for (const auto &q : p)
        magnitude = (std::max)({magnitude, std::abs(q.x), std::abs(q.y), std::abs(q.z)});
    int scale = 0;
    if (magnitude > 0) {
        const int exponent = std::ilogb(magnitude);
        if (exponent == FP_ILOGB0 || exponent == FP_ILOGBNAN)
            return false;
        scale = -exponent;
    }
    const auto scaled = [&](double value, double &out) noexcept {
        out = std::scalbn(value, scale);
        return std::isfinite(out) && (value == 0.0 || out != 0.0);
    };
    std::array<Vec3, 4> q{};
    for (std::size_t i = 0; i < 4; ++i)
        if (!scaled(p[i].x, q[i].x) || !scaled(p[i].y, q[i].y) || !scaled(p[i].z, q[i].z))
            return false;
    const auto ax = difference(q[1].x, q[0].x), ay = difference(q[1].y, q[0].y),
               az = difference(q[1].z, q[0].z);
    const auto bx = difference(q[2].x, q[0].x), by = difference(q[2].y, q[0].y),
               bz = difference(q[2].z, q[0].z);
    const auto cx = difference(q[3].x, q[0].x), cy = difference(q[3].y, q[0].y),
               cz = difference(q[3].z, q[0].z);
    const auto yz = crossTerm(by, cz, bz, cy), zx = crossTerm(bz, cx, bx, cz),
               xy = crossTerm(bx, cy, by, cx);
    result = sumExpansion(sumExpansion(productExpansion(ax, yz), productExpansion(ay, zx)),
                          productExpansion(az, xy));
    if (result.terms.size() > kT5SixVolumeExpansionTermsPerTet)
        return false;
    const int undo = -3 * scale;
    for (double &value : result.terms) {
        const double original = value;
        value = std::scalbn(value, undo);
        if (!std::isfinite(value) || (original != 0.0 && value == 0.0))
            return false;
    }
    work += result.terms.size();
    return true;
}
bool volumeInterval(const Expansion &expansion, double &lower, double &upper) noexcept {
    lower = 0;
    upper = 0;
    for (double value : expansion.terms) {
        lower = std::nextafter(lower + value, -std::numeric_limits<double>::infinity());
        upper = std::nextafter(upper + value, std::numeric_limits<double>::infinity());
    }
    lower = std::nextafter(lower / 6.0, -std::numeric_limits<double>::infinity());
    upper = std::nextafter(upper / 6.0, std::numeric_limits<double>::infinity());
    return std::isfinite(lower) && std::isfinite(upper);
}
std::uint64_t expansionDigest(const Expansion &expansion) noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    for (double v : expansion.terms) {
        auto bits = std::bit_cast<std::uint64_t>(v);
        for (int i = 0; i < 8; ++i) {
            h ^= static_cast<std::uint8_t>(bits >> (8 * i));
            h *= 1099511628211ULL;
        }
    }
    return h;
}

std::size_t meshAccounted(const mesh::MeshCapacitySnapshot &s, bool &ok) noexcept {
    if (!s.complete) {
        ok = false;
        return (std::numeric_limits<std::size_t>::max)();
    }
    std::size_t total = 0;
    ok = addBytes(total, s.nodeCoordinateCapacity, sizeof(Vec3)) &&
         addBytes(total, s.nodeClassificationCapacity, sizeof(mesh::Classification)) &&
         addBytes(total, s.cellCapacity, sizeof(mesh::Cell)) &&
         addBytes(total, s.connectivityCapacity, sizeof(Index)) &&
         addBytes(total, s.dedupBucketCount, sizeof(void *));
    const std::size_t entryWidth =
        sizeof(std::uint64_t) + sizeof(std::vector<Index>) + 2 * sizeof(void *);
    ok = ok && addBytes(total, s.dedupEntryCount, entryWidth) &&
         addBytes(total, s.dedupIndexCapacity, sizeof(Index));
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

template <class T> bool vectorBytes(std::size_t &total, const std::vector<T> &v) noexcept {
    return addBytes(total, v.capacity(), sizeof(T));
}

template <class K, class V, class C, class A>
bool mapBytes(std::size_t &total, const std::map<K, V, C, A> &values) noexcept {
    return addBytes(total, values.size(),
                    sizeof(typename std::map<K, V, C, A>::value_type) + 3U * sizeof(void *));
}

bool chainBytes(std::size_t &total, const SurfaceChain &chain) noexcept {
    return vectorBytes(total, chain.uses) && vectorBytes(total, chain.segmentEntityIds);
}

bool surfaceFaceBytes(std::size_t &total, const SurfaceFaceResult &face) noexcept {
    bool ok = checkedAdd(total, sizeof(SurfaceFaceResult), total) &&
              vectorBytes(total, face.chartCertificate.certifiedRegionBoundary) &&
              vectorBytes(total, face.localNodes) && vectorBytes(total, face.localParams) &&
              vectorBytes(total, face.localMappingParams) &&
              vectorBytes(total, face.localHandles) && vectorBytes(total, face.localCells) &&
              mapBytes(total, face.existingGlobalIndex) &&
              vectorBytes(total, face.publishedGlobalByLocalNode) &&
              vectorBytes(total, face.publishedGlobalCellByLocalCell) &&
              vectorBytes(total, face.boundaryLoops) && vectorBytes(total, face.normalizedChains) &&
              mapBytes(total, face.segmentOrigins) && vectorBytes(total, face.persistentEdges) &&
              vectorBytes(total, face.protectedRecovery.processed) &&
              vectorBytes(total, face.protectedRecovery.unattempted) &&
              vectorBytes(total, face.localRecovery.processed) &&
              vectorBytes(total, face.localRecovery.unattempted);
    for (const auto &loop : face.boundaryLoops)
        ok = ok && vectorBytes(total, loop);
    for (const auto &chain : face.normalizedChains)
        ok = ok && chainBytes(total, chain);
    if (face.localWorkspace)
        for (const auto bytes : {face.localWorkspace->nodes().stats().totalBytes,
                                 face.localWorkspace->cells().stats().totalBytes,
                                 face.localWorkspace->edges().stats().totalBytes,
                                 face.localWorkspace->faces().stats().totalBytes})
            ok = ok && checkedAdd(total, bytes, total);
    return ok;
}

std::size_t surfaceOutcomeBytes(const SurfaceSetOutcome &outcome, bool &ok) noexcept {
    std::size_t total = sizeof(SurfaceSetOutcome);
    ok = true;
    ok = vectorBytes(total, outcome.processed) && vectorBytes(total, outcome.unattempted);
    if (outcome.retainedArtifact) {
        ok = ok && checkedAdd(total, sizeof(SurfaceSetResult), total) &&
             vectorBytes(total, outcome.retainedArtifact->faces) &&
             mapBytes(total, outcome.retainedArtifact->sharedLedger);
        for (const auto &face : outcome.retainedArtifact->faces)
            if (face)
                ok = ok && surfaceFaceBytes(total, *face);
    }
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

std::size_t workspaceBytes(const TetAuthoritativeResult &a, bool &ok) noexcept {
    if (!a.workspace) {
        ok = false;
        return 0;
    }
    std::size_t total = 0;
    for (auto v :
         {a.workspace->nodes().stats().totalBytes, a.workspace->cells().stats().totalBytes,
          a.workspace->edges().stats().totalBytes, a.workspace->faces().stats().totalBytes})
        ok = ok && checkedAdd(total, v, total);
    const auto adjacency = a.workspace->adjacency().capacitySnapshot();
    if (!adjacency.complete) {
        ok = false;
        return 0;
    }
    const std::size_t entryWidth =
        sizeof(std::int32_t) + sizeof(std::uint32_t) + sizeof(std::vector<topo::FaceHandle>) +
        sizeof(std::vector<topo::CellAdjacency::Neighbor>) + 2 * sizeof(void *);
    ok = ok && addBytes(total, adjacency.entryBucketCount, sizeof(void *)) &&
         addBytes(total, adjacency.entryCount, entryWidth) &&
         addBytes(total, adjacency.entryFaceCapacity, sizeof(topo::FaceHandle)) &&
         addBytes(total, adjacency.entryNeighborCapacity, sizeof(topo::CellAdjacency::Neighbor)) &&
         addBytes(total, adjacency.boundaryCapacity, sizeof(topo::CellAdjacency::BoundaryFace)) &&
         addBytes(total, adjacency.interiorCapacity, sizeof(topo::CellAdjacency::InteriorFace)) &&
         addBytes(total, adjacency.nonManifoldCapacity,
                  sizeof(topo::CellAdjacency::NonManifoldFace)) &&
         addBytes(total, adjacency.nonManifoldCellCapacity, sizeof(topo::CellHandle));
    ok = ok && vectorBytes(total, a.nodeBindings) && vectorBytes(total, a.localCellByCandidate) &&
         vectorBytes(total, a.selectedCellBindings) &&
         vectorBytes(total, a.existingGlobalByLocalNodeSlot);
    return total;
}

bool constraintStateBytes(std::size_t &total, const bnd::TetFacetSetState &state) noexcept {
    return addBytes(total, state.edgeBindings.size(),
                    sizeof(decltype(state.edgeBindings)::value_type) + 3 * sizeof(void *)) &&
           addBytes(total, state.bindings.size(),
                    sizeof(decltype(state.bindings)::value_type) + 3 * sizeof(void *)) &&
           vectorBytes(total, state.canonicalEdgeSchedule) &&
           vectorBytes(total, state.canonicalSchedule);
}

std::size_t boundaryResultBytes(const BoundaryConformingV3Result &result, bool &ok) noexcept {
    std::size_t total = sizeof(BoundaryConformingV3Result);
    ok = true;
    if (result.finalAuthoritative) {
        const auto bytes = workspaceBytes(*result.finalAuthoritative, ok);
        ok = ok && checkedAdd(total, bytes, total);
    }
    ok = ok && constraintStateBytes(total, result.constraintState) &&
         vectorBytes(total, result.globalByAuthoritativeNode) &&
         vectorBytes(total, result.globalByAuthoritativeCell);
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

std::size_t boundaryOutcomeBytes(const BoundaryConformingV3Outcome &outcome, bool &ok) noexcept {
    std::size_t total = sizeof(BoundaryConformingV3Outcome);
    ok = true;
    if (outcome.surfaceOutcome) {
        bool surfaceOk = true;
        const auto bytes = surfaceOutcomeBytes(*outcome.surfaceOutcome, surfaceOk);
        ok = ok && surfaceOk && checkedAdd(total, bytes, total);
    }
    if (outcome.report)
        ok = ok && vectorBytes(total, outcome.report->processed) &&
             vectorBytes(total, outcome.report->unattempted);
    if (outcome.retainedResult) {
        bool resultOk = true;
        const auto bytes = boundaryResultBytes(*outcome.retainedResult, resultOk);
        ok = ok && resultOk && checkedAdd(total, bytes, total);
    }
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

std::size_t resultBytes(const T5Result &r, bool &ok) noexcept {
    std::size_t total = sizeof(T5Result);
    ok = true;
    if (r.authoritative)
        ok = checkedAdd(total, workspaceBytes(*r.authoritative, ok), total) && ok;
    ok = ok && constraintStateBytes(total, r.constraintState);
    ok = ok && vectorBytes(total, r.regions) && vectorBytes(total, r.materials) &&
         vectorBytes(total, r.interfaceNodes) && vectorBytes(total, r.interfaceEdges) &&
         vectorBytes(total, r.interfaceFacets) &&
         vectorBytes(total, r.interfaceNodeIncidentFacets) &&
         vectorBytes(total, r.interfaceEdgeIncidentArtifactFacets) &&
         vectorBytes(total, r.regionInterfaceFacets) && vectorBytes(total, r.regionTets) &&
         vectorBytes(total, r.materialRegions) && vectorBytes(total, r.sourceNodeMappings) &&
         vectorBytes(total, r.sourceSurfaceCellMappings) &&
         vectorBytes(total, r.sourceTetMappings) &&
         vectorBytes(total, r.exactVolumeExpansionTerms) &&
         vectorBytes(total, r.authorityCoreBytes);
    if (r.artifact) {
        ok = ok && vectorBytes(total, r.artifact->facets) &&
             vectorBytes(total, r.artifact->edges) && vectorBytes(total, r.artifact->regionTags);
        for (const auto &f : r.artifact->facets)
            ok = ok && vectorBytes(total, f.owners);
        for (const auto &e : r.artifact->edges)
            ok = ok && vectorBytes(total, e.incidentFacets);
    }
    if (r.privateMesh) {
        bool meshOk = true;
        const auto bytes = meshAccounted(r.privateMesh->capacitySnapshot(), meshOk);
        ok = ok && meshOk && checkedAdd(total, bytes, total);
    }
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

std::size_t retainedOutcomeBytes(const T5Outcome &out, bool &ok) noexcept {
    std::size_t total = 0;
    ok = true;
    if (out.retainedResult) {
        const auto bytes = resultBytes(*out.retainedResult, ok);
        ok = ok && checkedAdd(total, bytes, total);
    }
    if (out.report) {
        ok = ok &&
             addBytes(total, out.report->processedStorage.capacity(), sizeof(T5ProcessedRow)) &&
             addBytes(total, out.report->unattemptedStorage.capacity(), sizeof(T5WorkKey));
    }
    if (out.sourceReceipt) {
        ok = ok && checkedAdd(total, sizeof(T5V3SourceReceipt), total);
        if (out.sourceReceipt->surfaceOutcome) {
            bool surfaceOk = true;
            const auto surfaceBytes =
                surfaceOutcomeBytes(*out.sourceReceipt->surfaceOutcome, surfaceOk);
            ok = ok && surfaceOk && checkedAdd(total, surfaceBytes, total);
        }
        if (out.sourceReceipt->sourceReport)
            ok = ok &&
                 addBytes(total, out.sourceReceipt->sourceReport->processed.capacity(),
                          sizeof(BoundaryConformingV3ReportRow)) &&
                 addBytes(total, out.sourceReceipt->sourceReport->unattempted.capacity(),
                          sizeof(BoundaryConformingV3WorkKey));
    }
    return ok ? total : (std::numeric_limits<std::size_t>::max)();
}

struct BudgetStop {
    T5Reason reason;
    T5Phase phase;
};

class AllocationFlagScope {
  public:
    AllocationFlagScope(bool &flag, bool activate) noexcept : flag_(flag), previous_(flag) {
        if (activate)
            flag_ = true;
    }
    ~AllocationFlagScope() { flag_ = previous_; }
    AllocationFlagScope(const AllocationFlagScope &) = delete;
    AllocationFlagScope &operator=(const AllocationFlagScope &) = delete;
    void activate() noexcept { flag_ = true; }
    void restore() noexcept { flag_ = previous_; }

  private:
    bool &flag_;
    bool previous_ = false;
};

bool workKeyLess(const T5WorkKey &a, const T5WorkKey &b) noexcept {
    return std::tie(a.kind, a.region, a.material, a.globalNode, a.globalEdge, a.facet,
                    a.surfaceFace, a.normalizedLocalCellOrdinal, a.sourceSlot) <
           std::tie(b.kind, b.region, b.material, b.globalNode, b.globalEdge, b.facet,
                    b.surfaceFace, b.normalizedLocalCellOrdinal, b.sourceSlot);
}

void digestWord(std::uint64_t &h, std::uint64_t word) noexcept {
    for (unsigned byte = 0; byte < 8; ++byte) {
        h ^= static_cast<std::uint8_t>(word >> (byte * 8U));
        h *= 1099511628211ULL;
    }
}
template <class Handle> void digestHandle(std::uint64_t &h, Handle value) noexcept {
    digestWord(h, static_cast<std::uint64_t>(value.slot));
    digestWord(h, static_cast<std::uint64_t>(value.gen));
    digestWord(h, static_cast<std::uint64_t>(value.pool));
}
void digestClass(std::uint64_t &h, mesh::Classification value) noexcept {
    digestWord(h, static_cast<std::uint64_t>(value.dim));
    digestWord(h, static_cast<std::uint64_t>(value.tag));
}
std::uint64_t mappingDigest(const std::vector<T5SourceNodeMappingRow> &rows) noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    digestWord(h, rows.size());
    for (const auto &row : rows) {
        digestWord(h, row.authoritativeSlot);
        digestHandle(h, row.authoritativeNode);
        digestWord(h, static_cast<std::uint64_t>(row.sourceGlobalNode));
        digestWord(h, static_cast<std::uint64_t>(row.publishedGlobalNode));
    }
    return h;
}
std::uint64_t mappingDigest(const std::vector<T5SourceSurfaceCellMappingRow> &rows) noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    digestWord(h, rows.size());
    for (const auto &row : rows) {
        digestWord(h, row.face);
        digestWord(h, row.normalizedLocalCellOrdinal);
        digestWord(h, static_cast<std::uint64_t>(row.sourceGlobalCell));
        digestWord(h, static_cast<std::uint64_t>(row.publishedGlobalCell));
        digestClass(h, row.classification);
    }
    return h;
}
std::uint64_t mappingDigest(const std::vector<T5SourceTetMappingRow> &rows) noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    digestWord(h, rows.size());
    for (const auto &row : rows) {
        digestWord(h, row.authoritativeSlot);
        digestHandle(h, row.authoritativeCell);
        digestWord(h, static_cast<std::uint64_t>(row.sourceGlobalCell));
        digestWord(h, static_cast<std::uint64_t>(row.publishedGlobalCell));
        digestClass(h, row.classification);
    }
    return h;
}

bool consumeAudit(T5Metrics &metrics, std::uint64_t limit, std::uint64_t amount = 1) noexcept {
    if (amount > limit - metrics.auditWorkUnits)
        return false;
    metrics.auditWorkUnits += amount;
    return true;
}

bool completeSourceWork(const BoundaryConformingV3ConstructionMetrics &construction,
                        const ArtifactPreflightMetrics &initialArtifact,
                        const ArtifactPreflightMetrics &rebuiltArtifact,
                        const TetRecoveryMetrics &recovery, std::uint64_t &total) noexcept {
    total = construction.totalWorkUnits;
    const auto add = [&](auto... values) noexcept {
        return (checkedAdd(total, static_cast<std::uint64_t>(values), total) && ...);
    };
    const auto addArtifact = [&](const ArtifactPreflightMetrics &m) noexcept {
        return add(m.edgeEntries, m.facetEntries, m.treeNodes, m.axisSelectionVisits,
                   m.sortComparisons, m.nodePairVisits, m.aabbTests, m.edgeEdgePairs,
                   m.edgeFacetPairs, m.facetFacetPairs, m.exactEdgeEdgeTests, m.exactEdgeFacetTests,
                   m.exactFacetFacetTests, m.legalSharedContacts, m.rejectedContacts,
                   m.allPairsFallbacks);
    };
    return addArtifact(initialArtifact) && addArtifact(rebuiltArtifact) &&
           add(recovery.contextNodes, recovery.contextCells, recovery.contextEdges,
               recovery.contextFaces, recovery.contextAuditVisits, recovery.indexBuildComparisons,
               recovery.edgeAabbQueries, recovery.facetAabbQueries, recovery.indexNodesVisited,
               recovery.returnedTetCandidates, recovery.segmentTetTests, recovery.targetTetTests,
               recovery.candidateSurvivorTests, recovery.segmentWalkCells,
               recovery.edgeCavityCommits, recovery.predictedEdgeMatches,
               recovery.completeStarQueries, recovery.closureInsertions, recovery.apexPairs,
               recovery.facetCommits, recovery.indexRemovals, recovery.indexInsertions,
               recovery.indexRebalances, recovery.routeAssessments, recovery.dependencyEdges,
               recovery.dependencyInvalidations, recovery.selectedRebuilds,
               recovery.batchPreflightAssessments, recovery.alreadyPresentBatchEntities,
               recovery.protectedImportFallbacks, recovery.protectedConflicts,
               recovery.conflictSccVisits, recovery.fullStageRescans,
               recovery.productionGlobalScanFallbacks);
}

[[maybe_unused]] std::array<Index, 2> canonicalEdge(Index a, Index b) noexcept {
    return a < b ? std::array<Index, 2>{a, b} : std::array<Index, 2>{b, a};
}
std::array<Index, 3> canonicalFace(std::array<Index, 3> value) noexcept {
    std::sort(value.begin(), value.end());
    return value;
}
std::array<Index, 4> canonicalTet(std::array<Index, 4> value) noexcept {
    std::sort(value.begin(), value.end());
    return value;
}

struct Resolvers {
    std::map<int, const T5RegionBinding *> regionByTag;
    std::map<RegionPersistentId, const T5RegionBinding *> regionById;
    std::map<Index, topo::NodeHandle> nodeByGlobal;
    std::map<topo::NodeHandle, Index> globalByNode;
    std::map<topo::CellHandle, const T5RegionBinding *> regionByCell;
    std::map<topo::FaceHandle, const bnd::VolumeFacet *> artifactByFace;
    std::map<topo::FaceHandle, bnd::VolumeFacetKey> keyByFace;
    std::map<topo::EdgeHandle, std::vector<topo::CellHandle>> edgeCells;
    std::map<topo::EdgeHandle, std::vector<topo::FaceHandle>> edgeFaces;
};

bool rebindConstraintState(bnd::TetFacetSetState &state, const bnd::ClosedTri3Artifact &artifact,
                           const std::vector<Index> &globalBySlot, topo::TopologyWorkspace &ws,
                           Resolvers &resolvers) {
    if (state.artifactStamp != artifact.surfaceSetStamp ||
        state.edgeBindings.size() != artifact.edges.size() ||
        state.bindings.size() != artifact.facets.size())
        return false;
    for (const auto node : ws.nodes().activeBySequence()) {
        if (node.slot < 0 || static_cast<std::size_t>(node.slot) >= globalBySlot.size())
            return false;
        const Index global = globalBySlot[static_cast<std::size_t>(node.slot)];
        if (global < 0 || !resolvers.nodeByGlobal.emplace(global, node).second)
            return false;
        resolvers.globalByNode.emplace(node, global);
    }
    std::set<bnd::VolumeArtifactEdgeKey> edgeKeys;
    for (const auto &edge : artifact.edges) {
        if (!edgeKeys.insert(edge.key).second)
            return false;
        auto binding = state.edgeBindings.find(edge.key);
        const auto a = resolvers.nodeByGlobal.find(edge.key.canonicalGlobalNodes[0]);
        const auto b = resolvers.nodeByGlobal.find(edge.key.canonicalGlobalNodes[1]);
        if (binding == state.edgeBindings.end() || a == resolvers.nodeByGlobal.end() ||
            b == resolvers.nodeByGlobal.end() ||
            !sameClass(binding->second.classification, edge.classification))
            return false;
        const auto current = ws.edges().find(ws.nodes(), a->second, b->second);
        if (!current)
            return false;
        const auto data = ws.edges().get(*current);
        if (!data.requiredConstraint || !sameClass(data.cls, edge.classification))
            return false;
        binding->second.edge = *current;
    }
    std::set<bnd::VolumeFacetKey> facetKeys;
    for (const auto &facet : artifact.facets) {
        if (!facetKeys.insert(facet.key).second)
            return false;
        auto binding = state.bindings.find(facet.key);
        if (binding == state.bindings.end() ||
            !sameClass(binding->second.classification, facet.classification))
            return false;
        std::array<topo::NodeHandle, 3> vertices{};
        for (std::size_t i = 0; i < 3; ++i) {
            const auto found = resolvers.nodeByGlobal.find(facet.directedGlobalNodes[i]);
            if (found == resolvers.nodeByGlobal.end())
                return false;
            vertices[i] = found->second;
        }
        const auto current = ws.faces().find(ws.nodes(), vertices);
        if (!current)
            return false;
        const auto data = ws.faces().get(*current);
        if (!sameClass(data.cls, facet.classification) ||
            data.incident.size() != facet.owners.size())
            return false;
        std::sort(vertices.begin(), vertices.end(), [&](const auto a, const auto b) {
            return ws.nodes().seqOf(a) < ws.nodes().seqOf(b);
        });
        binding->second.canonicalVertices = vertices;
        binding->second.recoveredFace = *current;
        if (!resolvers.artifactByFace.emplace(*current, &facet).second)
            return false;
        resolvers.keyByFace.emplace(*current, facet.key);
    }
    if (edgeKeys.size() != state.edgeBindings.size() || facetKeys.size() != state.bindings.size())
        return false;
    state.workspaceStamp = topo::WorkspaceStamp::of(ws);
    return true;
}

bool buildTopologyResolvers(topo::TopologyWorkspace &ws, Resolvers &r, T5Metrics &metrics) {
    for (const auto cell : ws.cells().activeBySequence()) {
        const auto data = ws.cells().get(cell);
        if (data.type != topo::ElementType::Tet4 || data.cls.dim != 3 || data.cls.tag < 0)
            return false;
        const auto region = r.regionByTag.find(data.cls.tag);
        if (region == r.regionByTag.end())
            return false;
        r.regionByCell.emplace(cell, region->second);
        ++metrics.cellClassificationVisits;
        const auto cellEdges = ws.edges().cellEdges(ws.cells(), cell);
        if (cellEdges.size() != 6)
            return false;
        for (const auto edge : cellEdges)
            r.edgeCells[edge].push_back(cell);
    }
    for (const auto face : ws.faces().activeBySequence()) {
        const auto data = ws.faces().get(face);
        if (data.cornerCount != 3 || data.incident.empty() || data.incident.size() > 2)
            return false;
        for (std::size_t e = 0; e < 3; ++e) {
            const auto edge =
                ws.edges().find(ws.nodes(), data.corners[e], data.corners[(e + 1) % 3]);
            if (!edge)
                return false;
            r.edgeFaces[*edge].push_back(face);
        }
        metrics.faceOwnerVisits += data.incident.size();
    }
    for (auto &[edge, cells] : r.edgeCells) {
        (void)edge;
        std::sort(cells.begin(), cells.end(), [&](const auto a, const auto b) {
            return ws.cells().seqOf(a) < ws.cells().seqOf(b);
        });
        cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
    }
    for (auto &[edge, faces] : r.edgeFaces) {
        (void)edge;
        std::sort(faces.begin(), faces.end(), [&](const auto a, const auto b) {
            return ws.faces().seqOf(a) < ws.faces().seqOf(b);
        });
        faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    }
    return true;
}

std::optional<topo::CellHandle> incidenceCell(const topo::FacePool::Incidence &incidence,
                                              const topo::TopologyWorkspace &ws) noexcept {
    const topo::CellHandle cell{incidence.cellSlot, incidence.cellGen, ws.cells().poolId()};
    if (!ws.cells().isCurrent(cell))
        return std::nullopt;
    return cell;
}

bool facetOwnerRows(const bnd::VolumeFacet &facet, const bnd::PersistentTetFacetBinding &binding,
                    const topo::TopologyWorkspace &ws, const Resolvers &resolvers,
                    std::map<int, std::pair<topo::CellHandle, int>> &owners) {
    const auto face = ws.faces().get(binding.recoveredFace);
    const std::array<topo::NodeHandle, 3> stored{face.corners[0], face.corners[1], face.corners[2]};
    const auto canonical = binding.canonicalVertices;
    std::array<topo::NodeHandle, 3> directed{};
    for (std::size_t i = 0; i < 3; ++i) {
        const auto found = resolvers.nodeByGlobal.find(facet.directedGlobalNodes[i]);
        if (found == resolvers.nodeByGlobal.end())
            return false;
        directed[i] = found->second;
    }
    const int storedParity = permutationParity(stored, canonical);
    const int directedParity = permutationParity(directed, canonical);
    if (storedParity == 0 || directedParity == 0)
        return false;
    for (const auto &incidence : face.incident) {
        const auto cell = incidenceCell(incidence, ws);
        if (!cell)
            return false;
        const int outward = incidence.sense * storedParity * directedParity;
        const auto region = resolvers.regionByCell.find(*cell);
        if (region == resolvers.regionByCell.end())
            return false;
        const int tag = region->second->regionTag;
        const auto declared =
            std::find_if(facet.owners.begin(), facet.owners.end(),
                         [&](const auto &owner) { return owner.regionTag == tag; });
        if (declared == facet.owners.end() || declared->directedSense != outward ||
            !owners.emplace(tag, std::make_pair(*cell, outward)).second)
            return false;
    }
    return owners.size() == facet.owners.size();
}

bool faceContainsEdge(const topo::FacePool::FaceData &face, topo::NodeHandle a,
                      topo::NodeHandle b) noexcept {
    bool haveA = false, haveB = false;
    for (std::size_t i = 0; i < face.cornerCount; ++i) {
        haveA = haveA || face.corners[i] == a;
        haveB = haveB || face.corners[i] == b;
    }
    return haveA && haveB;
}

bool auditInterfaceEdgeFan(topo::EdgeHandle edge, const topo::TopologyWorkspace &ws,
                           const Resolvers &r, const std::set<topo::FaceHandle> &interfaceFaces,
                           const std::set<topo::FaceHandle> &exteriorFaces, T5Metrics &metrics,
                           std::uint64_t maxWork) {
    const auto edgeData = ws.edges().get(edge);
    const auto starIt = r.edgeCells.find(edge);
    const auto facesIt = r.edgeFaces.find(edge);
    if (starIt == r.edgeCells.end() || facesIt == r.edgeFaces.end() || starIt->second.empty())
        return false;
    const auto &cells = starIt->second;
    std::map<topo::CellHandle, std::set<topo::CellHandle>> graph;
    for (const auto cell : cells)
        graph.emplace(cell, std::set<topo::CellHandle>{});
    std::size_t transitions = 0, boundaryCount = 0;
    std::optional<std::pair<int, int>> regionPair;
    std::set<topo::FaceHandle> witnessedArtifactFaces;
    for (const auto face : facesIt->second) {
        if (!consumeAudit(metrics, maxWork))
            return false;
        const auto data = ws.faces().get(face);
        if (!faceContainsEdge(data, edgeData.a, edgeData.b))
            return false;
        if (data.incident.size() == 1) {
            ++boundaryCount;
            if (exteriorFaces.count(face) == 0 || interfaceFaces.count(face) != 0)
                return false;
            witnessedArtifactFaces.insert(face);
            const auto cell = incidenceCell(data.incident[0], ws);
            const auto artifact = r.artifactByFace.find(face);
            if (!cell || artifact == r.artifactByFace.end() ||
                artifact->second->owners.size() != 1 ||
                r.regionByCell.at(*cell)->regionTag != artifact->second->owners[0].regionTag)
                return false;
            continue;
        }
        if (data.incident.size() != 2)
            return false;
        const auto first = incidenceCell(data.incident[0], ws);
        const auto second = incidenceCell(data.incident[1], ws);
        if (!first || !second || graph.count(*first) == 0 || graph.count(*second) == 0)
            return false;
        graph[*first].insert(*second);
        graph[*second].insert(*first);
        const int a = r.regionByCell.at(*first)->regionTag;
        const int b = r.regionByCell.at(*second)->regionTag;
        if (a == b) {
            if (interfaceFaces.count(face) != 0 || exteriorFaces.count(face) != 0)
                return false;
        } else {
            ++transitions;
            if (interfaceFaces.count(face) == 0 || exteriorFaces.count(face) != 0)
                return false;
            witnessedArtifactFaces.insert(face);
            auto pair = std::minmax(a, b);
            const std::pair<int, int> ordered{pair.first, pair.second};
            if (regionPair && *regionPair != ordered)
                return false;
            regionPair = ordered;
            const auto artifact = r.artifactByFace.find(face);
            if (artifact == r.artifactByFace.end() || artifact->second->owners.size() != 2)
                return false;
            std::array<int, 2> declared{artifact->second->owners[0].regionTag,
                                        artifact->second->owners[1].regionTag};
            std::sort(declared.begin(), declared.end());
            if (declared != std::array<int, 2>{ordered.first, ordered.second})
                return false;
        }
    }
    // The Tet-edge link must be exactly one path or one cycle.
    std::size_t ends = 0;
    for (const auto &[cell, adjacent] : graph) {
        (void)cell;
        if (adjacent.size() > 2)
            return false;
        if (adjacent.size() == 1)
            ++ends;
        else if (adjacent.empty() && graph.size() != 1)
            return false;
    }
    std::set<topo::CellHandle> visited;
    std::queue<topo::CellHandle> pending;
    pending.push(graph.begin()->first);
    visited.insert(graph.begin()->first);
    while (!pending.empty()) {
        const auto cell = pending.front();
        pending.pop();
        for (const auto next : graph.at(cell))
            if (visited.insert(next).second)
                pending.push(next);
    }
    if (visited.size() != graph.size())
        return false;
    const bool cycle = graph.size() > 1 && ends == 0;
    const bool path = (graph.size() == 1 && ends == 0) || ends == 2;
    if (!cycle && !path)
        return false;
    if (transitions == 0) {
        if (!path || boundaryCount != 2)
            return false;
        const int tag = r.regionByCell.at(graph.begin()->first)->regionTag;
        for (const auto &[cell, adjacent] : graph) {
            (void)adjacent;
            if (r.regionByCell.at(cell)->regionTag != tag)
                return false;
        }
    } else if (cycle) {
        if (boundaryCount != 0 || transitions != 2 || !regionPair)
            return false;
    } else {
        if (boundaryCount != 2 || transitions != 1 || !regionPair)
            return false;
    }
    ++metrics.regionAdjacencyEdges;
    return true;
}

bool checkedRequirements(const T5Request &request, std::size_t regions, std::size_t materials,
                         std::size_t interfaceNodes, std::size_t interfaceEdges,
                         std::size_t interfaceFacets, std::size_t sourceNodes,
                         std::size_t sourceSurfaceCells, std::size_t sourceTets,
                         std::size_t incidentRefs, std::size_t reportRows, std::size_t volumeTerms,
                         BudgetStop &stop) noexcept {
    const auto reportLimit = [&](std::size_t need, std::size_t cap) {
        if (need <= cap)
            return true;
        stop = {T5Reason::ReportBudgetReached, T5Phase::ResultAssembly};
        return false;
    };
    if (!reportLimit(regions, request.maxRegionRows) ||
        !reportLimit(materials, request.maxMaterialRows) ||
        !reportLimit(interfaceNodes, request.maxInterfaceNodeRows) ||
        !reportLimit(interfaceEdges, request.maxInterfaceEdgeRows) ||
        !reportLimit(interfaceFacets, request.maxInterfaceFacetRows) ||
        !reportLimit(sourceNodes, request.maxSourceNodeMappingRows) ||
        !reportLimit(sourceSurfaceCells, request.maxSourceSurfaceCellMappingRows) ||
        !reportLimit(sourceTets, request.maxSourceTetMappingRows) ||
        !reportLimit(incidentRefs, request.maxIncidentReferences) ||
        !reportLimit(reportRows, request.maxReportRows))
        return false;
    if (volumeTerms > request.maxVolumeExpansionTerms) {
        stop = {T5Reason::ResultBudgetReached, T5Phase::RegionAuthority};
        return false;
    }
    return true;
}

void fillKnownUnattempted(T5Report &report, const std::vector<T5RegionBinding> &regions) {
    std::set<MaterialPersistentId> materials;
    for (const auto &row : regions)
        materials.insert(row.materialPersistentId);
    report.unattemptedStorage.clear();
    for (const auto &row : regions) {
        T5WorkKey key;
        key.kind = T5WorkKind::Region;
        key.region = row.regionPersistentId;
        report.unattemptedStorage.push_back(key);
    }
    for (const auto material : materials) {
        T5WorkKey key;
        key.kind = T5WorkKind::Material;
        key.material = material;
        report.unattemptedStorage.push_back(key);
    }
    std::sort(report.unattemptedStorage.begin(), report.unattemptedStorage.end(), workKeyLess);
    report.unattemptedUsed = report.unattemptedStorage.size();
    report.reportComplete = true;
}

bool validRegionBindings(const T5Request &request, std::vector<T5RegionBinding> &normalized,
                         T5Reason &reason) {
    normalized.assign(request.regionBindings.begin(), request.regionBindings.end());
    for (const auto &row : normalized)
        if (row.regionTag < 0 || row.regionPersistentId == 0 || row.materialPersistentId == 0 ||
            row.volumeClassification.dim != 3 || row.volumeClassification.tag != row.regionTag) {
            reason = T5Reason::InvalidRegionBinding;
            return false;
        }
    std::sort(normalized.begin(), normalized.end(), [](const auto &a, const auto &b) {
        return std::tie(a.regionPersistentId, a.regionTag, a.materialPersistentId) <
               std::tie(b.regionPersistentId, b.regionTag, b.materialPersistentId);
    });
    std::set<int> tags;
    std::set<RegionPersistentId> identities;
    for (const auto &row : normalized) {
        if (!tags.insert(row.regionTag).second) {
            reason = T5Reason::DuplicateRegionTag;
            return false;
        }
        if (!identities.insert(row.regionPersistentId).second) {
            reason = T5Reason::DuplicateRegionIdentity;
            return false;
        }
    }
    if (normalized.empty()) {
        reason = T5Reason::InvalidRequest;
        return false;
    }
    reason = T5Reason::None;
    return true;
}

template <class Row> bool reserveExact(std::vector<Row> &vector, std::size_t count) {
    if (count > vector.max_size())
        return false;
    vector.reserve(count);
    return vector.capacity() >= count;
}

} // namespace

#ifdef PAACS_TOPO_TESTING
namespace t5_detail {
ExactSixVolumeProbe exactSixVolumeForTesting(const std::array<Vec3, 4> &points) noexcept {
    ExactSixVolumeProbe probe;
    try {
        Expansion expansion;
        if (!exactSixVolume(points, expansion, probe.workUnits))
            return probe;
        probe.termCount = expansion.terms.size();
        probe.exactDigest = expansionDigest(expansion);
        probe.complete = volumeInterval(expansion, probe.volumeLower, probe.volumeUpper);
    } catch (...) {
        probe.complete = false;
    }
    return probe;
}
} // namespace t5_detail
#endif

T5Outcome
meshModelMultiRegionTet4(const geo::Model &model, const sizing::MeshSpec &spec,
                         const T5Request &request, mesh::Mesh &destination,
                         const VolumeBoundaryOwnershipProvider &ownershipProvider) noexcept {
    T5Outcome out;
    T5Phase activePhase = T5Phase::Preflight;
    AllocationFlagScope routeAllocations(t5_detail::routeAllocationActive(), true);
    AllocationFlagScope publicationAllocations(t5_detail::publicationPreparationActive(), false);
    try {
        out.report = std::make_unique<T5Report>();
    } catch (const std::bad_alloc &) {
        fail(out, T5Status::Failed, T5Reason::AllocationFailed, T5Phase::Preflight);
        return out;
    }
    try {
        std::vector<T5RegionBinding> normalizedRegions;
        T5Reason preflightReason = T5Reason::None;
        if (!validRegionBindings(request, normalizedRegions, preflightReason) ||
            request.sourceRequest.expectedSurfaceSetStamp.persistentId == 0 ||
            request.sourceRequest.expectedOwnershipProviderStamp.persistentId == 0) {
            fail(out, T5Status::Rejected,
                 preflightReason == T5Reason::None ? T5Reason::InvalidRequest : preflightReason,
                 T5Phase::Preflight);
            return out;
        }
        if (std::fegetround() != FE_TONEAREST) {
            fail(out, T5Status::Unsupported, T5Reason::FloatingPointEnvironmentUnsupported,
                 T5Phase::Preflight);
            return out;
        }
        std::set<MaterialPersistentId> inputMaterials;
        for (const auto &row : normalizedRegions)
            inputMaterials.insert(row.materialPersistentId);
        std::size_t inputRows = 0;
        if (!checkedAdd(normalizedRegions.size(), inputMaterials.size(), inputRows) ||
            normalizedRegions.size() > request.maxRegionRows ||
            inputMaterials.size() > request.maxMaterialRows || inputRows > request.maxReportRows) {
            fail(out, T5Status::Incomplete, T5Reason::ReportBudgetReached, T5Phase::Preflight);
            return out;
        }
        out.report->unattemptedStorage.reserve(inputRows);
        fillKnownUnattempted(*out.report, normalizedRegions);
        out.reportComplete = true;

        const auto entryProviderStamp = ownershipProvider.stamp();
        if (entryProviderStamp != request.sourceRequest.expectedOwnershipProviderStamp) {
            fail(out, T5Status::Rejected, T5Reason::StaleOwnershipProvider, T5Phase::Preflight);
            return out;
        }

        activePhase = T5Phase::SourceConstruction;
        if (ownershipProvider.stamp() != entryProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider,
                 T5Phase::SourceConstruction);
            return out;
        }
        auto sourceMesh = std::make_unique<mesh::Mesh>();
        out.sourceOutcome = std::make_unique<BoundaryConformingV3Outcome>();
        *out.sourceOutcome = meshModelBoundaryConforming(model, spec, request.sourceRequest,
                                                         *sourceMesh, ownershipProvider);
        auto &source = *out.sourceOutcome;
        if (ownershipProvider.stamp() != entryProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider,
                 T5Phase::SourceValidation);
            return out;
        }
        if (!source.published()) {
            out.status = sourceStatus(source.status);
            out.failure.reason = sourceReason(source.status);
            out.failure.phase = T5Phase::SourceConstruction;
            out.failure.sourceStatus = source.status;
            out.failure.sourceFailure = source.failure;
            out.reportComplete = out.report->reportComplete;
            return out;
        }

        activePhase = T5Phase::SourceValidation;
        if (!source.retainedResult || !source.report || !source.reportComplete ||
            !source.report->reportComplete || !source.surfaceOutcome ||
            source.surfaceOutcome->setStatus != SurfaceSetStatus::Published ||
            !source.surfaceOutcome->reportComplete || !source.surfaceOutcome->retainedArtifact) {
            fail(out, T5Status::Failed,
                 !source.reportComplete || !source.report ||
                         (source.report && !source.report->reportComplete)
                     ? T5Reason::SourceReportIncomplete
                     : T5Reason::MalformedSourceOutcome,
                 T5Phase::SourceValidation);
            return out;
        }
        auto &sourceResult = *source.retainedResult;
        if (!sourceResult.finalAuthoritative || !sourceResult.finalAuthoritative->workspace ||
            sourceResult.constraintState.artifactStamp !=
                request.sourceRequest.expectedSurfaceSetStamp ||
            sourceResult.globalByAuthoritativeNode.empty() ||
            sourceResult.globalByAuthoritativeCell.size() !=
                sourceResult.finalAuthoritative->selectedCellBindings.size()) {
            fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                 T5Phase::SourceValidation);
            return out;
        }
        if (source.surfaceOutcome->retainedArtifact->setStamp !=
            request.sourceRequest.expectedSurfaceSetStamp) {
            fail(out, T5Status::Rejected, T5Reason::StaleSurfaceArtifact,
                 T5Phase::SourceValidation);
            return out;
        }

        activePhase = T5Phase::ArtifactRevalidation;
        if (ownershipProvider.stamp() != entryProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider,
                 T5Phase::ArtifactRevalidation);
            return out;
        }
        auto rebuilt = buildClosedTri3Artifact(
            *source.surfaceOutcome->retainedArtifact, *sourceMesh,
            request.sourceRequest.expectedSurfaceSetStamp,
            request.sourceRequest.expectedOwnershipProviderStamp, ownershipProvider);
        if (!rebuilt.artifact) {
            out.status = sourceStatus(rebuilt.status);
            out.failure.reason =
                rebuilt.failure.reason == BoundaryConformingV3Reason::AllocationFailed
                    ? T5Reason::AllocationFailed
                    : T5Reason::ArtifactStateMismatch;
            out.failure.phase = T5Phase::ArtifactRevalidation;
            out.failure.sourceStatus = rebuilt.status;
            out.failure.sourceFailure = rebuilt.failure;
            return out;
        }
        if (ownershipProvider.stamp() != entryProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider,
                 T5Phase::ArtifactRevalidation);
            return out;
        }
        auto artifact = std::move(rebuilt.artifact);
        std::set<int> artifactTags(artifact->regionTags.begin(), artifact->regionTags.end());
        std::set<int> bindingTags;
        for (const auto &row : normalizedRegions)
            bindingTags.insert(row.regionTag);
        if (artifactTags != bindingTags) {
            bool unknown = false;
            for (const auto tag : artifactTags)
                unknown = unknown || !bindingTags.count(tag);
            fail(out, T5Status::Rejected,
                 unknown ? T5Reason::UnknownArtifactRegion : T5Reason::UnusedRegionBinding,
                 T5Phase::RegionAuthority);
            return out;
        }

        T5Metrics routeMetrics;
        if (!completeSourceWork(sourceResult.constructionMetrics,
                                sourceResult.artifactPreflightMetrics, rebuilt.metrics,
                                sourceResult.recoveryMetrics, routeMetrics.sourceWorkUnits)) {
            fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::SourceValidation);
            return out;
        }
        routeMetrics.nestedV3PeakAccountedBytes =
            sourceResult.memoryMetrics.peakAccountedOverlapBytes;
        Resolvers resolvers;
        for (const auto &row : normalizedRegions) {
            resolvers.regionByTag.emplace(row.regionTag, &row);
            resolvers.regionById.emplace(row.regionPersistentId, &row);
        }
        bnd::TetFacetSetState rebound = sourceResult.constraintState;
        auto &authoritative = *sourceResult.finalAuthoritative;
        auto &ws = *authoritative.workspace;
        if (!rebindConstraintState(rebound, *artifact, sourceResult.globalByAuthoritativeNode, ws,
                                   resolvers) ||
            !buildTopologyResolvers(ws, resolvers, routeMetrics)) {
            fail(out, T5Status::Failed, T5Reason::ArtifactStateMismatch,
                 T5Phase::ArtifactRevalidation);
            return out;
        }

        std::map<std::array<Index, 3>, Index> sourceTriByKey;
        for (Index cell = 0; cell < sourceMesh->cellCount(); ++cell) {
            if (sourceMesh->cell(cell).type != mesh::CellType::Tri3)
                continue;
            const auto nodes = sourceMesh->cellNodes(cell);
            if (!sourceTriByKey.emplace(canonicalFace({nodes[0], nodes[1], nodes[2]}), cell)
                     .second) {
                fail(out, T5Status::Failed, T5Reason::MalformedSurfaceOutcome,
                     T5Phase::SourceValidation);
                return out;
            }
        }
        std::vector<SurfacePlan> surfacePlans;
        surfacePlans.reserve(artifact->facets.size());
        std::map<bnd::VolumeFacetKey, const bnd::VolumeFacet *> facetByKey;
        std::map<bnd::VolumeArtifactEdgeKey, const bnd::VolumeArtifactEdge *> artifactEdgeByKey;
        std::set<bnd::VolumeFacetKey> interfaceFacetKeys;
        std::set<Index> interfaceGlobalNodes;
        for (const auto &facet : artifact->facets) {
            if (!facetByKey.emplace(facet.key, &facet).second) {
                fail(out, T5Status::Failed, T5Reason::ArtifactStateMismatch,
                     T5Phase::ArtifactRevalidation);
                return out;
            }
            const auto sourceCell = sourceTriByKey.find(canonicalFace(facet.directedGlobalNodes));
            if (sourceCell == sourceTriByKey.end() ||
                !sameClass(sourceMesh->cell(sourceCell->second).onEntity, facet.classification)) {
                fail(out, T5Status::Failed, T5Reason::MalformedSurfaceOutcome,
                     T5Phase::SourceValidation);
                return out;
            }
            surfacePlans.push_back(
                {facet.key, facet.directedGlobalNodes, facet.classification, sourceCell->second});
            if (facet.owners.size() == 2) {
                interfaceFacetKeys.insert(facet.key);
                interfaceGlobalNodes.insert(facet.directedGlobalNodes.begin(),
                                            facet.directedGlobalNodes.end());
            }
        }
        if (sourceTriByKey.size() != surfacePlans.size()) {
            fail(out, T5Status::Failed, T5Reason::MalformedSurfaceOutcome,
                 T5Phase::SourceValidation);
            return out;
        }
        std::vector<const bnd::VolumeArtifactEdge *> interfaceArtifactEdges;
        std::size_t edgeFacetReferences = 0;
        for (const auto &edge : artifact->edges) {
            if (!artifactEdgeByKey.emplace(edge.key, &edge).second) {
                fail(out, T5Status::Failed, T5Reason::ArtifactStateMismatch,
                     T5Phase::ArtifactRevalidation);
                return out;
            }
            bool isInterface = false;
            for (const auto key : edge.incidentFacets) {
                const auto facet = facetByKey.find(key);
                if (facet == facetByKey.end()) {
                    fail(out, T5Status::Failed, T5Reason::ArtifactStateMismatch,
                         T5Phase::ArtifactRevalidation);
                    return out;
                }
                isInterface = isInterface || facet->second->owners.size() == 2;
            }
            if (isInterface) {
                interfaceArtifactEdges.push_back(&edge);
                if (!checkedAdd(edgeFacetReferences, edge.incidentFacets.size(),
                                edgeFacetReferences)) {
                    fail(out, T5Status::Incomplete, T5Reason::ReportBudgetReached,
                         T5Phase::ResultAssembly);
                    return out;
                }
            }
        }

        std::vector<TetPlan> tetPlans;
        tetPlans.reserve(authoritative.selectedCellBindings.size());
        for (std::size_t ordinal = 0; ordinal < authoritative.selectedCellBindings.size();
             ++ordinal) {
            const auto cell = authoritative.selectedCellBindings[ordinal].local;
            if (!ws.cells().isCurrent(cell)) {
                fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                     T5Phase::SourceValidation);
                return out;
            }
            const auto data = ws.cells().get(cell);
            const auto region = resolvers.regionByTag.find(data.cls.tag);
            if (data.type != topo::ElementType::Tet4 || data.count != 4 ||
                region == resolvers.regionByTag.end() ||
                !sameClass(data.cls, region->second->volumeClassification)) {
                fail(out, T5Status::Failed, T5Reason::RegionCellMismatch, T5Phase::RegionAuthority);
                return out;
            }
            std::array<Index, 4> nodes{};
            for (std::size_t i = 0; i < 4; ++i) {
                const auto global = resolvers.globalByNode.find(data.nodes[i]);
                if (global == resolvers.globalByNode.end() || global->second < 0 ||
                    global->second >= sourceMesh->nodeCount()) {
                    fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                         T5Phase::SourceValidation);
                    return out;
                }
                nodes[i] = global->second;
            }
            const Index sourceCell = sourceResult.globalByAuthoritativeCell[ordinal];
            if (sourceCell < 0 || sourceCell >= sourceMesh->cellCount() ||
                sourceMesh->cell(sourceCell).type != mesh::CellType::Tet4 ||
                !sameClass(sourceMesh->cell(sourceCell).onEntity, data.cls)) {
                fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                     T5Phase::SourceValidation);
                return out;
            }
            const Index *sourceNodes = sourceMesh->cellNodes(sourceCell);
            if (std::array<Index, 4>{sourceNodes[0], sourceNodes[1], sourceNodes[2],
                                     sourceNodes[3]} != nodes) {
                fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                     T5Phase::SourceValidation);
                return out;
            }
            tetPlans.push_back(
                {cell, nodes, data.cls, sourceCell, region->second->regionPersistentId});
        }
        std::sort(tetPlans.begin(), tetPlans.end(), [](const auto &a, const auto &b) {
            const auto ac = canonicalTet(a.nodes), bc = canonicalTet(b.nodes);
            return std::tie(a.region, ac, a.nodes) < std::tie(b.region, bc, b.nodes);
        });

        const std::size_t regionCount = normalizedRegions.size();
        const std::size_t materialCount = inputMaterials.size();
        const std::size_t interfaceNodeCount = interfaceGlobalNodes.size();
        const std::size_t interfaceEdgeCount = interfaceArtifactEdges.size();
        const std::size_t interfaceFacetCount = interfaceFacetKeys.size();
        const std::size_t sourceNodeCount = resolvers.globalByNode.size();
        const std::size_t sourceSurfaceCount = surfacePlans.size();
        const std::size_t sourceTetCount = tetPlans.size();
        std::size_t reportRows = 0, publicationCells = 0, publicationConnectivity = 0,
                    incidentReferences = 0, volumeTerms = 0;
        bool arithmeticOk =
            checkedAdd(regionCount, materialCount, reportRows) &&
            checkedAdd(reportRows, interfaceNodeCount, reportRows) &&
            checkedAdd(reportRows, interfaceEdgeCount, reportRows) &&
            checkedAdd(reportRows, interfaceFacetCount, reportRows) &&
            checkedAdd(reportRows, sourceNodeCount, reportRows) &&
            checkedAdd(reportRows, sourceSurfaceCount, reportRows) &&
            checkedAdd(reportRows, sourceTetCount, reportRows) &&
            checkedAdd(surfacePlans.size(), tetPlans.size(), publicationCells) &&
            checkedMul(surfacePlans.size(), std::size_t{3}, publicationConnectivity);
        std::size_t tetConnectivity = 0, threeFacets = 0, twoFacets = 0;
        arithmeticOk =
            arithmeticOk && checkedMul(tetPlans.size(), std::size_t{4}, tetConnectivity) &&
            checkedAdd(publicationConnectivity, tetConnectivity, publicationConnectivity) &&
            checkedMul(interfaceFacetCount, std::size_t{3}, threeFacets) &&
            checkedMul(interfaceFacetCount, std::size_t{2}, twoFacets) &&
            checkedAdd(threeFacets, edgeFacetReferences, incidentReferences) &&
            checkedAdd(incidentReferences, twoFacets, incidentReferences) &&
            checkedAdd(incidentReferences, tetPlans.size(), incidentReferences) &&
            checkedAdd(incidentReferences, regionCount, incidentReferences) &&
            checkedMul(sourceTetCount, kT5SixVolumeExpansionTermsPerTet, volumeTerms) &&
            checkedMul(volumeTerms, std::size_t{2}, volumeTerms);
        if (!arithmeticOk) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::ResultAssembly);
            return out;
        }
        BudgetStop budget{};
        if (!checkedRequirements(request, regionCount, materialCount, interfaceNodeCount,
                                 interfaceEdgeCount, interfaceFacetCount, sourceNodeCount,
                                 sourceSurfaceCount, sourceTetCount, incidentReferences, reportRows,
                                 volumeTerms, budget)) {
            out.report->reportComplete = false;
            out.reportComplete = false;
            fail(out, T5Status::Incomplete, budget.reason, budget.phase);
            return out;
        }
        activePhase = T5Phase::ResultAssembly;
        if (!reserveExact(out.report->processedStorage, reportRows) ||
            !reserveExact(out.report->unattemptedStorage, reportRows)) {
            out.report->reportComplete = false;
            out.reportComplete = false;
            fail(out, T5Status::Incomplete, T5Reason::ReportBudgetReached, T5Phase::ResultAssembly);
            return out;
        }
        std::vector<T5WorkKey> expectedKeys;
        expectedKeys.reserve(reportRows);
        for (const auto &row : normalizedRegions) {
            T5WorkKey key;
            key.kind = T5WorkKind::Region;
            key.region = row.regionPersistentId;
            expectedKeys.push_back(key);
        }
        for (const auto material : inputMaterials) {
            T5WorkKey key;
            key.kind = T5WorkKind::Material;
            key.material = material;
            expectedKeys.push_back(key);
        }
        for (const auto global : interfaceGlobalNodes) {
            T5WorkKey key;
            key.kind = T5WorkKind::InterfaceNode;
            key.globalNode = global;
            expectedKeys.push_back(key);
        }
        for (const auto *edge : interfaceArtifactEdges) {
            T5WorkKey key;
            key.kind = T5WorkKind::InterfaceEdge;
            key.globalEdge = edge->key.canonicalGlobalNodes;
            expectedKeys.push_back(key);
        }
        for (const auto keyValue : interfaceFacetKeys) {
            T5WorkKey key;
            key.kind = T5WorkKind::InterfaceFacet;
            key.facet = keyValue;
            expectedKeys.push_back(key);
        }
        for (const auto &[node, global] : resolvers.globalByNode) {
            T5WorkKey key;
            key.kind = T5WorkKind::SourceNodeMapping;
            key.globalNode = global;
            key.sourceSlot = static_cast<std::uint64_t>(node.slot);
            expectedKeys.push_back(key);
        }
        for (const auto &plan : surfacePlans) {
            T5WorkKey key;
            key.kind = T5WorkKind::SourceSurfaceCellMapping;
            key.surfaceFace = plan.key.facePersistentId;
            key.normalizedLocalCellOrdinal = plan.key.normalizedLocalCellOrdinal;
            expectedKeys.push_back(key);
        }
        for (const auto &plan : tetPlans) {
            T5WorkKey key;
            key.kind = T5WorkKind::SourceTetMapping;
            key.sourceSlot = static_cast<std::uint64_t>(plan.cell.slot);
            expectedKeys.push_back(key);
        }
        std::sort(expectedKeys.begin(), expectedKeys.end(), workKeyLess);
        if (expectedKeys.size() != reportRows ||
            std::adjacent_find(expectedKeys.begin(), expectedKeys.end(),
                               [](const auto &a, const auto &b) {
                                   return !workKeyLess(a, b) && !workKeyLess(b, a);
                               }) != expectedKeys.end()) {
            fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::ResultAssembly);
            return out;
        }
        out.report->unattemptedStorage = expectedKeys;
        out.report->unattemptedUsed = expectedKeys.size();
        out.report->reportComplete = true;
        out.reportComplete = true;

        auto working = std::make_unique<T5Result>();
        auto &result = *working;
        result.metrics = routeMetrics;
        result.metrics.regionRows = regionCount;
        result.regions.reserve(regionCount);
        result.materials.reserve(materialCount);
        result.interfaceNodes.reserve(interfaceNodeCount);
        result.interfaceEdges.reserve(interfaceEdgeCount);
        result.interfaceFacets.reserve(interfaceFacetCount);
        result.interfaceNodeIncidentFacets.reserve(threeFacets);
        result.interfaceEdgeIncidentArtifactFacets.reserve(edgeFacetReferences);
        result.regionInterfaceFacets.reserve(twoFacets);
        result.regionTets.reserve(tetPlans.size());
        result.materialRegions.reserve(regionCount);
        result.sourceNodeMappings.reserve(sourceNodeCount);
        result.sourceSurfaceCellMappings.reserve(sourceSurfaceCount);
        result.sourceTetMappings.reserve(sourceTetCount);
        result.exactVolumeExpansionTerms.reserve(volumeTerms);

        std::map<Index, std::vector<bnd::VolumeFacetKey>> nodeIncident;
        std::set<topo::FaceHandle> interfaceFaces, exteriorFaces;
        std::map<bnd::VolumeFacetKey, std::size_t> surfacePlanByKey;
        for (std::size_t i = 0; i < surfacePlans.size(); ++i)
            surfacePlanByKey.emplace(surfacePlans[i].key, i);

        activePhase = T5Phase::InterfaceAudit;
        if (ownershipProvider.stamp() != entryProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider, T5Phase::InterfaceAudit);
            return out;
        }
        for (const auto &facet : artifact->facets) {
            if (!consumeAudit(result.metrics, request.maxAuditWork)) {
                fail(out, T5Status::Incomplete, T5Reason::WorkBudgetReached,
                     T5Phase::InterfaceAudit);
                return out;
            }
            const auto stateBinding = rebound.bindings.find(facet.key);
            if (stateBinding == rebound.bindings.end()) {
                fail(out, T5Status::Failed, T5Reason::ArtifactStateMismatch,
                     T5Phase::ArtifactRevalidation);
                return out;
            }
            std::map<int, std::pair<topo::CellHandle, int>> owners;
            if (!facetOwnerRows(facet, stateBinding->second, ws, resolvers, owners)) {
                ++result.audit.wrongOwnerSense;
                fail(out, T5Status::Failed, T5Reason::InterfaceOwnerMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            if (facet.owners.size() == 1) {
                exteriorFaces.insert(stateBinding->second.recoveredFace);
                ++result.audit.prescribedExterior;
                ++result.audit.matchedExterior;
                continue;
            }
            if (facet.owners.size() != 2) {
                fail(out, T5Status::Rejected, T5Reason::InvalidInterfaceOwnership,
                     T5Phase::InterfaceAudit);
                return out;
            }
            interfaceFaces.insert(stateBinding->second.recoveredFace);
            ++result.audit.prescribedInterface;
            const auto &a = *resolvers.regionByTag.at(facet.owners[0].regionTag);
            const auto &b = *resolvers.regionByTag.at(facet.owners[1].regionTag);
            const T5RegionBinding *first = &a;
            const T5RegionBinding *second = &b;
            if (second->regionPersistentId < first->regionPersistentId)
                std::swap(first, second);
            const auto firstOwner = owners.find(first->regionTag);
            const auto secondOwner = owners.find(second->regionTag);
            if (firstOwner == owners.end() || secondOwner == owners.end() ||
                firstOwner->second.second != -secondOwner->second.second) {
                ++result.audit.wrongOwnerSense;
                fail(out, T5Status::Failed, T5Reason::InterfaceOwnerMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            const auto firstCell = ws.cells().get(firstOwner->second.first);
            const auto secondCell = ws.cells().get(secondOwner->second.first);
            if (!sameClass(firstCell.cls, first->volumeClassification) ||
                !sameClass(secondCell.cls, second->volumeClassification)) {
                ++result.audit.wrongRegionClassification;
                fail(out, T5Status::Failed, T5Reason::RegionClassificationMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            T5InterfaceFacetRow row;
            row.facet = facet.key;
            row.authoritativeFace = stateBinding->second.recoveredFace;
            row.classification = facet.classification;
            row.firstRegion = first->regionPersistentId;
            row.secondRegion = second->regionPersistentId;
            row.firstMaterial = first->materialPersistentId;
            row.secondMaterial = second->materialPersistentId;
            row.firstOutwardSense = firstOwner->second.second;
            row.secondOutwardSense = secondOwner->second.second;
            row.firstOwner = firstOwner->second.first;
            row.secondOwner = secondOwner->second.first;
            row.canonicalPublishedNodes = facet.directedGlobalNodes;
            if (row.firstOutwardSense < 0)
                std::swap(row.canonicalPublishedNodes[1], row.canonicalPublishedNodes[2]);
            const auto plan = surfacePlanByKey.find(facet.key);
            if (plan == surfacePlanByKey.end()) {
                fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::InterfaceAudit);
                return out;
            }
            surfacePlans[plan->second].nodes = row.canonicalPublishedNodes;
            result.interfaceFacets.push_back(row);
            for (const auto node : facet.directedGlobalNodes)
                nodeIncident[node].push_back(facet.key);
            ++result.audit.matchedInterface;
            ++result.metrics.interfaceFacetInsertions;
        }

        // One global face pass proves every actual exterior and every
        // inter-region transition is named by exactly one artifact row.
        std::set<topo::FaceHandle> actualArtifactFaces;
        for (const auto face : ws.faces().activeBySequence()) {
            if (!consumeAudit(result.metrics, request.maxAuditWork)) {
                fail(out, T5Status::Incomplete, T5Reason::WorkBudgetReached,
                     T5Phase::InterfaceAudit);
                return out;
            }
            const auto data = ws.faces().get(face);
            if (data.incident.empty() || data.incident.size() > 2) {
                ++result.audit.undeclaredRegionFaces;
                fail(out, T5Status::Unsupported, T5Reason::UnsupportedRegionJunction,
                     T5Phase::ArtifactRevalidation);
                return out;
            }
            std::array<int, 2> tags{-1, -1};
            for (std::size_t i = 0; i < data.incident.size(); ++i) {
                const auto owner = incidenceCell(data.incident[i], ws);
                if (!owner || resolvers.regionByCell.count(*owner) == 0) {
                    fail(out, T5Status::Failed, T5Reason::RegionCellMismatch,
                         T5Phase::InterfaceAudit);
                    return out;
                }
                tags[i] = resolvers.regionByCell.at(*owner)->regionTag;
            }
            const auto artifactFace = resolvers.artifactByFace.find(face);
            if (data.incident.size() == 1) {
                if (artifactFace == resolvers.artifactByFace.end() ||
                    artifactFace->second->owners.size() != 1 ||
                    artifactFace->second->owners[0].regionTag != tags[0]) {
                    ++result.audit.unexpectedExterior;
                    fail(out, T5Status::Failed, T5Reason::BoundaryBijectionFailed,
                         T5Phase::BoundaryAudit);
                    return out;
                }
                actualArtifactFaces.insert(face);
            } else if (tags[0] != tags[1]) {
                if (artifactFace == resolvers.artifactByFace.end() ||
                    artifactFace->second->owners.size() != 2) {
                    ++result.audit.unexpectedInterface;
                    fail(out, T5Status::Failed, T5Reason::BoundaryBijectionFailed,
                         T5Phase::BoundaryAudit);
                    return out;
                }
                std::array<int, 2> declared{artifactFace->second->owners[0].regionTag,
                                            artifactFace->second->owners[1].regionTag};
                std::sort(declared.begin(), declared.end());
                std::sort(tags.begin(), tags.end());
                if (declared != tags) {
                    fail(out, T5Status::Failed, T5Reason::InterfaceOwnerMismatch,
                         T5Phase::InterfaceAudit);
                    return out;
                }
                actualArtifactFaces.insert(face);
            } else if (artifactFace != resolvers.artifactByFace.end()) {
                ++result.audit.undeclaredRegionFaces;
                fail(out, T5Status::Failed, T5Reason::InterfaceOwnerMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
        }
        if (actualArtifactFaces.size() != artifact->facets.size()) {
            result.audit.orphanExterior =
                result.audit.prescribedExterior - result.audit.matchedExterior;
            result.audit.missingInterface =
                result.audit.prescribedInterface - result.audit.matchedInterface;
            fail(out, T5Status::Failed, T5Reason::BoundaryBijectionFailed, T5Phase::BoundaryAudit);
            return out;
        }

        for (const auto global : interfaceGlobalNodes) {
            const auto local = resolvers.nodeByGlobal.find(global);
            if (local == resolvers.nodeByGlobal.end() || global < 0 ||
                global >= sourceMesh->nodeCount()) {
                fail(out, T5Status::Failed, T5Reason::InterfaceIdentityMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            auto &incident = nodeIncident[global];
            std::sort(incident.begin(), incident.end());
            incident.erase(std::unique(incident.begin(), incident.end()), incident.end());
            T5InterfaceNodeRow row;
            row.globalNode = global;
            row.authoritativeNode = local->second;
            row.classification = ws.nodes().get(local->second).cls;
            row.incidentFacetOffset = result.interfaceNodeIncidentFacets.size();
            row.incidentFacetCount = incident.size();
            result.interfaceNodeIncidentFacets.insert(result.interfaceNodeIncidentFacets.end(),
                                                      incident.begin(), incident.end());
            result.interfaceNodes.push_back(row);
            ++result.metrics.interfaceNodeInsertions;
        }

        for (const auto *artifactEdge : interfaceArtifactEdges) {
            const auto binding = rebound.edgeBindings.find(artifactEdge->key);
            if (binding == rebound.edgeBindings.end() ||
                !ws.edges().isCurrent(binding->second.edge)) {
                fail(out, T5Status::Failed, T5Reason::InterfaceIdentityMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            const auto edge = ws.edges().get(binding->second.edge);
            if (!edge.requiredConstraint || !sameClass(edge.cls, artifactEdge->classification)) {
                fail(out, T5Status::Failed, T5Reason::InterfaceClassificationMismatch,
                     T5Phase::InterfaceAudit);
                return out;
            }
            if (!auditInterfaceEdgeFan(binding->second.edge, ws, resolvers, interfaceFaces,
                                       exteriorFaces, result.metrics, request.maxAuditWork)) {
                fail(out, T5Status::Unsupported, T5Reason::UnsupportedRegionJunction,
                     T5Phase::ArtifactRevalidation);
                return out;
            }
            T5InterfaceEdgeRow row;
            row.canonicalGlobalNodes = artifactEdge->key.canonicalGlobalNodes;
            row.authoritativeEdge = binding->second.edge;
            row.classification = artifactEdge->classification;
            row.required = true;
            row.incidentArtifactFacetOffset = result.interfaceEdgeIncidentArtifactFacets.size();
            row.incidentArtifactFacetCount = artifactEdge->incidentFacets.size();
            result.interfaceEdgeIncidentArtifactFacets.insert(
                result.interfaceEdgeIncidentArtifactFacets.end(),
                artifactEdge->incidentFacets.begin(), artifactEdge->incidentFacets.end());
            result.interfaceEdges.push_back(row);
            ++result.metrics.interfaceEdgeInsertions;
        }
        result.audit.interfaceNodes = result.interfaceNodes.size();
        result.audit.interfaceEdges = result.interfaceEdges.size();
        result.audit.regionCount = regionCount;
        result.audit.materialCount = materialCount;

        activePhase = T5Phase::RegionAuthority;
        std::map<RegionPersistentId, Expansion> volumeByRegion;
        std::map<RegionPersistentId, std::vector<topo::CellHandle>> tetsByRegion;
        for (const auto &region : normalizedRegions)
            volumeByRegion[region.regionPersistentId].terms.push_back(0.0);
        for (const auto &plan : tetPlans) {
            if (!consumeAudit(result.metrics, request.maxAuditWork)) {
                fail(out, T5Status::Incomplete, T5Reason::WorkBudgetReached,
                     T5Phase::RegionAuthority);
                return out;
            }
            std::array<Vec3, 4> coordinates{};
            for (std::size_t i = 0; i < 4; ++i)
                coordinates[i] = sourceMesh->node(plan.nodes[i]);
            if (predicates::orient3d(coordinates[0], coordinates[1], coordinates[2],
                                     coordinates[3]) != predicates::Sign::Positive) {
                ++result.audit.invertedTets;
                fail(out, T5Status::Failed, T5Reason::RegionCellMismatch, T5Phase::RegionAuthority);
                return out;
            }
            Expansion tetVolume;
            std::uint64_t exactWork = 0;
            if (!exactSixVolume(coordinates, tetVolume, exactWork) ||
                tetVolume.terms.size() > kT5SixVolumeExpansionTermsPerTet) {
                fail(out, T5Status::Unsupported, T5Reason::VolumeRangeUnrepresentable,
                     T5Phase::RegionAuthority);
                return out;
            }
            if (!consumeAudit(result.metrics, request.maxAuditWork, exactWork)) {
                fail(out, T5Status::Incomplete, T5Reason::WorkBudgetReached,
                     T5Phase::RegionAuthority);
                return out;
            }
            result.metrics.volumeExpansionOperations += exactWork;
            auto &aggregate = volumeByRegion.at(plan.region);
            aggregate = sumExpansion(aggregate, tetVolume);
            if (aggregate.terms.size() > kT5SixVolumeExpansionTermsPerTet * tetPlans.size()) {
                fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::RegionAuthority);
                return out;
            }
            tetsByRegion[plan.region].push_back(plan.cell);
        }
        std::map<RegionPersistentId, std::vector<bnd::VolumeFacetKey>> facetsByRegion;
        for (const auto &row : result.interfaceFacets) {
            facetsByRegion[row.firstRegion].push_back(row.facet);
            facetsByRegion[row.secondRegion].push_back(row.facet);
        }
        for (auto &[region, facets] : facetsByRegion) {
            (void)region;
            std::sort(facets.begin(), facets.end());
        }
        std::map<MaterialPersistentId, Expansion> volumeByMaterial;
        std::map<MaterialPersistentId, std::vector<RegionPersistentId>> regionsByMaterial;
        for (const auto &binding : normalizedRegions) {
            T5RegionRow row;
            row.authority = binding;
            const auto &tets = tetsByRegion[binding.regionPersistentId];
            const auto &facets = facetsByRegion[binding.regionPersistentId];
            row.tetCount = tets.size();
            row.interfaceFacetCount = facets.size();
            row.tetReferenceOffset = result.regionTets.size();
            row.tetReferenceCount = tets.size();
            result.regionTets.insert(result.regionTets.end(), tets.begin(), tets.end());
            row.interfaceReferenceOffset = result.regionInterfaceFacets.size();
            row.interfaceReferenceCount = facets.size();
            result.regionInterfaceFacets.insert(result.regionInterfaceFacets.end(), facets.begin(),
                                                facets.end());
            // Discard the temporary pairwise append performed during facet
            // routing; region spans above are the sole authoritative order.
            const auto &expansion = volumeByRegion.at(binding.regionPersistentId);
            row.sixVolumeExpansionOffset = result.exactVolumeExpansionTerms.size();
            row.sixVolumeExpansionCount = expansion.terms.size();
            result.exactVolumeExpansionTerms.insert(result.exactVolumeExpansionTerms.end(),
                                                    expansion.terms.begin(), expansion.terms.end());
            if (!volumeInterval(expansion, row.volumeLower, row.volumeUpper)) {
                fail(out, T5Status::Unsupported, T5Reason::VolumeRangeUnrepresentable,
                     T5Phase::RegionAuthority);
                return out;
            }
            row.exactVolumeDigest = expansionDigest(expansion);
            result.regions.push_back(row);
            auto &material = volumeByMaterial[binding.materialPersistentId];
            if (material.terms.empty())
                material.terms.push_back(0.0);
            material = sumExpansion(material, expansion);
            regionsByMaterial[binding.materialPersistentId].push_back(binding.regionPersistentId);
        }
        for (const auto materialId : inputMaterials) {
            auto &regionIds = regionsByMaterial[materialId];
            std::sort(regionIds.begin(), regionIds.end());
            const auto &expansion = volumeByMaterial.at(materialId);
            T5MaterialRow row;
            row.materialPersistentId = materialId;
            row.regionReferenceOffset = result.materialRegions.size();
            row.regionReferenceCount = regionIds.size();
            result.materialRegions.insert(result.materialRegions.end(), regionIds.begin(),
                                          regionIds.end());
            for (const auto regionId : regionIds)
                row.tetCount += tetsByRegion.at(regionId).size();
            row.sixVolumeExpansionOffset = result.exactVolumeExpansionTerms.size();
            row.sixVolumeExpansionCount = expansion.terms.size();
            result.exactVolumeExpansionTerms.insert(result.exactVolumeExpansionTerms.end(),
                                                    expansion.terms.begin(), expansion.terms.end());
            if (!volumeInterval(expansion, row.volumeLower, row.volumeUpper)) {
                fail(out, T5Status::Unsupported, T5Reason::VolumeRangeUnrepresentable,
                     T5Phase::RegionAuthority);
                return out;
            }
            row.exactVolumeDigest = expansionDigest(expansion);
            result.materials.push_back(row);
        }
        if (result.exactVolumeExpansionTerms.size() > volumeTerms) {
            fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::RegionAuthority);
            return out;
        }

        // Predict every final publication index and retain separate surface
        // and authoritative-Tet provenance before consuming the v3 mesh.
        for (const auto node : ws.nodes().activeBySequence()) {
            const auto global = resolvers.globalByNode.at(node);
            result.sourceNodeMappings.push_back(
                {static_cast<std::uint64_t>(node.slot), node, global, global});
        }
        std::sort(result.sourceNodeMappings.begin(), result.sourceNodeMappings.end(),
                  [](const auto &a, const auto &b) {
                      return std::tie(a.sourceGlobalNode, a.authoritativeSlot) <
                             std::tie(b.sourceGlobalNode, b.authoritativeSlot);
                  });
        for (std::size_t i = 0; i < surfacePlans.size(); ++i) {
            const auto &plan = surfacePlans[i];
            result.sourceSurfaceCellMappings.push_back(
                {plan.key.facePersistentId,
                 static_cast<std::uint64_t>(plan.key.normalizedLocalCellOrdinal), plan.source,
                 static_cast<Index>(i), plan.classification});
        }
        for (std::size_t i = 0; i < tetPlans.size(); ++i) {
            const auto &plan = tetPlans[i];
            result.sourceTetMappings.push_back(
                {static_cast<std::uint64_t>(plan.cell.slot), plan.cell, plan.source,
                 static_cast<Index>(surfacePlans.size() + i), plan.classification});
        }
        if (result.sourceNodeMappings.size() != sourceNodeCount ||
            result.sourceSurfaceCellMappings.size() != sourceSurfaceCount ||
            result.sourceTetMappings.size() != sourceTetCount) {
            fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::ResultAssembly);
            return out;
        }
        std::vector<NodePlan> nodePlans;
        nodePlans.reserve(static_cast<std::size_t>(sourceMesh->nodeCount()));
        for (Index node = 0; node < sourceMesh->nodeCount(); ++node)
            nodePlans.push_back({sourceMesh->node(node), sourceMesh->nodeClass(node), node});
        if (nodePlans.size() != sourceNodeCount) {
            fail(out, T5Status::Failed, T5Reason::MalformedSourceOutcome,
                 T5Phase::SourceValidation);
            return out;
        }
        result.constraintState.swap(rebound);
        result.artifact = std::move(artifact);

        bool sourceMeshAccountedOk = true, sourceOutcomeAccountedOk = true;
        const auto sourceMeshBytes =
            meshAccounted(sourceMesh->capacitySnapshot(), sourceMeshAccountedOk);
        const auto sourceOutcomeAccounted = boundaryOutcomeBytes(source, sourceOutcomeAccountedOk);
        std::size_t preHandoffBytes = 0;
        bool preHandoffOk =
            sourceMeshAccountedOk && sourceOutcomeAccountedOk &&
            checkedAdd(sourceMeshBytes, sourceOutcomeAccounted, preHandoffBytes) &&
            addBytes(preHandoffBytes, out.report->processedStorage.capacity(),
                     sizeof(T5ProcessedRow)) &&
            addBytes(preHandoffBytes, out.report->unattemptedStorage.capacity(), sizeof(T5WorkKey));
        bool partialResultOk = true;
        const auto partialResultBytes = resultBytes(result, partialResultOk);
        preHandoffOk = preHandoffOk && partialResultOk &&
                       checkedAdd(preHandoffBytes, partialResultBytes, preHandoffBytes);
        if (!preHandoffOk || preHandoffBytes > request.maxPostSourceAccountedBytes) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::ResultAssembly);
            return out;
        }
        result.metrics.postSourceAccountedHighWaterBytes = preHandoffBytes;
        result.metrics.postSourceAccountedHighWaterPhase = T5Phase::ResultAssembly;

        // The final provider checkpoint follows every provider-dependent
        // artifact, owner-side and fan audit and immediately precedes the
        // consuming handoff.
        if (ownershipProvider.stamp() != entryProviderStamp ||
            entryProviderStamp != request.sourceRequest.expectedOwnershipProviderStamp) {
            fail(out, T5Status::Failed, T5Reason::StaleOwnershipProvider, T5Phase::InterfaceAudit);
            return out;
        }

        auto receipt = std::make_unique<T5V3SourceReceipt>();
        receipt->surfaceSetStamp = request.sourceRequest.expectedSurfaceSetStamp;
        receipt->providerStamp = entryProviderStamp;
        receipt->surfaceOutcome = std::move(source.surfaceOutcome);
        receipt->sourceReport = std::move(source.report);
        receipt->artifactMetrics = sourceResult.artifactPreflightMetrics;
        receipt->constructionMetrics = sourceResult.constructionMetrics;
        receipt->recoveryMetrics = sourceResult.recoveryMetrics;
        receipt->extractionMetrics = sourceResult.extractionMetrics;
        receipt->memoryMetrics = sourceResult.memoryMetrics;
        receipt->sourceAudit = sourceResult.boundaryAudit;
        receipt->sourceStrictHash = mesh::meshHash(*sourceMesh);
        receipt->sourceSemanticHash = mesh::semanticHashEx(*sourceMesh);
        receipt->sourceNodeMappingDigest = mappingDigest(result.sourceNodeMappings);
        receipt->sourceSurfaceCellMappingDigest = mappingDigest(result.sourceSurfaceCellMappings);
        receipt->sourceTetMappingDigest = mappingDigest(result.sourceTetMappings);
        receipt->sourceNodeCount = sourceNodeCount;
        receipt->sourceSurfaceCellCount = sourceSurfaceCount;
        receipt->sourceTetCount = sourceTetCount;

        result.authoritative = std::move(sourceResult.finalAuthoritative);
        source.retainedResult.reset();
        sourceMesh.reset();
        out.sourceReceipt = std::move(receipt);
        out.sourceOutcome.reset();
        out.retainedResult = std::move(working);
        auto &retained = *out.retainedResult;

        activePhase = T5Phase::MeshBuild;
        const auto publicationNodes = static_cast<std::size_t>(nodePlans.size());
        if (publicationNodes > request.maxPublicationNodes ||
            publicationCells > request.maxPublicationCells ||
            publicationConnectivity > request.maxPublicationConnectivityEntries ||
            publicationNodes >= static_cast<std::size_t>(kInvalidIndex) ||
            publicationCells >= static_cast<std::size_t>(kInvalidIndex) ||
            publicationConnectivity >= static_cast<std::size_t>(kInvalidIndex)) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::MeshBuild);
            return out;
        }
        publicationAllocations.activate();
        retained.privateMesh = std::make_unique<mesh::Mesh>();
        retained.privateMesh->reserveAppend(nodePlans.size(), publicationCells,
                                            publicationConnectivity);
        const auto preparedSnapshot = retained.privateMesh->capacitySnapshot();
        bool preparedMeshOk = true;
        (void)meshAccounted(preparedSnapshot, preparedMeshOk);
        if (!preparedMeshOk || preparedSnapshot.dedupEntryCount != 0 ||
            preparedSnapshot.dedupIndexCapacity != 0) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::MeshBuild);
            return out;
        }
        bool retainedBytesOk = true;
        auto retainedBytes = retainedOutcomeBytes(out, retainedBytesOk);
        if (!retainedBytesOk || retainedBytes > request.maxPostSourceAccountedBytes) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::MeshBuild);
            return out;
        }
        if (retainedBytes > retained.metrics.postSourceAccountedHighWaterBytes) {
            retained.metrics.postSourceAccountedHighWaterBytes = retainedBytes;
            retained.metrics.postSourceAccountedHighWaterPhase = T5Phase::MeshBuild;
        }
        for (const auto &node : nodePlans) {
            const Index published =
                retained.privateMesh->addNode(node.coordinate, node.classification);
            if (published != node.source) {
                fail(out, T5Status::Failed, T5Reason::MeshBuildFailed, T5Phase::MeshBuild);
                return out;
            }
        }
        for (std::size_t i = 0; i < surfacePlans.size(); ++i) {
            const auto &plan = surfacePlans[i];
            const Index published = retained.privateMesh->addCell(
                mesh::CellType::Tri3, std::vector<Index>(plan.nodes.begin(), plan.nodes.end()),
                plan.classification);
            if (published != static_cast<Index>(i)) {
                fail(out, T5Status::Failed, T5Reason::MeshBuildFailed, T5Phase::MeshBuild);
                return out;
            }
        }
        for (std::size_t i = 0; i < tetPlans.size(); ++i) {
            const auto &plan = tetPlans[i];
            const Index published = retained.privateMesh->addCell(
                mesh::CellType::Tet4, std::vector<Index>(plan.nodes.begin(), plan.nodes.end()),
                plan.classification);
            if (published != static_cast<Index>(surfacePlans.size() + i)) {
                fail(out, T5Status::Failed, T5Reason::MeshBuildFailed, T5Phase::MeshBuild);
                return out;
            }
        }
        const auto completedSnapshot = retained.privateMesh->capacitySnapshot();
        if (!completedSnapshot.complete || completedSnapshot.dedupEntryCount != 0 ||
            completedSnapshot.dedupIndexCapacity != 0) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached, T5Phase::MeshBuild);
            return out;
        }

        for (auto &row : retained.interfaceFacets) {
            const auto plan = surfacePlanByKey.find(row.facet);
            if (plan == surfacePlanByKey.end() ||
                plan->second >= retained.sourceSurfaceCellMappings.size()) {
                fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::MeshBuild);
                return out;
            }
            const auto &mapped = retained.sourceSurfaceCellMappings[plan->second];
            if (mapped.face != row.facet.facePersistentId ||
                mapped.normalizedLocalCellOrdinal != row.facet.normalizedLocalCellOrdinal) {
                fail(out, T5Status::Failed, T5Reason::InternalError, T5Phase::MeshBuild);
                return out;
            }
            ++retained.metrics.interfaceResolverLookups;
            row.publishedTri3Cell = mapped.publishedGlobalCell;
        }

        activePhase = T5Phase::BoundaryAudit;
        mesh::AuditOptions auditOptions;
        auditOptions.boundaryMode = mesh::BoundaryComparisonMode::Required;
        auditOptions.boundaryOrphanBudget = 0;
        auditOptions.boundaryUnexpectedBudget = 0;
        auditOptions.duplicateFacetsAreErrors = true;
        auditOptions.unusedNodesAreErrors = true;
        const auto flatAudit = mesh::auditMesh(*retained.privateMesh, auditOptions);
        retained.audit.unusedNodes = flatAudit.unusedNodes;
        retained.audit.invertedTets = flatAudit.nInverted;
        // MeshAudit's exterior-only mode deliberately counts each exact
        // two-sided interface Tri3 as one orphan. The authoritative audit
        // above proves the interface bijection and requires literal zero in
        // its own missing/unexpected/duplicate fields.
        if (flatAudit.boundaryOrphans != static_cast<Index>(retained.audit.prescribedInterface) ||
            flatAudit.boundaryUnexpected != 0 || flatAudit.boundaryDuplicates != 0 ||
            flatAudit.nonManifoldFaces != 0 || flatAudit.nInverted != 0 ||
            flatAudit.nNonFinite != 0 || flatAudit.unusedNodes != 0 ||
            retained.audit.prescribedExterior != retained.audit.matchedExterior ||
            retained.audit.prescribedInterface != retained.audit.matchedInterface ||
            retained.audit.orphanExterior != 0 || retained.audit.unexpectedExterior != 0 ||
            retained.audit.missingInterface != 0 || retained.audit.unexpectedInterface != 0 ||
            retained.audit.duplicateInterface != 0 || retained.audit.duplicateInterfaceNodes != 0 ||
            retained.audit.duplicateInterfaceEdges != 0 || retained.audit.wrongOwnerSense != 0 ||
            retained.audit.wrongRegionClassification != 0 ||
            retained.audit.undeclaredRegionFaces != 0) {
            fail(out, T5Status::Failed, T5Reason::BoundaryBijectionFailed, T5Phase::BoundaryAudit);
            return out;
        }
        retained.audit.passed = true;
        retained.strictHash = mesh::meshHash(*retained.privateMesh);
        retained.semanticHash = mesh::semanticHashEx(*retained.privateMesh);
        if (!retained.semanticHash.canonical) {
            fail(out, T5Status::Failed, T5Reason::BoundaryBijectionFailed, T5Phase::BoundaryAudit);
            return out;
        }
        retained.metrics.publicationRows = publicationCells;
        retained.metrics.perInterfaceGlobalTetScans = 0;
        retained.metrics.perRegionGlobalFaceScans = 0;
        retained.metrics.contextRebuilds = 0;
        retained.metrics.publicationTailAllocations = 0;
        {
            std::uint64_t total = 0;
            const auto add = [&](auto... values) noexcept {
                return (checkedAdd(total, static_cast<std::uint64_t>(values), total) && ...);
            };
            if (!add(retained.metrics.sourceWorkUnits, retained.metrics.regionRows,
                     retained.metrics.cellClassificationVisits, retained.metrics.faceOwnerVisits,
                     retained.metrics.interfaceNodeInsertions,
                     retained.metrics.interfaceEdgeInsertions,
                     retained.metrics.interfaceFacetInsertions,
                     retained.metrics.interfaceResolverLookups,
                     retained.metrics.regionAdjacencyEdges,
                     retained.metrics.volumeExpansionOperations, retained.metrics.publicationRows,
                     retained.metrics.auditWorkUnits, retained.metrics.perInterfaceGlobalTetScans,
                     retained.metrics.perRegionGlobalFaceScans, retained.metrics.contextRebuilds,
                     retained.metrics.publicationTailAllocations)) {
                fail(out, T5Status::Incomplete, T5Reason::WorkBudgetReached,
                     T5Phase::ResultAssembly);
                return out;
            }
            retained.metrics.totalWorkUnits = total;
        }
        retained.metrics.runComplete = true;

        out.report->processedStorage.resize(expectedKeys.size());
        for (std::size_t i = 0; i < expectedKeys.size(); ++i) {
            out.report->processedStorage[i].key = expectedKeys[i];
            out.report->processedStorage[i].status = T5WorkStatus::Published;
            out.report->processedStorage[i].canonicalOrdinal = i;
        }
        out.report->processedUsed = expectedKeys.size();
        out.report->unattemptedStorage.clear();
        out.report->unattemptedUsed = 0;
        out.report->reportComplete = true;
        out.reportComplete = true;

        activePhase = T5Phase::AuthoritySerialization;
        const auto authorityReason =
            t5_detail::buildAuthorityCore(retained, request.sourceRequest.expectedSurfaceSetStamp,
                                          entryProviderStamp, request.maxAuthorityCoreBytes);
        if (authorityReason != T5Reason::None) {
            fail(out,
                 authorityReason == T5Reason::AuthorityBudgetReached ? T5Status::Incomplete
                                                                     : T5Status::Failed,
                 authorityReason, T5Phase::AuthoritySerialization);
            return out;
        }
        retainedBytesOk = true;
        retainedBytes = retainedOutcomeBytes(out, retainedBytesOk);
        if (!retainedBytesOk || retainedBytes > request.maxPostSourceAccountedBytes) {
            fail(out, T5Status::Incomplete, T5Reason::ResultBudgetReached,
                 T5Phase::AuthoritySerialization);
            return out;
        }
        if (retainedBytes > retained.metrics.postSourceAccountedHighWaterBytes) {
            retained.metrics.postSourceAccountedHighWaterBytes = retainedBytes;
            retained.metrics.postSourceAccountedHighWaterPhase = T5Phase::AuthoritySerialization;
        }
        // The first serialization establishes the authority vector's actual
        // capacity. Re-encode once after recording that capacity high-water
        // so retained metrics and serialized metrics are identical. Payload
        // size is schema-fixed; this second pass cannot create a memory-policy
        // fixed-point loop.
        const auto finalAuthorityReason =
            t5_detail::buildAuthorityCore(retained, request.sourceRequest.expectedSurfaceSetStamp,
                                          entryProviderStamp, request.maxAuthorityCoreBytes);
        if (finalAuthorityReason != T5Reason::None) {
            fail(out,
                 finalAuthorityReason == T5Reason::AuthorityBudgetReached ? T5Status::Incomplete
                                                                          : T5Status::Failed,
                 finalAuthorityReason, T5Phase::AuthoritySerialization);
            return out;
        }
        retained.publicationReady = true;

        activePhase = T5Phase::Publication;
        publicationAllocations.restore();
        t5_detail::publicationTailActive() = true;
        destination.swap(*retained.privateMesh);
        retained.privateMesh.reset();
        t5_detail::publicationTailActive() = false;
        out.status = T5Status::Published;
        out.failure = {};
        out.destinationUnchanged = false;
        out.reportComplete = true;
        return out;
    } catch (const std::bad_alloc &) {
        t5_detail::publicationTailActive() = false;
        fail(out, T5Status::Failed, T5Reason::AllocationFailed, activePhase);
        return out;
    } catch (...) {
        t5_detail::publicationTailActive() = false;
        fail(out, T5Status::Failed, T5Reason::InternalError, activePhase);
        return out;
    }
}

} // namespace paacs::algo
