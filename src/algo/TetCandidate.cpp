// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/TetCandidate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <utility>

#include "paacs/predicates/Predicates.hpp"
#include "TetCandidateDetail.hpp"
#if defined(_MSC_VER)
#pragma warning(push)
// The frozen TransactionAccess test assertion consumes `pred` only when
// PAACS_TOPO_TESTING is enabled. Production CK-3 includes the workspace here
// without enabling test seams; suppress that pre-existing header-only C4100
// locally rather than editing frozen topology.
#pragma warning(disable : 4100)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
// The same frozen PAACS_TOPO_TESTING-only assertion leaves `pred` unused in
// ordinary builds. Keep the warning containment local to this additive TU.
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "paacs/topo/Compaction.hpp"
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace paacs::algo {

TetAuthoritativeResult::TetAuthoritativeResult() = default;
TetAuthoritativeResult::~TetAuthoritativeResult() = default;
TetAuthoritativeResult::TetAuthoritativeResult(
    TetAuthoritativeResult&&) noexcept = default;
TetAuthoritativeResult& TetAuthoritativeResult::operator=(
    TetAuthoritativeResult&&) noexcept = default;

namespace {

using predicates::Sign;
using Budget = detail::TetConstructionBudget;
using detail::constructionStep;

constexpr int kTetFace[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

using FaceKey = std::array<Index, 3>;
using CellKey = std::array<Index, 4>;

enum class ProtectedBoundaryImportPolicy : std::uint8_t {
    RequireComplete,
    ExactPresentOnly
};

bool finite(const Vec3& point) noexcept {
    return std::isfinite(point.x) && std::isfinite(point.y) &&
           std::isfinite(point.z);
}

bool pointLess(const Vec3& a, const Vec3& b) noexcept {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
}

bool pointEqual(const Vec3& a, const Vec3& b) noexcept {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

FaceKey faceKey(std::array<Index, 3> value) noexcept {
    std::sort(value.begin(), value.end());
    return value;
}

CellKey cellKey(const TetCandidateCell& cell) noexcept {
    CellKey key = cell.vertices;
    std::sort(key.begin(), key.end());
    return key;
}

bool reversedCycle(const std::array<Index, 3>& a,
                   const std::array<Index, 3>& b) noexcept {
    for (int start = 0; start < 3; ++start) {
        if (b[static_cast<std::size_t>(start)] != a[0]) continue;
        bool reverse = true;
        for (int i = 0; i < 3; ++i)
            reverse = reverse &&
                b[static_cast<std::size_t>((start - i + 3) % 3)] ==
                    a[static_cast<std::size_t>(i)];
        if (reverse) return true;
    }
    return false;
}

std::array<Index, 3> directedFace(const TetCandidateCell& cell,
                                  int localFace) noexcept {
    return {cell.vertices[static_cast<std::size_t>(kTetFace[localFace][0])],
            cell.vertices[static_cast<std::size_t>(kTetFace[localFace][1])],
            cell.vertices[static_cast<std::size_t>(kTetFace[localFace][2])]};
}

std::uint64_t vectorBytes(std::size_t capacity, std::size_t element) noexcept {
    constexpr auto max = (std::numeric_limits<std::uint64_t>::max)();
    if (element != 0 && capacity > max / element) return max;
    return static_cast<std::uint64_t>(capacity * element);
}

void saturatingAdd(std::uint64_t& value, std::uint64_t addend) noexcept {
    const auto max = (std::numeric_limits<std::uint64_t>::max)();
    value = addend > max - value ? max : value + addend;
}

bool checkedMetricAdd(std::uint64_t& value,
                      std::uint64_t addend = 1) noexcept {
    const auto max = (std::numeric_limits<std::uint64_t>::max)();
    if (addend > max - value) return false;
    value += addend;
    return true;
}

bool checkedMetricAddSize(std::uint64_t& value,
                          std::size_t addend) noexcept {
    if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
        if (addend >
            static_cast<std::size_t>(
                (std::numeric_limits<std::uint64_t>::max)()))
            return false;
    }
    return checkedMetricAdd(value,static_cast<std::uint64_t>(addend));
}

bool checkedMetricProductSize(std::uint64_t& value,std::size_t a,
                              std::size_t b) noexcept {
    if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
        const auto max=static_cast<std::size_t>(
            (std::numeric_limits<std::uint64_t>::max)());
        if (a>max||b>max) return false;
    }
    const auto aa=static_cast<std::uint64_t>(a);
    const auto bb=static_cast<std::uint64_t>(b);
    if (aa!=0&&bb>(std::numeric_limits<std::uint64_t>::max)()/aa)
        return false;
    return checkedMetricAdd(value,aa*bb);
}

struct CountingIndexLess {
    std::uint64_t* comparisons=nullptr;
    bool* complete=nullptr;

    bool operator()(Index a,Index b) const noexcept {
        if (!checkedMetricAdd(*comparisons)) *complete=false;
        return a<b;
    }
};

struct FaceUse {
    FaceKey key;
    Index cell = kInvalidIndex;
    std::uint8_t localFace = 0;
};

void appendFaceUses(std::vector<FaceUse>& uses, const TetCandidateCell& cell, Index owner) {
    for (std::uint8_t side=0;side<4;++side)
        uses.push_back({faceKey(directedFace(cell,side)),owner,side});
}

void sortFaceUses(std::vector<FaceUse>& uses, Budget* budget = nullptr) {
    std::sort(uses.begin(),uses.end(),[budget](const FaceUse& a,const FaceUse& b) {
        constructionStep(budget);
        return std::tie(a.key,a.cell,a.localFace)<std::tie(b.key,b.cell,b.localFace);
    });
}

std::size_t nextFace(const std::vector<FaceUse>& uses,std::size_t first, Budget* budget=nullptr) {
    auto end=first+1;
    while (end<uses.size() && uses[end].key==uses[first].key) {
        constructionStep(budget); ++end;
    }
    return end;
}

TetCandidateValidation invalid(TetCandidateReason reason,
                               Index entity = kInvalidIndex) noexcept {
    return TetCandidateValidation{false, reason, entity};
}

}  // namespace

std::uint64_t tetCandidateAccountedBytes(
    const TetCandidate& candidate) noexcept {
    std::uint64_t bytes = sizeof(TetCandidate);
    saturatingAdd(bytes,
                  vectorBytes(candidate.points.capacity(), sizeof(Vec3)));
    saturatingAdd(bytes, vectorBytes(candidate.pointToInput.capacity(),
                                     sizeof(Index)));
    saturatingAdd(bytes, vectorBytes(candidate.inputToPoint.capacity(),
                                     sizeof(Index)));
    saturatingAdd(bytes, vectorBytes(candidate.cells.capacity(),
                                     sizeof(TetCandidateCell)));
    saturatingAdd(bytes, vectorBytes(candidate.boundaryFaces.capacity(),
                                     sizeof(TetCandidateBoundaryFace)));
    return bytes;
}

static TetCandidateValidation validateTetCandidateImpl(
    const TetCandidate& candidate, std::uint64_t& hullTests,
    const std::vector<FaceUse>* constructedUses=nullptr, Budget* budget=nullptr) noexcept {
    try {
        constructionStep(budget,0);
        const std::size_t pointCount = candidate.points.size();
        const std::size_t cellCount = candidate.cells.size();
        if (pointCount < 4 || cellCount == 0)
            return invalid(TetCandidateReason::NonVolumetricInput);

        Vec3 lo=candidate.points.front(), hi=lo;
        for (std::size_t i = 0; i < pointCount; ++i) {
            constructionStep(budget);
            if (!finite(candidate.points[i]))
                return invalid(TetCandidateReason::NonfiniteInput,
                               static_cast<Index>(i));
            if (i > 0 && !pointLess(candidate.points[i - 1],
                                     candidate.points[i]))
                return invalid(TetCandidateReason::NonCanonicalOrdering,
                               static_cast<Index>(i));
            const auto& p=candidate.points[i];
            lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};
            hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};
        }
        if (candidate.pointToInput.size() != pointCount ||
            candidate.inputToPoint.empty())
            return invalid(TetCandidateReason::InvalidPointMapping);
        constructionStep(budget,pointCount);
        std::vector<bool> mapped(pointCount, false);
        for (std::size_t i = 0; i < candidate.inputToPoint.size(); ++i) {
            constructionStep(budget);
            const Index point = candidate.inputToPoint[i];
            if (point < 0 || static_cast<std::size_t>(point) >= pointCount)
                return invalid(TetCandidateReason::InvalidPointMapping,
                               static_cast<Index>(i));
            mapped[static_cast<std::size_t>(point)] = true;
        }
        for (std::size_t point = 0; point < pointCount; ++point) {
            constructionStep(budget);
            const Index input = candidate.pointToInput[point];
            if (input < 0 ||
                static_cast<std::size_t>(input) >=
                    candidate.inputToPoint.size() ||
                candidate.inputToPoint[static_cast<std::size_t>(input)] !=
                    static_cast<Index>(point) || !mapped[point])
                return invalid(TetCandidateReason::InvalidPointMapping,
                               static_cast<Index>(point));
        }

        constructionStep(budget,pointCount);
        constructionStep(budget,cellCount);
        std::vector<bool> used(pointCount, false);
        constructionStep(budget,0);
        std::vector<bool> sourceOrdinals(cellCount, false);
        std::vector<FaceUse> validationUses;
        if (!constructedUses && cellCount<=validationUses.max_size()/4)
            validationUses.reserve(cellCount*4);
        CellKey previous{};
        bool havePrevious = false;
        for (std::size_t ci = 0; ci < cellCount; ++ci) {
            constructionStep(budget);
            const auto& cell = candidate.cells[ci];
            std::array<Index, 4> distinct = cell.vertices;
            for (const Index vertex : distinct) {
                constructionStep(budget);
                if (vertex < 0 ||
                    static_cast<std::size_t>(vertex) >= pointCount)
                    return invalid(TetCandidateReason::InvalidCellIndex,
                                   static_cast<Index>(ci));
                used[static_cast<std::size_t>(vertex)] = true;
            }
            std::sort(distinct.begin(), distinct.end());
            if (std::adjacent_find(distinct.begin(), distinct.end()) !=
                distinct.end())
                return invalid(TetCandidateReason::RepeatedCellVertex,
                               static_cast<Index>(ci));
            if (predicates::orient3d(
                    candidate.points[static_cast<std::size_t>(cell.vertices[0])],
                    candidate.points[static_cast<std::size_t>(cell.vertices[1])],
                    candidate.points[static_cast<std::size_t>(cell.vertices[2])],
                    candidate.points[static_cast<std::size_t>(cell.vertices[3])])
                != Sign::Positive)
                return invalid(TetCandidateReason::NonPositiveOrientation,
                               static_cast<Index>(ci));
            const CellKey key = cellKey(cell);
            if (havePrevious && !(previous < key))
                return invalid(previous == key ? TetCandidateReason::DuplicateCell
                                               : TetCandidateReason::NonCanonicalOrdering,
                               static_cast<Index>(ci));
            previous = key;
            havePrevious = true;
            if (cell.sourceOrdinal < 0 ||
                static_cast<std::size_t>(cell.sourceOrdinal) >= cellCount ||
                sourceOrdinals[static_cast<std::size_t>(cell.sourceOrdinal)])
                return invalid(TetCandidateReason::NonCanonicalOrdering,
                               static_cast<Index>(ci));
            sourceOrdinals[static_cast<std::size_t>(cell.sourceOrdinal)] = true;
            if (!constructedUses) appendFaceUses(validationUses,cell,static_cast<Index>(ci));
        }
        for (std::size_t i = 0; i < used.size(); ++i) {
            constructionStep(budget);
            if (!used[i])
                return invalid(TetCandidateReason::MissingRequiredPoint,
                               static_cast<Index>(i));
        }

        if (!constructedUses) sortFaceUses(validationUses,budget);
        const auto& uses=constructedUses?*constructedUses:validationUses;
        std::size_t boundaryCount=0;
        std::uint64_t interiorFaces = 0;
        for (std::size_t first=0;first<uses.size();) {
            constructionStep(budget);
            const auto end=nextFace(uses,first,budget);
            const auto& a=uses[first];
            const auto& cellA=candidate.cells[static_cast<std::size_t>(a.cell)];
            if (end-first>2)
                return invalid(TetCandidateReason::NonManifoldFace,a.cell);
            if (end-first==1) {
                if (cellA.neighbours[a.localFace]!=kInvalidIndex)
                    return invalid(TetCandidateReason::AdjacencyMismatch,
                                   a.cell);
                if (boundaryCount>=candidate.boundaryFaces.size())
                    return invalid(TetCandidateReason::BoundaryExtractionMismatch);
                const auto& actual=candidate.boundaryFaces[boundaryCount];
                if (actual.vertices!=directedFace(cellA,a.localFace) ||
                    actual.owner!=a.cell || actual.ownerLocalFace!=a.localFace)
                    return invalid(TetCandidateReason::BoundaryExtractionMismatch,
                                   static_cast<Index>(boundaryCount));
                ++boundaryCount;
            } else {
                ++interiorFaces;
                const auto& b=uses[first+1];
                const auto& cellB=candidate.cells[static_cast<std::size_t>(b.cell)];
                if (cellA.neighbours[a.localFace]!=b.cell ||
                    cellB.neighbours[b.localFace]!=a.cell ||
                    !reversedCycle(directedFace(cellA,a.localFace),directedFace(cellB,b.localFace)))
                    return invalid(TetCandidateReason::AdjacencyMismatch,
                                   a.cell);
            }
            first=end;
        }
        if (boundaryCount!=candidate.boundaryFaces.size())
            return invalid(TetCandidateReason::BoundaryExtractionMismatch);

        // The Tet dual must be connected.
        constructionStep(budget,cellCount);
        std::vector<bool> reached(cellCount, false);
        constructionStep(budget,0);
        std::queue<Index> queue;
        reached[0] = true;
        constructionStep(budget,0);
        queue.push(0);
        std::size_t reachedCount = 0;
        while (!queue.empty()) {
            constructionStep(budget);
            const Index cell = queue.front();
            queue.pop();
            ++reachedCount;
            for (const Index neighbour :
                 candidate.cells[static_cast<std::size_t>(cell)].neighbours) {
                constructionStep(budget);
                if (neighbour == kInvalidIndex) continue;
                if (neighbour < 0 ||
                    static_cast<std::size_t>(neighbour) >= cellCount)
                    return invalid(TetCandidateReason::AdjacencyMismatch,
                                   cell);
                const auto neighbourIndex =
                    static_cast<std::size_t>(neighbour);
                if (!reached[neighbourIndex]) {
                    reached[neighbourIndex] = true;
                    constructionStep(budget,0);
                    queue.push(neighbour);
                }
            }
        }
        if (reachedCount != cellCount)
            return invalid(TetCandidateReason::DisconnectedComplex);

        // Boundary is one closed orientable sphere and every directed face is
        // a supporting face of the exact convex hull.
        using Edge = std::array<Index, 2>;
        std::map<Edge, std::vector<std::pair<Index, Index>>> boundaryEdges;
        std::set<Index> boundaryVertices;
        const auto outside=[&](Vec3 a, Vec3 b, Vec3 c, Vec3 point) {
            constructionStep(budget);
            ++hullTests;
            return predicates::orient3d(a,b,c,point)==Sign::Positive;
        };
        for (std::size_t fi = 0; fi < candidate.boundaryFaces.size(); ++fi) {
            constructionStep(budget);
            const auto& face = candidate.boundaryFaces[fi].vertices;
            for (const Index vertex : face) {
                constructionStep(budget); boundaryVertices.insert(vertex);
            }
            for (int edge = 0; edge < 3; ++edge) {
                constructionStep(budget);
                const Index a = face[static_cast<std::size_t>(edge)];
                const Index b = face[static_cast<std::size_t>((edge + 1) % 3)];
                Edge key{std::min(a, b), std::max(a, b)};
                constructionStep(budget,0);
                auto& owners=boundaryEdges[key];
                constructionStep(budget,0);
                owners.push_back({static_cast<Index>(fi), a < b ? 1 : -1});
            }
            const auto a=candidate.points[static_cast<std::size_t>(face[0])];
            const auto b=candidate.points[static_cast<std::size_t>(face[1])];
            const auto c=candidate.points[static_cast<std::size_t>(face[2])];
            // The determinant is affine in its fourth point. A supporting
            // half-space that contains the whole box contains every input.
            // If a synthetic corner is outside, test the actual point set.
            bool boxInside=true;
            for (unsigned corner=0; corner<8 && boxInside; ++corner) {
                constructionStep(budget);
                const Vec3 p{corner&1 ? hi.x : lo.x, corner&2 ? hi.y : lo.y,
                             corner&4 ? hi.z : lo.z};
                boxInside=!outside(a,b,c,p);
            }
            if (boxInside) continue;
            for (std::size_t pi = 0; pi < pointCount; ++pi) {
                constructionStep(budget);
                const Index point = static_cast<Index>(pi);
                if (point == face[0] || point == face[1] || point == face[2])
                    continue;
                if (outside(a,b,c,candidate.points[pi]))
                    return invalid(TetCandidateReason::BoundaryTopologyInvalid,
                                   static_cast<Index>(fi));
            }
        }
        constructionStep(budget,candidate.boundaryFaces.size());
        std::vector<std::vector<Index>> boundaryAdj(
            candidate.boundaryFaces.size());
        for (const auto& [edge, owners] : boundaryEdges) {
            constructionStep(budget);
            (void)edge;
            if (owners.size() != 2 || owners[0].second == owners[1].second)
                return invalid(TetCandidateReason::BoundaryTopologyInvalid);
            constructionStep(budget,0);
            boundaryAdj[static_cast<std::size_t>(owners[0].first)].push_back(
                owners[1].first);
            constructionStep(budget,0);
            boundaryAdj[static_cast<std::size_t>(owners[1].first)].push_back(
                owners[0].first);
        }
        const std::int64_t euler =
            static_cast<std::int64_t>(boundaryVertices.size()) -
            static_cast<std::int64_t>(boundaryEdges.size()) +
            static_cast<std::int64_t>(candidate.boundaryFaces.size());
        if (euler != 2)
            return invalid(TetCandidateReason::BoundaryTopologyInvalid);
        constructionStep(budget,candidate.boundaryFaces.size());
        std::vector<bool> boundaryReached(candidate.boundaryFaces.size(), false);
        constructionStep(budget,0);
        std::queue<Index> boundaryQueue;
        boundaryReached[0] = true;
        constructionStep(budget,0);
        boundaryQueue.push(0);
        std::size_t boundaryReachedCount = 0;
        while (!boundaryQueue.empty()) {
            constructionStep(budget);
            const Index face = boundaryQueue.front();
            boundaryQueue.pop();
            ++boundaryReachedCount;
            for (const Index neighbour :
                 boundaryAdj[static_cast<std::size_t>(face)]) {
                constructionStep(budget);
                const auto neighbourIndex =
                    static_cast<std::size_t>(neighbour);
                if (!boundaryReached[neighbourIndex]) {
                    boundaryReached[neighbourIndex] = true;
                    constructionStep(budget,0);
                    boundaryQueue.push(neighbour);
                }
            }
        }
        if (boundaryReachedCount != candidate.boundaryFaces.size())
            return invalid(TetCandidateReason::BoundaryTopologyInvalid);

        (void)interiorFaces;
        return TetCandidateValidation{true, TetCandidateReason::None,
                                      kInvalidIndex};
    } catch (const detail::TetConstructionStopped& stop) {
        return invalid(stop.cancelled ? TetCandidateReason::Cancelled : TetCandidateReason::BudgetExceeded);
    } catch (const std::bad_alloc&) {
        return invalid(TetCandidateReason::AllocationFailed);
    } catch (...) {
        return invalid(TetCandidateReason::ConstructionFailed);
    }
}

TetCandidateValidation validateTetCandidate(const TetCandidate& candidate) noexcept {
    std::uint64_t hullTests=0;
    return validateTetCandidateImpl(candidate,hullTests);
}

static TetCandidateBuildOutcome buildTetCandidateAttempt(
    const std::vector<Vec3>& input, bool exactCoordinates, Budget* budget=nullptr) noexcept {
    // Raw output and fixed candidate/validation envelopes are local to an
    // attempt. Traversal work and the high-water mark remain cumulative.
    struct MemoryScope {
        Budget* budget; std::uint64_t before;
        ~MemoryScope() { if (budget) budget->bytes=before; }
    } memoryScope{budget,budget ? budget->bytes : 0};
    TetCandidateBuildOutcome outcome;
    outcome.metrics.inputPoints = input.size();
    try {
        constructionStep(budget,0);
        for (std::size_t i = 0; i < input.size(); ++i) {
            constructionStep(budget);
            if (!finite(input[i])) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::NonfiniteInput;
                return outcome;
            }
        }

        auto& construction=outcome.metrics.construction;
        std::vector<Tet> raw=budget ?
            tet_candidate_detail::delaunay3DBounded(input,&construction,*budget,exactCoordinates) : exactCoordinates?
            tet_candidate_detail::delaunay3DExactInputCoordinates(input,&construction):delaunay3D(input,&construction);
        if (raw.empty()) {
            outcome.status = TetCandidateStatus::Rejected;
            outcome.reason = TetCandidateReason::NonVolumetricInput;
            return outcome;
        }

        if (budget) {
            budget->cells(raw.size());
            // Canonical vectors, all face uses, validation bitsets/queues,
            // boundary maps/sets and per-edge owner/adjacency vectors, with
            // conservative node and geometric-growth slack. Raw output was
            // admitted by Delaunay and remains charged during this overlap.
            budget->admit(sizeof(TetCandidate)+8192U);
            budget->admit(Budget::product(input.size(),256U));
            budget->admit(Budget::product(raw.size(),4096U));
            // Reserve slot-initialization work before vector/bitset growth.
            budget->step(Budget::product(input.size(),4U));
            budget->step(Budget::product(raw.size(),16U));
        }
        constructionStep(budget,0);
        auto candidate = std::make_unique<TetCandidate>();
        candidate->metrics.inputPoints = input.size();
        candidate->metrics.construction = construction;
        candidate->metrics.originalCoordinateRetryUsed=construction.originalCoordinateRetryUsed;
        if(construction.originalCoordinateRetryUsed&&!exactCoordinates)
            candidate->metrics.initialCoordinateFailure=TetCandidateReason::NonPositiveOrientation;
        constructionStep(budget,0);
        candidate->inputToPoint.resize(input.size(), kInvalidIndex);

        constructionStep(budget,0);
        std::vector<Index> order(input.size());
        for (std::size_t i = 0; i < input.size(); ++i) {
            constructionStep(budget); order[i] = static_cast<Index>(i);
        }
        std::sort(order.begin(), order.end(), [&](Index a, Index b) {
            constructionStep(budget);
            const Vec3& pa = input[static_cast<std::size_t>(a)];
            const Vec3& pb = input[static_cast<std::size_t>(b)];
            if (pointLess(pa, pb)) return true;
            if (pointLess(pb, pa)) return false;
            return a < b;
        });
        for (const Index inputIndex : order) {
            constructionStep(budget);
            const Vec3& point = input[static_cast<std::size_t>(inputIndex)];
            if (candidate->points.empty() ||
                !pointEqual(candidate->points.back(), point)) {
                constructionStep(budget,0);
                candidate->points.push_back(point);
                constructionStep(budget,0);
                candidate->pointToInput.push_back(inputIndex);
            }
            candidate->inputToPoint[static_cast<std::size_t>(inputIndex)] =
                static_cast<Index>(candidate->points.size() - 1);
        }

        constructionStep(budget,0);
        candidate->cells.reserve(raw.size());
        for (std::size_t source = 0; source < raw.size(); ++source) {
            constructionStep(budget);
            TetCandidateCell cell;
            for (int i = 0; i < 4; ++i) {
                constructionStep(budget);
                const Index inputIndex = raw[source].v[static_cast<std::size_t>(i)];
                if (inputIndex < 0 ||
                    static_cast<std::size_t>(inputIndex) >= input.size()) {
                    outcome.reason = TetCandidateReason::InvalidCellIndex;
                    return outcome;
                }
                cell.vertices[static_cast<std::size_t>(i)] =
                    candidate->inputToPoint[static_cast<std::size_t>(inputIndex)];
            }
            const Sign orientation = predicates::orient3d(
                candidate->points[static_cast<std::size_t>(cell.vertices[0])],
                candidate->points[static_cast<std::size_t>(cell.vertices[1])],
                candidate->points[static_cast<std::size_t>(cell.vertices[2])],
                candidate->points[static_cast<std::size_t>(cell.vertices[3])]);
            if (orientation == Sign::Negative)
                std::swap(cell.vertices[2], cell.vertices[3]);
            else if (orientation == Sign::Zero) {
                outcome.reason = TetCandidateReason::NonPositiveOrientation;
                return outcome;
            }
            cell.sourceOrdinal = static_cast<Index>(source);
            candidate->cells.push_back(cell);
        }
        std::sort(candidate->cells.begin(), candidate->cells.end(),
                  [budget](const TetCandidateCell& a, const TetCandidateCell& b) {
                      constructionStep(budget);
                      const CellKey ak = cellKey(a);
                      const CellKey bk = cellKey(b);
                      if (ak != bk) return ak < bk;
                      return a.vertices < b.vertices;
                  });

        std::vector<FaceUse> uses;
        constructionStep(budget,0);
        if (candidate->cells.size()<=uses.max_size()/4)
            uses.reserve(candidate->cells.size()*4);
        for (std::size_t i=0;i<candidate->cells.size();++i) {
            constructionStep(budget);
            appendFaceUses(uses,candidate->cells[i],static_cast<Index>(i));
        }
        sortFaceUses(uses,budget);
        for (std::size_t first=0;first<uses.size();) {
            constructionStep(budget);
            const auto end=nextFace(uses,first,budget);
            const auto& a=uses[first];
            auto& cellA=candidate->cells[static_cast<std::size_t>(a.cell)];
            if (end-first==1) {
                constructionStep(budget,0);
                candidate->boundaryFaces.push_back(
                    {directedFace(cellA,a.localFace),a.cell,a.localFace});
            } else if (end-first==2) {
                const auto& b=uses[first+1];
                cellA.neighbours[a.localFace]=b.cell;
                candidate->cells[static_cast<std::size_t>(b.cell)].neighbours[b.localFace]=a.cell;
                ++candidate->metrics.interiorFaces;
            } else {
                outcome.reason = TetCandidateReason::NonManifoldFace;
                return outcome;
            }
            first=end;
        }

        candidate->metrics.uniquePoints = candidate->points.size();
        candidate->metrics.candidateCells = candidate->cells.size();
        candidate->metrics.boundaryFaces = candidate->boundaryFaces.size();
        candidate->metrics.validationFaceUses = candidate->cells.size() * 4;
        candidate->metrics.candidateAccountedBytes =
            tetCandidateAccountedBytes(*candidate);
        candidate->metrics.constructionCandidateConservativeBytes =
            candidate->metrics.candidateAccountedBytes + construction.arenaBytes +
            static_cast<std::uint64_t>(raw.capacity() * sizeof(Tet));

        const TetCandidateValidation validation =
            validateTetCandidateImpl(*candidate,candidate->metrics.validationHullTests,&uses,budget);
        outcome.metrics=candidate->metrics;
        if (!validation.valid) {
            outcome.reason = validation.reason;
            return outcome;
        }
        constructionStep(budget,0);
        candidate->metrics.identityComplete = true;
        candidate->metrics.adjacencyComplete = true;
        candidate->metrics.boundaryComplete = true;
        candidate->metrics.runComplete = true;
        outcome.metrics = candidate->metrics;
        outcome.candidate = std::move(candidate);
        outcome.status = TetCandidateStatus::Built;
        outcome.reason = TetCandidateReason::None;
        return outcome;
    } catch (const detail::TetConstructionStopped& stop) {
        outcome.metrics.originalCoordinateRetryUsed=outcome.metrics.construction.originalCoordinateRetryUsed;
        outcome.reason=stop.cancelled ? TetCandidateReason::Cancelled : TetCandidateReason::BudgetExceeded;
        return outcome;
    } catch (const std::bad_alloc&) {
        outcome.metrics.originalCoordinateRetryUsed=outcome.metrics.construction.originalCoordinateRetryUsed;
        outcome.reason = TetCandidateReason::AllocationFailed;
        return outcome;
    } catch (...) {
        outcome.metrics.originalCoordinateRetryUsed=outcome.metrics.construction.originalCoordinateRetryUsed;
        outcome.reason = TetCandidateReason::ConstructionFailed;
        return outcome;
    }
}

static TetCandidateBuildOutcome buildTetCandidateImpl(
    const std::vector<Vec3>& input,Budget* budget) noexcept {
    auto initial=buildTetCandidateAttempt(input,false,budget);
    if(initial.built())return initial;
    bool retry=false;
    switch(initial.reason){
        case TetCandidateReason::NonVolumetricInput:
        case TetCandidateReason::MissingRequiredPoint:
        case TetCandidateReason::NonPositiveOrientation:
        case TetCandidateReason::NonManifoldFace:
        case TetCandidateReason::AdjacencyMismatch:
        case TetCandidateReason::DisconnectedComplex:
        case TetCandidateReason::BoundaryExtractionMismatch:
        case TetCandidateReason::BoundaryTopologyInvalid:
        case TetCandidateReason::ConstructionFailed:retry=true;break;
        default:break;
    }
    if(!retry)return initial;
    auto exact=buildTetCandidateAttempt(input,true,budget);
    tet_candidate_detail::retainDiscardedDelaunayAttempt(exact.metrics.construction,initial.metrics.construction);
    exact.metrics.originalCoordinateRetryUsed=true;
    exact.metrics.initialCoordinateFailure=initial.reason;
    exact.metrics.discardedCandidateWork=initial.metrics.discardedCandidateWork;
    for(const auto visits:{initial.metrics.inputPoints,initial.metrics.uniquePoints,
                          initial.metrics.candidateCells,initial.metrics.interiorFaces,
                          initial.metrics.boundaryFaces,initial.metrics.validationFaceUses,
                          initial.metrics.validationHullTests})
        saturatingAdd(exact.metrics.discardedCandidateWork,visits);
    if(exact.candidate)exact.candidate->metrics=exact.metrics;
    return exact;
}

TetCandidateBuildOutcome buildTetCandidate(const std::vector<Vec3>& input) noexcept {
    return buildTetCandidateImpl(input,nullptr);
}

TetCandidateBuildOutcome buildTetCandidate(
    const std::vector<Vec3>& input,const TetCandidateBuildLimits& limits) noexcept {
    Budget budget{limits.maxWork,limits.maxCells,limits.maxTemporaryBytes,limits.cancellation};
    TetCandidateBuildOutcome outcome;
    outcome.metrics.inputPoints=input.size();
    try {
        budget.step(0);
        // Existing exact-insphere scratch uses seven retained double vectors.
        // A term has at most 2*24*192+8 components; the accumulator has four
        // such terms. Four MiB covers these vectors, their growth overlap,
        // fixed predicate arrays and the private construction workspace.
        budget.admit(4U*1024U*1024U);
        outcome=buildTetCandidateImpl(input,&budget);
    } catch (const detail::TetConstructionStopped& stop) {
        outcome.reason=stop.cancelled ? TetCandidateReason::Cancelled : TetCandidateReason::BudgetExceeded;
    }
    outcome.metrics.boundedWork=budget.work;
    outcome.metrics.boundedPeakTemporaryBytes=budget.peakBytes;
    if (outcome.candidate) outcome.candidate->metrics=outcome.metrics;
    return outcome;
}

static TetHandoffOutcome importTetCandidateImpl(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification,
    std::span<const TetImportProtectedEdge> protectedEdges,
    std::span<const TetImportProtectedFace> protectedFaces,
    ProtectedBoundaryImportPolicy protectedPolicy) noexcept {
    TetHandoffOutcome outcome;
    try {
        TetHandoffMetrics partialImportMetrics{};
        bool partialMetricComplete=true;
        const auto failPartialMetricOverflow=[&]() noexcept {
            outcome.status=TetCandidateStatus::Failed;
            outcome.reason=TetCandidateReason::ConstructionFailed;
        };
        const auto validation = validateTetCandidate(candidate);
        if (!validation.valid) {
            outcome.reason = validation.reason;
            return outcome;
        }
        const std::size_t pointCount=candidate.points.size();
        const std::size_t cellCount=candidate.cells.size();
        const std::size_t boundaryCount=candidate.boundaryFaces.size();
        std::uint64_t faceUseCount=0;
        std::uint64_t boundaryCount64=0;
        const bool validationMetricsComplete=
            checkedMetricAdd(
                partialImportMetrics.importCandidateValidationCalls) &&
            checkedMetricProductSize(
                partialImportMetrics.importCandidateValidationPointRows,
                pointCount,3) &&
            checkedMetricAddSize(
                partialImportMetrics.importCandidateValidationPointRows,
                candidate.inputToPoint.size()) &&
            checkedMetricAddSize(
                partialImportMetrics.importCandidateValidationCellRows,
                cellCount) &&
            checkedMetricProductSize(faceUseCount,cellCount,4) &&
            checkedMetricAdd(
                partialImportMetrics.importCandidateValidationFaceUses,
                faceUseCount) &&
            checkedMetricAddSize(boundaryCount64,boundaryCount) &&
            checkedMetricProductSize(
                partialImportMetrics.importCandidateValidationBoundaryRows,
                boundaryCount,3) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationDualNeighbourVisits,
                cellCount,4) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationBoundaryEdgeUses,
                boundaryCount,3) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationBoundaryNeighbourVisits,
                boundaryCount,3) &&
            checkedMetricProductSize(
                partialImportMetrics.importCandidateValidationHullTests,
                boundaryCount,pointCount-3) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationSlotsInitialized,
                pointCount,2) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationSlotsInitialized,
                cellCount,2) &&
            checkedMetricProductSize(
                partialImportMetrics.
                    importCandidateValidationSlotsInitialized,
                boundaryCount,2);
        if (!validationMetricsComplete ||
            faceUseCount>(std::numeric_limits<std::uint64_t>::max)()-
                             boundaryCount64 ||
            ((faceUseCount+boundaryCount64)&UINT64_C(1))!=0 ||
            !checkedMetricAdd(
                partialImportMetrics.importCandidateValidationFaceRows,
                (faceUseCount+boundaryCount64)/2)) {
            failPartialMetricOverflow();
            return outcome;
        }
        if (existingGlobalByPoint.size() != candidate.points.size() ||
            expectedClassByPoint.size() != candidate.points.size()) {
            outcome.status = TetCandidateStatus::Rejected;
            outcome.reason = TetCandidateReason::InvalidGlobalMapping;
            return outcome;
        }
        if (cellClassification.dim != 3 || cellClassification.tag < 0) {
            outcome.status = TetCandidateStatus::Rejected;
            outcome.reason = TetCandidateReason::InvalidClassification;
            return outcome;
        }
        if (selectedCandidateCells.empty()) {
            outcome.status = TetCandidateStatus::Rejected;
            outcome.reason = TetCandidateReason::InvalidImportSelection;
            return outcome;
        }

        std::vector<Index> selected(selectedCandidateCells.begin(),
                                    selectedCandidateCells.end());
        if (!checkedMetricAddSize(
                partialImportMetrics.importSelectionRowsCopied,
                selected.size()) ||
            !checkedMetricAddSize(
                partialImportMetrics.importSelectionPointSlotsInitialized,
                candidate.points.size())) {
            failPartialMetricOverflow();
            return outcome;
        }
        std::sort(selected.begin(),selected.end(),
                  [&](Index a,Index b) noexcept {
            if (!checkedMetricAdd(
                    partialImportMetrics.importSelectionSortComparisons))
                partialMetricComplete=false;
            return a<b;
        });
        if (!partialMetricComplete) {
            failPartialMetricOverflow();
            return outcome;
        }
        const auto duplicateSelection=std::adjacent_find(
            selected.begin(),selected.end(),
            [&](Index a,Index b) noexcept {
                if (!checkedMetricAdd(
                        partialImportMetrics.
                            importSelectionDuplicateComparisons))
                    partialMetricComplete=false;
                return a==b;
            });
        if (!partialMetricComplete) {
            failPartialMetricOverflow();
            return outcome;
        }
        if (duplicateSelection != selected.end()) {
            outcome.status = TetCandidateStatus::Rejected;
            outcome.reason = TetCandidateReason::InvalidImportSelection;
            return outcome;
        }
        std::vector<bool> used(candidate.points.size(), false);
        for (const Index cell : selected) {
            if (!checkedMetricAdd(
                    partialImportMetrics.importSelectionCellVisits)) {
                failPartialMetricOverflow();
                return outcome;
            }
            if (cell < 0 ||
                static_cast<std::size_t>(cell) >= candidate.cells.size()) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidImportSelection;
                return outcome;
            }
            for (const Index point :
                 candidate.cells[static_cast<std::size_t>(cell)].vertices) {
                if (!checkedMetricAdd(
                        partialImportMetrics.importSelectionVertexVisits)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                used[static_cast<std::size_t>(point)] = true;
            }
        }

        std::set<Index,CountingIndexLess> seenGlobal(
            CountingIndexLess{
                &partialImportMetrics.importGlobalMapComparisons,
                &partialMetricComplete});
        for (std::size_t point = 0; point < used.size(); ++point) {
            if (!used[point]) continue;
            if (!checkedMetricAdd(
                    partialImportMetrics.importGlobalMapRows)) {
                failPartialMetricOverflow();
                return outcome;
            }
            const Index global = existingGlobalByPoint[point];
            if (global < 0 || global >= globalMesh.nodeCount()) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            const auto [globalIt,inserted]=seenGlobal.insert(global);
            (void)globalIt;
            if (!partialMetricComplete) {
                failPartialMetricOverflow();
                return outcome;
            }
            if (!inserted) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            if (!checkedMetricAdd(
                    partialImportMetrics.importGlobalMapInsertions)) {
                failPartialMetricOverflow();
                return outcome;
            }
            const Vec3& actualPoint = globalMesh.node(global);
            const Vec3& expectedPoint = candidate.points[point];
            const mesh::Classification actualClass =
                globalMesh.nodeClass(global);
            const mesh::Classification expectedClass =
                expectedClassByPoint[point];
            if (actualClass.dim < 0 || actualClass.dim > 3 ||
                actualClass.tag < 0 || expectedClass.dim < 0 ||
                expectedClass.dim > 3 || expectedClass.tag < 0) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidClassification;
                return outcome;
            }
            if (actualPoint.x != expectedPoint.x ||
                actualPoint.y != expectedPoint.y ||
                actualPoint.z != expectedPoint.z ||
                actualClass.dim != expectedClass.dim ||
                actualClass.tag != expectedClass.tag) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
        }

        // The partial v4 import has a deliberately stronger preflight than
        // the established all-or-nothing entry: every authoritative row is
        // canonicalized and validated before any private pool is built or
        // any metadata is applied. This keeps an invalid later row from
        // being mistaken for a merely absent constraint.
        std::vector<TetImportProtectedEdge> partialProtectedEdges;
        std::vector<TetImportProtectedFace> partialProtectedFaces;
        std::uint64_t partialPreflightMappingBytes=0;
        if (protectedPolicy ==
            ProtectedBoundaryImportPolicy::ExactPresentOnly) {
            partialProtectedEdges.assign(protectedEdges.begin(),
                                         protectedEdges.end());
            partialProtectedFaces.assign(protectedFaces.begin(),
                                         protectedFaces.end());
            if (!checkedMetricAddSize(
                    partialImportMetrics.partialImportRowsCopied,
                    partialProtectedEdges.size()) ||
                !checkedMetricAddSize(
                    partialImportMetrics.partialImportRowsCopied,
                    partialProtectedFaces.size())) {
                failPartialMetricOverflow();
                return outcome;
            }
            for (auto& row : partialProtectedEdges) {
                if (!checkedMetricAdd(
                        partialImportMetrics.partialImportEdgeValidationRows)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                std::sort(row.existingGlobalNodes.begin(),
                          row.existingGlobalNodes.end(),
                          [&](Index a,Index b) noexcept {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportSortComparisons))
                        partialMetricComplete=false;
                    return a<b;
                });
                if (!partialMetricComplete) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto& key = row.existingGlobalNodes;
                if (key[0] < 0 || key[1] >= globalMesh.nodeCount() ||
                    key[0] == key[1]) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                    return outcome;
                }
                if ((row.classification.dim != 1 &&
                     row.classification.dim != 2) ||
                    row.classification.tag < 0) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidClassification;
                    return outcome;
                }
            }
            std::sort(partialProtectedEdges.begin(),
                      partialProtectedEdges.end(),
                      [&](const auto& a,const auto& b) noexcept {
                          if (!checkedMetricAdd(
                                  partialImportMetrics.
                                      partialImportSortComparisons))
                              partialMetricComplete=false;
                          return a.existingGlobalNodes <
                                 b.existingGlobalNodes;
                      });
            if (!partialMetricComplete) {
                failPartialMetricOverflow();
                return outcome;
            }
            const auto duplicateEdge=std::adjacent_find(
                    partialProtectedEdges.begin(),
                    partialProtectedEdges.end(),
                    [&](const auto& a,const auto& b) noexcept {
                        if (!checkedMetricAdd(
                                partialImportMetrics.
                                    partialImportDuplicateKeyComparisons))
                            partialMetricComplete=false;
                        return a.existingGlobalNodes ==
                               b.existingGlobalNodes;
                    });
            if (!partialMetricComplete) {
                failPartialMetricOverflow();
                return outcome;
            }
            if (duplicateEdge != partialProtectedEdges.end()) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }

            for (auto& row : partialProtectedFaces) {
                if (!checkedMetricAdd(
                        partialImportMetrics.partialImportFaceValidationRows)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                std::sort(row.existingGlobalNodes.begin(),
                          row.existingGlobalNodes.end(),
                          [&](Index a,Index b) noexcept {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportSortComparisons))
                        partialMetricComplete=false;
                    return a<b;
                });
                if (!partialMetricComplete) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto& key = row.existingGlobalNodes;
                if (key[0] < 0 || key[2] >= globalMesh.nodeCount() ||
                    key[0] == key[1] || key[1] == key[2]) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                    return outcome;
                }
                if (row.classification.dim != 2 ||
                    row.classification.tag < 0) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidClassification;
                    return outcome;
                }
            }
            std::sort(partialProtectedFaces.begin(),
                      partialProtectedFaces.end(),
                      [&](const auto& a,const auto& b) noexcept {
                          if (!checkedMetricAdd(
                                  partialImportMetrics.
                                      partialImportSortComparisons))
                              partialMetricComplete=false;
                          return a.existingGlobalNodes <
                                 b.existingGlobalNodes;
                      });
            if (!partialMetricComplete) {
                failPartialMetricOverflow();
                return outcome;
            }
            const auto duplicateFace=std::adjacent_find(
                    partialProtectedFaces.begin(),
                    partialProtectedFaces.end(),
                    [&](const auto& a,const auto& b) noexcept {
                        if (!checkedMetricAdd(
                                partialImportMetrics.
                                    partialImportDuplicateKeyComparisons))
                            partialMetricComplete=false;
                        return a.existingGlobalNodes ==
                               b.existingGlobalNodes;
                    });
            if (!partialMetricComplete) {
                failPartialMetricOverflow();
                return outcome;
            }
            if (duplicateFace != partialProtectedFaces.end()) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }

            // Each row identity must resolve injectively to the supplied
            // candidate/global mapping even when the selected Tet subset does
            // not currently materialize that endpoint. Missing topology is
            // recoverable; missing or ambiguous authority is not.
            const auto nodeCount =
                static_cast<std::size_t>(globalMesh.nodeCount());
            constexpr auto kNotReferenced =
                (std::numeric_limits<std::size_t>::max)();
            constexpr auto kUnresolved = kNotReferenced - 1;
            std::vector<std::size_t> pointByProtectedGlobal(
                nodeCount,kNotReferenced);
            if (!checkedMetricAddSize(
                    partialImportMetrics.partialImportGlobalSlotsInitialized,
                    pointByProtectedGlobal.size())) {
                failPartialMetricOverflow();
                return outcome;
            }
            partialPreflightMappingBytes=vectorBytes(
                pointByProtectedGlobal.capacity(),sizeof(std::size_t));
            for (const auto& row : partialProtectedEdges) {
                for (const Index global : row.existingGlobalNodes) {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportProtectedEndpointVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    pointByProtectedGlobal[static_cast<std::size_t>(global)] =
                        kUnresolved;
                }
            }
            for (const auto& row : partialProtectedFaces) {
                for (const Index global : row.existingGlobalNodes) {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportProtectedEndpointVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    pointByProtectedGlobal[static_cast<std::size_t>(global)] =
                        kUnresolved;
                }
            }
            for (std::size_t point=0;
                 point<existingGlobalByPoint.size();++point) {
                if (!checkedMetricAdd(
                        partialImportMetrics.
                            partialImportCandidateMappingVisits)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const Index global=existingGlobalByPoint[point];
                if (global<0||global>=globalMesh.nodeCount()) continue;
                auto& mapped=pointByProtectedGlobal[
                    static_cast<std::size_t>(global)];
                if (mapped==kNotReferenced) continue;
                if (mapped!=kUnresolved) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                    return outcome;
                }
                if (!checkedMetricAdd(
                        partialImportMetrics.
                            partialImportAuthorityPointValidations)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto actualClass=globalMesh.nodeClass(global);
                const auto expectedClass=expectedClassByPoint[point];
                const auto& actualPoint=globalMesh.node(global);
                const auto& expectedPoint=candidate.points[point];
                if (actualClass.dim<0 || actualClass.dim>3 ||
                    actualClass.tag<0 || expectedClass.dim<0 ||
                    expectedClass.dim>3 || expectedClass.tag<0) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason =
                        TetCandidateReason::InvalidClassification;
                    return outcome;
                }
                if (actualPoint.x!=expectedPoint.x ||
                    actualPoint.y!=expectedPoint.y ||
                    actualPoint.z!=expectedPoint.z ||
                    actualClass.dim!=expectedClass.dim ||
                    actualClass.tag!=expectedClass.tag) {
                    outcome.status = TetCandidateStatus::Rejected;
                    outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                    return outcome;
                }
                mapped=point;
            }
            for (const auto& row : partialProtectedEdges) {
                for (const Index global : row.existingGlobalNodes) {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportProtectedEndpointVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    if (pointByProtectedGlobal[
                            static_cast<std::size_t>(global)]==kUnresolved) {
                        outcome.status = TetCandidateStatus::Rejected;
                        outcome.reason =
                            TetCandidateReason::InvalidGlobalMapping;
                        return outcome;
                    }
                }
            }
            for (const auto& row : partialProtectedFaces) {
                for (const Index global : row.existingGlobalNodes) {
                    if (!checkedMetricAdd(
                            partialImportMetrics.
                                partialImportProtectedEndpointVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    if (pointByProtectedGlobal[
                            static_cast<std::size_t>(global)]==kUnresolved) {
                        outcome.status = TetCandidateStatus::Rejected;
                        outcome.reason =
                            TetCandidateReason::InvalidGlobalMapping;
                        return outcome;
                    }
                }
            }
        }

        auto result = std::make_unique<TetAuthoritativeResult>();
        result->metrics=partialImportMetrics;
        result->metrics.selectedCells = selected.size();
        result->metrics.candidateAccountedBytes =
            tetCandidateAccountedBytes(candidate);
        topo::NodePool nodes;
        std::vector<topo::NodeHandle> localByPoint(candidate.points.size());
        std::vector<topo::NodeHandle> localByExistingGlobal(
            static_cast<std::size_t>(globalMesh.nodeCount()));
        if (!checkedMetricAddSize(
                result->metrics.importLocalPointSlotsInitialized,
                localByPoint.size()) ||
            !checkedMetricAddSize(
                result->metrics.importLocalGlobalSlotsInitialized,
                localByExistingGlobal.size())) {
            failPartialMetricOverflow();
            return outcome;
        }
        result->nodeBindings.reserve(candidate.points.size());
        for (std::size_t point = 0; point < used.size(); ++point) {
            if (!used[point]) continue;
            const Index global = existingGlobalByPoint[point];
            const auto local = nodes.add(candidate.points[point],
                                         globalMesh.nodeClass(global));
            localByPoint[point] = local;
            localByExistingGlobal[static_cast<std::size_t>(global)] = local;
            result->nodeBindings.push_back(
                {static_cast<Index>(point), existingGlobalByPoint[point],
                 local});
        }

        topo::CellPool cells;
        result->localCellByCandidate.resize(candidate.cells.size());
        if (!checkedMetricAddSize(
                result->metrics.importLocalCellSlotsInitialized,
                result->localCellByCandidate.size())) {
            failPartialMetricOverflow();
            return outcome;
        }
        result->selectedCellBindings.reserve(selected.size());
        for (const Index candidateCell : selected) {
            const auto& source =
                candidate.cells[static_cast<std::size_t>(candidateCell)];
            const std::array<topo::NodeHandle, 4> connectivity{
                localByPoint[static_cast<std::size_t>(source.vertices[0])],
                localByPoint[static_cast<std::size_t>(source.vertices[1])],
                localByPoint[static_cast<std::size_t>(source.vertices[2])],
                localByPoint[static_cast<std::size_t>(source.vertices[3])]};
            const auto local = cells.add(
                topo::ElementType::Tet4,
                std::span<const topo::NodeHandle>(connectivity.data(),
                                                  connectivity.size()),
                cellClassification, nodes);
            result->localCellByCandidate[static_cast<std::size_t>(candidateCell)] =
                local;
            result->selectedCellBindings.push_back({candidateCell, local});
        }
        topo::EdgePoolCompleteBuildMetrics edgeBuild;
        auto edges = topo::EdgePool::buildComplete(
            nodes,cells,&edgeBuild);
        topo::FacePoolCompleteBuildMetrics faceBuild;
        auto faces = topo::FacePool::buildComplete(
            nodes,cells,edges,&faceBuild);

        std::vector<std::pair<topo::EdgeHandle,mesh::Classification>>
            resolvedProtectedEdges;
        std::vector<std::pair<topo::FaceHandle,mesh::Classification>>
            resolvedProtectedFaces;
        resolvedProtectedEdges.reserve(protectedEdges.size());
        resolvedProtectedFaces.reserve(protectedFaces.size());
        std::set<std::array<Index,2>> seenProtectedEdges;
        std::set<std::array<Index,3>> seenProtectedFaces;
        if (protectedPolicy ==
            ProtectedBoundaryImportPolicy::ExactPresentOnly) {
            for (const auto& row : partialProtectedEdges) {
                const auto& key=row.existingGlobalNodes;
                const auto a=localByExistingGlobal[
                    static_cast<std::size_t>(key[0])];
                const auto b=localByExistingGlobal[
                    static_cast<std::size_t>(key[1])];
                if (!checkedMetricAdd(
                        result->metrics.partialImportResolutionNodeVisits)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const bool aCurrent=nodes.isCurrent(a);
                bool bCurrent=false;
                if (aCurrent) {
                    if (!checkedMetricAdd(
                            result->metrics.
                                partialImportResolutionNodeVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    bCurrent=nodes.isCurrent(b);
                }
                if (!aCurrent||!bCurrent) continue;
                if (!checkedMetricAdd(
                        result->metrics.partialImportEdgeResolutionLookups)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto edge=edges.find(nodes,a,b);
                if (!edge) continue;
                if (!checkedMetricAdd(
                        result->metrics.partialImportEdgeDataReads)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto data=edges.get(*edge);
                if (data.useCount==0) continue;
                resolvedProtectedEdges.push_back(
                    {*edge,row.classification});
            }
            for (const auto& row : partialProtectedFaces) {
                const auto& key=row.existingGlobalNodes;
                const std::array<topo::NodeHandle,3> local{
                    localByExistingGlobal[static_cast<std::size_t>(key[0])],
                    localByExistingGlobal[static_cast<std::size_t>(key[1])],
                    localByExistingGlobal[static_cast<std::size_t>(key[2])]};
                bool allCurrent=true;
                for (const auto node:local) {
                    if (!allCurrent) break;
                    if (!checkedMetricAdd(
                            result->metrics.
                                partialImportResolutionNodeVisits)) {
                        failPartialMetricOverflow();
                        return outcome;
                    }
                    allCurrent=nodes.isCurrent(node);
                }
                if (!allCurrent) continue;
                if (!checkedMetricAdd(
                        result->metrics.partialImportFaceResolutionLookups)) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                const auto face=faces.find(nodes,local);
                if (!face) continue;
                const auto data=faces.get(*face);
                if (!checkedMetricAdd(
                        result->metrics.partialImportFaceDataCopies) ||
                    !checkedMetricAddSize(
                        result->metrics.
                            partialImportFaceIncidenceRowsCopied,
                        data.incident.size())) {
                    failPartialMetricOverflow();
                    return outcome;
                }
                if (data.shape!=topo::FaceShape::Tri||
                    data.cornerCount!=3||data.incident.size()!=2||
                    data.incident[0].sense*data.incident[1].sense!=-1)
                    continue;
                resolvedProtectedFaces.push_back(
                    {*face,row.classification});
            }
        } else {
        for (const auto& row : protectedEdges) {
            auto key = row.existingGlobalNodes;
            std::sort(key.begin(),key.end());
            if (key[0] < 0 || key[1] >= globalMesh.nodeCount() ||
                key[0] == key[1]) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            if ((row.classification.dim != 1 &&
                 row.classification.dim != 2) ||
                row.classification.tag < 0) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidClassification;
                return outcome;
            }
            if (!seenProtectedEdges.insert(key).second) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            const auto a = localByExistingGlobal[
                static_cast<std::size_t>(key[0])];
            const auto b = localByExistingGlobal[
                static_cast<std::size_t>(key[1])];
            if (!nodes.isCurrent(a) || !nodes.isCurrent(b)) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            const auto edge = edges.find(nodes,a,b);
            if (!edge) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidImportSelection;
                return outcome;
            }
            resolvedProtectedEdges.push_back({*edge,row.classification});
        }
        for (const auto& row : protectedFaces) {
            auto key = row.existingGlobalNodes;
            std::sort(key.begin(),key.end());
            if (key[0] < 0 || key[2] >= globalMesh.nodeCount() ||
                key[0] == key[1] || key[1] == key[2]) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            if (row.classification.dim != 2 ||
                row.classification.tag < 0) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidClassification;
                return outcome;
            }
            if (!seenProtectedFaces.insert(key).second) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            const std::array<topo::NodeHandle,3> local{
                localByExistingGlobal[static_cast<std::size_t>(key[0])],
                localByExistingGlobal[static_cast<std::size_t>(key[1])],
                localByExistingGlobal[static_cast<std::size_t>(key[2])]};
            if (!nodes.isCurrent(local[0]) || !nodes.isCurrent(local[1]) ||
                !nodes.isCurrent(local[2])) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidGlobalMapping;
                return outcome;
            }
            const auto face = faces.find(nodes,local);
            if (!face) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::InvalidImportSelection;
                return outcome;
            }
            const auto data = faces.get(*face);
            if (data.incident.size() != 2 ||
                data.incident[0].sense*data.incident[1].sense != -1) {
                outcome.status = TetCandidateStatus::Rejected;
                outcome.reason = TetCandidateReason::NonManifoldFace;
                return outcome;
            }
            resolvedProtectedFaces.push_back({*face,row.classification});
        }
        }
        for (const auto& [edge,classification] : resolvedProtectedEdges) {
            edges.setClassification(edge,classification);
            edges.setRequiredConstraint(edge,true);
        }
        for (const auto& [face,classification] : resolvedProtectedFaces)
            faces.setClassification(face,classification);
        result->workspace = std::make_unique<topo::TopologyWorkspace>(
            std::move(nodes), std::move(cells), std::move(edges),
            std::move(faces));

        std::size_t maxSlot = 0;
        for (const auto& binding : result->nodeBindings)
            maxSlot = std::max(maxSlot,
                               static_cast<std::size_t>(binding.local.slot));
        if (!result->nodeBindings.empty())
            result->existingGlobalByLocalNodeSlot.assign(
                maxSlot + 1, kInvalidIndex);
        if (!checkedMetricAddSize(
                result->metrics.importPersistentNodeSlotsInitialized,
                result->existingGlobalByLocalNodeSlot.size())) {
            failPartialMetricOverflow();
            return outcome;
        }
        for (const auto& binding : result->nodeBindings)
            result->existingGlobalByLocalNodeSlot[
                static_cast<std::size_t>(binding.local.slot)] =
                    binding.existingGlobalIndex;

        const auto memory =
            topo::measureWorkspaceMemory(*result->workspace);
        result->metrics.importedNodes = result->nodeBindings.size();
        result->metrics.importedCells = result->selectedCellBindings.size();
        result->metrics.importEdgeValidationCellVisits =
            edgeBuild.validationCellVisits;
        result->metrics.importEdgeBuildCellVisits = edgeBuild.buildCellVisits;
        result->metrics.importEdgeUses = edgeBuild.edgeUses;
        result->metrics.importUniqueEdges = edgeBuild.uniqueEdges;
        result->metrics.importFaceValidationCellVisits =
            faceBuild.validationCellVisits;
        result->metrics.importFaceBuildCellVisits = faceBuild.buildCellVisits;
        result->metrics.importFaceUses = faceBuild.faceUses;
        result->metrics.importUniqueFaces = faceBuild.uniqueFaces;
        result->metrics.importEdgeAudits = edgeBuild.finalEdgeAudits
            + faceBuild.suppliedEdgeAudits;
        result->metrics.importFaceAudits = faceBuild.finalFaceAudits;
        result->metrics.importPerCellGlobalAudits =
            edgeBuild.perCellGlobalAudits + faceBuild.perCellGlobalAudits;
        result->metrics.authoritativeModeledStorageBytes = memory.totalBytes;
        std::uint64_t scratch =
            vectorBytes(selected.capacity(), sizeof(Index));
        saturatingAdd(scratch,
                      vectorBytes(used.capacity(), sizeof(bool)));
        saturatingAdd(scratch,
                      vectorBytes(localByPoint.capacity(),
                                  sizeof(topo::NodeHandle)));
        saturatingAdd(scratch,
                      vectorBytes(localByExistingGlobal.capacity(),
                                  sizeof(topo::NodeHandle)));
        saturatingAdd(scratch,
                      vectorBytes(resolvedProtectedEdges.capacity(),
                                  sizeof(decltype(resolvedProtectedEdges)::
                                      value_type)));
        saturatingAdd(scratch,
                      vectorBytes(resolvedProtectedFaces.capacity(),
                                  sizeof(decltype(resolvedProtectedFaces)::
                                      value_type)));
        saturatingAdd(scratch,
                      vectorBytes(partialProtectedEdges.capacity(),
                                  sizeof(TetImportProtectedEdge)));
        saturatingAdd(scratch,
                      vectorBytes(partialProtectedFaces.capacity(),
                                  sizeof(TetImportProtectedFace)));
        saturatingAdd(scratch,partialPreflightMappingBytes);
        result->metrics.importScratchModeledPayloadBytes = scratch;
        result->metrics.modeledOverlapPayloadBytes =
            result->metrics.candidateAccountedBytes;
        saturatingAdd(result->metrics.modeledOverlapPayloadBytes,
                      result->metrics.authoritativeModeledStorageBytes);
        saturatingAdd(result->metrics.modeledOverlapPayloadBytes,
                      result->metrics.importScratchModeledPayloadBytes);
        // At return the caller-owned candidate and the new authoritative
        // workspace coexist. No third complete Tet connectivity is retained
        // by the handoff result; publication is derived through local handles.
        result->metrics.completeTetRepresentationsAtOverlap = 2;
        result->metrics.completeTetRepresentationsBeforePublication = 2;
        result->metrics.workspaceValid =
            result->workspace->adjacencyReady() &&
            result->workspace->state() == topo::WorkspaceState::Ready;
        result->metrics.workAccountingComplete=true;
        result->metrics.runComplete = result->metrics.workspaceValid &&
            result->metrics.workAccountingComplete;
        if (!result->metrics.runComplete) {
            outcome.reason = TetCandidateReason::WorkspaceImportFailed;
            return outcome;
        }
        outcome.status = TetCandidateStatus::Built;
        outcome.reason = TetCandidateReason::None;
        outcome.result = std::move(result);
        return outcome;
    } catch (const std::bad_alloc&) {
        outcome.reason = TetCandidateReason::AllocationFailed;
        return outcome;
    } catch (...) {
        outcome.reason = TetCandidateReason::WorkspaceImportFailed;
        return outcome;
    }
}

TetHandoffOutcome importTetCandidate(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification) noexcept {
    return importTetCandidateImpl(
        candidate,selectedCandidateCells,existingGlobalByPoint,
        expectedClassByPoint,globalMesh,cellClassification,{},{},
        ProtectedBoundaryImportPolicy::RequireComplete);
}

TetHandoffOutcome importTetCandidateWithProtectedBoundary(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification,
    std::span<const TetImportProtectedEdge> protectedEdges,
    std::span<const TetImportProtectedFace> protectedFaces) noexcept {
    return importTetCandidateImpl(
        candidate,selectedCandidateCells,existingGlobalByPoint,
        expectedClassByPoint,globalMesh,cellClassification,
        protectedEdges,protectedFaces,
        ProtectedBoundaryImportPolicy::RequireComplete);
}

TetHandoffOutcome tet_candidate_detail::
importTetCandidateWithPartialProtectedBoundary(
    const TetCandidate& candidate,
    std::span<const Index> selectedCandidateCells,
    std::span<const Index> existingGlobalByPoint,
    std::span<const mesh::Classification> expectedClassByPoint,
    const mesh::Mesh& globalMesh,
    mesh::Classification cellClassification,
    std::span<const TetImportProtectedEdge> protectedEdges,
    std::span<const TetImportProtectedFace> protectedFaces) noexcept {
    return importTetCandidateImpl(
        candidate,selectedCandidateCells,existingGlobalByPoint,
        expectedClassByPoint,globalMesh,cellClassification,
        protectedEdges,protectedFaces,
        ProtectedBoundaryImportPolicy::ExactPresentOnly);
}

}  // namespace paacs::algo
