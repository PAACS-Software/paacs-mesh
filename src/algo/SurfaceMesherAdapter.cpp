// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4100)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "paacs/algo/SurfaceMesherAdapter.hpp"

#include "T2CSurfaceCache.hpp"
#include "PreparedSurfaceCandidate.hpp"
#include "ProductionSurfaceMesherDetail.hpp"
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace paacs::algo {

bool validSurfaceNodeIdentity(const SurfaceNodeIdentity& identity) noexcept {
    const std::array<Index,3> noFacet{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    const auto canonicalFiniteBits = [](std::uint64_t bits) noexcept {
        const double value=std::bit_cast<double>(bits);
        return std::isfinite(value) &&
            (value!=0.0 || bits==std::bit_cast<std::uint64_t>(0.0));
    };
    const bool noParent=identity.parentArtifact==ArtifactStamp{};
    const bool noOrdinal=identity.parentOrdinal==0 &&
                         identity.parentSubordinal==0;
    const bool noParameters=identity.parameterUBits==0 &&
                            identity.parameterVBits==0;
    const bool noDyadic=identity.dyadicNumerator==0 &&
                         identity.dyadicLevel==0;
    const bool canonicalFacet=
        identity.parentFacetGlobalNodes[0]>=0 &&
        identity.parentFacetGlobalNodes[0]<
            identity.parentFacetGlobalNodes[1] &&
        identity.parentFacetGlobalNodes[1]<
            identity.parentFacetGlobalNodes[2];

    switch(identity.domain) {
        case SurfaceNodeIdentityDomain::OriginalDiscretization:
            return noParent && identity.parentFacetGlobalNodes==noFacet &&
                   noOrdinal && noParameters && noDyadic;
        case SurfaceNodeIdentityDomain::RefinedCurve: {
            if(identity.parentArtifact.persistentId==0 ||
               identity.parentFacetGlobalNodes!=noFacet ||
               identity.parentSubordinal!=0 ||
               !canonicalFiniteBits(identity.parameterUBits) ||
               identity.parameterVBits!=0 || identity.dyadicLevel==0 ||
               identity.dyadicLevel>63 ||
               (identity.dyadicNumerator&1U)==0)
                return false;
            const std::uint64_t denominator=
                std::uint64_t{1}<<identity.dyadicLevel;
            return identity.dyadicNumerator<denominator;
        }
        case SurfaceNodeIdentityDomain::RefinedFaceEdge:
        case SurfaceNodeIdentityDomain::RefinedFacetInterior:
            return identity.parentArtifact.persistentId!=0 &&
                   canonicalFacet && noOrdinal &&
                   canonicalFiniteBits(identity.parameterUBits) &&
                   canonicalFiniteBits(identity.parameterVBits) && noDyadic;
        case SurfaceNodeIdentityDomain::InterfacePrefill:
            return identity.parentArtifact.persistentId!=0 &&
                   identity.parentFacetGlobalNodes==noFacet &&
                   identity.parentOrdinal!=0 &&
                   identity.parentSubordinal!=0 &&
                   canonicalFiniteBits(identity.parameterUBits) &&
                   canonicalFiniteBits(identity.parameterVBits) && noDyadic;
    }
    return false;
}

namespace {

bool sameClass(const mesh::Classification& a,
               const mesh::Classification& b) noexcept {
    return a.dim == b.dim && a.tag == b.tag;
}
bool samePoint(const Vec3& a, const Vec3& b) noexcept {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
bool sameParam(const Vec2& a, const Vec2& b) noexcept {
    return a.u == b.u && a.v == b.v;
}
bool finite(const Vec2& p) noexcept {
    return std::isfinite(p.u) && std::isfinite(p.v);
}
bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

bool identityMatchesFaceParameter(const SurfaceNodeIdentity& identity,
                                  const Vec2& parameter) noexcept {
    if(identity.domain==SurfaceNodeIdentityDomain::RefinedCurve)return true;
    const auto bits=[](double value) noexcept {
        return std::bit_cast<std::uint64_t>(value==0.0?0.0:value);
    };
    return identity.parameterUBits==bits(parameter.u) &&
           identity.parameterVBits==bits(parameter.v);
}

SurfaceFailure failure(SurfaceReason reason) noexcept {
    SurfaceFailure f;
    f.reason = reason;
    return f;
}
SurfaceFailure allocation(SurfaceAllocationPhase phase) noexcept {
    SurfaceFailure f = failure(SurfaceReason::AllocationFailed);
    f.allocation = phase;
    return f;
}
SurfaceFailure candidateFailure(SurfaceCandidateDetail detail) noexcept {
    SurfaceFailure f = failure(SurfaceReason::CandidateValidationFailed);
    f.candidate = detail;
    return f;
}
SurfaceFailure domainFailure(SurfaceDomainDetail detail) noexcept {
    SurfaceFailure f = failure(SurfaceReason::DomainExtractionFailed);
    f.domain = detail;
    return f;
}

using EdgeIndex = std::array<std::size_t, 2>;
EdgeIndex edgeIndex(std::size_t a, std::size_t b) noexcept {
    return a < b ? EdgeIndex{a, b} : EdgeIndex{b, a};
}
using HandleEdge = std::array<topo::NodeHandle, 2>;
HandleEdge handleEdge(topo::NodeHandle a, topo::NodeHandle b) noexcept {
    return a < b ? HandleEdge{a, b} : HandleEdge{b, a};
}

struct ChainRecord {
    SurfaceChain source;
    SurfaceChain chain;
    std::optional<std::size_t> loopOrdinal;
    std::vector<std::size_t> originalOrdinals;
};
struct FacePlan {
    const SurfaceFaceRequest* request = nullptr;
    std::map<FaceLocalVertexId, const SurfaceVertexUse*> vertices;
    std::vector<ChainRecord> chains;
    std::vector<std::vector<FaceLocalVertexId>> loops;
    std::vector<bnd::ConstraintSegment> segments;
    std::map<bnd::ConstraintEntityId, SurfaceSegmentOrigin> origins;
    std::set<bnd::ConstraintEntityId> protectedIds;
};
struct SetPlan {
    std::vector<FacePlan> faces;
    std::map<SharedCurveNodeKey, SharedCurveNodeBinding> ledger;
};

std::vector<SurfaceChainVertex> canonicalVertices(const SurfaceChain& chain,
                                                   bool& reversed,
                                                   std::size_t& rotation) {
    std::vector<SurfaceChainVertex> base = chain.uses;
    if (chain.closed && base.size() > 1 &&
        base.front().faceLocalVertexId == base.back().faceLocalVertexId &&
        base.front().sharedKey == base.back().sharedKey)
        base.pop_back();
    reversed = false;
    rotation = 0;
    auto lessUse = [](const SurfaceChainVertex& a,
                      const SurfaceChainVertex& b) {
        if (a.sharedKey && b.sharedKey) {
            if (*a.sharedKey < *b.sharedKey) return true;
            if (*b.sharedKey < *a.sharedKey) return false;
        } else if (a.sharedKey.has_value() != b.sharedKey.has_value()) {
            return a.sharedKey.has_value();
        }
        return a.faceLocalVertexId < b.faceLocalVertexId;
    };
    auto sequenceLess = [&](const std::vector<SurfaceChainVertex>& a,
                            const std::vector<SurfaceChainVertex>& b) {
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(),
                                            b.end(), lessUse);
    };
    if (!chain.closed) {
        auto rev = base;
        std::reverse(rev.begin(), rev.end());
        if (sequenceLess(rev, base)) {
            base.swap(rev);
            reversed = true;
        }
        return base;
    }
    std::vector<SurfaceChainVertex> best;
    for (int direction = 0; direction < 2; ++direction) {
        auto seq = base;
        if (direction != 0) std::reverse(seq.begin(), seq.end());
        for (std::size_t r = 0; r < seq.size(); ++r) {
            auto candidate = seq;
            std::rotate(candidate.begin(), candidate.begin() +
                            static_cast<std::ptrdiff_t>(r), candidate.end());
            if (best.empty() || sequenceLess(candidate, best)) {
                best = std::move(candidate);
                reversed = direction != 0;
                rotation = r;
            }
        }
    }
    return best;
}

bool validPolicy(const SurfacePhysicalPolicy& p) noexcept {
    return std::isfinite(p.minPhysicalDoubleArea) &&
           std::isfinite(p.maxVertexMappingError) &&
           std::isfinite(p.maxMappingDeviation) &&
           std::isfinite(p.minNormalAlignment) &&
           std::isfinite(p.maxArtifactPositionError) &&
           p.minPhysicalDoubleArea > 0.0 && p.maxVertexMappingError >= 0.0 &&
           p.maxMappingDeviation >= 0.0 && p.minNormalAlignment > 0.0 &&
           p.minNormalAlignment <= 1.0 && p.maxArtifactPositionError >= 0.0;
}

predicates::Sign ringWinding(const std::vector<FaceLocalVertexId>& ring,
                             const FacePlan& plan);
SurfaceFailure rebuildDerivedConstraints(FacePlan& plan);
SurfaceFailure normalizeFacePlan(FacePlan& plan,
                                 const ChartCertificate& certificate);

bool pointOnSegment(const Vec2& a, const Vec2& b, const Vec2& p) noexcept {
    if (predicates::orient2d(a, b, p) != predicates::Sign::Zero) return false;
    return p.u >= std::min(a.u, b.u) && p.u <= std::max(a.u, b.u) &&
           p.v >= std::min(a.v, b.v) && p.v <= std::max(a.v, b.v);
}

bool pointInClosure(const std::vector<Vec2>& ring, const Vec2& p) noexcept {
    bool inside = false;
    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
        const Vec2& a = ring[j];
        const Vec2& b = ring[i];
        if (pointOnSegment(a, b, p)) return true;
        const bool crosses = (a.v > p.v) != (b.v > p.v);
        if (!crosses) continue;
        const auto sign = predicates::orient2d(a, b, p);
        if ((b.v > a.v && sign == predicates::Sign::Positive) ||
            (b.v < a.v && sign == predicates::Sign::Negative))
            inside = !inside;
    }
    return inside;
}

bool segmentsConflict(const Vec2& a, const Vec2& b,
                      const Vec2& c, const Vec2& d,
                      bool sharedEndpoint) noexcept {
    const auto o1 = predicates::orient2d(a, b, c);
    const auto o2 = predicates::orient2d(a, b, d);
    const auto o3 = predicates::orient2d(c, d, a);
    const auto o4 = predicates::orient2d(c, d, b);
    if (o1 == predicates::Sign::Zero && pointOnSegment(a, b, c) &&
        !(sharedEndpoint && (sameParam(c, a) || sameParam(c, b)))) return true;
    if (o2 == predicates::Sign::Zero && pointOnSegment(a, b, d) &&
        !(sharedEndpoint && (sameParam(d, a) || sameParam(d, b)))) return true;
    if (o3 == predicates::Sign::Zero && pointOnSegment(c, d, a) &&
        !(sharedEndpoint && (sameParam(a, c) || sameParam(a, d)))) return true;
    if (o4 == predicates::Sign::Zero && pointOnSegment(c, d, b) &&
        !(sharedEndpoint && (sameParam(b, c) || sameParam(b, d)))) return true;
    return o1 != predicates::Sign::Zero && o2 != predicates::Sign::Zero &&
           o3 != predicates::Sign::Zero && o4 != predicates::Sign::Zero &&
           o1 != o2 && o3 != o4;
}

SurfaceFailure validateEmbedding(const FacePlan& plan) {
    struct SegmentGeom {
        bnd::ConstraintEntityId id;
        FaceLocalVertexId a, b;
        Vec2 pa, pb;
    };
    std::vector<SegmentGeom> geometry;
    geometry.reserve(plan.segments.size());
    for (const auto& s : plan.segments) {
        const auto pa = plan.vertices.at(s.a)->parameter.p;
        const auto pb = plan.vertices.at(s.b)->parameter.p;
        if (sameParam(pa, pb))
            return failure(SurfaceReason::InvalidConstraintEmbedding);
        geometry.push_back({s.id, s.a, s.b, pa, pb});
    }
    for (std::size_t i = 0; i < geometry.size(); ++i)
        for (std::size_t j = i + 1; j < geometry.size(); ++j) {
            const bool shared = geometry[i].a == geometry[j].a ||
                geometry[i].a == geometry[j].b ||
                geometry[i].b == geometry[j].a ||
                geometry[i].b == geometry[j].b;
            if (segmentsConflict(geometry[i].pa, geometry[i].pb,
                                 geometry[j].pa, geometry[j].pb, shared))
                return failure(SurfaceReason::InvalidConstraintEmbedding);
        }
    std::vector<std::vector<Vec2>> loops;
    for (const auto& ids : plan.loops) {
        std::vector<Vec2> ring;
        for (const auto id : ids) ring.push_back(plan.vertices.at(id)->parameter.p);
        if (ringWinding(ids, plan) == predicates::Sign::Zero)
            return failure(SurfaceReason::InvalidLoopTopology);
        loops.push_back(std::move(ring));
    }
    if (loops.empty()) return failure(SurfaceReason::InvalidLoopTopology);
    for (std::size_t h = 1; h < loops.size(); ++h) {
        if (!pointInClosure(loops[0], loops[h][0]))
            return failure(SurfaceReason::InvalidLoopTopology);
        for (std::size_t k = 1; k < loops.size(); ++k)
            if (h != k && pointInClosure(loops[k], loops[h][0]))
                return failure(SurfaceReason::InvalidLoopTopology);
    }
    for (const auto& rec : plan.chains)
        if (rec.chain.role == SurfaceChainRole::EmbeddedConstraint)
            for (const auto& use : rec.chain.uses) {
                const Vec2 p = plan.vertices.at(use.faceLocalVertexId)->parameter.p;
                if (!pointInClosure(loops[0], p))
                    return failure(SurfaceReason::InvalidConstraintEmbedding);
                for (std::size_t h = 1; h < loops.size(); ++h)
                    if (pointInClosure(loops[h], p))
                        return failure(SurfaceReason::InvalidConstraintEmbedding);
            }
    return {};
}

SurfaceFailure rebuildDerivedConstraints(FacePlan& plan) {
    plan.segments.clear();
    plan.origins.clear();
    plan.protectedIds.clear();
    for (auto& rec : plan.chains) {
        const auto& source = rec.source;
        const std::size_t sourceSegments = source.segmentEntityIds.size();
        const std::size_t normalizedSegments =
            rec.chain.closed ? rec.chain.uses.size()
                             : rec.chain.uses.size() - 1;
        if (sourceSegments != normalizedSegments)
            return failure(SurfaceReason::InvalidSegmentIdentity);
        rec.originalOrdinals.assign(normalizedSegments, 0);
        rec.chain.segmentEntityIds.assign(normalizedSegments, 0);
        for (std::size_t ni = 0; ni < normalizedSegments; ++ni) {
            const std::size_t nn = rec.chain.closed
                                       ? (ni + 1) % rec.chain.uses.size()
                                       : ni + 1;
            const auto& na = rec.chain.uses[ni];
            const auto& nb = rec.chain.uses[nn];
            std::size_t matched = sourceSegments;
            for (std::size_t oi = 0; oi < sourceSegments; ++oi) {
                const std::size_t on = source.closed
                                           ? (oi + 1) % source.uses.size()
                                           : oi + 1;
                const auto& oa = source.uses[oi];
                const auto& ob = source.uses[on];
                const bool same =
                    oa.faceLocalVertexId == na.faceLocalVertexId &&
                    ob.faceLocalVertexId == nb.faceLocalVertexId &&
                    oa.sharedKey == na.sharedKey && ob.sharedKey == nb.sharedKey;
                const bool reverse =
                    oa.faceLocalVertexId == nb.faceLocalVertexId &&
                    ob.faceLocalVertexId == na.faceLocalVertexId &&
                    oa.sharedKey == nb.sharedKey && ob.sharedKey == na.sharedKey;
                if (same || reverse) {
                    if (matched != sourceSegments)
                        return failure(SurfaceReason::InvalidSegmentIdentity);
                    matched = oi;
                }
            }
            if (matched == sourceSegments)
                return failure(SurfaceReason::InvalidSegmentIdentity);
            const auto entity = source.segmentEntityIds[matched];
            rec.originalOrdinals[ni] = matched;
            rec.chain.segmentEntityIds[ni] = entity;
            SurfaceSegmentOrigin origin;
            origin.bndEntityId = entity;
            origin.chainPersistentId = source.chainPersistentId;
            origin.originalSegmentOrdinal = matched;
            origin.normalizedSegmentOrdinal = ni;
            origin.endpoints = {na, nb};
            origin.classification = source.classification;
            if (!plan.origins.emplace(entity, origin).second)
                return failure(SurfaceReason::InvalidSegmentIdentity);
            plan.segments.push_back(
                {entity, na.faceLocalVertexId, nb.faceLocalVertexId,
                 plan.request->chart, source.classification});
            if (source.owningCurve ||
                source.role != SurfaceChainRole::EmbeddedConstraint)
                plan.protectedIds.insert(entity);
        }
    }
    std::sort(plan.segments.begin(), plan.segments.end(),
              [](const auto& a, const auto& b) { return a.id < b.id; });
    return {};
}

SurfaceFailure validateAndPlanFace(const SurfaceFaceRequest& f,
                                   FacePlan& plan) {
    plan.request = &f;
    if (f.faceClass.dim != 2 || f.faceClass.tag < 0)
        return failure(SurfaceReason::InvalidClassification);
    if (f.desiredNormalSign != 1 && f.desiredNormalSign != -1)
        return failure(SurfaceReason::InvalidPhysicalPolicy);
    if (!validPolicy(f.physicalPolicy))
        return failure(SurfaceReason::InvalidPhysicalPolicy);
    if (f.outer.role != SurfaceChainRole::OuterBoundary ||
        f.outer.pieces.empty())
        return failure(SurfaceReason::InvalidLoopTopology);
    for (const auto& h : f.holes)
        if (h.role != SurfaceChainRole::HoleBoundary || h.pieces.empty())
            return failure(SurfaceReason::InvalidLoopTopology);

    for (const auto& v : f.vertices) {
        if (v.parameter.chart != f.chart || !finite(v.parameter.p))
            return failure(SurfaceReason::NonfiniteParam);
        if (v.classification.dim < 0 || v.classification.dim > 2 ||
            v.classification.tag < 0)
            return failure(SurfaceReason::InvalidClassification);
        if (!plan.vertices.emplace(v.faceLocalVertexId, &v).second)
            return failure(SurfaceReason::InvalidConstraintEmbedding);
    }
    std::set<SurfaceNodeIdentity> declaredRefinedIdentities;
    for(const auto& refined:f.refinedNodes) {
        if(refined.parameter.chart!=f.chart||!finite(refined.parameter.p))
            return failure(SurfaceReason::NonfiniteParam);
        if(refined.classification.dim<0||refined.classification.dim>2||
           refined.classification.tag<0)
            return failure(SurfaceReason::InvalidClassification);
        if(!validSurfaceNodeIdentity(refined.identity)||
           refined.identity.domain==
               SurfaceNodeIdentityDomain::OriginalDiscretization||
           !declaredRefinedIdentities.insert(refined.identity).second)
            return failure(SurfaceReason::SharedIdentityConflict);
        if(refined.identity.domain==SurfaceNodeIdentityDomain::RefinedCurve) {
            if(refined.classification.dim!=1)
                return failure(SurfaceReason::InvalidClassification);
        } else {
            if(refined.identity.parentArtifact!=f.faceStamp||
               !identityMatchesFaceParameter(refined.identity,
                                             refined.parameter.p))
                return failure(SurfaceReason::SharedIdentityConflict);
            if(!sameClass(refined.classification,f.faceClass))
                return failure(SurfaceReason::InvalidClassification);
        }
    }
    std::set<FaceLocalVertexId> referenced;
    std::set<SurfaceChainId> chainIds;
    std::set<bnd::ConstraintEntityId> entityIds;
    std::set<SurfaceNodeIdentity> referencedRefinedCurveIdentities;

    auto addChain = [&](const SurfaceChain& input,
                        SurfaceChainRole expectedRole,
                        std::optional<std::size_t> loopOrdinal)
        -> SurfaceFailure {
        if (input.role != expectedRole)
            return failure(SurfaceReason::InvalidLoopTopology);
        if (input.classification.dim != 1 || input.classification.tag < 0)
            return failure(SurfaceReason::InvalidClassification);
        if (!chainIds.insert(input.chainPersistentId).second)
            return failure(SurfaceReason::InvalidSegmentIdentity);
        std::size_t n = input.uses.size();
        if (input.closed && n > 1 &&
            input.uses.front().faceLocalVertexId ==
                input.uses.back().faceLocalVertexId &&
            input.uses.front().sharedKey == input.uses.back().sharedKey)
            --n;
        const std::size_t expectedSegments = input.closed ? n : n - 1;
        if (n < (input.closed ? 3u : 2u) ||
            input.segmentEntityIds.size() != expectedSegments)
            return failure(SurfaceReason::InvalidSegmentIdentity);
        if (expectedRole != SurfaceChainRole::EmbeddedConstraint &&
            !input.owningCurve)
            return failure(SurfaceReason::InvalidConstraintEmbedding);
        for (const auto& use : input.uses) {
            if (plan.vertices.count(use.faceLocalVertexId) == 0)
                return failure(SurfaceReason::InvalidConstraintEmbedding);
            referenced.insert(use.faceLocalVertexId);
            if (input.owningCurve) {
                if (!use.sharedKey ||
                    use.sharedKey->curvePersistentId !=
                        input.owningCurve->persistentId ||
                    use.sharedKey->curveRevision != input.owningCurve->revision)
                    return failure(SurfaceReason::InvalidConstraintEmbedding);
                const auto& identity=use.sharedKey->identity;
                if(!validSurfaceNodeIdentity(identity))
                    return failure(SurfaceReason::SharedIdentityConflict);
                if(identity.domain!=
                       SurfaceNodeIdentityDomain::OriginalDiscretization) {
                    if(identity.domain!=SurfaceNodeIdentityDomain::RefinedCurve||
                       use.sharedKey->nodeOrdinal!=0||
                       identity.parentArtifact!=*input.owningCurve)
                        return failure(SurfaceReason::SharedIdentityConflict);
                    referencedRefinedCurveIdentities.insert(identity);
                }
            } else if (use.sharedKey) {
                return failure(SurfaceReason::InvalidConstraintEmbedding);
            }
        }
        for (const auto id : input.segmentEntityIds)
            if (!entityIds.insert(id).second)
                return failure(SurfaceReason::InvalidSegmentIdentity);

        bool reversed = false;
        std::size_t rotation = 0;
        ChainRecord rec;
        rec.source = input;
        if (rec.source.closed && rec.source.uses.size() > 1 &&
            rec.source.uses.front().faceLocalVertexId ==
                rec.source.uses.back().faceLocalVertexId &&
            rec.source.uses.front().sharedKey ==
                rec.source.uses.back().sharedKey)
            rec.source.uses.pop_back();
        rec.chain = rec.source;
        rec.loopOrdinal = loopOrdinal;
        if (expectedRole == SurfaceChainRole::EmbeddedConstraint)
            rec.chain.uses = canonicalVertices(rec.source, reversed, rotation);
        plan.chains.push_back(std::move(rec));
        return {};
    };

    auto addLoop = [&](const SurfaceLoop& loop,
                       SurfaceChainRole role) -> SurfaceFailure {
        std::vector<FaceLocalVertexId> ring;
        const std::size_t loopOrdinal = plan.loops.size();
        for (const auto& chain : loop.pieces) {
            SurfaceFailure e = addChain(chain, role, loopOrdinal);
            if (e.reason != SurfaceReason::None) return e;
            const auto& uses = plan.chains.back().chain.uses;
            if (chain.closed && loop.pieces.size() != 1)
                return failure(SurfaceReason::InvalidLoopTopology);
            if (!ring.empty() && ring.back() != uses.front().faceLocalVertexId)
                return failure(SurfaceReason::InvalidLoopTopology);
            const std::size_t begin = ring.empty() ? 0 : 1;
            for (std::size_t i = begin; i < uses.size(); ++i)
                ring.push_back(uses[i].faceLocalVertexId);
        }
        const bool closed = loop.pieces.front().closed ||
            (ring.size() > 1 && ring.front() == ring.back());
        if (ring.size() > 1 && ring.front() == ring.back()) ring.pop_back();
        if (ring.size() < 3 || !closed)
            return failure(SurfaceReason::InvalidLoopTopology);
        std::set<FaceLocalVertexId> unique(ring.begin(), ring.end());
        if (unique.size() != ring.size())
            return failure(SurfaceReason::InvalidLoopTopology);
        plan.loops.push_back(std::move(ring));
        return {};
    };

    SurfaceFailure e = addLoop(f.outer, SurfaceChainRole::OuterBoundary);
    if (e.reason != SurfaceReason::None) return e;
    for (const auto& hole : f.holes) {
        e = addLoop(hole, SurfaceChainRole::HoleBoundary);
        if (e.reason != SurfaceReason::None) return e;
    }
    for (const auto& embedded : f.embedded) {
        e = addChain(embedded, SurfaceChainRole::EmbeddedConstraint,
                     std::nullopt);
        if (e.reason != SurfaceReason::None) return e;
    }
    if (referenced.size() != plan.vertices.size())
        return failure(SurfaceReason::InvalidConstraintEmbedding);
    for(const auto& refined:f.refinedNodes)
        if(refined.identity.domain==SurfaceNodeIdentityDomain::RefinedCurve&&
           referencedRefinedCurveIdentities.count(refined.identity)==0)
            return failure(SurfaceReason::SharedIdentityConflict);
    e = rebuildDerivedConstraints(plan);
    if (e.reason != SurfaceReason::None) return e;
    return validateEmbedding(plan);
}

SurfaceFailure buildSetPlan(const SurfaceSetRequest& request,
                            const mesh::Mesh& destination,
                            SetPlan& plan) {
    if (request.setStamp.persistentId == 0)
        return failure(SurfaceReason::StaleFaceArtifact);
    std::set<FacePersistentId> faceIds;
    std::vector<const SurfaceFaceRequest*> ordered;
    ordered.reserve(request.faces.size());
    for (const auto& f : request.faces) {
        if (!faceIds.insert(f.faceStamp.persistentId).second)
            return failure(SurfaceReason::DuplicateFaceIdentity);
        ordered.push_back(&f);
    }
    std::sort(ordered.begin(), ordered.end(), [](const auto* a, const auto* b) {
        return a->faceStamp.persistentId < b->faceStamp.persistentId;
    });

    std::map<ArtifactStamp, mesh::Classification> curveClasses;
    std::set<SharedCurveNodeKey> referencedKeys;
    std::map<SharedCurveNodeKey,
             std::vector<std::pair<FacePersistentId, FaceLocalVertexId>>> uses;
    for (const auto* f : ordered) {
        FacePlan fp;
        SurfaceFailure e = validateAndPlanFace(*f, fp);
        if (e.reason != SurfaceReason::None) return e;
        for (const auto& rec : fp.chains) {
            if (rec.chain.owningCurve) {
                const auto [it, inserted] = curveClasses.emplace(
                    *rec.chain.owningCurve, rec.chain.classification);
                if (!inserted && !sameClass(it->second,
                                            rec.chain.classification))
                    return failure(SurfaceReason::SharedClassificationConflict);
            }
            for (const auto& use : rec.chain.uses)
                if (use.sharedKey) {
                    referencedKeys.insert(*use.sharedKey);
                    uses[*use.sharedKey].push_back(
                        {f->faceStamp.persistentId, use.faceLocalVertexId});
                }
        }
        plan.faces.push_back(std::move(fp));
    }

    for (const auto& binding : request.sharedBindings) {
        if(!validSurfaceNodeIdentity(binding.key.identity))
            return failure(SurfaceReason::SharedIdentityConflict);
        if(binding.key.identity.domain!=
               SurfaceNodeIdentityDomain::OriginalDiscretization&&
           (binding.key.identity.domain!=
                SurfaceNodeIdentityDomain::RefinedCurve||
            binding.key.nodeOrdinal!=0||
            binding.key.identity.parentArtifact.persistentId!=
                binding.key.curvePersistentId||
            binding.key.identity.parentArtifact.revision!=
                binding.key.curveRevision))
            return failure(SurfaceReason::SharedIdentityConflict);
        if (!plan.ledger.emplace(binding.key, binding).second)
            return failure(SurfaceReason::DuplicateSharedBinding);
    }
    for (const auto& key : referencedKeys)
        if (plan.ledger.count(key) == 0)
            return failure(SurfaceReason::MissingSharedBinding);
    std::map<SurfaceNodeIdentity,
             std::pair<Index,mesh::Classification>> refinedBindings;
    std::map<Index,SurfaceNodeIdentity> refinedIdentityByGlobal;
    for (const auto& [key, binding] : plan.ledger) {
        if (referencedKeys.count(key) == 0)
            return failure(SurfaceReason::UnusedSharedBinding);
        if (binding.existingGlobalIndex < 0 ||
            binding.existingGlobalIndex >= destination.nodeCount())
            return failure(SurfaceReason::UnknownSharedGlobalIndex);
        if(key.identity.domain==SurfaceNodeIdentityDomain::RefinedCurve) {
            if(binding.classification.dim!=1||binding.classification.tag<0)
                return failure(SurfaceReason::InvalidClassification);
            const auto [identityIt,identityInserted]=refinedBindings.emplace(
                key.identity,std::make_pair(binding.existingGlobalIndex,
                                            binding.classification));
            if(!identityInserted&&
               (identityIt->second.first!=binding.existingGlobalIndex||
                !sameClass(identityIt->second.second,binding.classification)))
                return failure(SurfaceReason::SharedIdentityConflict);
            const auto [globalIt,globalInserted]=
                refinedIdentityByGlobal.emplace(binding.existingGlobalIndex,
                                                key.identity);
            if(!globalInserted&&globalIt->second!=key.identity)
                return failure(SurfaceReason::SharedIdentityConflict);
        }
    }
    for (const auto& [key, keyUses] : uses) {
        const auto& binding = plan.ledger.at(key);
        if (!sameClass(binding.classification,
                       destination.nodeClass(binding.existingGlobalIndex)))
            return failure(SurfaceReason::SharedClassificationConflict);
        for (const auto& [faceId, localId] : keyUses) {
            const auto fit = std::find_if(
                plan.faces.begin(), plan.faces.end(), [&](const FacePlan& fp) {
                    return fp.request->faceStamp.persistentId == faceId;
                });
            if (fit == plan.faces.end() ||
                !sameClass(fit->vertices.at(localId)->classification,
                           binding.classification))
                return failure(SurfaceReason::SharedClassificationConflict);
        }
    }

    std::map<Index, std::vector<SharedCurveNodeKey>> reverse;
    for (const auto& [key, binding] : plan.ledger)
        reverse[binding.existingGlobalIndex].push_back(key);
    for (const auto& [global, keys] : reverse) {
        (void)global;
        if (keys.size() < 2) continue;
        const auto& cls = plan.ledger.at(keys.front()).classification;
        if (cls.dim != 0)
            return failure(SurfaceReason::SharedIdentityConflict);
        for (const auto& key : keys)
            if (!sameClass(plan.ledger.at(key).classification, cls))
                return failure(SurfaceReason::SharedClassificationConflict);
        std::map<SharedCurveNodeKey, std::set<SharedCurveNodeKey>> witnessGraph;
        for (const auto& fp : plan.faces) {
            std::map<FaceLocalVertexId, std::set<SharedCurveNodeKey>> byVertex;
            for (const auto& rec : fp.chains)
                for (const auto& use : rec.chain.uses)
                    if (use.sharedKey &&
                        plan.ledger.at(*use.sharedKey).existingGlobalIndex ==
                            global)
                        byVertex[use.faceLocalVertexId].insert(*use.sharedKey);
            for (const auto& [local, localKeys] : byVertex) {
                (void)local;
                for (const auto& a : localKeys)
                    for (const auto& b : localKeys)
                        witnessGraph[a].insert(b);
            }
        }
        std::set<SharedCurveNodeKey> witnessed;
        std::vector<SharedCurveNodeKey> keyStack{keys.front()};
        while (!keyStack.empty()) {
            const auto key = keyStack.back(); keyStack.pop_back();
            if (!witnessed.insert(key).second) continue;
            for (const auto& nextKey : witnessGraph[key])
                keyStack.push_back(nextKey);
        }
        if (witnessed.size() != keys.size())
            return failure(SurfaceReason::SharedIdentityConflict);
    }
    return {};
}

bool cyclicEdgesEqual(const std::vector<std::size_t>& ring,
                      const std::set<EdgeIndex>& actual) {
    if (ring.size() < 3) return false;
    std::set<EdgeIndex> expected;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        if (ring[i] == ring[(i + 1) % ring.size()]) return false;
        expected.insert(edgeIndex(ring[i], ring[(i + 1) % ring.size()]));
    }
    return expected == actual;
}

SurfaceFailure validateCandidate(const SurfaceFaceRequest& request,
                                 const FacePlan& plan,
                                 const SurfaceCandidate& candidate,
                                 const ChartCertificate& certificate) {
    if (!(candidate.providerStamp ==
          request.expectedCandidateProviderStamp))
        return candidateFailure(SurfaceCandidateDetail::ProviderStampChanged);
    if (!(candidate.faceStamp == request.faceStamp) ||
        candidate.chart != request.chart || candidate.nodes.empty() ||
        candidate.cells.empty())
        return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    if (candidate.declaredVertexMap.size() != plan.vertices.size())
        return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    struct RefinedExpectation {
        bnd::ParamPoint parameter{};
        mesh::Classification classification{};
    };
    std::map<SurfaceNodeIdentity,RefinedExpectation> expectedRefined;
    std::map<FaceLocalVertexId,SurfaceNodeIdentity> expectedRefinedVertex;
    for (const auto& expected : request.refinedNodes) {
        if (expected.identity.domain ==
                SurfaceNodeIdentityDomain::OriginalDiscretization ||
            expected.parameter.chart != request.chart ||
            !expectedRefined.emplace(
                expected.identity,
                RefinedExpectation{expected.parameter,
                                   expected.classification}).second)
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    }
    for (const auto& rec : plan.chains)
        for (const auto& use : rec.chain.uses)
            if (use.sharedKey && use.sharedKey->identity.domain !=
                                     SurfaceNodeIdentityDomain::
                                         OriginalDiscretization) {
                const auto vertex=plan.vertices.find(use.faceLocalVertexId);
                if(vertex==plan.vertices.end())
                    return candidateFailure(
                        SurfaceCandidateDetail::MalformedCandidate);
                const RefinedExpectation expectation{
                    vertex->second->parameter,
                    vertex->second->classification};
                const auto [it,inserted]=expectedRefined.emplace(
                    use.sharedKey->identity,expectation);
                if(!inserted&&
                   (!sameParam(it->second.parameter.p,
                               expectation.parameter.p)||
                    it->second.parameter.chart!=expectation.parameter.chart||
                    !sameClass(it->second.classification,
                               expectation.classification)))
                    return candidateFailure(
                        SurfaceCandidateDetail::MalformedCandidate);
                const auto [vertexIt,vertexInserted]=
                    expectedRefinedVertex.emplace(
                        use.faceLocalVertexId,use.sharedKey->identity);
                if(!vertexInserted&&vertexIt->second!=use.sharedKey->identity)
                    return candidateFailure(
                        SurfaceCandidateDetail::MalformedCandidate);
            }
    std::set<SurfaceNodeIdentity> observedRefined;
    std::set<std::size_t> mapped;
    for (const auto& [id, vertex] : plan.vertices) {
        const auto it = candidate.declaredVertexMap.find(id);
        if (it == candidate.declaredVertexMap.end() ||
            it->second >= candidate.nodes.size() ||
            !mapped.insert(it->second).second)
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        const auto& node = candidate.nodes[it->second];
        if (node.frameOnly || node.parameter.chart != request.chart ||
            !sameParam(node.parameter.p, vertex->parameter.p) ||
            !sameClass(node.node.cls, vertex->classification))
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        const auto refinedVertex=expectedRefinedVertex.find(id);
        if((refinedVertex==expectedRefinedVertex.end()&&
            node.persistentIdentity.has_value())||
           (refinedVertex!=expectedRefinedVertex.end()&&
            (!node.persistentIdentity||
             *node.persistentIdentity!=refinedVertex->second)))
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    }
    for (std::size_t i = 0; i < candidate.nodes.size(); ++i) {
        const auto& n = candidate.nodes[i];
        if (!finite(n.node.coord) || !finite(n.parameter.p) ||
            n.parameter.chart != request.chart)
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        if(n.persistentIdentity) {
            if(n.frameOnly||!validSurfaceNodeIdentity(*n.persistentIdentity))
                return candidateFailure(
                    SurfaceCandidateDetail::MalformedCandidate);
            const auto expected=expectedRefined.find(*n.persistentIdentity);
            if(expected==expectedRefined.end()||
               !observedRefined.insert(*n.persistentIdentity).second||
               !sameParam(expected->second.parameter.p,n.parameter.p)||
               expected->second.parameter.chart!=n.parameter.chart||
               !sameClass(expected->second.classification,n.node.cls))
                return candidateFailure(
                    SurfaceCandidateDetail::MalformedCandidate);
        }
        if (mapped.count(i) == 0) {
            if (!sameClass(n.node.cls, request.faceClass))
                return candidateFailure(
                    SurfaceCandidateDetail::MalformedCandidate);
            for (const auto& segment : plan.segments) {
                const auto& a = plan.vertices.at(segment.a)->parameter.p;
                const auto& b = plan.vertices.at(segment.b)->parameter.p;
                if (pointOnSegment(a, b, n.parameter.p))
                    return candidateFailure(
                        SurfaceCandidateDetail::MalformedCandidate);
            }
        }
    }
    if(observedRefined.size()!=expectedRefined.size())
        return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    const auto& region = certificate.certifiedRegionBoundary;
    for (std::size_t i = 0; i < region.size(); ++i)
        for (std::size_t j = i + 1; j < region.size(); ++j) {
            const std::size_t in = (i + 1) % region.size();
            const std::size_t jn = (j + 1) % region.size();
            if (i == j || in == j || jn == i) continue;
            if (segmentsConflict(region[i], region[in], region[j], region[jn],
                                 false))
                return candidateFailure(
                    SurfaceCandidateDetail::MalformedCandidate);
        }
    for (const auto& node : candidate.nodes)
        if (!pointInClosure(region, node.parameter.p))
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);

    std::vector<Vec2> frame;
    frame.reserve(candidate.frameBoundary.size());
    for (const auto i : candidate.frameBoundary) {
        if (i >= candidate.nodes.size())
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        frame.push_back(candidate.nodes[i].parameter.p);
    }
    for (const auto& [id, vertex] : plan.vertices) {
        (void)id;
        if (!pointInClosure(frame, vertex->parameter.p))
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        for (std::size_t i = 0; i < frame.size(); ++i)
            if (pointOnSegment(frame[i], frame[(i + 1) % frame.size()],
                               vertex->parameter.p))
                return candidateFailure(
                    SurfaceCandidateDetail::MalformedCandidate);
    }

    std::map<EdgeIndex, std::vector<std::size_t>> owners;
    std::set<std::array<std::size_t, 3>> triangles;
    predicates::Sign winding = predicates::Sign::Zero;
    std::vector<std::vector<std::size_t>> neighbors(candidate.cells.size());
    for (std::size_t ci = 0; ci < candidate.cells.size(); ++ci) {
        const auto& c = candidate.cells[ci];
        if (!sameClass(c.classification, request.faceClass) ||
            c.nodes[0] >= candidate.nodes.size() ||
            c.nodes[1] >= candidate.nodes.size() ||
            c.nodes[2] >= candidate.nodes.size() ||
            c.nodes[0] == c.nodes[1] || c.nodes[1] == c.nodes[2] ||
            c.nodes[2] == c.nodes[0])
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        const auto sign = predicates::orient2d(
            candidate.nodes[c.nodes[0]].parameter.p,
            candidate.nodes[c.nodes[1]].parameter.p,
            candidate.nodes[c.nodes[2]].parameter.p);
        if (sign == predicates::Sign::Zero ||
            (winding != predicates::Sign::Zero && sign != winding))
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        winding = sign;
        auto key = c.nodes;
        std::sort(key.begin(), key.end());
        if (!triangles.insert(key).second)
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
        for (std::size_t i = 0; i < 3; ++i) {
            auto& edgeOwners = owners[edgeIndex(c.nodes[i], c.nodes[(i + 1) % 3])];
            edgeOwners.push_back(ci);
            if (edgeOwners.size() > 2)
                return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
        }
    }
    std::set<EdgeIndex> boundary;
    for (const auto& [edge, edgeOwners] : owners) {
        if (edgeOwners.size() == 1) boundary.insert(edge);
        else if (edgeOwners.size() == 2) {
            neighbors[edgeOwners[0]].push_back(edgeOwners[1]);
            neighbors[edgeOwners[1]].push_back(edgeOwners[0]);
        } else {
            return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
        }
    }
    if (!cyclicEdgesEqual(candidate.frameBoundary, boundary))
        return candidateFailure(SurfaceCandidateDetail::PretrimmedDomain);
    for (const auto i : candidate.frameBoundary)
        if (i >= candidate.nodes.size() || !candidate.nodes[i].frameOnly)
            return candidateFailure(SurfaceCandidateDetail::MalformedCandidate);
    const std::int64_t chi = static_cast<std::int64_t>(candidate.nodes.size()) -
        static_cast<std::int64_t>(owners.size()) +
        static_cast<std::int64_t>(candidate.cells.size());
    if (chi != 1)
        return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
    std::set<std::size_t> reached;
    std::vector<std::size_t> stack{0};
    while (!stack.empty()) {
        const auto c = stack.back(); stack.pop_back();
        if (!reached.insert(c).second) continue;
        for (const auto n : neighbors[c]) stack.push_back(n);
    }
    if (reached.size() != candidate.cells.size())
        return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
    std::vector<std::vector<std::array<std::size_t, 2>>> links(
        candidate.nodes.size());
    for (const auto& c : candidate.cells)
        for (std::size_t i = 0; i < 3; ++i)
            links[c.nodes[i]].push_back(
                {c.nodes[(i + 1) % 3], c.nodes[(i + 2) % 3]});
    for (std::size_t v = 0; v < links.size(); ++v) {
        std::map<std::size_t, int> degree;
        std::map<std::size_t, std::vector<std::size_t>> graph;
        for (const auto& pair : links[v]) {
            ++degree[pair[0]]; ++degree[pair[1]];
            graph[pair[0]].push_back(pair[1]);
            graph[pair[1]].push_back(pair[0]);
        }
        if (degree.empty())
            return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
        int endpoints = 0;
        bool badDegree = false;
        for (const auto& [node, d] : degree) {
            (void)node;
            if (d == 1) ++endpoints;
            else if (d != 2) badDegree = true;
        }
        const bool boundaryNode =
            std::find(candidate.frameBoundary.begin(),
                      candidate.frameBoundary.end(), v) !=
            candidate.frameBoundary.end();
        if (badDegree || (boundaryNode ? endpoints != 2 : endpoints != 0))
            return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
        std::set<std::size_t> linkReached;
        std::vector<std::size_t> linkStack{degree.begin()->first};
        while (!linkStack.empty()) {
            const auto n = linkStack.back(); linkStack.pop_back();
            if (!linkReached.insert(n).second) continue;
            for (const auto x : graph[n]) linkStack.push_back(x);
        }
        if (linkReached.size() != degree.size())
            return candidateFailure(SurfaceCandidateDetail::DisconnectedFrame);
    }
    if (candidate.exteriorSeedCell >= candidate.cells.size())
        return candidateFailure(SurfaceCandidateDetail::MissingExteriorSeed);
    bool seedOnFrame = false;
    const auto& seed = candidate.cells[candidate.exteriorSeedCell];
    for (std::size_t i = 0; i < 3; ++i)
        if (boundary.count(edgeIndex(seed.nodes[i], seed.nodes[(i + 1) % 3])))
            seedOnFrame = true;
    if (!seedOnFrame)
        return candidateFailure(SurfaceCandidateDetail::MissingExteriorSeed);
    return {};
}

bool sameCertificate(const ChartCertificate& a,
                     const ChartCertificate& b) noexcept {
    if (!(a.providerStamp == b.providerStamp) || a.chart != b.chart ||
        !(a.certificateStamp == b.certificateStamp) ||
        a.chartPositiveToNormalSign != b.chartPositiveToNormalSign ||
        a.regular != b.regular || a.periodic != b.periodic ||
        a.seam != b.seam ||
        a.certifiedRegionBoundary.size() != b.certifiedRegionBoundary.size())
        return false;
    for (std::size_t i = 0; i < a.certifiedRegionBoundary.size(); ++i)
        if (!sameParam(a.certifiedRegionBoundary[i],
                       b.certifiedRegionBoundary[i]))
            return false;
    return true;
}

SurfaceFailure observeGeometryCertificate(
    const SurfaceFaceRequest& request,
    const SurfaceGeometryProvider& geometry,
    ChartCertificate& observed) noexcept {
    try {
        const auto stamp = geometry.stamp();
        if (!(stamp == request.expectedGeometryStamp))
            return failure(SurfaceReason::StaleGeometryProvider);
        auto cert = geometry.chartCertificate(request.chart);
        if (!(cert.providerStamp == stamp) || cert.chart != request.chart ||
            !(cert.certificateStamp ==
              request.expectedChartCertificateStamp))
            return failure(SurfaceReason::StaleGeometryProvider);
        if (!cert.regular)
            return failure(SurfaceReason::UnsupportedSingularChart);
        if (cert.periodic)
            return failure(SurfaceReason::UnsupportedPeriodicChart);
        if (cert.seam)
            return failure(SurfaceReason::UnsupportedSeamChart);
        if ((cert.chartPositiveToNormalSign != 1 &&
             cert.chartPositiveToNormalSign != -1) ||
            cert.certifiedRegionBoundary.size() < 3)
            return failure(SurfaceReason::UnsupportedSingularChart);
        for (const auto& p : cert.certifiedRegionBoundary)
            if (!finite(p)) return failure(SurfaceReason::NonfiniteParam);
        observed = std::move(cert);
        return {};
    } catch (const bnd::GeometryProviderDomainError&) {
        return failure(SurfaceReason::ProviderDomainRejected);
    } catch (const std::bad_alloc&) {
        return allocation(SurfaceAllocationPhase::GeometryProvider);
    } catch (...) {
        return failure(SurfaceReason::GeometryProviderError);
    }
}

SurfaceFailure geometryCheckpoint(const SurfaceFaceRequest& request,
                                  const SurfaceGeometryProvider& geometry,
                                  std::optional<ChartCertificate>& retained,
                                  bool first) noexcept {
    ChartCertificate observed;
    SurfaceFailure result =
        observeGeometryCertificate(request, geometry, observed);
    if (result.reason != SurfaceReason::None) return result;
    if (first) retained.emplace(std::move(observed));
    else if (!retained || !sameCertificate(*retained, observed))
        return failure(SurfaceReason::StaleGeometryProvider);
    return {};
}

SurfaceFailure verifyGeometryCheckpoint(
    const SurfaceFaceRequest& request,
    const SurfaceGeometryProvider& geometry,
    const ChartCertificate& retained) noexcept {
    ChartCertificate observed;
    SurfaceFailure result =
        observeGeometryCertificate(request, geometry, observed);
    if (result.reason != SurfaceReason::None) return result;
    return sameCertificate(retained, observed)
               ? SurfaceFailure{}
               : failure(SurfaceReason::StaleGeometryProvider);
}

const SurfaceChain* chainById(const FacePlan& plan, SurfaceChainId id) {
    for (const auto& rec : plan.chains)
        if (rec.chain.chainPersistentId == id) return &rec.chain;
    return nullptr;
}

SurfaceFailure translateBnd(const bnd::SteinerResult& r) noexcept {
    SurfaceFailure f;
    f.bndReject = r.reject;
    f.bndRefuse = r.refuse;
    f.bndFail = r.fail;
    if (r.status == bnd::SteinerStatus::Rejected) {
        f.reason = SurfaceReason::RecoveryStopped;
    } else if (r.status == bnd::SteinerStatus::Refused) {
        f.reason = r.refuse == bnd::SteinerRefuse::MinimumFeatureReached
                       ? SurfaceReason::MinimumFeatureReached
                       : r.refuse == bnd::SteinerRefuse::SplitBudgetExceeded
                             ? SurfaceReason::SplitBudgetExceeded
                             : SurfaceReason::RecoveryStopped;
    } else {
        f.reason = SurfaceReason::RecoveryStopped;
    }
    return f;
}

SurfaceFaceStatus statusForFailure(const SurfaceFailure& f) noexcept {
    switch (f.reason) {
        case SurfaceReason::UnsupportedPeriodicChart:
        case SurfaceReason::UnsupportedSeamChart:
        case SurfaceReason::UnsupportedSingularChart:
        case SurfaceReason::BoundarySteinerProhibited:
        case SurfaceReason::CurveArtifactSubdivisionRequired:
        case SurfaceReason::MinimumFeatureReached:
        case SurfaceReason::SplitBudgetExceeded:
            return SurfaceFaceStatus::Refused;
        case SurfaceReason::InvalidLoopTopology:
        case SurfaceReason::InvalidConstraintEmbedding:
        case SurfaceReason::InvalidSegmentIdentity:
        case SurfaceReason::InvalidPhysicalPolicy:
        case SurfaceReason::StaleGeometryProvider:
        case SurfaceReason::StaleCandidateProvider:
        case SurfaceReason::StaleFaceArtifact:
        case SurfaceReason::StaleCurveArtifact:
        case SurfaceReason::DuplicateFaceIdentity:
        case SurfaceReason::DuplicateSharedBinding:
        case SurfaceReason::MissingSharedBinding:
        case SurfaceReason::UnusedSharedBinding:
        case SurfaceReason::UnknownSharedGlobalIndex:
        case SurfaceReason::SharedClassificationConflict:
        case SurfaceReason::SharedIdentityConflict:
        case SurfaceReason::SharedMappedPositionMismatch:
        case SurfaceReason::InvalidClassification:
        case SurfaceReason::NonfiniteParam:
        case SurfaceReason::NonfiniteMappedPosition:
        case SurfaceReason::NonfiniteMappedNormal:
        case SurfaceReason::ProviderDomainRejected:
            return SurfaceFaceStatus::Rejected;
        default: return SurfaceFaceStatus::Failed;
    }
}

SurfaceSetStatus setStatusFor(SurfaceFaceStatus status) noexcept {
    if (status == SurfaceFaceStatus::Rejected) return SurfaceSetStatus::Rejected;
    if (status == SurfaceFaceStatus::Refused) return SurfaceSetStatus::Refused;
    return SurfaceSetStatus::Failed;
}

struct FaceBuild {
    SurfaceFailure error{};
    std::unique_ptr<SurfaceFaceResult> result;
};

std::map<HandleEdge, std::vector<topo::CellHandle>> edgeOwners(
    const topo::TopologyWorkspace& ws) {
    std::map<HandleEdge, std::vector<topo::CellHandle>> owners;
    for (const auto& ch : ws.cells().activeBySequence()) {
        const auto cell = ws.cells().get(ch);
        for (std::size_t i = 0; i < 3; ++i)
            owners[handleEdge(cell.nodes[i], cell.nodes[(i + 1) % 3])]
                .push_back(ch);
    }
    return owners;
}

topo::NodeHandle thirdNode(const topo::CellPool::CellSnapshot& c,
                           const HandleEdge& edge) {
    for (std::size_t i = 0; i < 3; ++i)
        if (!(c.nodes[i] == edge[0]) && !(c.nodes[i] == edge[1]))
            return c.nodes[i];
    return {};
}

predicates::Sign ringWinding(
    const std::vector<FaceLocalVertexId>& ring,
    const FacePlan& plan) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < ring.size(); ++i) {
        const auto& p = plan.vertices.at(ring[i])->parameter.p;
        const auto& q = plan.vertices.at(ring[best])->parameter.p;
        if (p.u < q.u || (p.u == q.u && p.v < q.v)) best = i;
    }
    const auto& a = plan.vertices.at(ring[(best + ring.size() - 1) %
                                          ring.size()])->parameter.p;
    const auto& b = plan.vertices.at(ring[best])->parameter.p;
    const auto& c = plan.vertices.at(ring[(best + 1) % ring.size()])->parameter.p;
    return predicates::orient2d(a, b, c);
}

std::optional<SharedCurveNodeKey> canonicalVertexKey(
    const FacePlan& plan, FaceLocalVertexId id) {
    std::optional<SharedCurveNodeKey> key;
    for (const auto& rec : plan.chains) {
        if (!rec.loopOrdinal) continue;
        for (const auto& use : rec.source.uses)
            if (use.faceLocalVertexId == id && use.sharedKey &&
                (!key || *use.sharedKey < *key))
                key = use.sharedKey;
    }
    return key;
}

bool vertexTokenLess(const FacePlan& plan, FaceLocalVertexId a,
                     FaceLocalVertexId b) {
    const auto ka = canonicalVertexKey(plan, a);
    const auto kb = canonicalVertexKey(plan, b);
    if (ka && kb) {
        if (*ka < *kb) return true;
        if (*kb < *ka) return false;
    } else if (ka.has_value() != kb.has_value()) {
        return ka.has_value();
    }
    return a < b;
}

bool ringSequenceLess(const FacePlan& plan,
                      const std::vector<FaceLocalVertexId>& a,
                      const std::vector<FaceLocalVertexId>& b) {
    return std::lexicographical_compare(
        a.begin(), a.end(), b.begin(), b.end(),
        [&](FaceLocalVertexId x, FaceLocalVertexId y) {
            return vertexTokenLess(plan, x, y);
        });
}

void canonicalRotateRing(FacePlan& plan,
                         std::vector<FaceLocalVertexId>& ring) {
    std::size_t best = 0;
    for (std::size_t start = 1; start < ring.size(); ++start) {
        bool decided = false;
        for (std::size_t offset = 0; offset < ring.size(); ++offset) {
            const auto candidate = ring[(start + offset) % ring.size()];
            const auto current = ring[(best + offset) % ring.size()];
            if (candidate == current) continue;
            if (vertexTokenLess(plan, candidate, current)) best = start;
            decided = true;
            break;
        }
        (void)decided;
    }
    std::rotate(ring.begin(),
                ring.begin() + static_cast<std::ptrdiff_t>(best), ring.end());
}

int directedEdgeInRing(const std::vector<FaceLocalVertexId>& ring,
                       FaceLocalVertexId a, FaceLocalVertexId b,
                       std::size_t* position = nullptr) noexcept {
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto j = (i + 1) % ring.size();
        if (ring[i] == a && ring[j] == b) {
            if (position) *position = i;
            return 1;
        }
        if (ring[i] == b && ring[j] == a) {
            if (position) *position = i;
            return -1;
        }
    }
    return 0;
}

SurfaceFailure normalizeFacePlan(FacePlan& plan,
                                 const ChartCertificate& certificate) {
    if (plan.loops.empty())
        return failure(SurfaceReason::InvalidLoopTopology);
    const int outerWinding = certificate.chartPositiveToNormalSign *
                             plan.request->desiredNormalSign;
    for (std::size_t i = 0; i < plan.loops.size(); ++i) {
        auto& ring = plan.loops[i];
        const auto winding = ringWinding(ring, plan);
        if (winding == predicates::Sign::Zero)
            return failure(SurfaceReason::InvalidLoopTopology);
        const int actual = winding == predicates::Sign::Positive ? 1 : -1;
        const int required = i == 0 ? outerWinding : -outerWinding;
        if (actual != required) std::reverse(ring.begin(), ring.end());
        canonicalRotateRing(plan, ring);
    }

    struct HoleRecord {
        std::size_t oldOrdinal = 0;
        std::vector<FaceLocalVertexId> ring;
    };
    std::vector<HoleRecord> holes;
    holes.reserve(plan.loops.size() - 1);
    for (std::size_t i = 1; i < plan.loops.size(); ++i)
        holes.push_back({i, std::move(plan.loops[i])});
    std::sort(holes.begin(), holes.end(), [&](const auto& a, const auto& b) {
        return ringSequenceLess(plan, a.ring, b.ring);
    });
    std::map<std::size_t, std::size_t> newLoopOrdinal;
    std::vector<std::vector<FaceLocalVertexId>> normalizedLoops;
    normalizedLoops.reserve(plan.loops.size());
    normalizedLoops.push_back(std::move(plan.loops[0]));
    newLoopOrdinal.emplace(0, 0);
    for (std::size_t i = 0; i < holes.size(); ++i) {
        newLoopOrdinal.emplace(holes[i].oldOrdinal, i + 1);
        normalizedLoops.push_back(std::move(holes[i].ring));
    }
    plan.loops = std::move(normalizedLoops);

    for (auto& rec : plan.chains) {
        if (!rec.loopOrdinal) {
            bool reversed = false;
            std::size_t rotation = 0;
            rec.chain = rec.source;
            rec.chain.uses = canonicalVertices(rec.source, reversed, rotation);
            continue;
        }
        rec.loopOrdinal = newLoopOrdinal.at(*rec.loopOrdinal);
        const auto& ring = plan.loops[*rec.loopOrdinal];
        rec.chain = rec.source;
        if (rec.source.closed) {
            if (rec.source.uses.size() != ring.size())
                return failure(SurfaceReason::InvalidLoopTopology);
            rec.chain.uses.clear();
            rec.chain.uses.reserve(ring.size());
            for (const auto id : ring) {
                const auto it = std::find_if(
                    rec.source.uses.begin(), rec.source.uses.end(),
                    [&](const auto& use) {
                        return use.faceLocalVertexId == id;
                    });
                if (it == rec.source.uses.end())
                    return failure(SurfaceReason::InvalidLoopTopology);
                rec.chain.uses.push_back(*it);
            }
        } else {
            const int direction = directedEdgeInRing(
                ring, rec.source.uses[0].faceLocalVertexId,
                rec.source.uses[1].faceLocalVertexId);
            if (direction == 0)
                return failure(SurfaceReason::InvalidLoopTopology);
            if (direction < 0)
                std::reverse(rec.chain.uses.begin(), rec.chain.uses.end());
            for (std::size_t i = 0; i + 1 < rec.chain.uses.size(); ++i)
                if (directedEdgeInRing(
                        ring, rec.chain.uses[i].faceLocalVertexId,
                        rec.chain.uses[i + 1].faceLocalVertexId) != 1)
                    return failure(SurfaceReason::InvalidLoopTopology);
        }
    }

    auto chainPosition = [&](const ChainRecord& rec) {
        std::size_t position = 0;
        if (rec.loopOrdinal && rec.chain.uses.size() >= 2)
            (void)directedEdgeInRing(
                plan.loops[*rec.loopOrdinal],
                rec.chain.uses[0].faceLocalVertexId,
                rec.chain.uses[1].faceLocalVertexId, &position);
        return position;
    };
    auto chainUsesLess = [](const SurfaceChain& a, const SurfaceChain& b) {
        auto lessUse = [](const SurfaceChainVertex& x,
                          const SurfaceChainVertex& y) {
            if (x.sharedKey && y.sharedKey) {
                if (*x.sharedKey < *y.sharedKey) return true;
                if (*y.sharedKey < *x.sharedKey) return false;
            } else if (x.sharedKey.has_value() != y.sharedKey.has_value()) {
                return x.sharedKey.has_value();
            }
            return x.faceLocalVertexId < y.faceLocalVertexId;
        };
        return std::lexicographical_compare(a.uses.begin(), a.uses.end(),
                                            b.uses.begin(), b.uses.end(),
                                            lessUse);
    };
    std::sort(plan.chains.begin(), plan.chains.end(),
              [&](const ChainRecord& a, const ChainRecord& b) {
        const int groupA = a.chain.role == SurfaceChainRole::OuterBoundary
                               ? 0
                               : a.chain.role == SurfaceChainRole::HoleBoundary
                                     ? 1
                                     : 2;
        const int groupB = b.chain.role == SurfaceChainRole::OuterBoundary
                               ? 0
                               : b.chain.role == SurfaceChainRole::HoleBoundary
                                     ? 1
                                     : 2;
        if (groupA != groupB) return groupA < groupB;
        if (a.loopOrdinal != b.loopOrdinal)
            return a.loopOrdinal < b.loopOrdinal;
        if (a.loopOrdinal) {
            const auto pa = chainPosition(a);
            const auto pb = chainPosition(b);
            if (pa != pb) return pa < pb;
        } else {
            if (chainUsesLess(a.chain, b.chain)) return true;
            if (chainUsesLess(b.chain, a.chain)) return false;
        }
        return a.chain.chainPersistentId < b.chain.chainPersistentId;
    });
    return rebuildDerivedConstraints(plan);
}

void floodCells(const topo::TopologyWorkspace& ws,
                const std::map<HandleEdge, std::vector<topo::CellHandle>>& owners,
                const std::set<HandleEdge>& blocked,
                topo::CellHandle seed,
                std::set<topo::CellHandle>& removed) {
    std::vector<topo::CellHandle> stack{seed};
    while (!stack.empty()) {
        const auto ch = stack.back(); stack.pop_back();
        if (!ws.cells().isCurrent(ch) || !removed.insert(ch).second) continue;
        const auto cell = ws.cells().get(ch);
        for (std::size_t i = 0; i < 3; ++i) {
            const auto edge = handleEdge(cell.nodes[i], cell.nodes[(i + 1) % 3]);
            if (blocked.count(edge)) continue;
            const auto it = owners.find(edge);
            if (it == owners.end()) continue;
            for (const auto neighbor : it->second)
                if (!(neighbor == ch)) stack.push_back(neighbor);
        }
    }
}

SurfaceFailure physicalAudit(const SurfaceFaceRequest& request,
                             const ChartCertificate& cert,
                             const SurfaceGeometryProvider& geometry,
                             const SurfaceFaceResult& result) noexcept {
    try {
        const auto& uv = result.localMappingParams;
        if (uv.size() != result.localParams.size())
            return failure(SurfaceReason::LocalAuditFailed);
        for (std::size_t i = 0; i < result.localNodes.size(); ++i) {
            const Vec3 mapped = geometry.pointAt(request.chart, uv[i]);
            if (!finite(mapped))
                return failure(SurfaceReason::NonfiniteMappedPosition);
            if ((mapped - result.localNodes[i].coord).norm() >
                request.physicalPolicy.maxVertexMappingError)
                return failure(SurfaceReason::PhysicalOrientationFailed);
        }
        const int requiredWinding =
            cert.chartPositiveToNormalSign * request.desiredNormalSign;
        for (const auto& cell : result.localCells) {
            const auto s = predicates::orient2d(uv[cell.nodes[0]],
                                                uv[cell.nodes[1]],
                                                uv[cell.nodes[2]]);
            const int w = s == predicates::Sign::Positive ? 1 : -1;
            if (s == predicates::Sign::Zero || w != requiredWinding)
                return failure(SurfaceReason::PhysicalOrientationFailed);
            const Vec3 p0 = result.localNodes[cell.nodes[0]].coord;
            const Vec3 p1 = result.localNodes[cell.nodes[1]].coord;
            const Vec3 p2 = result.localNodes[cell.nodes[2]].coord;
            const Vec3 cross = (p1 - p0).cross(p2 - p0);
            const double area = cross.norm();
            if (!std::isfinite(area) ||
                area < request.physicalPolicy.minPhysicalDoubleArea)
                return failure(SurfaceReason::PhysicalOrientationFailed);
            const std::array<Vec2, 7> samples{
                uv[cell.nodes[0]], uv[cell.nodes[1]], uv[cell.nodes[2]],
                (uv[cell.nodes[0]] + uv[cell.nodes[1]]) * 0.5,
                (uv[cell.nodes[1]] + uv[cell.nodes[2]]) * 0.5,
                (uv[cell.nodes[2]] + uv[cell.nodes[0]]) * 0.5,
                (uv[cell.nodes[0]] + uv[cell.nodes[1]] + uv[cell.nodes[2]]) *
                    (1.0 / 3.0)};
            const std::array<Vec3, 7> linear{
                p0, p1, p2,
                (p0 + p1) * 0.5, (p1 + p2) * 0.5, (p2 + p0) * 0.5,
                (p0 + p1 + p2) * (1.0 / 3.0)};
            for (std::size_t i = 0; i < samples.size(); ++i) {
                const Vec3 normal = geometry.normalAt(request.chart, samples[i]);
                const Vec3 mapped = geometry.pointAt(request.chart, samples[i]);
                const double normalNorm = normal.norm();
                if (!finite(normal) || !std::isfinite(normalNorm) ||
                    normalNorm == 0.0)
                    return failure(SurfaceReason::NonfiniteMappedNormal);
                if (!finite(mapped))
                    return failure(SurfaceReason::NonfiniteMappedPosition);
                const double align = request.desiredNormalSign *
                    cross.dot(normal) / (area * normalNorm);
                if (align < request.physicalPolicy.minNormalAlignment ||
                    (mapped - linear[i]).norm() >
                        request.physicalPolicy.maxMappingDeviation)
                    return failure(SurfaceReason::PhysicalOrientationFailed);
            }
        }
        return {};
    } catch (const bnd::GeometryProviderDomainError&) {
        return failure(SurfaceReason::ProviderDomainRejected);
    } catch (const std::bad_alloc&) {
        return allocation(SurfaceAllocationPhase::Audit);
    } catch (...) {
        return failure(SurfaceReason::GeometryProviderError);
    }
}


struct PreparedFaceBuildState {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::NodeHandle> candidateNodes;
    std::vector<topo::CellHandle> candidateCells;
    bnd::ConstraintSet constraints;
    bnd::SteinerSegmentsReport protectedReport;
    bnd::SteinerSegmentsReport localReport;
};

bool exactNodeData(const topo::NodePool::NodeData& a,
                   const topo::NodePool::NodeData& b) noexcept {
    return a.coord.x == b.coord.x && a.coord.y == b.coord.y &&
           a.coord.z == b.coord.z && sameClass(a.cls, b.cls);
}

const bnd::ConstraintSegment* segmentById(
        const std::vector<bnd::ConstraintSegment>& segments,
        bnd::ConstraintEntityId id) noexcept {
    for (const auto& segment : segments)
        if (segment.id == id) return &segment;
    return nullptr;
}

std::optional<PreparedFaceBuildState> takePreparedFaceBuild(
        const FacePlan& plan, const SurfaceCandidate& candidate) {
    if (!candidate.preparedRecovery ||
        candidate.preparedRecovery.use_count() != 1)
        return std::nullopt;
    auto& prepared = *candidate.preparedRecovery;
    const auto& request = *plan.request;
    if (!prepared.finalized || !prepared.workspace ||
        !(prepared.providerStamp == candidate.providerStamp) ||
        !(prepared.faceStamp == candidate.faceStamp) ||
        prepared.chart != candidate.chart ||
        !(prepared.geometryStamp == request.expectedGeometryStamp) ||
        !prepared.workspaceStamp.matches(*prepared.workspace) ||
        prepared.frameBoundary != candidate.frameBoundary ||
        prepared.exteriorSeedCell != candidate.exteriorSeedCell ||
        prepared.candidateNodes.size() != candidate.nodes.size() ||
        prepared.candidateCells.size() != candidate.cells.size() ||
        prepared.constraints.segments.size() != plan.segments.size() ||
        !prepared.recovery.unattempted.empty() ||
        prepared.recovery.processed.size() != plan.segments.size())
        return std::nullopt;

    const auto activeNodes = prepared.workspace->nodes().activeBySequence();
    const auto activeCells = prepared.workspace->cells().activeBySequence();
    if (activeNodes != prepared.candidateNodes ||
        activeCells != prepared.candidateCells ||
        candidate.exteriorSeedCell >= prepared.candidateCells.size())
        return std::nullopt;

    for (std::size_t i = 0; i < candidate.nodes.size(); ++i) {
        const auto handle = prepared.candidateNodes[i];
        if (!prepared.workspace->nodes().isCurrent(handle) ||
            !exactNodeData(prepared.workspace->nodes().get(handle),
                           candidate.nodes[i].node))
            return std::nullopt;
        const auto param = prepared.constraints.params.find(handle);
        const auto auth = prepared.constraints.authParams.find(handle);
        if (param == prepared.constraints.params.end() ||
            param->second.chart != candidate.nodes[i].parameter.chart ||
            param->second.p.u != candidate.nodes[i].parameter.p.u ||
            param->second.p.v != candidate.nodes[i].parameter.p.v ||
            auth == prepared.constraints.authParams.end() ||
            auth->second.affine ||
            auth->second.chart != candidate.nodes[i].parameter.chart ||
            auth->second.storedUv.u != candidate.nodes[i].parameter.p.u ||
            auth->second.storedUv.v != candidate.nodes[i].parameter.p.v)
            return std::nullopt;
    }
    for (const auto& [id, index] : candidate.declaredVertexMap) {
        if (index >= prepared.candidateNodes.size()) return std::nullopt;
        const auto it = prepared.constraints.vertices.find(id);
        if (it == prepared.constraints.vertices.end() ||
            it->second != prepared.candidateNodes[index])
            return std::nullopt;
    }
    if (prepared.constraints.vertices.size() !=
        candidate.declaredVertexMap.size())
        return std::nullopt;

    for (std::size_t i = 0; i < candidate.cells.size(); ++i) {
        const auto& candidateCell = candidate.cells[i];
        if (candidateCell.nodes[0] >= prepared.candidateNodes.size() ||
            candidateCell.nodes[1] >= prepared.candidateNodes.size() ||
            candidateCell.nodes[2] >= prepared.candidateNodes.size())
            return std::nullopt;
        const auto cell = prepared.workspace->cells().get(
            prepared.candidateCells[i]);
        if (cell.type != topo::ElementType::Tri3 ||
            !sameClass(cell.cls, candidateCell.classification) ||
            cell.nodes[0] != prepared.candidateNodes[candidateCell.nodes[0]] ||
            cell.nodes[1] != prepared.candidateNodes[candidateCell.nodes[1]] ||
            cell.nodes[2] != prepared.candidateNodes[candidateCell.nodes[2]])
            return std::nullopt;
    }

    for (const auto& segment : plan.segments) {
        const auto* stored = segmentById(prepared.constraints.segments,
                                         segment.id);
        if (!stored || stored->chart != segment.chart ||
            !sameClass(stored->cls, segment.cls) ||
            !((stored->a == segment.a && stored->b == segment.b) ||
              (stored->a == segment.b && stored->b == segment.a)))
            return std::nullopt;
    }
    for (const auto& row : prepared.recovery.processed) {
        const auto* segment = segmentById(plan.segments, row.entity);
        if (!segment ||
            (row.result.status != bnd::SteinerStatus::AlreadyPresent &&
             row.result.status != bnd::SteinerStatus::Recovered) ||
            !row.result.createdParams.empty() ||
            row.result.recoveredChain.size() < 2)
            return std::nullopt;
        for (const auto handle : row.result.recoveredChain)
            if (!prepared.workspace->nodes().isCurrent(handle))
                return std::nullopt;
        for (std::size_t i = 0;
             i + 1 < row.result.recoveredChain.size(); ++i) {
            const auto edge = prepared.workspace->edges().find(
                prepared.workspace->nodes(), row.result.recoveredChain[i],
                row.result.recoveredChain[i + 1]);
            if (!edge) return std::nullopt;
            const auto data = prepared.workspace->edges().get(*edge);
            if (!data.requiredConstraint ||
                !sameClass(data.cls, segment->cls))
                return std::nullopt;
        }
    }

    PreparedFaceBuildState out;
    const std::size_t protectedCount = static_cast<std::size_t>(std::count_if(
        plan.segments.begin(), plan.segments.end(), [&](const auto& segment) {
            return plan.protectedIds.count(segment.id) != 0;
        }));
    out.protectedReport.processed.reserve(protectedCount);
    out.localReport.processed.reserve(plan.segments.size() - protectedCount);
    for (auto& row : prepared.recovery.processed) {
        if (plan.protectedIds.count(row.entity) != 0)
            out.protectedReport.processed.push_back(std::move(row));
        else
            out.localReport.processed.push_back(std::move(row));
    }
    out.workspace = std::move(prepared.workspace);
    out.candidateNodes = std::move(prepared.candidateNodes);
    out.candidateCells = std::move(prepared.candidateCells);
    out.constraints = std::move(prepared.constraints);
    out.constraints.stamp = bnd::OwnershipStamp::of(*out.workspace);
    return out;
}

FaceBuild buildFace(const FacePlan& plan, const SurfaceCandidate& candidate,
                    const mesh::Mesh& destination,
                    const std::map<SharedCurveNodeKey,
                                   SharedCurveNodeBinding>& ledger,
                    const SurfaceGeometryProvider& geometry,
                    const bnd::SteinerPolicy& localPolicy,
                    const ChartCertificate& certificate) {
    const auto& request = *plan.request;
    FaceBuild built;
    SurfaceAllocationPhase activePhase = SurfaceAllocationPhase::Candidate;
    try {
        std::unique_ptr<topo::TopologyWorkspace> workspace;
        std::vector<topo::NodeHandle> candidateNodes;
        std::vector<topo::CellHandle> candidateCells;
        bnd::ConstraintSet constraints;
        bnd::SteinerSegmentsReport protectedReport;
        bnd::SteinerSegmentsReport localReport;
        bnd::ConstraintSet authority;
        bool preparedRecoveryConsumed = false;

        auto prepared = takePreparedFaceBuild(plan, candidate);
        if (prepared) {
            activePhase = SurfaceAllocationPhase::Recovery;
            workspace = std::move(prepared->workspace);
            candidateNodes = std::move(prepared->candidateNodes);
            candidateCells = std::move(prepared->candidateCells);
            constraints = std::move(prepared->constraints);
            protectedReport = std::move(prepared->protectedReport);
            localReport = std::move(prepared->localReport);
            authority = std::move(constraints);
            built.error = verifyGeometryCheckpoint(
                request, geometry, certificate);
            if (built.error.reason != SurfaceReason::None) return built;
            preparedRecoveryConsumed = true;
        } else {
            topo::NodePool nodes;
            candidateNodes.reserve(candidate.nodes.size());
            for (const auto& n : candidate.nodes)
                candidateNodes.push_back(nodes.add(n.node.coord, n.node.cls));
            topo::CellPool cells;
            candidateCells.reserve(candidate.cells.size());
            for (const auto& c : candidate.cells) {
                const std::array<topo::NodeHandle, 3> cn{
                    candidateNodes[c.nodes[0]], candidateNodes[c.nodes[1]],
                    candidateNodes[c.nodes[2]]};
                candidateCells.push_back(cells.add(
                    topo::ElementType::Tri3,
                    std::span<const topo::NodeHandle>(cn.data(), cn.size()),
                    c.classification, nodes));
            }
            topo::EdgePool edges;
            for (const auto ch : cells.activeBySequence())
                edges.insertCellEdges(nodes, cells, ch);
            workspace = std::make_unique<topo::TopologyWorkspace>(
                std::move(nodes), std::move(cells), std::move(edges),
                topo::FacePool{});

            activePhase = SurfaceAllocationPhase::Recovery;
            for (std::size_t i = 0; i < candidate.nodes.size(); ++i) {
                constraints.params.emplace(candidateNodes[i],
                                           candidate.nodes[i].parameter);
                bnd::AuthoritativeParamRecord auth;
                auth.affine = false;
                auth.storedUv = candidate.nodes[i].parameter.p;
                auth.chart = request.chart;
                constraints.authParams.emplace(candidateNodes[i], auth);
            }
            for (const auto& [id, index] : candidate.declaredVertexMap)
                constraints.vertices.emplace(id, candidateNodes[index]);
            constraints.segments = plan.segments;
            constraints.stamp = bnd::OwnershipStamp::of(*workspace);

            built.error = verifyGeometryCheckpoint(
                request, geometry, certificate);
            if (built.error.reason != SurfaceReason::None) return built;

            bnd::SessionStatus sessionStatus = bnd::SessionStatus::Ok;
            auto session = bnd::RecoverySession::open(
                *workspace, constraints, sessionStatus);
            if (!session) {
                built.error =
                    sessionStatus == bnd::SessionStatus::IndexBuildFailed
                        ? allocation(SurfaceAllocationPhase::Recovery)
                        : failure(SurfaceReason::RecoveryStopped);
                return built;
            }
            std::set<topo::NodeHandle> sharedDeclared;
            for (const auto& rec : plan.chains)
                if (rec.chain.owningCurve)
                    for (const auto& use : rec.chain.uses)
                        sharedDeclared.insert(
                            constraints.vertices.at(use.faceLocalVertexId));
            for (const auto& segment : plan.segments) {
                if (plan.protectedIds.count(segment.id) == 0) continue;
                std::vector<topo::NodeHandle> chain;
                const auto route = bnd::routeEntity(
                    *workspace, session->input(), session->index(), segment,
                    chain);
                const auto& origin = plan.origins.at(segment.id);
                const auto* owner = chainById(
                    plan, origin.chainPersistentId);
                bool declared = !chain.empty();
                for (const auto h : chain)
                    if (sharedDeclared.count(h) == 0) declared = false;
                if ((route != bnd::Routing::AlreadyPresent &&
                     route != bnd::Routing::MarkOnly) || !declared) {
                    built.error = failure(
                        owner && owner->role ==
                                     SurfaceChainRole::EmbeddedConstraint
                            ? SurfaceReason::CurveArtifactSubdivisionRequired
                            : SurfaceReason::BoundarySteinerProhibited);
                    return built;
                }
            }

            bnd::ConstraintSet protectedSet = constraints;
            protectedSet.segments.erase(
                std::remove_if(
                    protectedSet.segments.begin(),
                    protectedSet.segments.end(), [&](const auto& segment) {
                        return plan.protectedIds.count(segment.id) == 0;
                    }),
                protectedSet.segments.end());
            if (!protectedSet.segments.empty()) {
                bnd::SteinerPolicy protectedPolicy = localPolicy;
                protectedPolicy.steinerEnabled = false;
                protectedReport = bnd::recoverSteinerSegments(
                    *workspace, protectedSet, geometry,
                    request.expectedGeometryStamp, protectedPolicy);
                for (const auto& row : protectedReport.processed)
                    if (row.result.status !=
                            bnd::SteinerStatus::AlreadyPresent &&
                        row.result.status != bnd::SteinerStatus::Recovered) {
                        built.error =
                            row.result.reject ==
                                    bnd::SteinerReject::PslgRejected
                                ? allocation(
                                      SurfaceAllocationPhase::Recovery)
                                : translateBnd(row.result);
                        return built;
                    }
                if (!protectedReport.unattempted.empty() ||
                    protectedReport.processed.size() !=
                        protectedSet.segments.size()) {
                    built.error = allocation(
                        SurfaceAllocationPhase::Recovery);
                    return built;
                }
            }

            bnd::ConstraintSet localSet = protectedSet;
            localSet.segments = constraints.segments;
            localSet.segments.erase(
                std::remove_if(
                    localSet.segments.begin(), localSet.segments.end(),
                    [&](const auto& segment) {
                        return plan.protectedIds.count(segment.id) != 0;
                    }),
                localSet.segments.end());
            localSet.stamp = bnd::OwnershipStamp::of(*workspace);
            if (!localSet.segments.empty()) {
                localReport = bnd::recoverSteinerSegments(
                    *workspace, localSet, geometry,
                    request.expectedGeometryStamp, localPolicy);
                for (const auto& row : localReport.processed)
                    if (row.result.status !=
                            bnd::SteinerStatus::AlreadyPresent &&
                        row.result.status != bnd::SteinerStatus::Recovered) {
                        built.error =
                            row.result.reject ==
                                    bnd::SteinerReject::PslgRejected
                                ? allocation(
                                      SurfaceAllocationPhase::Recovery)
                                : translateBnd(row.result);
                        return built;
                    }
                if (!localReport.unattempted.empty() ||
                    localReport.processed.size() !=
                        localSet.segments.size()) {
                    built.error = allocation(
                        SurfaceAllocationPhase::Recovery);
                    return built;
                }
            }
            authority = localSet.segments.empty()
                            ? std::move(protectedSet)
                            : std::move(localSet);
        }

        std::map<topo::NodeHandle,
                 std::optional<SurfaceNodeIdentity>> candidateIdentity;
        for (std::size_t i = 0; i < candidateNodes.size(); ++i)
            candidateIdentity.emplace(candidateNodes[i],
                                      candidate.nodes[i].persistentIdentity);

        activePhase = SurfaceAllocationPhase::DomainExtraction;
        if (workspace->state() == topo::WorkspaceState::CommittedAdjacencyStale)
            workspace->rebuildAdjacency();
        const auto owners = edgeOwners(*workspace);
        std::set<HandleEdge> blocked;
        std::map<bnd::ConstraintEntityId, const bnd::SteinerResult*> recovered;
        for (const auto& row : protectedReport.processed)
            recovered[row.entity] = &row.result;
        for (const auto& row : localReport.processed)
            recovered[row.entity] = &row.result;
        for (const auto& segment : plan.segments) {
            const auto rr = recovered.find(segment.id);
            if (rr == recovered.end()) {
                built.error = failure(SurfaceReason::RecoveryStopped);
                return built;
            }
            const auto& origin = plan.origins.at(segment.id);
            const auto* owner = chainById(plan, origin.chainPersistentId);
            if (!owner) {
                built.error = failure(SurfaceReason::InternalError);
                return built;
            }
            for (std::size_t i = 0;
                 i + 1 < rr->second->recoveredChain.size(); ++i) {
                const auto edge = handleEdge(rr->second->recoveredChain[i],
                                             rr->second->recoveredChain[i + 1]);
                if (owner->role == SurfaceChainRole::OuterBoundary ||
                    owner->role == SurfaceChainRole::HoleBoundary) {
                    const auto own = owners.find(edge);
                    if (own == owners.end() || own->second.size() != 2) {
                        built.error = domainFailure(
                            SurfaceDomainDetail::OneSidedDomainBoundary);
                        return built;
                    }
                    const auto pa = authority.params.at(edge[0]).p;
                    const auto pb = authority.params.at(edge[1]).p;
                    const auto p0 = authority.params.at(
                        thirdNode(workspace->cells().get(own->second[0]), edge)).p;
                    const auto p1 = authority.params.at(
                        thirdNode(workspace->cells().get(own->second[1]), edge)).p;
                    const auto s0 = predicates::orient2d(pa, pb, p0);
                    const auto s1 = predicates::orient2d(pa, pb, p1);
                    if (s0 == predicates::Sign::Zero ||
                        s1 == predicates::Sign::Zero || s0 == s1) {
                        built.error = domainFailure(
                            SurfaceDomainDetail::OneSidedDomainBoundary);
                        return built;
                    }
                    blocked.insert(edge);
                }
            }
        }
        if (candidate.exteriorSeedCell >= candidateCells.size() ||
            !workspace->cells().isCurrent(
                candidateCells[candidate.exteriorSeedCell])) {
            built.error = domainFailure(SurfaceDomainDetail::FloodFailed);
            return built;
        }
        std::set<topo::CellHandle> removed;
        floodCells(*workspace, owners, blocked,
                   candidateCells[candidate.exteriorSeedCell], removed);

        for (std::size_t hi = 1; hi < plan.loops.size(); ++hi) {
            const auto& ring = plan.loops[hi];
            const auto winding = ringWinding(ring, plan);
            if (winding == predicates::Sign::Zero) {
                built.error = failure(SurfaceReason::InvalidLoopTopology);
                return built;
            }
            const auto a = authority.vertices.at(ring[0]);
            const auto b = authority.vertices.at(ring[1]);
            const auto edge = handleEdge(a, b);
            const auto own = owners.find(edge);
            if (own == owners.end() || own->second.size() != 2) {
                built.error = domainFailure(
                    SurfaceDomainDetail::OneSidedDomainBoundary);
                return built;
            }
            std::optional<topo::CellHandle> holeSeed;
            const auto pa = authority.params.at(a).p;
            const auto pb = authority.params.at(b).p;
            for (const auto ch : own->second) {
                const auto third = thirdNode(workspace->cells().get(ch), edge);
                if (predicates::orient2d(pa, pb,
                        authority.params.at(third).p) == winding)
                    holeSeed = ch;
            }
            if (!holeSeed || removed.count(*holeSeed)) {
                built.error = domainFailure(SurfaceDomainDetail::FloodFailed);
                return built;
            }
            floodCells(*workspace, owners, blocked, *holeSeed, removed);
        }

        std::vector<topo::CellHandle> retainedCells;
        std::set<topo::NodeHandle> retainedNodes;
        for (const auto ch : workspace->cells().activeBySequence())
            if (removed.count(ch) == 0) {
                retainedCells.push_back(ch);
                const auto c = workspace->cells().get(ch);
                for (std::size_t i = 0; i < 3; ++i)
                    retainedNodes.insert(c.nodes[i]);
            }
        if (retainedCells.empty()) {
            built.error = domainFailure(SurfaceDomainDetail::FloodFailed);
            return built;
        }
        bool noFrameSurvivors = true;
        for (std::size_t i = 0; i < candidate.nodes.size(); ++i)
            if (candidate.nodes[i].frameOnly &&
                retainedNodes.count(candidateNodes[i]) != 0)
                noFrameSurvivors = false;

        topo::NodePool keptNodes;
        std::map<topo::NodeHandle, topo::NodeHandle> oldToNew;
        std::vector<topo::NodeHandle> oldNodeOrder;
        for (const auto old : workspace->nodes().activeBySequence())
            if (retainedNodes.count(old)) {
                const auto d = workspace->nodes().get(old);
                oldToNew.emplace(old, keptNodes.add(d.coord, d.cls));
                oldNodeOrder.push_back(old);
            }
        topo::CellPool keptCells;
        for (const auto oldCell : retainedCells) {
            const auto c = workspace->cells().get(oldCell);
            const std::array<topo::NodeHandle, 3> cn{
                oldToNew.at(c.nodes[0]), oldToNew.at(c.nodes[1]),
                oldToNew.at(c.nodes[2])};
            keptCells.add(topo::ElementType::Tri3,
                          std::span<const topo::NodeHandle>(cn.data(), cn.size()),
                          c.cls, keptNodes);
        }
        topo::EdgePool keptEdges;
        for (const auto ch : keptCells.activeBySequence())
            keptEdges.insertCellEdges(keptNodes, keptCells, ch);
        for (const auto oldEdge : workspace->edges().activeBySequence()) {
            const auto d = workspace->edges().get(oldEdge);
            if (!d.requiredConstraint || oldToNew.count(d.a) == 0 ||
                oldToNew.count(d.b) == 0)
                continue;
            const auto newEdge = keptEdges.find(
                keptNodes, oldToNew.at(d.a), oldToNew.at(d.b));
            if (!newEdge) {
                built.error = domainFailure(
                    SurfaceDomainDetail::RetainedRebuildFailed);
                return built;
            }
            keptEdges.setRequiredConstraint(*newEdge, true);
            keptEdges.setClassification(*newEdge, d.cls);
        }
        auto kept = std::make_unique<topo::TopologyWorkspace>(
            std::move(keptNodes), std::move(keptCells), std::move(keptEdges),
            topo::FacePool{});

        auto result = std::make_unique<SurfaceFaceResult>();
        result->faceStamp = request.faceStamp;
        result->chart = request.chart;
        result->faceClass = request.faceClass;
        result->desiredNormalSign = request.desiredNormalSign;
        result->physicalPolicy = request.physicalPolicy;
        result->expectedGeometryStamp = request.expectedGeometryStamp;
        result->chartCertificate = certificate;
        result->normalizedChains.reserve(plan.chains.size());
        for (const auto& rec : plan.chains)
            result->normalizedChains.push_back(rec.chain);
        result->segmentOrigins = plan.origins;
        result->protectedRecovery = std::move(protectedReport);
        result->localRecovery = std::move(localReport);
        result->preparedRecoveryConsumed = preparedRecoveryConsumed;
        result->localWorkspace = std::move(kept);

        std::map<topo::NodeHandle, std::size_t> localIndex;
        for (std::size_t i = 0; i < oldNodeOrder.size(); ++i) {
            const auto old = oldNodeOrder[i];
            const auto current = oldToNew.at(old);
            localIndex.emplace(current, i);
            result->localHandles.push_back(current);
            result->localNodes.push_back(result->localWorkspace->nodes().get(current));
            const auto identity = candidateIdentity.find(old);
            if(identity!=candidateIdentity.end()&&identity->second)
                result->refinedIdentityByLocalNode.emplace(
                    i,*identity->second);
            const auto p = authority.params.at(old).p;
            result->localMappingParams.push_back(p);
            const auto auth = authority.authParams.find(old);
            if (auth != authority.authParams.end() && auth->second.affine)
                result->localParams.push_back(bnd::AffineParam{
                    auth->second.lo, auth->second.hi, auth->second.num,
                    auth->second.level});
            else
                result->localParams.push_back(bnd::StoredParam{p});
        }
        std::set<SurfaceNodeIdentity> candidateRefinedIdentities;
        for(const auto& node:candidate.nodes)
            if(node.persistentIdentity)
                candidateRefinedIdentities.insert(*node.persistentIdentity);
        std::set<SurfaceNodeIdentity> retainedRefinedIdentities;
        for(const auto& [local,identity]:result->refinedIdentityByLocalNode) {
            (void)local;
            retainedRefinedIdentities.insert(identity);
        }
        if(candidateRefinedIdentities!=retainedRefinedIdentities) {
            built.error=failure(SurfaceReason::SharedIdentityConflict);
            return built;
        }
        std::size_t ci = 0;
        for (const auto ch : result->localWorkspace->cells().activeBySequence()) {
            const auto c = result->localWorkspace->cells().get(ch);
            SurfaceFaceCell out;
            out.classification = c.cls;
            for (std::size_t i = 0; i < 3; ++i)
                out.nodes[i] = localIndex.at(c.nodes[i]);
            result->localCells.push_back(out);
            ++ci;
        }
        (void)ci;
        for (const auto& ring : plan.loops) {
            std::vector<std::size_t> out;
            out.reserve(ring.size());
            for (const auto id : ring) {
                const auto old = authority.vertices.at(id);
                if (oldToNew.count(old) == 0) {
                    built.error = domainFailure(
                        SurfaceDomainDetail::RetainedRebuildFailed);
                    return built;
                }
                out.push_back(localIndex.at(oldToNew.at(old)));
            }
            result->boundaryLoops.push_back(std::move(out));
        }
        for (const auto& rec : plan.chains)
            for (const auto& use : rec.chain.uses)
                if (use.sharedKey) {
                    const auto old = authority.vertices.at(
                        use.faceLocalVertexId);
                    if (oldToNew.count(old) == 0) {
                        built.error = domainFailure(
                            SurfaceDomainDetail::RetainedRebuildFailed);
                        return built;
                    }
                    const auto local = localIndex.at(oldToNew.at(old));
                    const auto global = ledger.at(*use.sharedKey).existingGlobalIndex;
                    const auto inserted = result->existingGlobalIndex.emplace(
                        local, global);
                    if (!inserted.second && inserted.first->second != global) {
                        built.error = failure(SurfaceReason::SharedIdentityConflict);
                        return built;
                    }
                    if ((result->localNodes[local].coord -
                         destination.node(global)).norm() >
                        request.physicalPolicy.maxArtifactPositionError) {
                        built.error = failure(
                            SurfaceReason::SharedMappedPositionMismatch);
                        return built;
                    }
                }

        for (const auto& segment : plan.segments) {
            const auto rr = recovered.at(segment.id);
            const auto& origin = plan.origins.at(segment.id);
            const auto* chain = chainById(plan, origin.chainPersistentId);
            for (std::size_t i = 0; i + 1 < rr->recoveredChain.size(); ++i) {
                if (oldToNew.count(rr->recoveredChain[i]) == 0 ||
                    oldToNew.count(rr->recoveredChain[i + 1]) == 0) {
                    built.error = failure(SurfaceReason::LocalAuditFailed);
                    return built;
                }
                PersistentConstraintEdge e;
                e.faceStamp = request.faceStamp;
                e.bndEntityId = segment.id;
                e.chainPersistentId = origin.chainPersistentId;
                e.originalSegmentOrdinal = origin.originalSegmentOrdinal;
                e.normalizedSegmentOrdinal = origin.normalizedSegmentOrdinal;
                e.role = chain->role;
                e.localEndpoints = {oldToNew.at(rr->recoveredChain[i]),
                                    oldToNew.at(rr->recoveredChain[i + 1])};
                const auto originalA = authority.vertices.at(
                    origin.endpoints[0].faceLocalVertexId);
                const auto originalB = authority.vertices.at(
                    origin.endpoints[1].faceLocalVertexId);
                const auto keyFor = [&](topo::NodeHandle old) {
                    if (old == originalA) return origin.endpoints[0].sharedKey;
                    if (old == originalB) return origin.endpoints[1].sharedKey;
                    return std::optional<SharedCurveNodeKey>{};
                };
                e.endpointSharedKeys = {keyFor(rr->recoveredChain[i]),
                                        keyFor(rr->recoveredChain[i + 1])};
                e.classification = origin.classification;
                result->persistentEdges.push_back(std::move(e));
            }
        }
        activePhase = SurfaceAllocationPhase::Audit;
        built.error = verifyGeometryCheckpoint(request, geometry, certificate);
        if (built.error.reason != SurfaceReason::None) return built;
        SurfaceFailure p = physicalAudit(request, certificate, geometry, *result);
        if (p.reason != SurfaceReason::None) {
            built.error = p;
            return built;
        }
        auto& audit = result->audit;
        const auto& rws = *result->localWorkspace;
        audit.workspaceReady = topo::workspaceTopologyReady(
            rws.nodes(), rws.cells(), rws.edges(), rws.faces()) &&
            rws.adjacencyReady();
        audit.nodeLinksValid =
            bnd::auditTri3Manifold(rws) == bnd::Tri3AuditResult::Ok;
        audit.noFrameSurvivors = noFrameSurvivors;
        audit.physicalOrientation = true;
        std::map<EdgeIndex, std::vector<std::size_t>> localOwners;
        std::map<EdgeIndex, std::array<std::size_t, 2>> firstDirection;
        std::set<std::array<std::size_t, 3>> triKeys;
        bool classes = true;
        for (std::size_t i = 0; i < result->localCells.size(); ++i) {
            const auto& c = result->localCells[i];
            if (!sameClass(c.classification, request.faceClass)) classes = false;
            auto key = c.nodes;
            std::sort(key.begin(), key.end());
            if (!triKeys.insert(key).second) audit.noDuplicateTriangles = false;
            for (std::size_t j = 0; j < 3; ++j) {
                const std::array<std::size_t, 2> directed{
                    c.nodes[j], c.nodes[(j + 1) % 3]};
                const auto edge = edgeIndex(directed[0], directed[1]);
                if (localOwners[edge].empty())
                    firstDirection.emplace(edge, directed);
                localOwners[edge].push_back(i);
            }
        }
        audit.noDuplicateTriangles =
            triKeys.size() == result->localCells.size();
        std::set<std::array<std::size_t, 2>> actualBoundary;
        std::set<std::array<std::size_t, 2>> expectedBoundary;
        bool usesValid = true;
        for (const auto& [edge, own] : localOwners) {
            if (own.size() == 1)
                actualBoundary.insert(firstDirection.at(edge));
            else if (own.size() != 2) usesValid = false;
        }
        for (const auto& ring : result->boundaryLoops)
            for (std::size_t i = 0; i < ring.size(); ++i)
                if (!expectedBoundary.insert(
                         {ring[i], ring[(i + 1) % ring.size()]}).second)
                    usesValid = false;
        audit.boundaryCycleCount = result->boundaryLoops.size();
        audit.eulerCharacteristic =
            static_cast<std::int64_t>(result->localNodes.size()) -
            static_cast<std::int64_t>(localOwners.size()) +
            static_cast<std::int64_t>(result->localCells.size());
        audit.boundaryCyclesExact =
            usesValid && actualBoundary == expectedBoundary &&
            audit.boundaryCycleCount == 1 + request.holes.size() &&
            audit.eulerCharacteristic ==
                1 - static_cast<std::int64_t>(request.holes.size());
        std::set<std::size_t> connected;
        std::vector<std::size_t> stack;
        if (!result->localCells.empty()) stack.push_back(0);
        while (!stack.empty()) {
            const auto c = stack.back(); stack.pop_back();
            if (!connected.insert(c).second) continue;
            const auto& cell = result->localCells[c];
            for (std::size_t j = 0; j < 3; ++j)
                for (const auto n : localOwners[
                         edgeIndex(cell.nodes[j], cell.nodes[(j + 1) % 3])])
                    if (n != c) stack.push_back(n);
        }
        audit.connectedGenusZero =
            connected.size() == result->localCells.size() &&
            audit.eulerCharacteristic ==
                1 - static_cast<std::int64_t>(request.holes.size());
        bool persistent = true;
        for (const auto& e : result->persistentEdges) {
            const auto eh = rws.edges().find(rws.nodes(), e.localEndpoints[0],
                                             e.localEndpoints[1]);
            if (!eh) { persistent = false; continue; }
            const auto d = rws.edges().get(*eh);
            if (!d.requiredConstraint || !sameClass(d.cls, e.classification))
                persistent = false;
            const auto a = localIndex.at(e.localEndpoints[0]);
            const auto b = localIndex.at(e.localEndpoints[1]);
            const auto own = localOwners.find(edgeIndex(a, b));
            if (own == localOwners.end() ||
                (e.role == SurfaceChainRole::EmbeddedConstraint &&
                 own->second.size() != 2))
                persistent = false;
        }
        audit.constraintsPersistent = persistent;
        audit.classificationsPreserved = classes;
        bool shared = true;
        for (const auto& [local, global] : result->existingGlobalIndex)
            if (local >= result->localNodes.size() || global < 0 ||
                global >= destination.nodeCount() ||
                !sameClass(result->localNodes[local].cls,
                           destination.nodeClass(global)) ||
                (result->localNodes[local].coord - destination.node(global)).norm() >
                    request.physicalPolicy.maxArtifactPositionError)
                shared = false;
        audit.sharedMappingsPreserved = shared;
        if (!audit.passed()) {
            built.error = failure(SurfaceReason::LocalAuditFailed);
            return built;
        }
        built.result = std::move(result);
        return built;
    } catch (const std::bad_alloc&) {
        built.error = allocation(activePhase);
        return built;
    } catch (const std::exception&) {
        built.error = failure(SurfaceReason::InternalError);
        return built;
    } catch (...) {
        built.error = failure(SurfaceReason::InternalError);
        return built;
    }
}

}  // namespace

bool SurfaceTopologyAudit::passed() const noexcept {
    return workspaceReady && connectedGenusZero && boundaryCyclesExact &&
           nodeLinksValid && constraintsPersistent &&
           classificationsPreserved && sharedMappingsPreserved &&
           noFrameSurvivors && noDuplicateTriangles && physicalOrientation;
}

SurfaceFailure publishSurfaceSet(mesh::Mesh& destination,
                                 SurfaceSetResult& artifact,
                                 const SurfaceGeometryProvider& geometry,
                                 bool& published) noexcept {
#if defined(PAACS_TOPO_TESTING)
    const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
        t2c_detail::SurfaceCacheAllocationZone::SetAudit);
#endif
    published = false;
    try {
        struct FacePublicationPlan {
            SurfaceFaceResult* face = nullptr;
            std::vector<Index> localToGlobal;
            std::vector<Index> cellToGlobal;
            std::vector<std::vector<Index>> connectivity;
            std::map<topo::NodeHandle, std::size_t> localIndex;
        };
        std::vector<FacePublicationPlan> plans;
        plans.reserve(artifact.faces.size());
        std::size_t newNodes = 0;
        std::size_t newCells = 0;
        std::size_t newConnectivity = 0;
        Index next = destination.nodeCount();
        Index nextCell = destination.cellCount();
        std::set<std::array<Index, 3>> triangleKeys;
        for (Index ci = 0; ci < destination.cellCount(); ++ci) {
            const auto& cell = destination.cell(ci);
            if (cell.type != mesh::CellType::Tri3) continue;
            const Index* conn = destination.cellNodes(ci);
            std::array<Index, 3> key{conn[0], conn[1], conn[2]};
            std::sort(key.begin(), key.end());
            triangleKeys.insert(key);
        }
        for (auto& facePtr : artifact.faces) {
            if (!facePtr || !facePtr->audit.passed() ||
                !facePtr->localWorkspace ||
                facePtr->localNodes.size() != facePtr->localHandles.size() ||
                facePtr->localNodes.size() != facePtr->localParams.size() ||
                facePtr->localNodes.size() !=
                    facePtr->localMappingParams.size())
                return failure(SurfaceReason::PublicationValidationFailed);
            auto& face = *facePtr;
            SurfaceFaceRequest checkpoint;
            checkpoint.chart = face.chart;
            checkpoint.expectedGeometryStamp = face.expectedGeometryStamp;
            checkpoint.expectedChartCertificateStamp =
                face.chartCertificate.certificateStamp;
            SurfaceFailure check = verifyGeometryCheckpoint(
                checkpoint, geometry, face.chartCertificate);
            if (check.reason != SurfaceReason::None) return check;
            FacePublicationPlan plan;
            plan.face = &face;
            plan.localToGlobal.resize(face.localNodes.size(), kInvalidIndex);
            plan.cellToGlobal.reserve(face.localCells.size());
            face.publishedGlobalByLocalNode.assign(face.localNodes.size(),
                                                   kInvalidIndex);
            face.publishedGlobalCellByLocalCell.assign(face.localCells.size(),
                                                       kInvalidIndex);
            for (std::size_t i = 0; i < face.localHandles.size(); ++i)
                plan.localIndex.emplace(face.localHandles[i], i);
            for (std::size_t i = 0; i < face.localNodes.size(); ++i) {
                const auto existing = face.existingGlobalIndex.find(i);
                if (existing != face.existingGlobalIndex.end()) {
                    if (existing->second < 0 ||
                        existing->second >= destination.nodeCount() ||
                        !sameClass(face.localNodes[i].cls,
                                   destination.nodeClass(existing->second)) ||
                        (face.localNodes[i].coord -
                         destination.node(existing->second)).norm() >
                            face.physicalPolicy.maxArtifactPositionError)
                        return failure(
                            SurfaceReason::PublicationValidationFailed);
                    plan.localToGlobal[i] = existing->second;
                } else {
                    if (next == std::numeric_limits<Index>::max())
                        return failure(
                            SurfaceReason::PublicationValidationFailed);
                    plan.localToGlobal[i] = next++;
                    ++newNodes;
                }
            }
            plan.connectivity.reserve(face.localCells.size());
            for (const auto& cell : face.localCells) {
                if (!sameClass(cell.classification, face.faceClass) ||
                    cell.nodes[0] >= plan.localToGlobal.size() ||
                    cell.nodes[1] >= plan.localToGlobal.size() ||
                    cell.nodes[2] >= plan.localToGlobal.size())
                    return failure(
                        SurfaceReason::PublicationValidationFailed);
                std::array<Index, 3> key{
                    plan.localToGlobal[cell.nodes[0]],
                    plan.localToGlobal[cell.nodes[1]],
                    plan.localToGlobal[cell.nodes[2]]};
                auto sorted = key;
                std::sort(sorted.begin(), sorted.end());
                if (!triangleKeys.insert(sorted).second)
                    return failure(
                        SurfaceReason::PublicationValidationFailed);
                plan.connectivity.push_back({key[0], key[1], key[2]});
                if (nextCell == std::numeric_limits<Index>::max())
                    return failure(
                        SurfaceReason::PublicationValidationFailed);
                plan.cellToGlobal.push_back(nextCell++);
                ++newCells;
                if (newConnectivity >
                    std::numeric_limits<std::size_t>::max() - 3)
                    return failure(
                        SurfaceReason::PublicationValidationFailed);
                newConnectivity += 3;
            }
            plans.push_back(std::move(plan));
        }
        artifact.publishedGlobalByRefinedIdentity.clear();
        std::map<Index,SurfaceNodeIdentity> refinedIdentityByGlobal;
        for(const auto& plan:plans)
            for(const auto& [local,identity]:
                    plan.face->refinedIdentityByLocalNode) {
                if(local>=plan.localToGlobal.size()||
                   !validSurfaceNodeIdentity(identity)||
                   identity.domain==
                       SurfaceNodeIdentityDomain::OriginalDiscretization)
                    return failure(SurfaceReason::PublicationValidationFailed);
                const Index global=plan.localToGlobal[local];
                const auto [it,inserted]=
                    artifact.publishedGlobalByRefinedIdentity.emplace(
                        identity,global);
                if(!inserted&&it->second!=global)
                    return failure(SurfaceReason::SharedIdentityConflict);
                const auto [globalIt,globalInserted]=
                    refinedIdentityByGlobal.emplace(global,identity);
                if(!globalInserted&&globalIt->second!=identity)
                    return failure(SurfaceReason::SharedIdentityConflict);
            }
        std::set<SurfaceNodeIdentity> sharedRefinedIdentities;
        for (const auto& [key, binding] : artifact.sharedLedger) {
            if (binding.existingGlobalIndex < 0 ||
                binding.existingGlobalIndex >= destination.nodeCount() ||
                !sameClass(binding.classification,
                           destination.nodeClass(binding.existingGlobalIndex)))
                return failure(SurfaceReason::PublicationValidationFailed);
            if(key.identity.domain!=
                   SurfaceNodeIdentityDomain::OriginalDiscretization) {
                if(key.identity.domain!=SurfaceNodeIdentityDomain::RefinedCurve)
                    return failure(SurfaceReason::SharedIdentityConflict);
                const auto identity=
                    artifact.publishedGlobalByRefinedIdentity.find(key.identity);
                if(identity==artifact.publishedGlobalByRefinedIdentity.end()||
                   identity->second!=binding.existingGlobalIndex)
                    return failure(SurfaceReason::SharedIdentityConflict);
                sharedRefinedIdentities.insert(key.identity);
            }
        }
        for(const auto& [identity,global]:
                artifact.publishedGlobalByRefinedIdentity) {
            (void)global;
            if(identity.domain==SurfaceNodeIdentityDomain::RefinedCurve&&
               sharedRefinedIdentities.count(identity)==0)
                return failure(SurfaceReason::SharedIdentityConflict);
        }
        for (const auto& [key, binding] : artifact.sharedLedger) {
            (void)binding;
            struct MappedUse { Vec3 point; double tolerance; };
            std::vector<MappedUse> mapped;
            for (const auto& facePlan : plans) {
                const auto& face = *facePlan.face;
                for (const auto& edge : face.persistentEdges)
                    for (std::size_t endpoint = 0; endpoint < 2; ++endpoint)
                        if (edge.endpointSharedKeys[endpoint] &&
                            *edge.endpointSharedKeys[endpoint] == key) {
                            const auto it = facePlan.localIndex.find(
                                edge.localEndpoints[endpoint]);
                            if (it == facePlan.localIndex.end())
                                return failure(
                                    SurfaceReason::PublicationValidationFailed);
                            const auto& node = face.localNodes[it->second];
                            bool exists = false;
                            for (const auto& prior : mapped)
                                if (samePoint(prior.point, node.coord) &&
                                    prior.tolerance ==
                                        face.physicalPolicy.maxArtifactPositionError)
                                    exists = true;
                            if (!exists)
                                mapped.push_back({
                                    node.coord,
                                    face.physicalPolicy.maxArtifactPositionError});
                        }
            }
            for (std::size_t i = 0; i < mapped.size(); ++i)
                for (std::size_t j = i + 1; j < mapped.size(); ++j)
                    if ((mapped[i].point - mapped[j].point).norm() >
                        std::min(mapped[i].tolerance, mapped[j].tolerance))
                        return failure(
                            SurfaceReason::SharedMappedPositionMismatch);
        }
        try {
            destination.reserveAppend(newNodes, newCells, newConnectivity);
        } catch (const std::bad_alloc&) {
            return allocation(SurfaceAllocationPhase::PublicationPreparation);
        } catch (const std::length_error&) {
            return failure(SurfaceReason::PublicationValidationFailed);
        }

        auto append = [&]() noexcept {
#if defined(PAACS_TOPO_TESTING)
            const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
                t2c_detail::SurfaceCacheAllocationZone::CommitTail);
#endif
            for (auto& plan : plans) {
                auto& face = *plan.face;
                for (std::size_t i = 0; i < face.localNodes.size(); ++i)
                    if (face.existingGlobalIndex.count(i) == 0) {
                        const Index added = destination.addNode(
                            face.localNodes[i].coord, face.localNodes[i].cls);
                        if (added != plan.localToGlobal[i]) std::terminate();
                    }
                for (std::size_t i = 0; i < face.localNodes.size(); ++i)
                    face.publishedGlobalByLocalNode[i] = plan.localToGlobal[i];
                for (std::size_t i = 0; i < face.localCells.size(); ++i) {
                    const Index added = destination.addCell(
                        mesh::CellType::Tri3, plan.connectivity[i],
                        face.localCells[i].classification);
                    if (added != plan.cellToGlobal[i]) std::terminate();
                    face.publishedGlobalCellByLocalCell[i] = added;
                }
                for (auto& edge : face.persistentEdges)
                    edge.publishedGlobalEndpoints = std::array<Index, 2>{
                        plan.localToGlobal[
                            plan.localIndex.at(edge.localEndpoints[0])],
                        plan.localToGlobal[
                            plan.localIndex.at(edge.localEndpoints[1])]};
            }
        };
        append();
        published = true;
        return {};
    } catch (const std::bad_alloc&) {
        return allocation(SurfaceAllocationPhase::PublicationPreparation);
    } catch (...) {
        return failure(SurfaceReason::PublicationValidationFailed);
    }
}

SurfaceSetOutcome buildAndPublishSurfaceSet(
    const SurfaceSetRequest& request, mesh::Mesh& destination,
    const SurfaceGeometryProvider& geometry,
    const SurfaceCandidateProvider& candidateProvider,
    const bnd::SteinerPolicy& localSteinerPolicy) noexcept {
    SurfaceSetOutcome out;
    try {
        out.processed.reserve(request.faces.size());
        out.unattempted.reserve(request.faces.size());
    } catch (...) {
        out.setStatus = SurfaceSetStatus::Failed;
        out.setFailure = allocation(SurfaceAllocationPhase::SetMapping);
        out.reportComplete = false;
        return out;
    }
    out.reportComplete = true;
    std::vector<FacePersistentId> ids;
    try {
        ids.reserve(request.faces.size());
        for (const auto& face : request.faces)
            ids.push_back(face.faceStamp.persistentId);
        std::sort(ids.begin(), ids.end());
        SetPlan plan;
        SurfaceFailure preflight = buildSetPlan(request, destination, plan);
        if (preflight.reason != SurfaceReason::None) {
            const auto status = statusForFailure(preflight);
            for (const auto id : ids)
                out.processed.emplace_back(
                    SurfaceFaceRow{id, status, preflight, nullptr});
            out.setStatus = setStatusFor(status);
            out.setFailure = preflight;
            return out;
        }
        CandidateProviderStamp providerStamp;
        try {
            providerStamp = candidateProvider.stamp();
        } catch (const std::bad_alloc&) {
            out.setFailure = allocation(SurfaceAllocationPhase::Candidate);
            for (const auto id : ids) out.unattempted.push_back(id);
            return out;
        } catch (...) {
            out.setFailure = failure(SurfaceReason::CandidateProviderError);
            for (const auto id : ids) out.unattempted.push_back(id);
            return out;
        }
        for (std::size_t fi = 0; fi < plan.faces.size(); ++fi) {
            auto& fp = plan.faces[fi];
            const auto& face = *fp.request;
            SurfaceFailure err;
            std::optional<ChartCertificate> certificate;
            err = geometryCheckpoint(face, geometry, certificate, true);
            if (err.reason == SurfaceReason::None)
                err = normalizeFacePlan(fp, *certificate);
            if (err.reason == SurfaceReason::None &&
                !(providerStamp == face.expectedCandidateProviderStamp))
                err = failure(SurfaceReason::StaleCandidateProvider);
            SurfaceCandidate candidate;
            if (err.reason == SurfaceReason::None) {
                try {
                    candidate = candidateProvider.build(face);
                } catch (const std::bad_alloc&) {
                    err = allocation(SurfaceAllocationPhase::Candidate);
                } catch (...) {
                    err = failure(SurfaceReason::CandidateProviderError);
                }
            }
            if (err.reason == SurfaceReason::None)
                err = geometryCheckpoint(face, geometry, certificate, false);
            if (err.reason == SurfaceReason::None)
                err = validateCandidate(face, fp, candidate, *certificate);
            FaceBuild built;
            if (err.reason == SurfaceReason::None)
                built = buildFace(fp, candidate, destination, plan.ledger,
                                  geometry, localSteinerPolicy, *certificate);
            if (err.reason == SurfaceReason::None) err = built.error;
            if (err.reason != SurfaceReason::None) {
                const auto status = statusForFailure(err);
                out.processed.emplace_back(SurfaceFaceRow{
                    face.faceStamp.persistentId, status, err, nullptr});
                for (std::size_t later = fi + 1; later < plan.faces.size();
                     ++later)
                    out.unattempted.push_back(
                        plan.faces[later].request->faceStamp.persistentId);
                out.setStatus = setStatusFor(status);
                out.setFailure = err;
                return out;
            }
            out.processed.emplace_back(SurfaceFaceRow{
                face.faceStamp.persistentId, SurfaceFaceStatus::Ready, {},
                std::move(built.result)});
        }

        auto artifact = std::make_unique<SurfaceSetResult>();
        artifact->setStamp = request.setStamp;
        artifact->faces.reserve(out.processed.size());
        artifact->sharedLedger = std::move(plan.ledger);
        for (auto& row : out.processed)
            artifact->faces.push_back(std::move(row.readyResult));
        bool published = false;
        SurfaceFailure publication =
            publishSurfaceSet(destination, *artifact, geometry, published);
        out.retainedArtifact = std::move(artifact);
        if (!published) {
            out.setStatus = setStatusFor(statusForFailure(publication));
            out.setFailure = publication;
            return out;
        }
        out.setStatus = SurfaceSetStatus::Published;
        out.setFailure = {};
        return out;
    } catch (const std::bad_alloc&) {
        out.setStatus = SurfaceSetStatus::Failed;
        out.setFailure = allocation(SurfaceAllocationPhase::SetMapping);
        for (const auto& face : request.faces) {
            const auto id = face.faceStamp.persistentId;
            bool seen = false;
            for (const auto& row : out.processed)
                if (row.facePersistentId == id) seen = true;
            for (const auto pending : out.unattempted)
                if (pending == id) seen = true;
            if (!seen) out.unattempted.push_back(id);
        }
        std::sort(out.unattempted.begin(), out.unattempted.end());
        return out;
    } catch (...) {
        out.setStatus = SurfaceSetStatus::Failed;
        out.setFailure = failure(SurfaceReason::InternalError);
        for (const auto& face : request.faces) {
            const auto id = face.faceStamp.persistentId;
            bool seen = false;
            for (const auto& row : out.processed)
                if (row.facePersistentId == id) seen = true;
            for (const auto pending : out.unattempted)
                if (pending == id) seen = true;
            if (!seen) out.unattempted.push_back(id);
        }
        std::sort(out.unattempted.begin(), out.unattempted.end());
        return out;
    }
}

SurfaceSetOutcome t2c_detail::buildAndPublishSurfaceSetCached(
    const SurfaceSetRequest& request,mesh::Mesh& destination,
    const SurfaceGeometryProvider& geometry,
    const SurfaceCandidateProvider& candidateProvider,
    const bnd::SteinerPolicy& localSteinerPolicy,
    const std::map<FacePersistentId,std::vector<std::uint8_t>>& authorities,
    SurfaceFaceCacheSession& cache,
    const SurfaceCacheSizingEvaluator& sizingEvaluator,
    const SurfaceCacheMetricsSource& metricsSource) noexcept {
    SurfaceSetOutcome out;
    auto* cacheReport=cache.diagnosticChannel();
    const auto appendProcessed=[&](FacePersistentId face,
            SurfaceCacheFaceDisposition disposition,
            SurfaceFaceStatus status,const SurfaceFailure& failure) noexcept {
        if(!cacheReport||
           cacheReport->processedUsed>=cacheReport->processedStorage.size()) {
            if(cacheReport)cacheReport->reportComplete=false;
            return false;
        }
        cacheReport->processedStorage[cacheReport->processedUsed++]={
            face,disposition,status,failure};
        return true;
    };
    const auto appendUnattempted=[&](FacePersistentId face) noexcept {
        if(!cacheReport||cacheReport->unattemptedUsed>=
                             cacheReport->unattemptedStorage.size()) {
            if(cacheReport)cacheReport->reportComplete=false;
            return false;
        }
        cacheReport->unattemptedStorage[cacheReport->unattemptedUsed++]=face;
        return true;
    };
    const auto finishCacheReport=[&](bool complete,
                                     FacePersistentId failedFace) noexcept {
        if(!cacheReport)return;
        cacheReport->reportComplete=complete&&
            cacheReport->processedUsed+cacheReport->unattemptedUsed==
                request.faces.size();
        cacheReport->hasFailedFace=failedFace!=0;
        cacheReport->failedFacePersistentId=failedFace;
        // The scalar channel and the fixed report envelope are one outcome.
        // Keep their completeness/count fields synchronized on every exit,
        // including failures raised internally by cache preparation helpers.
        const auto diagnostic=cache.diagnostic();
        if(diagnostic.status!=SurfaceFaceCacheStatus::None)
            cache.setDiagnostic(diagnostic.status,diagnostic.reason,
                diagnostic.phase,cacheReport->reportComplete,
                failedFace!=0 ? failedFace : diagnostic.failedFace,
                cacheReport->processedUsed,cacheReport->unattemptedUsed);
    };
    const auto appendRemainingProcessed=[&](std::size_t first,
            const std::vector<FacePersistentId>& ordered,
            const SurfaceFailure& failure) noexcept {
        const auto status=statusForFailure(failure);
        bool complete=true;
        for(std::size_t i=first;i<ordered.size();++i)
            complete=appendProcessed(ordered[i],
                SurfaceCacheFaceDisposition::Stopped,status,failure)&&complete;
        return complete;
    };
    std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>> nextSlots;
    struct ReplacementScratchGuard {
        SurfaceFaceCacheSession* cache = nullptr;
        std::vector<std::shared_ptr<const CachedSurfaceFaceArtifact>>* rows =
            nullptr;
        ~ReplacementScratchGuard() noexcept {
            if(cache&&rows)
                cache->recycleReplacementScratch(std::move(*rows));
        }
    } replacementScratchGuard{&cache,&nextSlots};
    std::vector<FacePersistentId> ids;
    auto finishMissing=[&]() noexcept {
        try {
            for(const auto id:ids) {
                bool seen=false;
                for(const auto& row:out.processed)
                    if(row.facePersistentId==id)seen=true;
                for(const auto pending:out.unattempted)
                    if(pending==id)seen=true;
                if(!seen)out.unattempted.push_back(id);
            }
            std::sort(out.unattempted.begin(),out.unattempted.end());
        } catch(...) { out.reportComplete=false; }
    };
    try {
#if defined(PAACS_TOPO_TESTING)
        const t2c_detail::SurfaceCacheAllocationZoneScope allocationZone(
            t2c_detail::SurfaceCacheAllocationZone::DiagnosticEnvelope);
#endif
        out.processed.reserve(request.faces.size());
        out.unattempted.reserve(request.faces.size());
        ids.reserve(request.faces.size());
        for(const auto& face:request.faces)
            ids.push_back(face.faceStamp.persistentId);
        std::sort(ids.begin(),ids.end());
        out.reportComplete=true;
    } catch(...) {
        out.setStatus=SurfaceSetStatus::Failed;
        out.setFailure=allocation(SurfaceAllocationPhase::SetMapping);
        out.reportComplete=false;
        // The v4 coordinator prepares this fixed envelope before entering
        // the cached adapter.  Even when one of the adapter's own dynamic
        // report vectors cannot be reserved, retain a complete, allocation-
        // free cache report: every face was stopped before processing.
        bool cacheReportComplete=cacheReport!=nullptr&&
            cacheReport->processedStorage.size()>=request.faces.size();
        if(cacheReport) {
            cacheReport->processedUsed=0;
            cacheReport->unattemptedUsed=0;
            if(cacheReportComplete) {
                for(const auto& face:request.faces) {
                    cacheReport->processedStorage[
                        cacheReport->processedUsed++]={
                            face.faceStamp.persistentId,
                            SurfaceCacheFaceDisposition::Stopped,
                            SurfaceFaceStatus::Failed,out.setFailure};
                }
                std::sort(cacheReport->processedStorage.begin(),
                    cacheReport->processedStorage.begin()+
                        static_cast<std::ptrdiff_t>(
                            cacheReport->processedUsed),
                    [](const SurfaceCacheFaceReportRow& a,
                       const SurfaceCacheFaceReportRow& b) {
                        return a.facePersistentId<b.facePersistentId;
                    });
            }
        }
        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::AllocationFailed,
            SurfaceFaceCachePhase::Envelope,cacheReportComplete,0,
            cacheReportComplete ? request.faces.size() : 0,0);
        finishCacheReport(cacheReportComplete,0);
        return out;
    }
    if(!cache.acquireReplacementScratch(request.faces.size(),nextSlots)) {
        const auto diagnostic=cache.diagnostic();
        out.setStatus=SurfaceSetStatus::Failed;
        out.setFailure=diagnostic.reason==
                SurfaceFaceCacheReason::AllocationFailed
            ? allocation(SurfaceAllocationPhase::SetMapping)
            : failure(SurfaceReason::InternalError);
        bool reportOk=true;
        for(const auto id:ids) {
            out.processed.emplace_back(SurfaceFaceRow{
                id,SurfaceFaceStatus::Failed,out.setFailure,nullptr});
            reportOk=appendProcessed(id,
                SurfaceCacheFaceDisposition::Stopped,
                SurfaceFaceStatus::Failed,out.setFailure)&&reportOk;
        }
        finishCacheReport(reportOk,0);
        return out;
    }
    try {
        SetPlan plan;
        SurfaceFailure preflight=buildSetPlan(request,destination,plan);
        if(preflight.reason!=SurfaceReason::None) {
            const auto status=statusForFailure(preflight);
            bool reportOk=true;
            for(const auto id:ids) {
                out.processed.emplace_back(
                    SurfaceFaceRow{id,status,preflight,nullptr});
                reportOk=appendProcessed(id,
                    SurfaceCacheFaceDisposition::Stopped,status,preflight)&&
                    reportOk;
            }
            out.setStatus=setStatusFor(status); out.setFailure=preflight;
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::FaceBuildStopped,
                SurfaceFaceCachePhase::Preparation,reportOk,0,ids.size(),0);
            finishCacheReport(reportOk,0);
            return out;
        }
        // Dependency and invariant-authority preflight happens while the
        // production request is assembled, before this adapter is entered.
        // If it has already stopped the attempt, materialize the complete
        // face report now and do not touch either provider.  Treating this as
        // the first dirty-face failure would lose the whole-cache taxonomy
        // and could perform work after a hard cache budget was exhausted.
        const auto cachePreflight=cache.diagnostic();
        if(cachePreflight.status!=SurfaceFaceCacheStatus::None) {
            const SurfaceFailure cacheFailure=
                cachePreflight.reason==
                    SurfaceFaceCacheReason::AllocationFailed
                ? allocation(SurfaceAllocationPhase::SetMapping)
                : failure(SurfaceReason::InternalError);
            const auto status=statusForFailure(cacheFailure);
            bool reportOk=true;
            for(const auto id:ids) {
                out.processed.emplace_back(
                    SurfaceFaceRow{id,status,cacheFailure,nullptr});
                reportOk=appendProcessed(id,
                    SurfaceCacheFaceDisposition::Stopped,status,
                    cacheFailure)&&reportOk;
            }
            out.setStatus=setStatusFor(status);
            out.setFailure=cacheFailure;
            finishCacheReport(reportOk,cachePreflight.failedFace);
            return out;
        }
        CandidateProviderStamp providerStamp;
        try {
#if defined(PAACS_TOPO_TESTING)
            t2c_detail::noteSurfaceCacheProviderCallForTesting(
                t2c_detail::SurfaceCacheProviderCallKind::CandidateStamp);
#endif
            providerStamp=candidateProvider.stamp();
        }
        catch(const std::bad_alloc&) {
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=allocation(SurfaceAllocationPhase::Candidate);
            out.unattempted=ids;
            const bool reportOk=appendRemainingProcessed(
                0,ids,out.setFailure);
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AllocationFailed,
                SurfaceFaceCachePhase::Preparation,reportOk,0,ids.size(),0);
            finishCacheReport(reportOk,0);
            return out;
        } catch(...) {
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=failure(SurfaceReason::CandidateProviderError);
            out.unattempted=ids;
            const bool reportOk=appendRemainingProcessed(
                0,ids,out.setFailure);
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AuthorityChanged,
                SurfaceFaceCachePhase::Preparation,reportOk,0,ids.size(),0);
            finishCacheReport(reportOk,0);
            return out;
        }
        for(std::size_t fi=0;fi<plan.faces.size();++fi) {
            auto& fp=plan.faces[fi];
            const auto& face=*fp.request;
            SurfaceFailure err;
            std::optional<ChartCertificate> certificate;
            err=geometryCheckpoint(face,geometry,certificate,true);
            if(err.reason==SurfaceReason::None)
                err=normalizeFacePlan(fp,*certificate);
            if(err.reason==SurfaceReason::None&&
               !(providerStamp==face.expectedCandidateProviderStamp))
                err=failure(SurfaceReason::StaleCandidateProvider);

            const auto authority=authorities.find(face.faceStamp.persistentId);
            if(err.reason==SurfaceReason::None&&authority==authorities.end())
                err=failure(SurfaceReason::InternalError);
            std::unique_ptr<SurfaceFaceResult> result;
            std::shared_ptr<const CachedSurfaceFaceArtifact> entry;
            if(err.reason==SurfaceReason::None)
                entry=cache.find(face.faceStamp.persistentId,authority->second);
            const bool cleanHit=static_cast<bool>(entry);
            if(err.reason==SurfaceReason::None&&
               !cache.noteLookup(face.faceStamp.persistentId,!cleanHit))
                err=failure(SurfaceReason::InternalError);
            if(entry) {
                if(!cache.verifyAndReplaySizing(*entry,sizingEvaluator))
                    err=failure(SurfaceReason::InternalError);
                if(err.reason==SurfaceReason::None) {
                    try {
                        result=rehydrateSurfaceFaceFromCache(
                            *entry,plan.ledger,destination);
                    } catch(const std::bad_alloc&) {
                        err=allocation(SurfaceAllocationPhase::SetMapping);
                        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::AllocationFailed,
                            SurfaceFaceCachePhase::Rehydration,true,
                            face.faceStamp.persistentId,fi,
                            plan.faces.size()-fi-1);
                    } catch(...) {
                        err=failure(SurfaceReason::InternalError);
                        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::Invalid,
                            SurfaceFaceCachePhase::Rehydration,true,
                            face.faceStamp.persistentId,fi,
                            plan.faces.size()-fi-1);
                    }
                }
                if(err.reason==SurfaceReason::None&&
                   (!auditRehydratedSurfaceFace(*result,destination))) {
                    err=failure(SurfaceReason::LocalAuditFailed);
                    cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                        SurfaceFaceCacheReason::Invalid,
                        SurfaceFaceCachePhase::SetAudit,true,
                        face.faceStamp.persistentId,fi,
                        plan.faces.size()-fi-1);
                }
                if(err.reason==SurfaceReason::None)
                    err=verifyGeometryCheckpoint(face,geometry,*certificate);
                if(err.reason==SurfaceReason::None&&
                   !cache.noteReused(*entry))
                    err=cache.diagnostic().reason==
                            SurfaceFaceCacheReason::AllocationFailed
                        ? allocation(SurfaceAllocationPhase::SetMapping)
                        : failure(SurfaceReason::InternalError);
            } else if(err.reason==SurfaceReason::None) {
                SurfaceCandidate candidate;
                try {
#if defined(PAACS_TOPO_TESTING)
                    t2c_detail::noteSurfaceCacheProviderCallForTesting(
                        t2c_detail::SurfaceCacheProviderCallKind::
                            CandidateBuild,face.faceStamp.persistentId);
#endif
                    candidate=candidateProvider.build(face);
                }
                catch(const std::bad_alloc&) {
                    err=allocation(SurfaceAllocationPhase::Candidate);
                } catch(...) {
                    err=failure(SurfaceReason::CandidateProviderError);
                }
                if(err.reason==SurfaceReason::None)
                    err=geometryCheckpoint(face,geometry,certificate,false);
                if(err.reason==SurfaceReason::None)
                    err=validateCandidate(face,fp,candidate,*certificate);
                FaceBuild built;
                if(err.reason==SurfaceReason::None)
                    built=buildFace(fp,candidate,destination,plan.ledger,
                                    geometry,localSteinerPolicy,*certificate);
                if(err.reason==SurfaceReason::None)err=built.error;
                if(err.reason==SurfaceReason::None) {
                    ProductionSurfaceFaceMetrics metric{};
                    metric.facePersistentId=face.faceStamp.persistentId;
                    (void)metricsSource.metricForFace(
                        face.faceStamp.persistentId,metric);
                    try {
                        entry=snapshotSurfaceFaceForCache(
                            *built.result,authority->second,
                            cache.observedSizing(face.faceStamp.persistentId),
                            metric,plan.ledger);
                        result=std::move(built.result);
                        if(!cache.noteBuilt(*entry))
                            err=cache.diagnostic().reason==
                                    SurfaceFaceCacheReason::AllocationFailed
                                ? allocation(
                                    SurfaceAllocationPhase::SetMapping)
                                : failure(SurfaceReason::InternalError);
                    } catch(const std::bad_alloc&) {
                        err=allocation(SurfaceAllocationPhase::SetMapping);
                        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::AllocationFailed,
                            SurfaceFaceCachePhase::Preparation,true,
                            face.faceStamp.persistentId,fi,
                            plan.faces.size()-fi-1);
                    } catch(...) {
                        err=failure(SurfaceReason::InternalError);
                        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                            SurfaceFaceCacheReason::Invalid,
                            SurfaceFaceCachePhase::Preparation,true,
                            face.faceStamp.persistentId,fi,
                            plan.faces.size()-fi-1);
                    }
                }
            }
            if(err.reason!=SurfaceReason::None) {
                const auto status=statusForFailure(err);
                out.processed.emplace_back(SurfaceFaceRow{
                    face.faceStamp.persistentId,status,err,nullptr});
                bool reportOk=appendProcessed(face.faceStamp.persistentId,
                    cleanHit ? SurfaceCacheFaceDisposition::CleanReplay
                             : SurfaceCacheFaceDisposition::DirtyBuild,
                    status,err);
                const bool replayFailure=cleanHit||
                    cache.diagnostic().phase==
                        SurfaceFaceCachePhase::SetAudit;
                for(std::size_t later=fi+1;later<plan.faces.size();++later)
                    out.unattempted.push_back(
                        plan.faces[later].request->faceStamp.persistentId);
                if(replayFailure) {
                    for(std::size_t later=fi+1;later<plan.faces.size();++later)
                        reportOk=appendProcessed(
                            plan.faces[later].request->faceStamp.persistentId,
                            SurfaceCacheFaceDisposition::Stopped,status,err)&&
                            reportOk;
                } else {
                    for(std::size_t later=fi+1;later<plan.faces.size();++later)
                        reportOk=appendUnattempted(
                            plan.faces[later].request->faceStamp.persistentId)&&
                            reportOk;
                }
                out.setStatus=setStatusFor(status); out.setFailure=err;
                if(cache.diagnostic().status==SurfaceFaceCacheStatus::None)
                    cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                        SurfaceFaceCacheReason::FaceBuildStopped,
                        SurfaceFaceCachePhase::FaceBuild,reportOk,
                        face.faceStamp.persistentId,
                        replayFailure ? plan.faces.size() : fi+1,
                        replayFailure ? 0 : plan.faces.size()-fi-1);
                finishCacheReport(reportOk,face.faceStamp.persistentId);
                return out;
            }
            {
#if defined(PAACS_TOPO_TESTING)
                const t2c_detail::SurfaceCacheAllocationZoneScope
                    mappingAllocationZone(
                        t2c_detail::SurfaceCacheAllocationZone::
                            LocalGlobalMapping);
#endif
                nextSlots.push_back(std::move(entry));
                out.processed.emplace_back(SurfaceFaceRow{
                    face.faceStamp.persistentId,SurfaceFaceStatus::Ready,{},
                    std::move(result)});
            }
            if(!appendProcessed(face.faceStamp.persistentId,
                    cleanHit ? SurfaceCacheFaceDisposition::CleanReplay
                             : SurfaceCacheFaceDisposition::DirtyBuild,
                    SurfaceFaceStatus::Ready,{})) {
                out.setStatus=SurfaceSetStatus::Failed;
                out.setFailure=failure(SurfaceReason::InternalError);
                cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                    SurfaceFaceCacheReason::Invalid,
                    SurfaceFaceCachePhase::Preparation,false,
                    face.faceStamp.persistentId,fi+1,0);
                finishCacheReport(false,face.faceStamp.persistentId);
                return out;
            }
        }

        if(!cache.verifyFinalSizing(nextSlots,sizingEvaluator)) {
            const auto failedFace=cache.diagnostic().failedFace;
            const SurfaceFailure sizingFailure=
                cache.diagnostic().reason==
                    SurfaceFaceCacheReason::AllocationFailed
                ? allocation(SurfaceAllocationPhase::SetMapping)
                : failure(SurfaceReason::InternalError);
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=sizingFailure;
            if(cacheReport) {
                for(std::size_t i=0;i<cacheReport->processedUsed;++i) {
                    auto& row=cacheReport->processedStorage[i];
                    if(failedFace!=0&&row.facePersistentId==failedFace) {
                        row.surfaceStatus=SurfaceFaceStatus::Failed;
                        row.originalSurfaceFailure=sizingFailure;
                    }
                }
            }
            finishCacheReport(true,failedFace);
            return out;
        }
        if(!cache.noteSetAudit(static_cast<std::uint64_t>(
                request.faces.size()+nextSlots.size()))) {
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=cache.diagnostic().reason==
                    SurfaceFaceCacheReason::AllocationFailed
                ? allocation(SurfaceAllocationPhase::SetMapping)
                : failure(SurfaceReason::InternalError);
            finishCacheReport(cacheReport&&
                cacheReport->processedUsed==request.faces.size(),0);
            return out;
        }
        if(!cache.prepareReplacement(std::move(nextSlots),ids)) {
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=cache.diagnostic().reason==
                    SurfaceFaceCacheReason::AllocationFailed
                ? allocation(SurfaceAllocationPhase::SetMapping)
                : failure(SurfaceReason::InternalError);
            finishCacheReport(cacheReport&&
                cacheReport->processedUsed==request.faces.size(),0);
            return out;
        }

        // The candidate-provider authority is checked only after the final
        // sizing replay, complete set audit, and cache replacement have all
        // been prepared.  It is therefore the last fallible authority read
        // before publication preparation.  A failure abandons every cache
        // placeholder before any publication object is constructed.
        CandidateProviderStamp finalProviderStamp{};
        try {
#if defined(PAACS_TOPO_TESTING)
            t2c_detail::noteSurfaceCacheProviderCallForTesting(
                t2c_detail::SurfaceCacheProviderCallKind::CandidateStamp);
#endif
            finalProviderStamp=candidateProvider.stamp();
        } catch(const std::bad_alloc&) {
            cache.abandonPrepared();
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=allocation(SurfaceAllocationPhase::Candidate);
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AllocationFailed,
                SurfaceFaceCachePhase::SetAudit,true,0,
                request.faces.size(),0);
            finishCacheReport(true,0);
            return out;
        } catch(...) {
            cache.abandonPrepared();
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=failure(SurfaceReason::CandidateProviderError);
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AuthorityChanged,
                SurfaceFaceCachePhase::SetAudit,true,0,
                request.faces.size(),0);
            finishCacheReport(true,0);
            return out;
        }
        if(!(finalProviderStamp==providerStamp)) {
            cache.abandonPrepared();
            out.setStatus=SurfaceSetStatus::Failed;
            out.setFailure=failure(SurfaceReason::StaleCandidateProvider);
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                SurfaceFaceCacheReason::AuthorityChanged,
                SurfaceFaceCachePhase::SetAudit,true,0,
                request.faces.size(),0);
            finishCacheReport(true,0);
            return out;
        }
        auto artifact=std::make_unique<SurfaceSetResult>();
        artifact->setStamp=request.setStamp;
        artifact->faces.reserve(out.processed.size());
        artifact->sharedLedger=std::move(plan.ledger);
        for(auto& row:out.processed)
            artifact->faces.push_back(std::move(row.readyResult));
        bool published=false;
        const auto publication=publishSurfaceSet(
            destination,*artifact,geometry,published);
        out.retainedArtifact=std::move(artifact);
        if(!published) {
            cache.abandonPrepared();
            out.setStatus=setStatusFor(statusForFailure(publication));
            out.setFailure=publication;
            cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
                publication.reason==SurfaceReason::AllocationFailed
                    ? SurfaceFaceCacheReason::AllocationFailed
                    : SurfaceFaceCacheReason::Invalid,
                SurfaceFaceCachePhase::SetAudit,true,0,ids.size(),0);
            finishCacheReport(true,0);
            return out;
        }
        cache.commitPrepared();
        cache.setDiagnostic(SurfaceFaceCacheStatus::Ready,
            SurfaceFaceCacheReason::None,SurfaceFaceCachePhase::Commit,
            true,0,ids.size(),0);
        finishCacheReport(true,0);
        out.setStatus=SurfaceSetStatus::Published; out.setFailure={};
        return out;
    } catch(const std::bad_alloc&) {
        cache.abandonPrepared();
        out.setStatus=SurfaceSetStatus::Failed;
        out.setFailure=allocation(SurfaceAllocationPhase::SetMapping);
        finishMissing();
        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::AllocationFailed,
            SurfaceFaceCachePhase::Preparation,out.reportComplete,0,
            out.processed.size(),out.unattempted.size());
        const bool reportOk=appendRemainingProcessed(
            cacheReport ? cacheReport->processedUsed : 0,ids,out.setFailure);
        finishCacheReport(reportOk,0);
    } catch(...) {
        cache.abandonPrepared();
        out.setStatus=SurfaceSetStatus::Failed;
        out.setFailure=failure(SurfaceReason::InternalError);
        finishMissing();
        cache.setDiagnostic(SurfaceFaceCacheStatus::Failed,
            SurfaceFaceCacheReason::Invalid,
            SurfaceFaceCachePhase::Preparation,out.reportComplete,0,
            out.processed.size(),out.unattempted.size());
        const bool reportOk=appendRemainingProcessed(
            cacheReport ? cacheReport->processedUsed : 0,ids,out.setFailure);
        finishCacheReport(reportOk,0);
    }
    return out;
}

namespace detail {

CandidateBoundaryRecoveryResult recoverCandidateBoundaryForCk2(
    const SurfaceFaceRequest& request, SurfaceCandidate& candidate,
    const SurfaceGeometryProvider& geometry,
    bool allowDirectConstraintInitialization) {
    detail::CandidateBoundaryRecoveryResult result;
    auto prepared = std::make_shared<detail::PreparedSurfaceCandidate>();
    prepared->providerStamp = candidate.providerStamp;
    prepared->faceStamp = candidate.faceStamp;
    prepared->chart = candidate.chart;
    prepared->geometryStamp = request.expectedGeometryStamp;

    topo::NodePool nodes;
    std::vector<topo::NodeHandle> handles;
    handles.reserve(candidate.nodes.size());
    for (const auto& node : candidate.nodes)
        handles.push_back(nodes.add(node.node.coord, node.node.cls));
    topo::CellPool cells;
    for (const auto& cell : candidate.cells) {
        const std::array<topo::NodeHandle, 3> cn{
            handles[cell.nodes[0]], handles[cell.nodes[1]],
            handles[cell.nodes[2]]};
        cells.add(topo::ElementType::Tri3,
                  std::span<const topo::NodeHandle>(cn.data(), cn.size()),
                  cell.classification, nodes);
    }
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        edges.insertCellEdges(nodes, cells, cell);

    bnd::ConstraintSet constraints;
    for (std::size_t i = 0; i < handles.size(); ++i) {
        constraints.params.emplace(handles[i], candidate.nodes[i].parameter);
        bnd::AuthoritativeParamRecord auth;
        auth.affine = false;
        auth.storedUv = candidate.nodes[i].parameter.p;
        auth.chart = request.chart;
        constraints.authParams.emplace(handles[i], auth);
    }
    for (const auto& [id, index] : candidate.declaredVertexMap)
        constraints.vertices.emplace(id, handles[index]);
    auto addChain = [&](const SurfaceChain& chain) {
        const std::size_t count = chain.closed ? chain.uses.size()
                                               : chain.uses.size() - 1;
        for (std::size_t i = 0; i < count; ++i)
            constraints.segments.push_back(
                {chain.segmentEntityIds[i],
                 chain.uses[i].faceLocalVertexId,
                 chain.uses[(i + 1) % chain.uses.size()].faceLocalVertexId,
                 request.chart, chain.classification});
    };
    for (const auto& chain : request.outer.pieces) addChain(chain);
    for (const auto& hole : request.holes)
        for (const auto& chain : hole.pieces) addChain(chain);
    for (const auto& chain : request.embedded)
        if (chain.owningCurve) addChain(chain);
    std::sort(constraints.segments.begin(), constraints.segments.end(),
              [](const bnd::ConstraintSegment& a,
                 const bnd::ConstraintSegment& b) { return a.id < b.id; });
    struct DirectSegmentRow {
        bnd::ConstraintEntityId entity = 0;
        topo::NodeHandle a{};
        topo::NodeHandle b{};
        topo::EdgeHandle edge{};
        mesh::Classification classification{};
        bool alreadyMarked = false;
    };
    std::vector<DirectSegmentRow> directRows;
    directRows.reserve(constraints.segments.size());
    bool directCandidate = allowDirectConstraintInitialization;
    if (directCandidate)
        for (const auto& segment : constraints.segments) {
            const auto ia = constraints.vertices.find(segment.a);
            const auto ib = constraints.vertices.find(segment.b);
            if (ia == constraints.vertices.end() ||
                ib == constraints.vertices.end()) {
                directCandidate = false;
                break;
            }
            const auto edge = edges.find(nodes, ia->second, ib->second);
            if (!edge) {
                directCandidate = false;
                break;
            }
            const auto data = edges.get(*edge);
            const mesh::Classification defaultClass{};
            const bool exactClass = sameClass(data.cls, segment.cls);
            const bool defaultClassified = sameClass(data.cls, defaultClass);
            if ((data.requiredConstraint && !exactClass) ||
                (!data.requiredConstraint && !defaultClassified &&
                 !exactClass)) {
                directCandidate = false;
                break;
            }
            directRows.push_back({segment.id, ia->second, ib->second, *edge,
                                  segment.cls,
                                  data.requiredConstraint && exactClass});
        }
    if (directCandidate) {
        for (const auto& row : directRows) {
            if (row.alreadyMarked) continue;
            edges.setClassification(row.edge, row.classification);
            edges.setRequiredConstraint(row.edge, true);
        }
    }

    auto workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes), std::move(cells), std::move(edges),
        topo::FacePool{});
    constraints.stamp = bnd::OwnershipStamp::of(*workspace);

    // BND's qualified five-outcome set wrapper deliberately folds a null
    // session index into PslgRejected. The v3 contract is stricter: every
    // allocation failure must remain allocation-attributed. Establish the
    // immutable input audits and complete index before entering that wrapper.
    if (!bnd::preflightConstraintSet(*workspace, constraints).empty() ||
        !bnd::preflightSourceEmbedding(*workspace, constraints).empty() ||
        !bnd::preflightPslg(constraints).empty()) {
        result.reject = bnd::SteinerReject::PslgRejected;
        return result;
    }

    bnd::SteinerSegmentsReport report;
    if (directCandidate) {
        // The candidate topology was constructed with every requested
        // segment as a literal edge.  Marking those edges while the transient
        // workspace is being assembled avoids opening a general recovery
        // session, building its adjacency index and preparing one mark-only
        // cavity transaction per segment.  No authoritative topology exists
        // yet; the fully stamped receipt is still revalidated by the adapter.
        report.processed.reserve(directRows.size());
        for (const auto& row : directRows) {
            bnd::SteinerResult recovered;
            recovered.status = row.alreadyMarked
                ? bnd::SteinerStatus::AlreadyPresent
                : bnd::SteinerStatus::Recovered;
            recovered.recoveredChain = {row.a, row.b};
            report.processed.push_back({row.entity, std::move(recovered)});
        }
    } else {
        auto preflightIndex = bnd::buildTriAdjacency(*workspace);
        if (!preflightIndex) {
            result.reject = bnd::SteinerReject::PslgRejected;
            return result;
        }
        preflightIndex.reset();

        bnd::SteinerPolicy policy;
        policy.steinerEnabled = false;
        report = bnd::recoverSteinerSegments(
            *workspace, constraints, geometry, request.expectedGeometryStamp,
            policy);
    }
    // A normal geometric refusal also leaves an unattempted suffix. Preserve
    // that first failing result before diagnosing an incomplete report.
    for (const auto& row : report.processed) {
        if (row.result.status != bnd::SteinerStatus::AlreadyPresent &&
            row.result.status != bnd::SteinerStatus::Recovered) {
            if (row.result.reject == bnd::SteinerReject::PslgRejected)
                result.fail = bnd::SteinerFail::CandidateAllocFailed;
            else
                result.reject = row.result.reject;
            result.refuse = row.result.refuse;
            if (result.fail == bnd::SteinerFail::None)
                result.fail = row.result.fail;
            return result;
        }
        if (!row.result.createdParams.empty()) {
            result.refuse = bnd::SteinerRefuse::SteinerDisabled;
            return result;
        }
        ++result.recoveredSegments;
    }
    if (!report.unattempted.empty() ||
        report.processed.size() != constraints.segments.size()) {
        result.fail = bnd::SteinerFail::CandidateAllocFailed;
        return result;
    }

    std::map<topo::NodeHandle, std::size_t> localIndex;
    for (std::size_t i = 0; i < handles.size(); ++i)
        localIndex.emplace(handles[i], i);
    const auto activeCells = workspace->cells().activeBySequence();
    std::vector<SurfaceCandidateCell> recoveredCells;
    recoveredCells.reserve(activeCells.size());
    for (const auto ch : activeCells) {
        const auto cell = workspace->cells().get(ch);
        recoveredCells.push_back(
            {{{localIndex.at(cell.nodes[0]), localIndex.at(cell.nodes[1]),
               localIndex.at(cell.nodes[2])}}, cell.cls});
    }
    candidate.cells.swap(recoveredCells);

    constraints.stamp = bnd::OwnershipStamp::of(*workspace);
    prepared->workspaceStamp = constraints.stamp;
    prepared->candidateNodes = handles;
    prepared->candidateCells = activeCells;
    prepared->constraints = std::move(constraints);
    prepared->recovery = std::move(report);
    prepared->workspace = std::move(workspace);
    candidate.preparedRecovery = std::move(prepared);
    result.complete = true;
    return result;
}

}  // namespace detail

#if defined(PAACS_TOPO_TESTING)
CandidateBoundaryRecoveryTestingOutcome recoverCandidateBoundaryForTesting(
    const SurfaceFaceRequest& request, SurfaceCandidate& candidate,
    const SurfaceGeometryProvider& geometry) {
    const auto recovered = detail::recoverCandidateBoundaryForCk2(
        request,candidate,geometry);
    return {recovered.complete,recovered.reject,recovered.refuse,
            recovered.fail};
}
#endif

}  // namespace paacs::algo
