// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/FacetConflictResolver.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <limits>
#include <new>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace paacs::algo {
namespace {

enum class WorkState : std::uint8_t { Ready, Budget, Overflow };

FacetConflictResolveResult failure(FacetConflictResolveStatus status,
                                   FacetConflictResolveReason reason,
                                   std::uint64_t workUnits = 0) noexcept {
    FacetConflictResolveResult result;
    result.status = status;
    result.reason = reason;
    result.workUnits = workUnits;
    return result;
}

WorkState consumeWork(std::uint64_t& workUnits,
                      std::uint64_t maxWorkUnits) noexcept {
    if (workUnits == std::numeric_limits<std::uint64_t>::max())
        return WorkState::Overflow;
    if (workUnits >= maxWorkUnits) return WorkState::Budget;
    ++workUnits;
    return WorkState::Ready;
}

FacetConflictResolveResult workFailure(WorkState state,
                                       std::uint64_t workUnits) noexcept {
    return failure(
        FacetConflictResolveStatus::Unsupported,
        state == WorkState::Budget
            ? FacetConflictResolveReason::WorkBudgetReached
            : FacetConflictResolveReason::WorkAccountingOverflow,
        workUnits);
}

bool validProviderStamp(const bnd::GeometryProviderStamp& stamp) noexcept {
    return stamp.identity != 0 && stamp.revision != 0;
}

bool validArtifactStamp(const ArtifactStamp& stamp) noexcept {
    return stamp.persistentId != 0 && stamp.revision != 0;
}

std::uint64_t canonicalParameterBits(double value) noexcept {
    return value == 0.0 ? UINT64_C(0) :
                          std::bit_cast<std::uint64_t>(value);
}

bool sameParameterBits(const Vec2& a, const Vec2& b) noexcept {
    return canonicalParameterBits(a.u) == canonicalParameterBits(b.u) &&
           canonicalParameterBits(a.v) == canonicalParameterBits(b.v);
}

bool compatibleAuthority(const FacetConflictResolveResult& a,
                         const FacetConflictResolveResult& b) noexcept {
    if (!(a.faceStamp == b.faceStamp) || a.faceTag != b.faceTag ||
        !(a.chartCertificate.providerStamp ==
          b.chartCertificate.providerStamp) ||
        !(a.chartCertificate.certificateStamp ==
          b.chartCertificate.certificateStamp))
        return false;
    return true;
}

bool canonicalMappingLess(const FacetConflictResolveResult& a,
                          const FacetConflictResolveResult& b) noexcept {
    return std::tie(a.chart,
                    a.chartCertificate.certificateStamp.persistentId,
                    a.chartCertificate.certificateStamp.revision) <
           std::tie(b.chart,
                    b.chartCertificate.certificateStamp.persistentId,
                    b.chartCertificate.certificateStamp.revision);
}

FacetConflictResolveResult resolveFace(
    const SurfaceFaceResult& face,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits,
    std::uint64_t& workUnits) noexcept {
    WorkState state = consumeWork(workUnits, maxWorkUnits);
    if (state != WorkState::Ready) return workFailure(state, workUnits);
    if (!validArtifactStamp(face.faceStamp) ||
        face.faceStamp.persistentId != key.facePersistentId)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::InvalidFaceIdentity,
                       workUnits);
    if (face.faceClass.dim != 2 || face.faceClass.tag < 0)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::InvalidFaceClassification,
                       workUnits);
    if (face.localNodes.size() != face.localMappingParams.size() ||
        face.localNodes.size() != face.publishedGlobalByLocalNode.size() ||
        face.localCells.size() !=
            face.publishedGlobalCellByLocalCell.size())
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::IncompleteFacePayload,
                       workUnits);

    state = consumeWork(workUnits, maxWorkUnits);
    if (state != WorkState::Ready) return workFailure(state, workUnits);
    if (!validProviderStamp(face.expectedGeometryStamp) ||
        !validProviderStamp(face.chartCertificate.providerStamp) ||
        !validArtifactStamp(face.chartCertificate.certificateStamp) ||
        !(face.expectedGeometryStamp ==
          face.chartCertificate.providerStamp) ||
        face.chartCertificate.chart != face.chart)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::InvalidChartCertificate,
                       workUnits);

    if (scratch.size() < face.localCells.size())
        return failure(FacetConflictResolveStatus::Unsupported,
                       FacetConflictResolveReason::ScratchCapacityReached,
                       workUnits);

    const std::size_t rowCount = face.localCells.size();
    for (std::size_t ci = 0; ci < rowCount; ++ci) {
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        const auto& cell = face.localCells[ci];
        if (cell.classification.dim != face.faceClass.dim ||
            cell.classification.tag != face.faceClass.tag)
            return failure(FacetConflictResolveStatus::Invalid,
                           FacetConflictResolveReason::InvalidLocalCell,
                           workUnits);
        std::array<Index, 3> published{};
        for (std::size_t vi = 0; vi < published.size(); ++vi) {
            state = consumeWork(workUnits, maxWorkUnits);
            if (state != WorkState::Ready)
                return workFailure(state, workUnits);
            const std::size_t local = cell.nodes[vi];
            if (local >= face.publishedGlobalByLocalNode.size())
                return failure(FacetConflictResolveStatus::Invalid,
                               FacetConflictResolveReason::InvalidLocalCell,
                               workUnits);
            published[vi] = face.publishedGlobalByLocalNode[local];
            if (published[vi] < 0)
                return failure(
                    FacetConflictResolveStatus::Invalid,
                    FacetConflictResolveReason::InvalidPublishedNode,
                    workUnits);
        }
        for (const auto pair :
             std::array<std::array<std::size_t, 2>, 3>{{
                 {{0, 1}}, {{1, 2}}, {{0, 1}}}}) {
            state = consumeWork(workUnits, maxWorkUnits);
            if (state != WorkState::Ready)
                return workFailure(state, workUnits);
            if (published[pair[1]] < published[pair[0]])
                std::swap(published[pair[0]], published[pair[1]]);
        }
        if (published[0] == published[1] || published[1] == published[2])
            return failure(FacetConflictResolveStatus::Invalid,
                           FacetConflictResolveReason::InvalidLocalCell,
                           workUnits);
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        scratch[ci] = {published, ci};
    }

    for (std::size_t i = 0; i < rowCount; ++i) {
        std::size_t smallest = i;
        for (std::size_t j = i + 1; j < rowCount; ++j) {
            state = consumeWork(workUnits, maxWorkUnits);
            if (state != WorkState::Ready)
                return workFailure(state, workUnits);
            if (scratch[j].canonicalGlobalNodes <
                scratch[smallest].canonicalGlobalNodes)
                smallest = j;
        }
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        std::swap(scratch[i], scratch[smallest]);
    }
    for (std::size_t i = 1; i < rowCount; ++i) {
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        if (scratch[i - 1].canonicalGlobalNodes ==
            scratch[i].canonicalGlobalNodes)
            return failure(
                FacetConflictResolveStatus::Unsupported,
                FacetConflictResolveReason::DuplicateCanonicalFacet,
                workUnits);
    }
    state = consumeWork(workUnits, maxWorkUnits);
    if (state != WorkState::Ready) return workFailure(state, workUnits);
    if (key.normalizedLocalCellOrdinal >= rowCount)
        return failure(FacetConflictResolveStatus::Unsupported,
                       FacetConflictResolveReason::FacetOrdinalOutOfRange,
                       workUnits);

    const auto& selected = scratch[key.normalizedLocalCellOrdinal];
    const auto& cell = face.localCells[selected.localCell];
    FacetConflictResolveResult result;
    result.faceStamp = face.faceStamp;
    result.faceTag = face.faceClass.tag;
    result.chart = face.chart;
    result.chartCertificate = {
        face.chartCertificate.providerStamp,
        face.chartCertificate.chart,
        face.chartCertificate.certificateStamp,
        face.chartCertificate.certifiedRegionBoundary,
        face.chartCertificate.chartPositiveToNormalSign,
        face.chartCertificate.regular,
        face.chartCertificate.periodic,
        face.chartCertificate.seam};
    for (std::size_t sourceIndex = 0; sourceIndex < cell.nodes.size();
         ++sourceIndex) {
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        const std::size_t local = cell.nodes[sourceIndex];
        if (local >= face.publishedGlobalByLocalNode.size() ||
            local >= face.localMappingParams.size())
            return failure(FacetConflictResolveStatus::Invalid,
                           FacetConflictResolveReason::InvalidLocalCell,
                           workUnits);
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        result.vertices[sourceIndex] = {
            face.publishedGlobalByLocalNode[local], face.chart,
            face.localMappingParams[local]};
    }
    std::array<FacetBarycentreVertex, 3> canonicalVertices{};
    for (std::size_t sourceIndex = 0; sourceIndex < result.vertices.size();
         ++sourceIndex) {
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        const auto id = result.vertices[sourceIndex].publishedGlobalNode;
        std::size_t position = selected.canonicalGlobalNodes.size();
        for (std::size_t candidate = 0;
             candidate < selected.canonicalGlobalNodes.size(); ++candidate) {
            state = consumeWork(workUnits, maxWorkUnits);
            if (state != WorkState::Ready)
                return workFailure(state, workUnits);
            if (selected.canonicalGlobalNodes[candidate] == id) {
                position = candidate;
                break;
            }
        }
        if (position == selected.canonicalGlobalNodes.size())
            return failure(FacetConflictResolveStatus::Invalid,
                           FacetConflictResolveReason::InvalidLocalCell,
                           workUnits);
        state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        canonicalVertices[position] = result.vertices[sourceIndex];
    }
    state = consumeWork(workUnits, maxWorkUnits);
    if (state != WorkState::Ready) return workFailure(state, workUnits);
    result.vertices = canonicalVertices;
    result.status = FacetConflictResolveStatus::Resolved;
    result.reason = FacetConflictResolveReason::None;
    result.workUnits = workUnits;
    return result;
}

FacetConflictResolveResult resolveImpl(
    const SurfaceSetResult& retained,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits) noexcept {
    std::uint64_t workUnits = 0;
    if (key.facePersistentId == 0)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::InvalidFaceIdentity);

    FacetConflictResolveResult selected;
    bool found = false;
    for (const auto& owner : retained.faces) {
        const WorkState state = consumeWork(workUnits, maxWorkUnits);
        if (state != WorkState::Ready) return workFailure(state, workUnits);
        if (!owner || owner->faceStamp.persistentId != key.facePersistentId)
            continue;
        FacetConflictResolveResult candidate = resolveFace(
            *owner, key, scratch, maxWorkUnits, workUnits);
        if (!candidate.resolved()) return candidate;
        if (!found) {
            const WorkState selectState = consumeWork(
                workUnits, maxWorkUnits);
            if (selectState != WorkState::Ready)
                return workFailure(selectState, workUnits);
            selected = std::move(candidate);
            found = true;
            continue;
        }
        WorkState compareState = consumeWork(workUnits, maxWorkUnits);
        if (compareState != WorkState::Ready)
            return workFailure(compareState, workUnits);
        if (!compatibleAuthority(selected, candidate))
            return failure(FacetConflictResolveStatus::Unsupported,
                           FacetConflictResolveReason::AmbiguousFace,
                           workUnits);
        for (std::size_t i = 0; i < selected.vertices.size(); ++i) {
            compareState = consumeWork(workUnits, maxWorkUnits);
            if (compareState != WorkState::Ready)
                return workFailure(compareState, workUnits);
            if (selected.vertices[i].publishedGlobalNode !=
                    candidate.vertices[i].publishedGlobalNode ||
                !sameParameterBits(selected.vertices[i].parameter,
                                   candidate.vertices[i].parameter))
                return failure(FacetConflictResolveStatus::Unsupported,
                               FacetConflictResolveReason::AmbiguousFace,
                               workUnits);
        }
        if (canonicalMappingLess(candidate, selected)) {
            compareState = consumeWork(workUnits, maxWorkUnits);
            if (compareState != WorkState::Ready)
                return workFailure(compareState, workUnits);
            selected = std::move(candidate);
        }
    }
    if (!found)
        return failure(FacetConflictResolveStatus::Unsupported,
                       FacetConflictResolveReason::FaceNotFound, workUnits);
    selected.workUnits = workUnits;
    return selected;
}

}  // namespace

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetResult& retained,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits) noexcept {
    try {
        return resolveImpl(retained, key, scratch, maxWorkUnits);
    } catch (const std::bad_alloc&) {
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::AllocationFailed);
    } catch (...) {
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::UnexpectedFailure);
    }
}

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetResult& retained,
    const bnd::VolumeFacetKey& key) noexcept {
    std::size_t required = 0;
    for (const auto& owner : retained.faces)
        if (owner && owner->faceStamp.persistentId == key.facePersistentId)
            required = (std::max)(required, owner->localCells.size());
    try {
        std::vector<FacetConflictRankedFacetScratch> scratch(required);
        return resolveFacetConflict(
            retained, key, scratch,
            std::numeric_limits<std::uint64_t>::max());
    } catch (const std::bad_alloc&) {
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::AllocationFailed);
    } catch (...) {
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::UnexpectedFailure);
    }
}

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetOutcome& outcome,
    const bnd::VolumeFacetKey& key,
    std::span<FacetConflictRankedFacetScratch> scratch,
    std::uint64_t maxWorkUnits) noexcept {
    if (!outcome.retainedArtifact)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::MissingRetainedArtifact);
    return resolveFacetConflict(
        *outcome.retainedArtifact, key, scratch, maxWorkUnits);
}

FacetConflictResolveResult resolveFacetConflict(
    const SurfaceSetOutcome& outcome,
    const bnd::VolumeFacetKey& key) noexcept {
    if (!outcome.retainedArtifact)
        return failure(FacetConflictResolveStatus::Invalid,
                       FacetConflictResolveReason::MissingRetainedArtifact);
    return resolveFacetConflict(*outcome.retainedArtifact, key);
}

}  // namespace paacs::algo
