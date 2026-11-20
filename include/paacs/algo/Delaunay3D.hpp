// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include "paacs/core/Types.hpp"

namespace paacs::algo {

struct Tet { std::array<Index, 4> v; };

struct Delaunay3DStats {
    std::uint64_t inputPoints = 0;
    std::uint64_t pointsInserted = 0;
    std::uint64_t pointsDuplicate = 0;
    std::uint64_t pointsSkipped = 0;
    std::uint64_t walkFallbackScans = 0;
    std::uint64_t cavityTrimResurrections = 0;
    std::uint64_t bootstrapDistinctTests = 0;
    std::uint64_t bootstrapOrientationTests = 0;
    std::uint64_t walkCellVisits = 0;
    std::uint64_t walkOrientationTests = 0;
    std::uint64_t conflictQueries = 0;
    std::uint64_t conflictOrientationTests = 0;
    std::uint64_t conflictInsphereTests = 0;
    std::uint64_t fallbackCellTests = 0;
    std::uint64_t cavityCellVisits = 0;
    std::uint64_t cavityNeighbourVisits = 0;
    std::uint64_t trimCellVisits = 0;
    std::uint64_t trimOrientationTests = 0;
    std::uint64_t shellFacesEmitted = 0;
    std::uint64_t shellEdgePairInsertions = 0;
    std::uint64_t shellEdgePairProbes = 0;
    std::uint64_t peakShellPairTableSlots = 0;
    std::uint64_t adjacencyAuditFaceVisits = 0;
    std::uint64_t peakCavityCells = 0;
    std::uint64_t peakShellFaces = 0;
    std::uint64_t outputCells = 0;
    bool auditPassed = false;
    bool runComplete = false;
    std::uint64_t arenaPeakLiveCells = 0;
    std::uint64_t arenaMaterializedSlots = 0;
    std::uint64_t arenaBytes = 0;
    std::uint64_t arenaHotBytes = 0;
    std::uint64_t arenaColdBytes = 0;
    std::uint64_t arenaTotalCellsCreated = 0;
    std::uint64_t outputBytes = 0;
    std::uint64_t scratchBytes = 0;
    std::uint64_t workspaceRetainedBytes = 0;
    std::uint64_t directShellTets = 0;
    std::uint64_t genericShellTets = 0;
    std::uint64_t reciprocalFaceScansAvoided = 0;
    std::uint64_t finiteCavityCertificates = 0;
    std::uint64_t finiteCavityTrimTestsAvoided = 0;
    bool duplicateTetDetected = false;
    bool usedCompactArena = false;
    bool usedStableRadixOrder = false;
    bool usedLinearFinalization = false;
    bool workspaceReused = false;
    bool compactEnvelopeFallback = false;
    bool collectedDetailedOperationCounts = false;
    double normalizationOrderSeconds = 0.0;
    double insertionSeconds = 0.0;
    double adjacencyAuditSeconds = 0.0;
    double finalizationSeconds = 0.0;
    double totalSeconds = 0.0;

    // V0.9.3 deterministic-parallel accounting.
    bool usedParallelDiscovery = false;
    bool parallelDegeneracyFallback = false;
    int requestedThreads = 1;
    std::uint32_t packetSize = 0;
    std::uint64_t speculativeCandidates = 0;
    std::uint64_t speculativeDiscoveries = 0;
    std::uint64_t parallelRounds = 0;
    std::uint64_t speculativeCandidatesReused = 0;
    std::uint64_t speculativeCandidatesRecomputed = 0;
    std::uint64_t speculativePredicateCalls = 0;
    std::uint64_t serialEquivalentPredicateCalls = 0;
    std::uint64_t speculativeCavityVisits = 0;
    std::uint64_t serialEquivalentCavityVisits = 0;
    std::uint64_t parallelHazardPredicateTests = 0;
    std::uint64_t parallelHazardBroadPhaseCulls = 0;
    std::uint64_t parallelConflictPairChecks = 0;
    std::uint64_t parallelWriteMarkConflicts = 0;
    std::uint64_t parallelCandidateAabbCulls = 0;
    std::uint64_t parallelReadyObservations = 0;
    std::uint64_t parallelHullObservations = 0;
    std::uint64_t parallelSkippedObservations = 0;
    std::uint64_t parallelHullPrepassCandidates = 0;
    std::uint64_t parallelHullPrepassOutside = 0;
    double parallelHullPrepassSeconds = 0.0;
    double parallelDiscoverySeconds = 0.0;
    double canonicalCommitSeconds = 0.0;
    // Original-coordinate positivity audit and a bounded exact-coordinate
    // retry. Accepted counts describe the last construction; discarded work
    // retains explicitly counted visits from rejected construction attempts.
    std::uint64_t originalOrientationAuditTests = 0;
    std::uint64_t discardedConstructionWork = 0;
    bool originalCoordinateRetryUsed = false;
};

enum class DelaunayParallelPolicy : std::uint8_t {
    SerialExact = 0,
    StrictSerialEquivalent = 1,
    UniqueDelaunayParallel = 2
};

struct Delaunay3DOptions {
    bool collectDetailedOperationCounts = false;
    int requestedThreads = 1;
    std::uint32_t packetSize = 128;
    // Backward-compatible switch: when true and parallelPolicy remains
    // SerialExact, UniqueDelaunayParallel is selected.
    bool deterministicParallel = false;
    DelaunayParallelPolicy parallelPolicy = DelaunayParallelPolicy::SerialExact;
};

class DelaunayWorkspace {
public:
    struct Impl;
    DelaunayWorkspace();
    ~DelaunayWorkspace();
    DelaunayWorkspace(const DelaunayWorkspace&) = delete;
    DelaunayWorkspace& operator=(const DelaunayWorkspace&) = delete;
    DelaunayWorkspace(DelaunayWorkspace&&) noexcept;
    DelaunayWorkspace& operator=(DelaunayWorkspace&&) noexcept;
    [[nodiscard]] std::uint64_t retainedBytes() const noexcept;
    void release() noexcept;
private:
    std::unique_ptr<Impl> impl_;
    friend std::vector<Tet> delaunay3D(const std::vector<Vec3>&,
                                       DelaunayWorkspace&,
                                       const Delaunay3DOptions&,
                                       Delaunay3DStats*);
    friend std::vector<Tet> delaunay3D(const std::vector<Vec3>&,
                                       DelaunayWorkspace&,
                                       Delaunay3DStats*);
};

std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            DelaunayWorkspace& workspace,
                            const Delaunay3DOptions& options,
                            Delaunay3DStats* stats = nullptr);
std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            DelaunayWorkspace& workspace,
                            Delaunay3DStats* stats = nullptr);
std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            const Delaunay3DOptions& options,
                            Delaunay3DStats* stats = nullptr);
std::vector<Tet> delaunay3D(const std::vector<Vec3>& points,
                            Delaunay3DStats* stats = nullptr);

} // namespace paacs::algo
