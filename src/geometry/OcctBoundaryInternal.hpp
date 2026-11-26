// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private retained native tessellation; never installed with the public API.
#include "paacs/geometry/OcctBoundary.hpp"
#include <limits>
#include <optional>
#include <span>
namespace paacs::geo::detail {
class OcctSpatialSizing;
enum class OcctBoundaryAuditPolicy { PerBoundary, Deferred };
enum class OcctBoundaryFacetPolicy { Immediate, DeferredCandidate };
// Applied during prepare() to the single owned compound. Native tags refer to
// this model instance; the caller's edge span need only outlive that call.
// All incident faces are prepared together so neighbouring bodies retain the
// same curves, face triangulations and node keys. Empty selection is a no-op.
struct OcctBoundaryConstructionOptions {
    std::span<const int> refinedEdges{};
    double localLinearDeflection = 0;
    unsigned maximumRefinedEdges = 16, maximumIncidentFaces = 64;
};
struct OcctBoundaryFacetAuditOutcome {
    bool completed = false, sourceModelUnmodified = false;
    OcctBoundaryReason reason = OcctBoundaryReason::InternalError;
    std::uint64_t work = 0, samples = 0, closestTrimmedFaceSamples = 0;
    double maximumFacetSampleError = 0;
    int maximumErrorFace = -1;
    Vec3 maximumErrorSample{}, maximumErrorClosest{};
};
struct OcctBoundaryAliasOptions {
    // Partition consumers select original exterior faces; shared cuts stay fixed.
    std::span<const int> faceTags{};
    unsigned maximumFlips=64;
    std::uint64_t maximumWorkingBytes=64ULL*1024*1024;
    unsigned maximumRebuiltCharts=0;
    // A positive cap enables bounded curve redistribution on crossing trim
    // charts of a single native body.
    double curveMetricSizeCap=0;
    unsigned curvaturePointsPerTurn=0;
    const OcctSpatialSizing* sizeField=nullptr;
};
// One bounded construction witness. Tags identify the entity/use at the
// actual refusal, not a nearby entity inferred afterwards. Strings are static
// literals from the relevant private builder; no geometry is retained here.
struct OcctBoundaryConstructionWitness {
    const char* stage="not_reported";
    const char* detail="not_reported";
    int edgeTag=-1,faceTag=-1;
    std::optional<std::uint64_t> useId;
    // Conservative construction classification, not an assertion that the CAD
    // or an alternate surface is valid. Unknown/resource failures remain false.
    bool geometricRefusal=false;
};
struct OcctBoundaryAliasReceipt {
    std::uint64_t flips=0,maximumWorkingPayloadEstimate=0,rebuiltCharts=0,reusedMidpoints=0,rephasedEdges=0;
    double projectionSeconds=0;
    std::vector<CadSurfaceSpacingReport> surfaceSpacing;
    OcctBoundaryConstructionWitness witness;
};
// Build one solid from the original CAD curve and face uses. A caller still
// needs the geometric shell certificate before using it for volume meshing.
OcctBoundaryOutcome buildSpacedOcctBoundary(const OcctModel&,int solidTag,
    const OcctBoundaryOptions&,double targetSize,std::uint64_t maximumWorkingBytes,
    mesh::Mesh&,OcctBoundaryAliasReceipt&,bool regularFront=false,bool improveQuality=false,
    unsigned curvaturePointsPerTurn=0,const OcctSpatialSizing* sizeField=nullptr) noexcept;
class OcctBoundarySession final {
public:
    ~OcctBoundarySession();
    // Retains model by reference. Caller must keep it alive and immutable.
    // All bodies are deep-copied/tessellated once, so actual shared faces keep
    // common triangulations and keys. Work admission is cumulative per session;
    // output limits concern each extracted body's boundary, not the assembly.
    static std::unique_ptr<OcctBoundarySession> prepare(const OcctModel&,
        const OcctBoundaryOptions&, OcctBoundaryOutcome& preparation,
        OcctBoundaryAuditPolicy = OcctBoundaryAuditPolicy::PerBoundary) noexcept;
    static std::unique_ptr<OcctBoundarySession> prepare(const OcctModel&,
        const OcctBoundaryOptions&, OcctBoundaryOutcome& preparation,
        OcctBoundaryAuditPolicy, OcctBoundaryFacetPolicy) noexcept;
    // Explicit private preparation policy, applied only to the owned copy.
    // Existing entrypoints retain the SDK's default analytic-surface policy.
    static std::unique_ptr<OcctBoundarySession> prepare(const OcctModel&,
        const OcctBoundaryOptions&, OcctBoundaryOutcome& preparation,
        OcctBoundaryAuditPolicy, OcctBoundaryFacetPolicy,
        bool controlDeflectionOnAllSurfaces) noexcept;
    static std::unique_ptr<OcctBoundarySession> prepare(const OcctModel&,
        const OcctBoundaryOptions&, OcctBoundaryOutcome& preparation,
        OcctBoundaryAuditPolicy, OcctBoundaryFacetPolicy,
        bool controlDeflectionOnAllSurfaces,
        const OcctBoundaryConstructionOptions*) noexcept;
    // Interleaved volume/assembly stages can tighten this call's admission to
    // their remaining shared budget without increasing the session ceiling.
    OcctBoundaryOutcome boundary(int solidTag, mesh::Mesh& destination,
        std::uint64_t maximumAdditionalWork = (std::numeric_limits<std::uint64_t>::max)()) noexcept;
    OcctBoundaryOutcome boundary(int solidTag, mesh::Mesh& destination,
        const OcctBoundaryAliasOptions&,OcctBoundaryAliasReceipt&,
        std::uint64_t maximumAdditionalWork = (std::numeric_limits<std::uint64_t>::max)()) noexcept;
    // Deferred extraction is private staging only: each boundary receipt leaves
    // sourceModelUnmodified false. The transactional consumer must call this
    // AFTER every extraction/validation and BEFORE publishing its final mesh.
    // It compares the whole source's BREP bytes with the preparation snapshot;
    // success is reason=None/sourceModelUnmodified=true, not a mesh publication.
    // Public buildOcctBoundary and default sessions always audit per boundary.
    OcctBoundaryOutcome auditSource(
        std::uint64_t maximumAdditionalWork = (std::numeric_limits<std::uint64_t>::max)()) noexcept;
    std::uint64_t work() const noexcept;
private:
    struct Impl;
    explicit OcctBoundarySession(std::unique_ptr<Impl>);
    std::unique_ptr<Impl> impl_;
};
// PRIVATE staging only: native nodes/curves, provenance, topology, orientation,
// source bytes and surface-volume checks remain mandatory. Four facet samples
// are postponed; this boundary must not be represented as a fully audited
// public surface or published volume. Consumers must require the audit below
// after feasibility and before their final transactional publication.
OcctBoundaryOutcome buildOcctBoundaryCandidate(const OcctModel&, int solidTag,
    const OcctBoundaryOptions&, mesh::Mesh& destination) noexcept;
OcctBoundaryFacetAuditOutcome auditOcctBoundaryFacets(const OcctModel&, int solidTag,
    const mesh::Mesh&, const OcctBoundaryOptions&,
    std::uint64_t maximumSnapshotBytes = 16ULL*1024*1024) noexcept;
} // namespace paacs::geo::detail
