// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T5: bounded additive multi-region Tet4 conformity over boundary-conforming-v3.

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

namespace paacs::algo {

inline constexpr const char* kT5RouteId = "multi-region-tet4-v6";
inline constexpr std::uint32_t kT5RouteVersion = 6;
inline constexpr std::size_t kT5SixVolumeExpansionTermsPerTet = 192;
using RegionPersistentId = std::uint64_t;
using MaterialPersistentId = std::uint64_t;

struct T5RegionBinding {
    int regionTag = -1;
    RegionPersistentId regionPersistentId = 0;
    MaterialPersistentId materialPersistentId = 0;
    mesh::Classification volumeClassification{};
};

struct T5Request {
    BoundaryConformingV3Request sourceRequest{};
    std::span<const T5RegionBinding> regionBindings{};
    std::size_t maxRegionRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxMaterialRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxInterfaceNodeRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxInterfaceEdgeRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxInterfaceFacetRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxSourceNodeMappingRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxSourceSurfaceCellMappingRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxSourceTetMappingRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxIncidentReferences = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxVolumeExpansionTerms = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxReportRows = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxPublicationNodes = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxPublicationCells = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxPublicationConnectivityEntries = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxAuthorityCoreBytes = (std::numeric_limits<std::size_t>::max)();
    std::size_t maxPostSourceAccountedBytes = (std::numeric_limits<std::size_t>::max)();
    std::uint64_t maxAuditWork = (std::numeric_limits<std::uint64_t>::max)();
};

enum class T5Status : std::uint8_t {
    Published, Incomplete, Rejected, Unsupported, Refused, Failed
};
enum class T5Phase : std::uint8_t {
    None, Preflight, SourceConstruction, SourceValidation,
    ArtifactRevalidation, RegionAuthority, InterfaceAudit, ResultAssembly,
    MeshBuild, BoundaryAudit, AuthoritySerialization, Publication
};
enum class T5Reason : std::uint8_t {
    None, InvalidRequest, InvalidRegionBinding, DuplicateRegionTag,
    DuplicateRegionIdentity, UnknownArtifactRegion, UnusedRegionBinding,
    RegionClassificationMismatch, StaleSurfaceArtifact,
    StaleOwnershipProvider, SourceRejected, SourceUnsupported, SourceRefused,
    SourceFailed, SourceReportIncomplete, MalformedSourceOutcome,
    MalformedSurfaceOutcome, ArtifactStateMismatch,
    InvalidInterfaceOwnership, UnsupportedRegionJunction,
    DisconnectedRegionTagReuse, InterfaceIdentityMismatch,
    InterfaceOwnerMismatch, InterfaceClassificationMismatch,
    DuplicateInterfacePublication, RegionCellMismatch,
    BoundaryBijectionFailed, ReportBudgetReached, ResultBudgetReached,
    WorkBudgetReached, AuthorityBudgetReached, VolumeRangeUnrepresentable,
    FloatingPointEnvironmentUnsupported, MeshBuildFailed,
    PublicationPreparationFailed, AllocationFailed, InternalError
};
enum class T5WorkKind : std::uint8_t {
    Region, Material, InterfaceNode, InterfaceEdge, InterfaceFacet,
    SourceNodeMapping, SourceSurfaceCellMapping, SourceTetMapping
};

struct T5WorkKey {
    T5WorkKind kind = T5WorkKind::Region;
    RegionPersistentId region = 0;
    MaterialPersistentId material = 0;
    Index globalNode = kInvalidIndex;
    std::array<Index,2> globalEdge{kInvalidIndex,kInvalidIndex};
    bnd::VolumeFacetKey facet{};
    FacePersistentId surfaceFace = 0;
    std::uint64_t normalizedLocalCellOrdinal = 0;
    std::uint64_t sourceSlot = 0;
};
struct T5Failure {
    T5Reason reason = T5Reason::None;
    T5Phase phase = T5Phase::None;
    T5WorkKey work{};
    BoundaryConformingV3Status sourceStatus = BoundaryConformingV3Status::Failed;
    BoundaryConformingV3Failure sourceFailure{};
};
enum class T5WorkStatus : std::uint8_t {
    Verified, Published, Rejected, Unsupported, Refused, Failed,
    TerminalIncomplete
};
struct T5ProcessedRow {
    T5WorkKey key{};
    T5WorkStatus status = T5WorkStatus::Failed;
    std::uint64_t canonicalOrdinal = 0;
    T5Failure failure{};
};

struct T5Metrics {
    std::uint64_t sourceWorkUnits = 0, regionRows = 0,
                  cellClassificationVisits = 0, faceOwnerVisits = 0,
                  interfaceNodeInsertions = 0,
                  interfaceEdgeInsertions = 0,
                  interfaceFacetInsertions = 0,
                  interfaceResolverLookups = 0, regionAdjacencyEdges = 0,
                  volumeExpansionOperations = 0, publicationRows = 0,
                  auditWorkUnits = 0, totalWorkUnits = 0,
                  perInterfaceGlobalTetScans = 0,
                  perRegionGlobalFaceScans = 0, contextRebuilds = 0,
                  publicationTailAllocations = 0;
    std::size_t nestedV3PeakAccountedBytes = 0,
                postSourceAccountedHighWaterBytes = 0;
    T5Phase postSourceAccountedHighWaterPhase = T5Phase::None;
    bool runComplete = false;
};
struct T5Audit {
    std::uint64_t prescribedExterior = 0, matchedExterior = 0,
                  orphanExterior = 0, unexpectedExterior = 0,
                  prescribedInterface = 0, matchedInterface = 0,
                  missingInterface = 0, unexpectedInterface = 0,
                  duplicateInterface = 0, interfaceNodes = 0,
                  duplicateInterfaceNodes = 0, interfaceEdges = 0,
                  duplicateInterfaceEdges = 0, wrongOwnerSense = 0,
                  wrongRegionClassification = 0, undeclaredRegionFaces = 0,
                  invertedTets = 0, unusedNodes = 0, regionCount = 0,
                  materialCount = 0;
    bool passed = false;
};

struct T5RegionRow {
    T5RegionBinding authority{};
    std::uint64_t tetCount = 0, interfaceFacetCount = 0,
                  tetReferenceOffset = 0, tetReferenceCount = 0,
                  interfaceReferenceOffset = 0,
                  interfaceReferenceCount = 0,
                  sixVolumeExpansionOffset = 0,
                  sixVolumeExpansionCount = 0;
    double volumeLower = 0, volumeUpper = 0;
    std::uint64_t exactVolumeDigest = 0;
};
struct T5MaterialRow {
    MaterialPersistentId materialPersistentId = 0;
    std::uint64_t regionReferenceOffset = 0, regionReferenceCount = 0,
                  tetCount = 0, sixVolumeExpansionOffset = 0,
                  sixVolumeExpansionCount = 0;
    double volumeLower = 0, volumeUpper = 0;
    std::uint64_t exactVolumeDigest = 0;
};
struct T5InterfaceNodeRow {
    Index globalNode = kInvalidIndex;
    topo::NodeHandle authoritativeNode{};
    mesh::Classification classification{};
    std::uint64_t incidentFacetOffset = 0, incidentFacetCount = 0;
};
struct T5InterfaceEdgeRow {
    std::array<Index,2> canonicalGlobalNodes{kInvalidIndex,kInvalidIndex};
    topo::EdgeHandle authoritativeEdge{};
    mesh::Classification classification{};
    bool required = false;
    std::uint64_t incidentArtifactFacetOffset = 0,
                  incidentArtifactFacetCount = 0;
};
struct T5InterfaceFacetRow {
    bnd::VolumeFacetKey facet{};
    topo::FaceHandle authoritativeFace{};
    std::array<Index,3> canonicalPublishedNodes{
        kInvalidIndex,kInvalidIndex,kInvalidIndex};
    mesh::Classification classification{};
    RegionPersistentId firstRegion = 0, secondRegion = 0;
    MaterialPersistentId firstMaterial = 0, secondMaterial = 0;
    int firstOutwardSense = 0, secondOutwardSense = 0;
    topo::CellHandle firstOwner{}, secondOwner{};
    Index publishedTri3Cell = kInvalidIndex;
};
struct T5SourceNodeMappingRow {
    std::uint64_t authoritativeSlot = 0;
    topo::NodeHandle authoritativeNode{};
    Index sourceGlobalNode = kInvalidIndex, publishedGlobalNode = kInvalidIndex;
};
struct T5SourceSurfaceCellMappingRow {
    FacePersistentId face = 0;
    std::uint64_t normalizedLocalCellOrdinal = 0;
    Index sourceGlobalCell = kInvalidIndex, publishedGlobalCell = kInvalidIndex;
    mesh::Classification classification{};
};
struct T5SourceTetMappingRow {
    std::uint64_t authoritativeSlot = 0;
    topo::CellHandle authoritativeCell{};
    Index sourceGlobalCell = kInvalidIndex, publishedGlobalCell = kInvalidIndex;
    mesh::Classification classification{};
};

struct T5Report {
    T5Report() = default;
    T5Report(const T5Report&) = delete;
    T5Report& operator=(const T5Report&) = delete;
    T5Report(T5Report&&) = delete;
    T5Report& operator=(T5Report&&) = delete;
    std::vector<T5ProcessedRow> processedStorage;
    std::vector<T5WorkKey> unattemptedStorage;
    std::uint64_t processedUsed = 0, unattemptedUsed = 0;
    bool reportComplete = false;
};
struct T5V3SourceReceipt {
    T5V3SourceReceipt() = default;
    T5V3SourceReceipt(const T5V3SourceReceipt&) = delete;
    T5V3SourceReceipt& operator=(const T5V3SourceReceipt&) = delete;
    T5V3SourceReceipt(T5V3SourceReceipt&&) = delete;
    T5V3SourceReceipt& operator=(T5V3SourceReceipt&&) = delete;
    const char* routeId = kBoundaryConformingV3RouteId;
    std::uint32_t routeVersion = kBoundaryConformingV3RouteVersion;
    ArtifactStamp surfaceSetStamp{};
    VolumeBoundaryProviderStamp providerStamp{};
    std::unique_ptr<SurfaceSetOutcome> surfaceOutcome;
    std::unique_ptr<BoundaryConformingV3ReportEnvelope> sourceReport;
    ArtifactPreflightMetrics artifactMetrics{};
    BoundaryConformingV3ConstructionMetrics constructionMetrics{};
    TetRecoveryMetrics recoveryMetrics{};
    VolumeExtractionMetrics extractionMetrics{};
    BoundaryConformingV3MemoryMetrics memoryMetrics{};
    BoundaryConformingV3Audit sourceAudit{};
    std::uint64_t sourceStrictHash = 0;
    mesh::SemanticHash sourceSemanticHash{};
    std::uint64_t sourceNodeMappingDigest = 0,
                  sourceSurfaceCellMappingDigest = 0,
                  sourceTetMappingDigest = 0, sourceNodeCount = 0,
                  sourceSurfaceCellCount = 0, sourceTetCount = 0;
};
struct T5Result {
    T5Result() = default;
    T5Result(const T5Result&) = delete;
    T5Result& operator=(const T5Result&) = delete;
    T5Result(T5Result&&) = delete;
    T5Result& operator=(T5Result&&) = delete;
    const char* routeId = kT5RouteId;
    std::uint32_t routeVersion = kT5RouteVersion;
    std::unique_ptr<TetAuthoritativeResult> authoritative;
    bnd::TetFacetSetState constraintState;
    std::unique_ptr<bnd::ClosedTri3Artifact> artifact;
    std::vector<T5RegionRow> regions;
    std::vector<T5MaterialRow> materials;
    std::vector<T5InterfaceNodeRow> interfaceNodes;
    std::vector<T5InterfaceEdgeRow> interfaceEdges;
    std::vector<T5InterfaceFacetRow> interfaceFacets;
    std::vector<bnd::VolumeFacetKey> interfaceNodeIncidentFacets;
    std::vector<bnd::VolumeFacetKey> interfaceEdgeIncidentArtifactFacets;
    std::vector<bnd::VolumeFacetKey> regionInterfaceFacets;
    std::vector<topo::CellHandle> regionTets;
    std::vector<RegionPersistentId> materialRegions;
    std::vector<T5SourceNodeMappingRow> sourceNodeMappings;
    std::vector<T5SourceSurfaceCellMappingRow> sourceSurfaceCellMappings;
    std::vector<T5SourceTetMappingRow> sourceTetMappings;
    std::vector<double> exactVolumeExpansionTerms;
    std::vector<std::uint8_t> authorityCoreBytes;
    std::array<std::uint8_t,32> strictAuthorityCoreDigest{};
    std::array<std::uint8_t,32> semanticAuthorityCoreDigest{};
    T5Metrics metrics{};
    T5Audit audit{};
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
    std::unique_ptr<mesh::Mesh> privateMesh;
    bool publicationReady = false;
};
struct T5Outcome {
    T5Outcome();
    ~T5Outcome();
    T5Outcome(const T5Outcome&) = delete;
    T5Outcome& operator=(const T5Outcome&) = delete;
    T5Outcome(T5Outcome&&) noexcept;
    T5Outcome& operator=(T5Outcome&&) noexcept;
    T5Status status = T5Status::Failed;
    T5Failure failure{};
    bool destinationUnchanged = true, reportComplete = false;
    std::unique_ptr<BoundaryConformingV3Outcome> sourceOutcome;
    std::unique_ptr<T5V3SourceReceipt> sourceReceipt;
    std::unique_ptr<T5Report> report;
    std::unique_ptr<T5Result> retainedResult;
    bool published() const noexcept { return status == T5Status::Published; }
};

enum class T5ExportFormat : std::uint8_t { VtkLegacy, Msh41 };
struct T5ExportBindingV1 {
    std::uint32_t schemaVersion = 1;
    T5ExportFormat format = T5ExportFormat::VtkLegacy;
    std::uint32_t exporterVersion = 1;
    std::array<std::uint8_t,32> exportedBytesSha256{},
                               strictAuthorityCoreDigest{},
                               semanticAuthorityCoreDigest{},
                               exportBindingDigest{};
};
enum class T5ExportReason : std::uint8_t {
    None, InvalidResult, MeshWriteFailed, CoreWriteFailed,
    BindingWriteFailed, ReadFailed, DigestMismatch, MalformedBinding
};
struct T5ExportOutcome {
    bool complete = false;
    T5ExportReason reason = T5ExportReason::InvalidResult;
    T5ExportBindingV1 binding{};
};
struct T5ExportReloadAudit {
    bool passed = false;
    T5ExportReason reason = T5ExportReason::ReadFailed;
    std::uint64_t nodeCount = 0, cellCount = 0;
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
    bool meshHashesMatch = false;
    bool authorityProjectionValid = false;
    bool authorityMappingsMatch = false;
};

namespace t5_detail {
inline bool& routeAllocationActive() noexcept {
    static thread_local bool active = false;
    return active;
}
inline bool& publicationPreparationActive() noexcept {
    static thread_local bool active = false;
    return active;
}
inline bool& publicationTailActive() noexcept {
    static thread_local bool active = false;
    return active;
}
std::array<std::uint8_t,32> sha256(
    std::span<const std::uint8_t> bytes) noexcept;
T5Reason buildAuthorityCore(T5Result& result, ArtifactStamp surfaceStamp,
                            VolumeBoundaryProviderStamp providerStamp,
                            std::size_t maxBytes) noexcept;
bool decodeAuthorityCoreForTesting(std::span<const std::uint8_t> bytes,
                                   std::uint64_t& regionCount,
                                   std::uint64_t& interfaceFacetCount) noexcept;
#ifdef PAACS_TOPO_TESTING
struct ExactSixVolumeProbe {
    bool complete = false;
    std::size_t termCount = 0;
    std::uint64_t workUnits = 0;
    std::uint64_t exactDigest = 0;
    double volumeLower = 0;
    double volumeUpper = 0;
};
ExactSixVolumeProbe exactSixVolumeForTesting(
    const std::array<Vec3,4>& points) noexcept;
#endif
}

T5Outcome meshModelMultiRegionTet4(
    const geo::Model& model, const sizing::MeshSpec& spec,
    const T5Request& request, mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider) noexcept;
T5ExportOutcome writeT5VtkExport(
    const mesh::Mesh& mesh, const T5Result& result,
    const std::string& meshPath, const std::string& corePath,
    const std::string& bindingPath) noexcept;
T5ExportOutcome writeT5MshExport(
    const mesh::Mesh& mesh, const T5Result& result,
    const std::string& meshPath, const std::string& corePath,
    const std::string& bindingPath) noexcept;
T5ExportReloadAudit reloadT5ExportBinding(
    const std::string& meshPath, const std::string& corePath,
    const std::string& bindingPath) noexcept;

static_assert(kT5SixVolumeExpansionTermsPerTet == 192);
static_assert(std::is_standard_layout_v<T5RegionBinding>);
static_assert(std::is_trivially_copyable_v<T5RegionBinding>);
static_assert(std::is_standard_layout_v<T5WorkKey>);
static_assert(std::is_trivially_copyable_v<T5WorkKey>);
static_assert(std::is_standard_layout_v<T5Failure>);
static_assert(std::is_trivially_copyable_v<T5Failure>);
static_assert(std::is_standard_layout_v<T5ProcessedRow>);
static_assert(std::is_trivially_copyable_v<T5ProcessedRow>);
static_assert(std::is_standard_layout_v<T5Metrics>);
static_assert(std::is_trivially_copyable_v<T5Metrics>);
static_assert(std::is_standard_layout_v<T5Audit>);
static_assert(std::is_trivially_copyable_v<T5Audit>);
static_assert(std::is_nothrow_move_constructible_v<T5Outcome>);
static_assert(!std::is_move_constructible_v<T5Result>);
static_assert(!std::is_move_constructible_v<T5Report>);
static_assert(!std::is_move_constructible_v<T5V3SourceReceipt>);
static_assert(noexcept(std::declval<mesh::Mesh&>().swap(
    std::declval<mesh::Mesh&>())));

} // namespace paacs::algo
