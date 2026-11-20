// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// T4-B: bounded Tet10 curving and complete elementwise Jacobian
// certification over the accepted T4-A Tet4 result.

#include "paacs/algo/Tet10JacobianCertificate.hpp"
#include "paacs/algo/TetQualityMesher.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

namespace paacs::algo {

inline constexpr const char* kT4BRouteId = "curved-tet10-v5";
inline constexpr std::uint32_t kT4BRouteVersion = 1;

enum class Tet10PlacementKind : std::uint8_t { Curve, Surface };

struct Tet10PlacementStamp {
    std::uint64_t persistentId = 0;
    std::uint64_t revision = 0;
    friend bool operator==(const Tet10PlacementStamp&,
                           const Tet10PlacementStamp&) = default;
};

struct Tet10SourceEdgeKey {
    std::array<RefinementVertexKey,2> vertices{};
    friend bool operator==(const Tet10SourceEdgeKey&,
                           const Tet10SourceEdgeKey&) = default;
    friend bool operator<(const Tet10SourceEdgeKey& a,
                          const Tet10SourceEdgeKey& b) noexcept {
        return a.vertices<b.vertices;
    }
};

struct Tet10PlacementRequest {
    std::array<RefinementVertexKey,2> endpoints{};
    std::array<Vec3,2> endpointCoordinates{};
    mesh::Classification classification{};
    Tet10PlacementKind kind = Tet10PlacementKind::Surface;
    std::span<const bnd::VolumeFacetKey> incidentFacets{};
};

enum class Tet10PlacementStatus : std::uint8_t {
    Placed,
    Rejected,
    Unsupported,
    Failed
};

struct Tet10PlacementResult {
    Tet10PlacementStatus status = Tet10PlacementStatus::Failed;
    std::uint32_t providerReason = 0;
    Tet10PlacementStamp observedStamp{};
    Vec3 coordinate{};
};

class Tet10PlacementProvider {
public:
    virtual ~Tet10PlacementProvider() = default;
    virtual Tet10PlacementStamp stamp() const = 0;
    virtual Tet10PlacementResult place(
        const Tet10PlacementRequest& request) const = 0;
};

struct T4BPolicy {
    double minimumCertifiedJacobianRatio = 0.0;
    std::uint32_t maxCertificateDepth = 0;
    std::uint64_t maxCertificateCells = 0;
    std::uint64_t maxCoefficientTests = 0;
    std::uint64_t maxWorkUnits = 0;
    std::uint64_t maxResultNodes = 0;
    std::uint64_t maxResultCells = 0;
    std::uint64_t maxResultConnectivity = 0;
    std::uint64_t maxPlacementRows = 0;
    std::uint64_t maxCertificateRows = 0;
    std::uint64_t maxUnattemptedRows = 0;
};

struct T4BRequest {
    T4ARequest source{};
    Tet10PlacementStamp expectedPlacementStamp{};
    T4BPolicy policy{};
};

enum class T4BStatus : std::uint8_t {
    Published,
    Incomplete,
    Rejected,
    Unsupported,
    Failed
};

enum class T4BPhase : std::uint8_t {
    ReportShell,
    RequestPreflight,
    SourceT4A,
    SourceValidation,
    EdgeExtraction,
    Placement,
    CertificatePreparation,
    JacobianCertification,
    FinalAudit,
    PublicationPreparation,
    Publication
};

enum class T4BReason : std::uint16_t {
    None,
    InvalidRequest,
    CapacityOverflow,
    SourceT4AStopped,
    InvalidAuthoritativeSource,
    ForeignPlacementProvider,
    StalePlacementProvider,
    PlacementProviderIdentityChanged,
    PlacementRejected,
    PlacementUnsupported,
    PlacementFailed,
    PlacementProviderThrew,
    NonfinitePlacement,
    UnsupportedSourceCell,
    ConnectivityOverflow,
    MissingEdgeBinding,
    DuplicateEdgeBinding,
    ForeignEdgeBinding,
    StaleEdgeBinding,
    ConflictingEdgeClassification,
    AmbiguousBoundaryOwnership,
    MidpointRangeUnrepresentable,
    ResultBudgetReached,
    ReportBudgetReached,
    JacobianCertificationBudgetReached,
    JacobianBelowThreshold,
    NonpositiveJacobian,
    CertificateRangeUnrepresentable,
    InvalidElevatedTopology,
    BoundaryBijectionFailed,
    PublicationPreparationFailed,
    AllocationFailed,
    InternalError
};

struct T4BFailure {
    T4BReason reason = T4BReason::None;
    T4BPhase phase = T4BPhase::RequestPreflight;
    Tet10SourceEdgeKey edge{};
    std::uint64_t cellOrdinal = 0;
};

enum class T4BPlacementDisposition : std::uint8_t {
    InteriorMidpoint,
    Placed,
    Unsupported,
    Rejected,
    Failed,
    TerminalFailure
};

struct T4BPlacementRow {
    Tet10SourceEdgeKey edge{};
    Index midsideIndex = kInvalidIndex;
    mesh::Classification classification{};
    Vec3 coordinate{};
    Tet10PlacementStatus providerStatus = Tet10PlacementStatus::Failed;
    std::uint32_t providerReason = 0;
    T4BPlacementDisposition disposition =
        T4BPlacementDisposition::TerminalFailure;
};

struct T4BCertificateRow {
    std::uint64_t sourceCellOrdinal = 0;
    Tet10CertificateStatus status = Tet10CertificateStatus::InvalidInput;
    Tet10CertificateWitness witness{};
    Tet10CertificateMetrics metrics{};
    double conservativeMinimumLowerBound =
        (std::numeric_limits<double>::infinity)();
};

enum class T4BBoundaryCompatibility : std::uint8_t {
    ProviderTrusted
};

struct T4BReport {
    T4BReport() = default;
    T4BReport(const T4BReport&) = delete;
    T4BReport& operator=(const T4BReport&) = delete;
    T4BReport(T4BReport&&) = delete;
    T4BReport& operator=(T4BReport&&) = delete;
    std::vector<T4BPlacementRow> placements;
    std::vector<T4BCertificateRow> certificates;
    std::vector<Tet10SourceEdgeKey> unattemptedEdges;
    std::vector<std::uint64_t> unattemptedCells;
    T4BBoundaryCompatibility boundaryCompatibility =
        T4BBoundaryCompatibility::ProviderTrusted;
};

struct T4BMetrics {
    std::uint64_t sourceCellsScanned = 0;
    std::uint64_t sourceEdgeUses = 0;
    std::uint64_t uniqueEdges = 0;
    std::uint64_t interiorEdges = 0;
    std::uint64_t curveEdges = 0;
    std::uint64_t surfaceEdges = 0;
    std::uint64_t providerStampCalls = 0;
    std::uint64_t providerCalls = 0;
    std::uint64_t coefficientTests = 0;
    std::uint64_t certificateSubcells = 0;
    std::uint64_t auditRows = 0;
    std::uint64_t workUnits = 0;
    std::uint64_t accountedHighWaterBytes = 0;
    std::uint64_t perEdgeGlobalScans = 0;
    std::uint64_t perCellGlobalScans = 0;
    std::uint64_t contextRebuilds = 0;
    std::uint64_t publicationTailAllocations = 0;
    std::uint32_t maximumCertificateDepth = 0;
    bool runComplete = false;
};

struct T4BSourceReceipt {
    T4BSourceReceipt() = default;
    T4BSourceReceipt(const T4BSourceReceipt&) = delete;
    T4BSourceReceipt& operator=(const T4BSourceReceipt&) = delete;
    T4BSourceReceipt(T4BSourceReceipt&&) = delete;
    T4BSourceReceipt& operator=(T4BSourceReceipt&&) = delete;
    std::unique_ptr<T3AOutcome> sourceOutcome;
    std::unique_ptr<T4AT3SourceReceipt> sourceReceipt;
    std::unique_ptr<T4AReport> sourceReport;
    T4AQualitySummary before{},after{};
    T4AMetrics metrics{};
    BoundaryConformingV3Audit boundaryAudit{};
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
};

struct T4BResult {
    T4BResult() = default;
    T4BResult(const T4BResult&) = delete;
    T4BResult& operator=(const T4BResult&) = delete;
    T4BResult(T4BResult&&) = delete;
    T4BResult& operator=(T4BResult&&) = delete;
    const char* routeId = kT4BRouteId;
    std::uint32_t routeVersion = kT4BRouteVersion;
    std::unique_ptr<T4AMutableArtifact> sourceArtifact;
    std::unique_ptr<mesh::Mesh> privateMesh;
    std::vector<Index> sourceNodeToOutput;
    std::vector<Index> sourceCellToOutput;
    std::vector<T4BPlacementRow> edgeToMidside;
    std::vector<T4BCertificateRow> certificates;
    T4BMetrics metrics{};
    mesh::AuditResult audit{};
    std::uint64_t strictHash = 0;
    mesh::SemanticHash semanticHash{};
    bool certificationComplete = false;
    bool publicationReady = false;
};

struct T4BOutcome {
    T4BOutcome() = default;
    T4BOutcome(const T4BOutcome&) = delete;
    T4BOutcome& operator=(const T4BOutcome&) = delete;
    T4BOutcome(T4BOutcome&&) noexcept = default;
    T4BOutcome& operator=(T4BOutcome&&) noexcept = default;
    T4BStatus status = T4BStatus::Failed;
    T4BFailure failure{};
    bool destinationUnchanged = true;
    bool reportComplete = false;
    std::unique_ptr<T4AOutcome> sourceOutcome;
    std::unique_ptr<T4BSourceReceipt> sourceReceipt;
    std::unique_ptr<T4BReport> report;
    std::unique_ptr<T4BResult> retainedResult;
    bool published() const noexcept { return status==T4BStatus::Published; }
};

T4BOutcome meshModelCurvedTet10(
    const geo::Model& model,const sizing::MeshSpec& spec,
    const T4BRequest& request,mesh::Mesh& destination,
    const VolumeBoundaryOwnershipProvider& ownershipProvider,
    const Tet10PlacementProvider& placementProvider) noexcept;

#if defined(PAACS_TOPO_TESTING)
namespace t4b_detail {
inline bool& routeActive() noexcept {
    static thread_local bool value=false;return value;
}
inline bool& publicationPreparationActive() noexcept {
    static thread_local bool value=false;return value;
}
inline bool& publicationTailActive() noexcept {
    static thread_local bool value=false;return value;
}

T4BReason validateSourceAuthorityForTesting(
    const mesh::Mesh& source,T4AMutableArtifact& artifact,
    const ArtifactStamp& expectedArtifact) noexcept;
}
#endif

static_assert(std::is_standard_layout_v<Tet10SourceEdgeKey>);
static_assert(std::is_trivially_copyable_v<Tet10SourceEdgeKey>);
static_assert(std::is_standard_layout_v<T4BPolicy>);
static_assert(std::is_trivially_copyable_v<T4BPolicy>);
static_assert(std::is_standard_layout_v<T4BPlacementRow>);
static_assert(std::is_trivially_copyable_v<T4BPlacementRow>);
static_assert(std::is_nothrow_move_constructible_v<T4BPlacementRow>);
static_assert(std::is_standard_layout_v<T4BCertificateRow>);
static_assert(std::is_trivially_copyable_v<T4BCertificateRow>);
static_assert(std::is_nothrow_move_constructible_v<T4BCertificateRow>);
static_assert(std::is_standard_layout_v<T4BMetrics>);
static_assert(std::is_trivially_copyable_v<T4BMetrics>);
static_assert(std::is_nothrow_move_constructible_v<T4BOutcome>);
static_assert(std::is_nothrow_move_assignable_v<T4BOutcome>);
static_assert(!std::is_move_constructible_v<T4BReport>);
static_assert(!std::is_move_constructible_v<T4BSourceReceipt>);
static_assert(!std::is_move_constructible_v<T4BResult>);

} // namespace paacs::algo
