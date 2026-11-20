// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paacs/algo/Meshers.hpp"
#include "paacs/algo/ProductionModelMesher.hpp"
#include "paacs/core/BuildInfo.hpp"
#include "paacs/mesh/MeshAudit.hpp"

#include "test_harness.hpp"

extern int gFailures;

namespace {
using namespace paacs;

struct MeshSnapshot {
    std::vector<double> coordinates;
    std::vector<std::int64_t> topology;
    friend bool operator==(const MeshSnapshot&, const MeshSnapshot&) = default;
};

MeshSnapshot snapshot(const mesh::Mesh& value) {
    MeshSnapshot out;
    out.topology.push_back(value.nodeCount());
    out.topology.push_back(value.cellCount());
    out.topology.push_back(
        static_cast<std::int64_t>(value.dedupTolerance() * 1.0e15));
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const auto& p = value.node(i);
        out.coordinates.insert(out.coordinates.end(), {p.x, p.y, p.z});
        out.topology.push_back(value.nodeClass(i).dim);
        out.topology.push_back(value.nodeClass(i).tag);
    }
    for (Index i = 0; i < value.cellCount(); ++i) {
        const auto& cell = value.cell(i);
        out.topology.push_back(static_cast<std::int64_t>(cell.type));
        out.topology.push_back(cell.onEntity.dim);
        out.topology.push_back(cell.onEntity.tag);
        for (int n = 0; n < mesh::nodesPerCell(cell.type); ++n)
            out.topology.push_back(value.cellNodes(i)[n]);
    }
    return out;
}

void seedDestination(mesh::Mesh& value) {
    value.setDedupTolerance(0.125);
    const Index a = value.addNode({-3, -2, -1}, {0, 91});
    const Index b = value.addNode({-2, -2, -1}, {0, 92});
    value.addCell(mesh::CellType::Line2, {a, b}, {1, 93});
}

sizing::MeshSpec constantSpec(double h, int threads = 1) {
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(h);
    spec.threads = threads;
    return spec;
}

enum class SizingAdversary {
    Nan,
    PositiveInfinity,
    Zero,
    Negative,
    Throwing,
    Allocation,
    FaceCenterZero,
    SolidCenterZero,
    ShieldZero,
    LatticeZero,
    TranslatedCenterTiny,
    Tiny
};

class AdversarialField final : public sizing::SizeField {
public:
    explicit AdversarialField(SizingAdversary kind) : kind_(kind) {}

    double eval(const Vec3& p) const override {
        switch (kind_) {
            case SizingAdversary::Nan:
                return (std::numeric_limits<double>::quiet_NaN)();
            case SizingAdversary::PositiveInfinity:
                return (std::numeric_limits<double>::infinity)();
            case SizingAdversary::Zero:
                return 0.0;
            case SizingAdversary::Negative:
                return -0.25;
            case SizingAdversary::Throwing:
                throw std::runtime_error("sizing evaluation adversary");
            case SizingAdversary::Allocation:
                throw std::bad_alloc{};
            case SizingAdversary::FaceCenterZero: {
                const int boundary = (p.x == 0.0 || p.x == 1.0 ? 1 : 0) +
                                     (p.y == 0.0 || p.y == 1.0 ? 1 : 0) +
                                     (p.z == 0.0 || p.z == 1.0 ? 1 : 0);
                const int center = (p.x == 0.5 ? 1 : 0) +
                                   (p.y == 0.5 ? 1 : 0) +
                                   (p.z == 0.5 ? 1 : 0);
                return boundary == 1 && center == 2 ? 0.0 : 0.3;
            }
            case SizingAdversary::SolidCenterZero:
                return p.x == 0.5 && p.y == 0.5 && p.z == 0.5
                           ? 0.0
                           : 0.3;
            case SizingAdversary::ShieldZero: {
                const double radius = std::hypot(p.x, p.y);
                return p.z > 0.0 && p.z < 1.0 &&
                               radius > 0.35 && radius < 0.49
                           ? 0.0
                           : 0.15;
            }
            case SizingAdversary::LatticeZero: {
                const bool interior = p.x > 0.0 && p.x < 1.0 &&
                                      p.y > 0.0 && p.y < 1.0 &&
                                      p.z > 0.0 && p.z < 1.0;
                const bool center =
                    p.x == 0.5 && p.y == 0.5 && p.z == 0.5;
                return interior && !center ? -0.3 : 0.3;
            }
            case SizingAdversary::TranslatedCenterTiny:
                return p.x == 1.0e9 + 0.5 &&
                               p.y == 1.0e9 + 0.5 &&
                               p.z == 1.0e9 + 0.5
                           ? 1.0e-8
                           : 0.3;
            case SizingAdversary::Tiny:
                return (std::numeric_limits<double>::denorm_min)();
        }
        return 0.3;
    }

private:
    SizingAdversary kind_;
};

sizing::MeshSpec adversarialSpec(SizingAdversary kind) {
    sizing::MeshSpec spec;
    spec.field = std::make_shared<AdversarialField>(kind);
    spec.threads = 1;
    return spec;
}

void expectStructuredSizingStop(
    const geo::Model& model, const sizing::MeshSpec& invalid,
    algo::ProductionModelReason reason, algo::ProductionModelPhase phase,
    double retrySize) {
    mesh::Mesh destination;
    seedDestination(destination);
    const auto before = snapshot(destination);
    const auto outcome =
        algo::meshModelConstrained(model, invalid, destination);
    CHECK(!outcome.published());
    CHECK(outcome.status == algo::ProductionModelStatus::Failed);
    CHECK(outcome.reason == reason);
    CHECK(outcome.phase == phase);
    CHECK(outcome.destinationUnchanged);
    CHECK(snapshot(destination) == before);
    if (phase == algo::ProductionModelPhase::SurfaceConstruction &&
        reason == algo::ProductionModelReason::SizingEvaluationFailed) {
        CHECK(outcome.surface.reason ==
              algo::ProductionSurfaceReason::SizingEvaluationFailed);
        CHECK(outcome.surface.phase == algo::ProductionSurfacePhase::Spacing);
    }
    const auto retry = algo::meshModelConstrained(
        model, constantSpec(retrySize), destination);
    CHECK(retry.published());
}

enum class FacePresentation { Natural, Rotated, Reversed };
enum class FaceFailure { None, Geometry, Candidate };

class PresentedFace final : public geo::GFace {
public:
    PresentedFace(const geo::GFace& source, FacePresentation presentation,
                  FaceFailure failure, int failingTag)
        : GFace(source.tag()), source_(source), presentation_(presentation),
          failure_(source.tag() == failingTag ? failure : FaceFailure::None) {}

    Vec3 pointAt(const Vec2& uv) const override {
        if (failure_ == FaceFailure::Geometry)
            throw std::runtime_error("constrained-v2 geometry adversary");
        if (failure_ == FaceFailure::Candidate) {
            Vec2 lo, hi;
            source_.paramBounds(lo, hi);
            const Vec2 center{(lo.u + hi.u) * 0.5, (lo.v + hi.v) * 0.5};
            const bool strictInterior = uv.u > lo.u && uv.u < hi.u &&
                                        uv.v > lo.v && uv.v < hi.v;
            const bool isCenter = uv.u == center.u && uv.v == center.v;
            if (strictInterior && !isCenter)
                throw std::runtime_error("constrained-v2 candidate adversary");
        }
        return source_.pointAt(uv);
    }
    Vec2 paramAt(const Vec3& point) const override {
        return source_.paramAt(point);
    }
    Vec3 normalAt(const Vec2& uv) const override {
        return source_.normalAt(uv);
    }
    void paramBounds(Vec2& lo, Vec2& hi) const override {
        source_.paramBounds(lo, hi);
    }
    std::vector<int> boundaryEdgeTags() const override {
        auto tags = source_.boundaryEdgeTags();
        if (presentation_ == FacePresentation::Rotated && tags.size() > 1)
            std::rotate(tags.begin(), tags.begin() + 1, tags.end());
        else if (presentation_ == FacePresentation::Reversed)
            std::reverse(tags.begin(), tags.end());
        return tags;
    }
    double curvatureAt(const Vec2& uv) const override {
        return source_.curvatureAt(uv);
    }
    Vec2 metricScale() const override { return source_.metricScale(); }
    bool containsParam(const Vec2& uv) const override {
        return source_.containsParam(uv);
    }

private:
    const geo::GFace& source_;
    FacePresentation presentation_;
    FaceFailure failure_;
};

class PresentedModel final : public geo::Model {
public:
    PresentedModel(const geo::Model& source, FacePresentation presentation,
                   FaceFailure failure = FaceFailure::None,
                   bool duplicateFirstEdge = false)
        : source_(source), duplicateFirstEdge_(duplicateFirstEdge) {
        const auto sourceFaces = source.faces();
        const int failingTag = sourceFaces.empty() ? -1 : sourceFaces.front()->tag();
        for (const auto* face : sourceFaces)
            faces_.push_back(std::make_unique<PresentedFace>(
                *face, presentation, failure, failingTag));
    }
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        auto out = source_.edges();
        if (duplicateFirstEdge_ && !out.empty()) out.push_back(out.front());
        return out;
    }
    std::vector<const geo::GFace*> faces() const override {
        std::vector<const geo::GFace*> out;
        for (const auto& face : faces_) out.push_back(face.get());
        return out;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override {
        return source_.edge(tag);
    }
    const geo::GFace* face(int tag) const override {
        for (const auto& face : faces_)
            if (face->tag() == tag) return face.get();
        return nullptr;
    }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }

private:
    const geo::Model& source_;
    bool duplicateFirstEdge_;
    std::vector<std::unique_ptr<PresentedFace>> faces_;
};

struct BaselineCase {
    const char* name;
    Vec3 lo;
    Vec3 hi;
    double h;
    bool cylinder;
    bool elevate;
    int threads;
    std::uint64_t semantic;
    std::uint64_t strict;
    Index nodes;
    Index cells;
    Index tet4;
    Index tet10;
    Index tri3;
    Index tri6;
    Index boundary;
};

// Semantic values are the independently reviewed constrained-v2 observations;
// exact counts are reproduced by the versioned recorder. Reference-environment
// strict values are filled by the canonical GCC qualification run.
const std::array<BaselineCase, 7> kV2Cases{{
    {"unit-box", {0,0,0}, {1,1,1}, 0.3, false, false, 1,
     0x76c447c6cca9e3d8ULL, 0x8345b73101487f94ULL,
     88,438,246,0,156,0,156},
    {"cylinder", {}, {}, 0.15, true, false, 1,
     0x21ada68f58654d55ULL, 0xb0cbadfc5f7bbfe7ULL,
     543,3011,2369,0,580,0,580},
    {"thin-box", {0,0,0}, {1,1,1e-3}, 0.2, false, false, 1,
     0x74cd3aefe1bcf639ULL, 0x25914c912d1cf4a0ULL,
     80,374,174,0,156,0,156},
    {"translated-box", {1e9,1e9,1e9}, {1e9+1,1e9+1,1e9+1}, 0.3,
     false, false, 1, 0xceacb1fb013150bfULL, 0xcc53e0ec29535420ULL,
     88,438,246,0,156,0,156},
    {"scaled-box", {0,0,0}, {1e6,1e6,1e6}, 3e5, false, false, 1,
     0x6552f4202c20d9a3ULL, 0xa687a5ec558bea5cULL,
     88,438,246,0,156,0,156},
    {"tet10-box", {0,0,0}, {1,1,1}, 0.35, false, true, 1,
     0xfa7cbcdd76d2326aULL, 0x5b5ee511dda5ea53ULL,
     358,321,0,177,0,108,108},
    {"threads4-cylinder", {}, {}, 0.15, true, false, 4,
     0x21ada68f58654d55ULL, 0xb0cbadfc5f7bbfe7ULL,
     543,3011,2369,0,580,0,580},
}};

void testConstrainedV2Baseline() {
    TEST("constrained_v2.baseline_manifest");
    auto kernel = geo::makeBuiltinKernel();
    for (const auto& expected : kV2Cases) {
        auto model = expected.cylinder
            ? kernel->makeCylinder({0,0,0}, 0.5, 1.0)
            : kernel->makeBox(expected.lo, expected.hi);
        auto spec = constantSpec(expected.h, expected.threads);
        if (expected.elevate) {
            spec.order = 2;
            spec.orderFactor = 1.0;
        }
        mesh::Mesh result;
        const auto outcome = algo::meshModelConstrained(*model, spec, result);
        CHECK(outcome.published());
        CHECK(outcome.tetHandoff.runComplete);
        CHECK(outcome.tetHandoff.workAccountingComplete);
        CHECK(outcome.tetHandoff.workspaceValid);
        CHECK(outcome.tetHandoff.candidateReleasedBeforePublication);
        CHECK(outcome.tetHandoff.filterRepresentationReleasedBeforePublication);
        CHECK(outcome.tetHandoff.completeTetRepresentationsAtOverlap == 2);
        CHECK(outcome.tetHandoff.completeTetRepresentationsBeforePublication ==
              1);
        CHECK(outcome.tetHandoff.importedCells ==
              static_cast<std::uint64_t>(expected.elevate
                                             ? expected.tet10
                                             : expected.tet4));
        CHECK(outcome.tetHandoff.selectedCells ==
              outcome.tetHandoff.importedCells);
        CHECK(outcome.tetHandoff.importCandidateValidationCalls>=1);
        CHECK(outcome.tetHandoff.importGlobalMapRows==
              outcome.tetHandoff.importedNodes);
        CHECK(outcome.tetHandoff.importGlobalMapInsertions==
              outcome.tetHandoff.importedNodes);
        CHECK(outcome.tetHandoff.importGlobalMapComparisons>0);
        CHECK(outcome.tetHandoff.importEdgeValidationCellVisits==
              outcome.tetHandoff.importedCells);
        CHECK(outcome.tetHandoff.importFaceValidationCellVisits==
              outcome.tetHandoff.importedCells);
        CHECK(outcome.tetHandoff.modeledOverlapPayloadBytes >=
              outcome.tetHandoff.candidateAccountedBytes);
        CHECK(outcome.tetHandoff.modeledOverlapPayloadBytes >=
              outcome.tetHandoff.authoritativeModeledStorageBytes);
        CHECK(outcome.tetHandoff
                  .constructionCandidateConservativeAccountedBytes >=
              outcome.tetHandoff.candidateAccountedBytes);
        CHECK(outcome.tetHandoff.filterConnectivityStageAccountedBytes > 0);
        CHECK(outcome.tetHandoff.callerMappingAccountedBytes > 0);
        CHECK(outcome.tetHandoff.modeledOverlapPayloadBytes >=
              outcome.tetHandoff.publishedPayloadAccountedBytes);
        CHECK(outcome.tetHandoff.publishedPayloadAccountedBytes > 0);
        CHECK(outcome.routeVersion == algo::kConstrainedV2RouteVersion);
        CHECK(outcome.audit.inverted == 0);
        CHECK(outcome.audit.invalidTet10 == 0);
        CHECK(outcome.audit.nonFinite == 0);
        CHECK(outcome.audit.malformed == 0);
        CHECK(outcome.audit.deferred == 0);
        CHECK(outcome.audit.orphan == 0);
        CHECK(outcome.audit.unexpected == 0);
        CHECK(outcome.audit.duplicate == 0);
        CHECK(outcome.audit.nonManifold == 0);
        CHECK(outcome.audit.prescribed == expected.boundary);
        CHECK(outcome.audit.matched == expected.boundary);
        CHECK(outcome.audit.connectedComponents == 1);
        mesh::AuditOptions options;
        options.boundaryMode = mesh::BoundaryComparisonMode::Required;
        options.expectedComponents = 1;
        options.duplicateFacetsAreErrors = true;
        const auto audit = mesh::auditMesh(result, options);
        const auto semantic = mesh::semanticHashEx(result);
        CHECK(semantic.canonical);
        CHECK(semantic.value == expected.semantic);
        CHECK(audit.nInverted == 0);
        CHECK(audit.nInvalidTet10 == 0);
        CHECK(audit.boundaryOrphans == 0);
        CHECK(audit.boundaryUnexpected == 0);
        CHECK(audit.boundaryPrescribed == audit.boundaryMatched);
        if (expected.nodes >= 0) {
            CHECK(audit.nNodes == expected.nodes);
            CHECK(result.cellCount() == expected.cells);
            CHECK(audit.nTet4 == expected.tet4);
            CHECK(audit.nTet10 == expected.tet10);
            CHECK(audit.nTri3 == expected.tri3);
            CHECK(audit.nTri6 == expected.tri6);
            CHECK(audit.boundaryPrescribed == expected.boundary);
        }
        if (core::isStrictBaselineReferenceEnv()) {
            CHECK(expected.strict != 0);
            CHECK(mesh::meshHash(result) == expected.strict);
        }
        std::printf("  [v2 %s] nodes=%lld cells=%lld tet4=%lld tet10=%lld "
                    "tri3=%lld tri6=%lld boundary=%lld strict=0x%016llx "
                    "semantic=0x%016llx\n", expected.name,
                    (long long)audit.nNodes, (long long)result.cellCount(),
                    (long long)audit.nTet4, (long long)audit.nTet10,
                    (long long)audit.nTri3, (long long)audit.nTri6,
                    (long long)audit.boundaryPrescribed,
                    (unsigned long long)mesh::meshHash(result),
                    (unsigned long long)semantic.value);
    }
}

void testBoundaryPresentations() {
    TEST("constrained_v2.rotated_reversed_boundary_presentations");
    auto kernel = geo::makeBuiltinKernel();
    auto source = kernel->makeBox({0,0,0}, {1,1,1});
    PresentedModel rotated(*source, FacePresentation::Rotated);
    PresentedModel reversed(*source, FacePresentation::Reversed);
    auto spec = constantSpec(0.3);
    mesh::Mesh naturalMesh, rotatedMesh, reversedMesh;
    const auto natural = algo::meshModelConstrained(*source, spec, naturalMesh);
    const auto rot = algo::meshModelConstrained(rotated, spec, rotatedMesh);
    const auto rev = algo::meshModelConstrained(reversed, spec, reversedMesh);
    CHECK(natural.published());
    CHECK(rot.published());
    CHECK(rev.published());
    CHECK(snapshot(naturalMesh) == snapshot(rotatedMesh));
    CHECK(snapshot(naturalMesh) == snapshot(reversedMesh));
}

void testStructuredFailures() {
    TEST("constrained_v2.structured_failure_invariance");
    auto kernel = geo::makeBuiltinKernel();
    auto source = kernel->makeBox({0,0,0}, {1,1,1});

    auto run = [&](const geo::Model& model, const sizing::MeshSpec& spec) {
        mesh::Mesh destination;
        seedDestination(destination);
        const auto before = snapshot(destination);
        auto outcome = algo::meshModelConstrained(model, spec, destination);
        CHECK(!outcome.published());
        CHECK(outcome.destinationUnchanged);
        CHECK(snapshot(destination) == before);
        return outcome;
    };

    auto invalid = constantSpec(0.3);
    invalid.order = 0;
    auto invalidResult = run(*source, invalid);
    CHECK(invalidResult.status == algo::ProductionModelStatus::Rejected);
    CHECK(invalidResult.reason ==
          algo::ProductionModelReason::InvalidSizingPolicy);
    CHECK(invalidResult.phase == algo::ProductionModelPhase::Preflight);

    PresentedModel ambiguous(*source, FacePresentation::Natural,
                             FaceFailure::None, true);
    auto ambiguousResult = run(ambiguous, constantSpec(0.3));
    CHECK(ambiguousResult.reason == algo::ProductionModelReason::SurfaceStopped);
    CHECK(ambiguousResult.phase ==
          algo::ProductionModelPhase::SurfaceConstruction);
    CHECK(ambiguousResult.surface.reason ==
          algo::ProductionSurfaceReason::AmbiguousBoundaryDiscretization);

    PresentedModel badGeometry(*source, FacePresentation::Natural,
                               FaceFailure::Geometry);
    auto geometryResult = run(badGeometry, constantSpec(0.3));
    CHECK(geometryResult.reason == algo::ProductionModelReason::SurfaceStopped);
    CHECK(geometryResult.surface.reason ==
          algo::ProductionSurfaceReason::GeometryProviderFailure);
    CHECK(geometryResult.surface.phase ==
          algo::ProductionSurfacePhase::BoundaryMapping);

    PresentedModel badCandidate(*source, FacePresentation::Natural,
                                FaceFailure::Candidate);
    auto candidateResult = run(badCandidate, constantSpec(0.3));
    CHECK(candidateResult.reason == algo::ProductionModelReason::SurfaceStopped);
    CHECK(candidateResult.surface.reason ==
          algo::ProductionSurfaceReason::CandidateConstructionFailure);
    CHECK(candidateResult.surface.phase ==
          algo::ProductionSurfacePhase::ScalableConstruction);
}

void testCheckedEdgeSizingEvaluations() {
    TEST("constrained_v2.checked_edge_sizing_evaluations");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0,0,0}, {1,1,1});
    const SizingAdversary invalidValues[]{
        SizingAdversary::Nan,
        SizingAdversary::PositiveInfinity,
        SizingAdversary::Zero,
        SizingAdversary::Negative,
        SizingAdversary::Throwing};
    for (const auto kind : invalidValues)
        expectStructuredSizingStop(
            *model, adversarialSpec(kind),
            algo::ProductionModelReason::SizingEvaluationFailed,
            algo::ProductionModelPhase::EdgeConstruction, 0.3);

    expectStructuredSizingStop(
        *model, adversarialSpec(SizingAdversary::Allocation),
        algo::ProductionModelReason::AllocationFailed,
        algo::ProductionModelPhase::EdgeConstruction, 0.3);

    auto overflow = constantSpec(
        (std::numeric_limits<double>::max)() * 0.5);
    overflow.order = 2;
    overflow.orderFactor = 4.0;
    expectStructuredSizingStop(
        *model, overflow,
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::EdgeConstruction, 0.3);

    // The base value itself is finite and positive. The segment integral's
    // rounded count is not representable as int after the historical 1e-12
    // floor, so the checked route must stop before the cast/allocation.
    expectStructuredSizingStop(
        *model, adversarialSpec(SizingAdversary::Tiny),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::EdgeConstruction, 0.3);
}

void testCheckedPhaseSizingEvaluations() {
    TEST("constrained_v2.checked_surface_volume_sizing_phases");
    auto kernel = geo::makeBuiltinKernel();
    auto box = kernel->makeBox({0,0,0}, {1,1,1});

    // Edge samples are valid; a face-interior center is the first invalid
    // actual evaluation and must retain surface phase attribution.
    expectStructuredSizingStop(
        *box, adversarialSpec(SizingAdversary::FaceCenterZero),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::SurfaceConstruction, 0.3);

    // Boundary and face evaluations are valid. The solid-centre reference
    // evaluation is invalid and must stop before PointGrid/lattice setup.
    expectStructuredSizingStop(
        *box, adversarialSpec(SizingAdversary::SolidCenterZero),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::VolumeConstruction, 0.3);

    // The volume reference is valid; an actual interior-lattice query fails.
    expectStructuredSizingStop(
        *box, adversarialSpec(SizingAdversary::LatticeZero),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::VolumeConstruction, 0.3);

    auto cylinder = kernel->makeCylinder({0,0,0}, 0.5, 1.0);
    // Surface points and the solid centre are valid; the inward shield sample
    // itself is the first invalid evaluation.
    expectStructuredSizingStop(
        *cylinder, adversarialSpec(SizingAdversary::ShieldZero),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::VolumeConstruction, 0.15);
}

void testCheckedDerivedProgress() {
    TEST("constrained_v2.checked_derived_spacing_progress");
    auto kernel = geo::makeBuiltinKernel();
    // At this translation the centre-only 1e-8 size is smaller than the
    // binary-double ULP. The field is valid at every query, but lo + h == lo.
    // The checked route rejects the non-advancing lattice increment instead
    // of entering a loop.
    constexpr double origin = 1.0e9;
    auto translated = kernel->makeBox(
        {origin, origin, origin},
        {origin + 1.0, origin + 1.0, origin + 1.0});
    expectStructuredSizingStop(
        *translated, adversarialSpec(SizingAdversary::TranslatedCenterTiny),
        algo::ProductionModelReason::SizingEvaluationFailed,
        algo::ProductionModelPhase::VolumeConstruction, 0.3);
}

void testThreadDeterminism() {
    TEST("constrained_v2.thread_determinism_1_4_auto");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeCylinder({0,0,0}, 0.5, 1.0);
    std::array<MeshSnapshot, 3> snapshots;
    std::array<algo::TetHandoffMetrics, 3> handoffMetrics;
    const int threads[3]{1, 4, 0};
    for (std::size_t i = 0; i < 3; ++i) {
        mesh::Mesh result;
        const auto outcome = algo::meshModelConstrained(
            *model, constantSpec(0.15, threads[i]), result);
        CHECK(outcome.published());
        snapshots[i] = snapshot(result);
        handoffMetrics[i] = outcome.tetHandoff;
    }
    CHECK(snapshots[0] == snapshots[1]);
    CHECK(snapshots[0] == snapshots[2]);
    CHECK(handoffMetrics[0] == handoffMetrics[1]);
    CHECK(handoffMetrics[0] == handoffMetrics[2]);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* constrainedV2Tests(std::size_t& count) {
    static const Entry tests[] = {
        {"constrained_v2.manifest", testConstrainedV2Baseline},
        {"constrained_v2.presentations", testBoundaryPresentations},
        {"constrained_v2.failures", testStructuredFailures},
        {"constrained_v2.sizing_edge", testCheckedEdgeSizingEvaluations},
        {"constrained_v2.sizing_phases", testCheckedPhaseSizingEvaluations},
        {"constrained_v2.sizing_progress", testCheckedDerivedProgress},
        {"constrained_v2.determinism", testThreadDeterminism},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
