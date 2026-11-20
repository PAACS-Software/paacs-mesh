// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <type_traits>

#include "paacs/algo/SurfaceMesherAdapter.hpp"

namespace paacs::algo {

// Plain-value, allocation-free preparation of the canonical parametric point
// used by the bounded facet-conflict refinement proposal.  The caller resolves
// the normalized facet ordinal and supplies exactly its three retained rows;
// this helper never scans a SurfaceFaceResult or calls geometry.
enum class FacetBarycentreStatus : std::uint8_t {
    Ready,
    InvalidFaceIdentity,
    InvalidVertexIdentity,
    DuplicateVertexIdentity,
    NonfiniteParameter,
    ChartMismatch,
    UnsupportedChart,
    InvalidChartBoundary,
    OutsideCertifiedRegion,
    DegenerateFacet,
    ArithmeticRange,
    InvalidCertificateAuthority,
    WorkBudgetReached,
    WorkAccountingOverflow
};

struct FacetBarycentreVertex {
    Index publishedGlobalNode = kInvalidIndex;
    bnd::ChartId chart = 0;
    Vec2 parameter{};
    friend bool operator==(const FacetBarycentreVertex& a,
                           const FacetBarycentreVertex& b) noexcept {
        return a.publishedGlobalNode==b.publishedGlobalNode&&
               a.chart==b.chart&&a.parameter.u==b.parameter.u&&
               a.parameter.v==b.parameter.v;
    }
};

struct FacetBarycentreRequest {
    ArtifactStamp faceStamp{};
    bnd::ChartId chart = 0;
    bnd::GeometryProviderStamp expectedProviderStamp{};
    ArtifactStamp expectedCertificateStamp{};
    std::array<FacetBarycentreVertex,3> vertices{};
};

// Allocation-free view of the retained immutable certificate.  The view is
// valid only while the retained SurfaceFaceResult that owns its boundary rows
// remains alive.
struct FacetBarycentreCertificateView {
    bnd::GeometryProviderStamp providerStamp{};
    bnd::ChartId chart = 0;
    ArtifactStamp certificateStamp{};
    std::span<const Vec2> certifiedRegionBoundary{};
    int chartPositiveToNormalSign = 1;
    bool regular = true;
    bool periodic = false;
    bool seam = false;
};

struct FacetBarycentreResult {
    FacetBarycentreStatus status = FacetBarycentreStatus::InvalidFaceIdentity;
    std::array<FacetBarycentreVertex,3> canonicalVertices{};
    Vec2 parameter{};
    SurfaceNodeIdentity identity{};
    std::uint64_t workUnits = 0;

    bool ready() const noexcept {
        return status == FacetBarycentreStatus::Ready;
    }
};

// Canonical order is ascending publishedGlobalNode.  Each component is
// evaluated in that order with three rounded multiplications by the binary64
// constant nearest 1/3 followed by two rounded additions. Negative zero is
// normalized before evaluation and in the result.  This translation unit is
// compiled with floating-point contraction disabled in every canonical build.
FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    const ChartCertificate& certificate) noexcept;

FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    const ChartCertificate& certificate,
    std::uint64_t maxWorkUnits) noexcept;

FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    FacetBarycentreCertificateView certificate,
    std::uint64_t maxWorkUnits) noexcept;

static_assert(std::is_standard_layout_v<FacetBarycentreVertex>);
static_assert(std::is_trivially_copyable_v<FacetBarycentreVertex>);
static_assert(std::is_standard_layout_v<FacetBarycentreRequest>);
static_assert(std::is_trivially_copyable_v<FacetBarycentreRequest>);
static_assert(std::is_standard_layout_v<FacetBarycentreResult>);
static_assert(std::is_trivially_copyable_v<FacetBarycentreResult>);

}  // namespace paacs::algo
