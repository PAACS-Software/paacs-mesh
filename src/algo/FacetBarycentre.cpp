// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/FacetBarycentre.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

#include "paacs/predicates/Predicates.hpp"

namespace paacs::algo {
namespace {

static_assert(std::numeric_limits<double>::is_iec559);

constexpr double kOneThird = 0x1.5555555555555p-2;

enum class WorkState : std::uint8_t { Ready, Budget, Overflow };
enum class ClosureState : std::uint8_t { Outside, Inside, WorkStopped };

struct WorkCounter {
    std::uint64_t used = 0;
    std::uint64_t limit = 0;
    WorkState state = WorkState::Ready;

    bool take() noexcept {
        if (used == std::numeric_limits<std::uint64_t>::max()) {
            state = WorkState::Overflow;
            return false;
        }
        if (used >= limit) {
            state = WorkState::Budget;
            return false;
        }
        ++used;
        return true;
    }
};

bool finite(const Vec2& value) noexcept {
    return std::isfinite(value.u) && std::isfinite(value.v);
}

bool validProviderStamp(const bnd::GeometryProviderStamp& stamp) noexcept {
    return stamp.identity != 0 && stamp.revision != 0;
}

bool validArtifactStamp(const ArtifactStamp& stamp) noexcept {
    return stamp.persistentId != 0 && stamp.revision != 0;
}

double canonicalZero(double value) noexcept {
    return value == 0.0 ? 0.0 : value;
}

ClosureState pointOnSegment(const Vec2& a, const Vec2& b,
                            const Vec2& point,
                            WorkCounter& work) noexcept {
    if (!work.take()) return ClosureState::WorkStopped;
    if (predicates::orient2d(a, b, point) != predicates::Sign::Zero)
        return ClosureState::Outside;
    if (!work.take()) return ClosureState::WorkStopped;
    const double minU = (std::min)(a.u, b.u);
    if (!work.take()) return ClosureState::WorkStopped;
    const double maxU = (std::max)(a.u, b.u);
    if (!work.take()) return ClosureState::WorkStopped;
    const double minV = (std::min)(a.v, b.v);
    if (!work.take()) return ClosureState::WorkStopped;
    const double maxV = (std::max)(a.v, b.v);
    if (!work.take()) return ClosureState::WorkStopped;
    const bool aboveU = point.u >= minU;
    if (!work.take()) return ClosureState::WorkStopped;
    const bool belowU = point.u <= maxU;
    if (!work.take()) return ClosureState::WorkStopped;
    const bool aboveV = point.v >= minV;
    if (!work.take()) return ClosureState::WorkStopped;
    const bool belowV = point.v <= maxV;
    const bool on = aboveU && belowU && aboveV && belowV;
    return on ? ClosureState::Inside : ClosureState::Outside;
}

ClosureState pointInClosure(std::span<const Vec2> ring,
                            const Vec2& point,
                            WorkCounter& work) noexcept {
    if (ring.size() < 3) return ClosureState::Outside;
    bool inside = false;
    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size();
         j = i++) {
        const auto& a = ring[j];
        const auto& b = ring[i];
        const ClosureState segment = pointOnSegment(a, b, point, work);
        if (segment == ClosureState::WorkStopped) return segment;
        if (segment == ClosureState::Inside) return segment;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool aAbove = a.v > point.v;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool bAbove = b.v > point.v;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool crosses = aAbove != bAbove;
        if (!crosses) continue;
        if (!work.take()) return ClosureState::WorkStopped;
        const auto sign = predicates::orient2d(a, b, point);
        if (!work.take()) return ClosureState::WorkStopped;
        const bool ascending = b.v > a.v;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool descending = b.v < a.v;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool positive = sign == predicates::Sign::Positive;
        if (!work.take()) return ClosureState::WorkStopped;
        const bool negative = sign == predicates::Sign::Negative;
        if ((ascending && positive) || (descending && negative))
            inside = !inside;
    }
    return inside ? ClosureState::Inside : ClosureState::Outside;
}

bool barycentreComponent(double a, double b, double c,
                         WorkCounter& work, double& result) noexcept {
    if (!work.take()) return false;
    a = canonicalZero(a);
    if (!work.take()) return false;
    b = canonicalZero(b);
    if (!work.take()) return false;
    c = canonicalZero(c);
    if (!work.take()) return false;
    const double scaledA = a * kOneThird;
    if (!work.take()) return false;
    const double scaledB = b * kOneThird;
    if (!work.take()) return false;
    const double sum01 = scaledA + scaledB;
    if (!work.take()) return false;
    const double scaledC = c * kOneThird;
    if (!work.take()) return false;
    const double sum = sum01 + scaledC;
    if (!work.take()) return false;
    result = canonicalZero(sum);
    return true;
}

void applyWorkFailure(FacetBarycentreResult& result,
                      const WorkCounter& work) noexcept {
    result.status = work.state == WorkState::Budget
        ? FacetBarycentreStatus::WorkBudgetReached
        : FacetBarycentreStatus::WorkAccountingOverflow;
    result.workUnits = work.used;
}

}  // namespace

FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    FacetBarycentreCertificateView certificate,
    std::uint64_t maxWorkUnits) noexcept {
    FacetBarycentreResult result;
    WorkCounter work{0, maxWorkUnits, WorkState::Ready};
    result.canonicalVertices = request.vertices;
    for (const auto pair :
         std::array<std::array<std::size_t, 2>, 3>{{
             {{0, 1}}, {{1, 2}}, {{0, 1}}}}) {
        if (!work.take()) {
            applyWorkFailure(result, work);
            return result;
        }
        if (result.canonicalVertices[pair[1]].publishedGlobalNode <
            result.canonicalVertices[pair[0]].publishedGlobalNode)
            std::swap(result.canonicalVertices[pair[0]],
                      result.canonicalVertices[pair[1]]);
    }

    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!validArtifactStamp(request.faceStamp)) return result;

    for (const auto& vertex : result.canonicalVertices) {
        if (!work.take()) {
            applyWorkFailure(result, work);
            return result;
        }
        if (vertex.publishedGlobalNode < 0) {
            result.status = FacetBarycentreStatus::InvalidVertexIdentity;
            result.workUnits = work.used;
            return result;
        }
        if (!finite(vertex.parameter)) {
            result.status = FacetBarycentreStatus::NonfiniteParameter;
            result.workUnits = work.used;
            return result;
        }
        if (vertex.chart != request.chart) {
            result.status = FacetBarycentreStatus::ChartMismatch;
            result.workUnits = work.used;
            return result;
        }
    }
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (result.canonicalVertices[0].publishedGlobalNode ==
            result.canonicalVertices[1].publishedGlobalNode ||
        result.canonicalVertices[1].publishedGlobalNode ==
            result.canonicalVertices[2].publishedGlobalNode) {
        result.status = FacetBarycentreStatus::DuplicateVertexIdentity;
        result.workUnits = work.used;
        return result;
    }

    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!validProviderStamp(request.expectedProviderStamp) ||
        !validArtifactStamp(request.expectedCertificateStamp) ||
        !validProviderStamp(certificate.providerStamp) ||
        !validArtifactStamp(certificate.certificateStamp) ||
        !(certificate.providerStamp == request.expectedProviderStamp) ||
        !(certificate.certificateStamp ==
          request.expectedCertificateStamp)) {
        result.status = FacetBarycentreStatus::InvalidCertificateAuthority;
        result.workUnits = work.used;
        return result;
    }
    if (certificate.chart != request.chart) {
        result.status = FacetBarycentreStatus::ChartMismatch;
        result.workUnits = work.used;
        return result;
    }
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!certificate.regular || certificate.periodic || certificate.seam) {
        result.status = FacetBarycentreStatus::UnsupportedChart;
        result.workUnits = work.used;
        return result;
    }
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if ((certificate.chartPositiveToNormalSign != 1 &&
         certificate.chartPositiveToNormalSign != -1) ||
        certificate.certifiedRegionBoundary.size() < 3) {
        result.status = FacetBarycentreStatus::InvalidChartBoundary;
        result.workUnits = work.used;
        return result;
    }
    for (const auto& point : certificate.certifiedRegionBoundary) {
        if (!work.take()) {
            applyWorkFailure(result, work);
            return result;
        }
        if (!finite(point)) {
            result.status = FacetBarycentreStatus::InvalidChartBoundary;
            result.workUnits = work.used;
            return result;
        }
    }
    for (const auto& vertex : result.canonicalVertices) {
        const ClosureState closure = pointInClosure(
            certificate.certifiedRegionBoundary, vertex.parameter, work);
        if (closure == ClosureState::WorkStopped) {
            applyWorkFailure(result, work);
            return result;
        }
        if (closure == ClosureState::Outside) {
            result.status = FacetBarycentreStatus::OutsideCertifiedRegion;
            result.workUnits = work.used;
            return result;
        }
    }
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (predicates::orient2d(result.canonicalVertices[0].parameter,
                             result.canonicalVertices[1].parameter,
                             result.canonicalVertices[2].parameter) ==
        predicates::Sign::Zero) {
        result.status = FacetBarycentreStatus::DegenerateFacet;
        result.workUnits = work.used;
        return result;
    }

    if (!barycentreComponent(result.canonicalVertices[0].parameter.u,
                             result.canonicalVertices[1].parameter.u,
                             result.canonicalVertices[2].parameter.u,
                             work, result.parameter.u) ||
        !barycentreComponent(result.canonicalVertices[0].parameter.v,
                             result.canonicalVertices[1].parameter.v,
                             result.canonicalVertices[2].parameter.v,
                             work, result.parameter.v)) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!finite(result.parameter)) {
        result.status = FacetBarycentreStatus::ArithmeticRange;
        result.workUnits = work.used;
        return result;
    }
    const ClosureState resultClosure = pointInClosure(
        certificate.certifiedRegionBoundary, result.parameter, work);
    if (resultClosure == ClosureState::WorkStopped) {
        applyWorkFailure(result, work);
        return result;
    }
    if (resultClosure == ClosureState::Outside) {
        result.status = FacetBarycentreStatus::OutsideCertifiedRegion;
        result.workUnits = work.used;
        return result;
    }

    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    result.identity.domain =
        SurfaceNodeIdentityDomain::RefinedFacetInterior;
    result.identity.parentArtifact = request.faceStamp;
    for (std::size_t i = 0; i < result.canonicalVertices.size(); ++i)
        result.identity.parentFacetGlobalNodes[i] =
            result.canonicalVertices[i].publishedGlobalNode;
    result.identity.parameterUBits =
        std::bit_cast<std::uint64_t>(result.parameter.u);
    result.identity.parameterVBits =
        std::bit_cast<std::uint64_t>(result.parameter.v);
    if (!work.take()) {
        applyWorkFailure(result, work);
        return result;
    }
    if (!validSurfaceNodeIdentity(result.identity)) {
        result.status = FacetBarycentreStatus::InvalidVertexIdentity;
        result.workUnits = work.used;
        return result;
    }
    result.status = FacetBarycentreStatus::Ready;
    result.workUnits = work.used;
    return result;
}

FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    const ChartCertificate& certificate,
    std::uint64_t maxWorkUnits) noexcept {
    return computeFacetBarycentre(
        request,
        FacetBarycentreCertificateView{
            certificate.providerStamp,
            certificate.chart,
            certificate.certificateStamp,
            certificate.certifiedRegionBoundary,
            certificate.chartPositiveToNormalSign,
            certificate.regular,
            certificate.periodic,
            certificate.seam},
        maxWorkUnits);
}

FacetBarycentreResult computeFacetBarycentre(
    const FacetBarycentreRequest& request,
    const ChartCertificate& certificate) noexcept {
    return computeFacetBarycentre(
        request, certificate,
        std::numeric_limits<std::uint64_t>::max());
}

}  // namespace paacs::algo
