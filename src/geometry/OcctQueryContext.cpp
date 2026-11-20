// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctQueryContext.hpp"
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>
#include <array>
#include <bit>
#include <functional>

namespace paacs::geo::detail {
namespace {
thread_local OcctQueryContext* current = nullptr;
using PointBits = std::array<std::uint64_t,3>;
PointBits pointBits(Vec3 point) {
    return {std::bit_cast<std::uint64_t>(point.x), std::bit_cast<std::uint64_t>(point.y),
            std::bit_cast<std::uint64_t>(point.z)};
}
}
struct OcctQueryContext::Impl {
    struct Projection {
        TopoDS_Face face;
        PointBits sample{};
        Vec3 closest{};
    };
    static constexpr std::size_t projectionSlots = 32'768;
    explicit Impl(bool memoize)
        : projections(memoize ? std::make_unique<Projection[]>(projectionSlots) : nullptr) {}

    Projection& projection(const TopoDS_Face& face, const PointBits& sample) {
        std::uint64_t hash = std::hash<TopoDS_Shape>{}(face) ^ std::uint64_t(face.Orientation());
        for (const auto coordinate : sample) {
            hash ^= coordinate;
            hash = std::rotl(hash, 25) * 0x9E3779B97F4A7C15ULL;
        }
        return projections[(hash ^ (hash >> 32)) & (projectionSlots - 1)];
    }

    occ::handle<Geom_Surface> surface;
    std::array<std::uint64_t,4> bounds{};
    GeomAPI_ProjectPointOnSurf projector;
    TopoDS_Face distanceTarget;
    BRepExtrema_DistShapeShape distance;
    std::unique_ptr<Projection[]> projections;
    OcctQueryCounters counters;
    bool initialized = false;
};
OcctQueryContext::OcctQueryContext(bool enabled, bool memoizeTrimmedFacePoints)
    : impl_(enabled ? std::make_unique<Impl>(memoizeTrimmedFacePoints) : nullptr), previous_(current) {
    current = enabled ? this : nullptr;
}
OcctQueryContext::~OcctQueryContext() { current = previous_; }
OcctQueryCounters OcctQueryContext::counters() const noexcept {
    return impl_ ? impl_->counters : OcctQueryCounters{};
}
std::uint64_t OcctQueryContext::maximumRetainedPayloadBytes(bool memoizeTrimmedFacePoints) noexcept {
    const auto projectionBytes = memoizeTrimmedFacePoints ? Impl::projectionSlots * sizeof(Impl::Projection) : 0;
    return sizeof(Impl) + sizeof(OcctQueryContext) + projectionBytes;
}
OcctQueryContext* OcctQueryContext::active() noexcept { return current; }
GeomAPI_ProjectPointOnSurf& OcctQueryContext::surfaceProjector(
    const occ::handle<Geom_Surface>& surface,double u0,double u1,double v0,double v1,
    const gp_Pnt& point) {
    auto& state=*impl_;++state.counters.surfaceProjectionQueries;
    const std::array<std::uint64_t,4> bounds{std::bit_cast<std::uint64_t>(u0),
        std::bit_cast<std::uint64_t>(u1),std::bit_cast<std::uint64_t>(v0),
        std::bit_cast<std::uint64_t>(v1)};
    if(!state.initialized||state.surface!=surface||state.bounds!=bounds) {
        state.initialized=false;++state.counters.surfaceProjectionInitializations;
        // Identical overload/default tolerance and Grad solver to fresh queries.
        state.projector.Init(surface,u0,u1,v0,v1);
        state.surface=surface;state.bounds=bounds;state.initialized=true;
    } else ++state.counters.surfaceProjectionReuses;
    state.projector.Perform(point);
    return state.projector;
}
BRepExtrema_DistShapeShape& OcctQueryContext::trimmedFaceDistance(
    const TopoDS_Face& face, const gp_Pnt& point) {
    auto& state = *impl_;
    try {
        if (!state.distanceTarget.IsEqual(face)) {
            state.distance.LoadS2(face);
            state.distanceTarget = face;
        }
        state.distance.LoadS1(BRepBuilderAPI_MakeVertex(point).Vertex());
        state.distance.Perform();
    } catch (...) {
        state.distanceTarget.Nullify();
        throw;
    }
    return state.distance;
}
bool OcctQueryContext::findTrimmedFacePoint(const TopoDS_Face& face, Vec3 sample, Vec3& closest) {
    auto& state = *impl_;
    ++state.counters.trimmedProjectionQueries;
    if (!state.projections) return false;
    const auto key = pointBits(sample);
    const auto& entry = state.projection(face, key);
    if (entry.face.IsNull() || !entry.face.IsEqual(face) || entry.sample != key) return false;
    closest = entry.closest;
    ++state.counters.trimmedProjectionMemoHits;
    return true;
}
void OcctQueryContext::rememberTrimmedFacePoint(const TopoDS_Face& face, Vec3 sample, Vec3 closest) {
    auto& state = *impl_;
    if (!state.projections) return;
    const auto key = pointBits(sample);
    state.projection(face, key) = {face, key, closest};
}
} // namespace paacs::geo::detail
