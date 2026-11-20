// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private per-call optional query state; not installed with the public API.
#include "paacs/core/Types.hpp"
#include <Standard_Handle.hxx>
#include <cstdint>
#include <memory>
class Geom_Surface;
class GeomAPI_ProjectPointOnSurf;
class BRepExtrema_DistShapeShape;
class TopoDS_Face;
class gp_Pnt;

namespace paacs::geo::detail {
struct OcctQueryCounters {
    std::uint64_t surfaceProjectionQueries = 0;
    std::uint64_t surfaceProjectionInitializations = 0;
    std::uint64_t surfaceProjectionReuses = 0;
    std::uint64_t trimmedProjectionQueries = 0;
    std::uint64_t trimmedProjectionMemoHits = 0;
};
// Scope must be created/destroyed on the querying thread. Nested scopes restore
// their caller's context. Disabled scopes use the unchanged fresh-query path.
// Retains one surface projector and one trimmed-face distance tool, released
// at scope exit. An enabled scope can also memoize exact trimmed-face results.
// Native solver allocations remain outside the logical payload.
class OcctQueryContext final {
public:
    explicit OcctQueryContext(bool reuseSurfaceProjection = false,
                              bool memoizeTrimmedFacePoints = false);
    ~OcctQueryContext();
    OcctQueryContext(const OcctQueryContext&) = delete;
    OcctQueryContext& operator=(const OcctQueryContext&) = delete;
    OcctQueryCounters counters() const noexcept;
    static std::uint64_t maximumRetainedPayloadBytes(bool memoizeTrimmedFacePoints = false) noexcept;
    static OcctQueryContext* active() noexcept;
    GeomAPI_ProjectPointOnSurf& surfaceProjector(const occ::handle<Geom_Surface>&,
        double u0, double u1, double v0, double v1, const gp_Pnt&);
    BRepExtrema_DistShapeShape& trimmedFaceDistance(const TopoDS_Face&, const gp_Pnt&);
    bool findTrimmedFacePoint(const TopoDS_Face&, Vec3 sample, Vec3& closest);
    void rememberTrimmedFacePoint(const TopoDS_Face&, Vec3 sample, Vec3 closest);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    OcctQueryContext* previous_ = nullptr;
};
} // namespace paacs::geo::detail
