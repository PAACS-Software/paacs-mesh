// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional-target access; never installed with the public API.
#include "paacs/geometry/OcctModel.hpp"
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
namespace paacs::geo::detail {
struct OcctModelAccess {
    static const TopoDS_Shape &shape(const OcctModel &model) noexcept;
    static TopoDS_Solid solid(const OcctModel &model, int tag);
    static TopoDS_Face face(const OcctModel &model, int tag);
    static TopoDS_Edge edge(const OcctModel &model, int tag);
    // On the same imported pcurve, fit the original 3D edge point only
    // between the neighboring planned samples. Leaves uv unchanged on refusal.
    static bool fitParameterOnUse(const OcctModel&, std::uint64_t useId,
                                  double previous, double at, double next, Vec2& uv);
    // Private derived-geometry snapshot. The caller owns the shape and must
    // cease mutating its shared OCCT handles after successful construction.
    // No repair, file transfer, coordinate conversion or identity welding.
    static NativeCadImportOutcome fromShape(const TopoDS_Shape& ownedShape,
                                            const NativeCadImportOptions& options = {}) noexcept;
};
} // namespace paacs::geo::detail
