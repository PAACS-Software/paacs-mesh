// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Included inside bindings.cpp's anonymous namespace after CadModel. Its
// public geometry include belongs at the top of that translation unit.
py::dict previewCad(const CadModel& holder, const py::object& deflection,
                    const py::object& maximumPoints, const py::object& maximumTriangles) {
    geo::OcctPreviewOptions options;
    options.linearDeflection = positiveNumber(deflection, "linear_deflection");
    options.maximumPoints = positiveInteger(maximumPoints, "max_points");
    options.maximumTriangles = positiveInteger(maximumTriangles, "max_triangles");
    if (options.maximumPoints > 250000 || options.maximumTriangles > 500000)
        throw py::value_error("Preview limits are 250000 points and 500000 triangles");
    geo::OcctPreview display;
    geo::OcctPreviewOutcome outcome;
    { py::gil_scoped_release release;
      std::lock_guard lock(nativeOperationMutex);
      outcome = geo::previewOcctCad(*holder.imported.model, options, display);
    }
    if (!outcome.completed)
        throw std::runtime_error(std::string("CAD display: ") + geo::occtPreviewReasonName(outcome.reason));
    py::dict result;
    py::list points;
    for (const auto& point : display.points) points.append(py::make_tuple(point.x, point.y, point.z));
    result["points"] = std::move(points);
    result["triangles"] = display.triangles;
    result["face_tags"] = display.faceTags;
    result["linear_deflection"] = options.linearDeflection;
    result["angular_deflection_radians"] = options.angularDeflection;
    return result;
}
