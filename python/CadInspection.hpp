// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Included inside bindings.cpp's anonymous namespace after CadModel.
// All includes belong at the top of bindings.cpp; this helper uses only the
// public geometry API and copies DTOs before constructing Python objects.
namespace cad_inspection_detail {
struct Budget {
    std::uint64_t maximum, used = 0;
    void consume(std::size_t count) {
        if (count > maximum - used)
            throw std::length_error("CAD inspection topology exceeds max_uses; request a smaller page");
        used += static_cast<std::uint64_t>(count);
    }
};
struct EntityDto {
    int tag = -1, dimension = -1;
    double tolerance = 0;
    std::vector<std::string> names;
    const char* nameStatus = "missing";
};
struct FaceDto {
    EntityDto entity;
    geo::OcctFaceProperties properties;
    std::optional<geo::OcctSurfaceDescriptor> descriptor;
    std::vector<geo::OcctFaceOwner> owners;
    std::vector<int> edgeTags, adjacentFaces;
    std::size_t wireCount = 0, holeCount = 0;
    std::optional<std::vector<geo::OcctWireUse>> wires;
};
struct BodyDto {
    EntityDto entity;
    geo::OcctSolidProperties properties;
    std::vector<geo::OcctFaceUse> faces;
};
struct EdgeDto {
    EntityDto entity;
    geo::CurveKind curve = geo::CurveKind::Unknown;
    int startVertex = -1, endVertex = -1;
    double length = 0, first = 0, last = 0;
    bool periodic = false;
    std::vector<int> faceTags;
};
struct VertexDto {
    EntityDto entity;
    Vec3 position;
};
struct Snapshot {
    std::vector<FaceDto> faces;
    std::vector<BodyDto> bodies;
    std::vector<EdgeDto> edges;
    std::vector<VertexDto> vertices;
};

const char* orientationName(geo::CadOrientation orientation) noexcept {
    switch (orientation) {
    case geo::CadOrientation::Forward: return "forward";
    case geo::CadOrientation::Reversed: return "reversed";
    case geo::CadOrientation::Internal: return "internal";
    case geo::CadOrientation::External: return "external";
    }
    return "unknown";
}
const char* curveName(geo::CurveKind kind) noexcept {
    switch (kind) {
    case geo::CurveKind::Unknown: return "unknown";
    case geo::CurveKind::Line: return "line";
    case geo::CurveKind::Circle: return "circle";
    case geo::CurveKind::Other: return "other";
    }
    return "unknown";
}
bool finitePoint(const Vec3& point) noexcept {
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
}
EntityDto entityDto(const geo::OcctModel& model, int dimension, int tag,
                    const geo::OcctSourceNameSummary& summary) {
    EntityDto result;
    result.tag = tag;
    result.dimension = dimension;
    result.tolerance = model.tolerance(dimension, tag);
    const auto names = model.sourceNames(dimension, tag);
    result.names.assign(names.begin(), names.end());
    result.nameStatus = !summary.requested ? "disabled" : !summary.complete ? "incomplete"
                      : names.empty() ? "missing" : names.size() == 1 ? "unique" : "ambiguous";
    return result;
}
void validateTags(const geo::OcctModel& model, const std::string& kind,
                  const std::vector<int>& tags) {
    if (kind != "faces" && kind != "bodies" && kind != "edges" && kind != "vertices")
        throw std::invalid_argument("kind must be faces, bodies, edges or vertices");
    if (tags.size() > 1000)
        throw std::invalid_argument("CAD inspection accepts at most 1000 entity tags per page");
    std::set<int> seen;
    for (const int tag : tags) {
        if (tag <= 0 || !seen.insert(tag).second)
            throw std::invalid_argument("CAD inspection tags must be positive and unique");
        const bool exists = kind == "faces" ? model.face(tag) != nullptr
                          : kind == "bodies" ? model.solid(tag) != nullptr
                          : kind == "edges" ? model.edge(tag) != nullptr
                                            : model.vertex(tag) != nullptr;
        if (!exists)
            throw std::invalid_argument("Unknown CAD " + kind + " tag " + std::to_string(tag));
    }
}
Snapshot collect(const geo::OcctModel& model, const std::string& kind,
                  const std::vector<int>& tags, bool topology, std::uint64_t maximumUses) {
    validateTags(model, kind, tags);
    Snapshot result;
    Budget budget{maximumUses};
    const auto summary = model.sourceNameSummary();
    if (kind == "faces") {
        result.faces.reserve(tags.size());
        for (const int tag : tags) {
            FaceDto record;
            record.entity = entityDto(model, 2, tag, summary);
            const auto properties = model.faceProperties(tag);
            if (!properties) throw std::runtime_error("Native face inspection omitted a valid tag");
            record.properties = *properties;
            if (const auto* descriptor = model.surfaceDescriptor(tag)) record.descriptor = *descriptor;
            const auto owners = model.faceOwners(tag);
            budget.consume(owners.size());
            record.owners.assign(owners.begin(), owners.end());
            const auto wires = model.boundaryUses(tag);
            record.wireCount = wires.size();
            if (topology) record.wires.emplace();
            std::set<int> edgeTags;
            for (const auto& wire : wires) {
                if (!wire.outer) ++record.holeCount;
                budget.consume(wire.edges.size());
                for (const auto& edge : wire.edges) edgeTags.insert(edge.edgeTag);
                if (record.wires) record.wires->push_back(wire);
            }
            record.edgeTags.assign(edgeTags.begin(), edgeTags.end());
            std::set<int> neighbors;
            for (const int edgeTag : record.edgeTags) {
                // Each edge's own-face lookup is already bounded by the
                // selected wire uses. Count every other face incidence,
                // including duplicates across different shared edges.
                for (const int neighbor : model.edgeFaces(edgeTag)) {
                    if (neighbor == tag) continue;
                    budget.consume(1);
                    neighbors.insert(neighbor);
                }
            }
            record.adjacentFaces.assign(neighbors.begin(), neighbors.end());
            result.faces.push_back(std::move(record));
        }
    } else if (kind == "bodies") {
        result.bodies.reserve(tags.size());
        for (const int tag : tags) {
            BodyDto record;
            record.entity = entityDto(model, 3, tag, summary);
            const auto properties = model.solidProperties(tag);
            if (!properties) throw std::runtime_error("Native body inspection omitted a valid tag");
            record.properties = *properties;
            const auto uses = model.faceUses(tag);
            budget.consume(uses.size());
            record.faces.assign(uses.begin(), uses.end());
            result.bodies.push_back(std::move(record));
        }
    } else if (kind == "edges") {
        result.edges.reserve(tags.size());
        for (const int tag : tags) {
            EdgeDto record;
            record.entity = entityDto(model, 1, tag, summary);
            const auto* edge = model.edge(tag);
            record.curve = edge->kind();
            record.length = edge->length();
            record.startVertex = edge->startVertexTag();
            record.endVertex = edge->endVertexTag();
            edge->parameterBounds(record.first, record.last);
            record.periodic = edge->periodic();
            const auto faces = model.edgeFaces(tag);
            budget.consume(faces.size());
            record.faceTags.assign(faces.begin(), faces.end());
            result.edges.push_back(std::move(record));
        }
    } else {
        result.vertices.reserve(tags.size());
        for (const int tag : tags) {
            VertexDto record;
            record.entity = entityDto(model, 0, tag, summary);
            record.position = model.vertex(tag)->position();
            result.vertices.push_back(std::move(record));
        }
    }
    return result;
}

py::object pointObject(const Vec3& point) {
    return finitePoint(point) ? py::cast(std::array{point.x, point.y, point.z}) : py::object(py::none());
}
py::object pointObject(const std::optional<Vec3>& point) {
    return point ? pointObject(*point) : py::object(py::none());
}
py::object boundsObject(const geo::OcctBounds& bounds) {
    if (!finitePoint(bounds.lower) || !finitePoint(bounds.upper)) return py::none();
    py::dict result;
    result["lower"] = pointObject(bounds.lower);
    result["upper"] = pointObject(bounds.upper);
    return result;
}
py::object boundsObject(const std::optional<geo::OcctBounds>& bounds) {
    return bounds ? boundsObject(*bounds) : py::object(py::none());
}
py::object scalarObject(const std::optional<double>& value) {
    return value ? finiteValue(*value) : py::object(py::none());
}
py::dict entityObject(const EntityDto& value) {
    py::dict result;
    result["tag"] = value.tag;
    result["dimension"] = value.dimension;
    result["tolerance"] = finiteValue(value.tolerance);
    result["names"] = value.names;
    result["name_status"] = value.nameStatus;
    result["name"] = std::string_view(value.nameStatus) == "unique"
        ? py::cast(value.names.front()) : py::object(py::none());
    return result;
}
py::dict summaryObject(const geo::OcctSourceNameSummary& value) {
    py::dict result, unmatched;
    result["scope"] = "direct_shape_representation_items";
    result["requested"] = value.requested;
    result["complete"] = value.complete;
    result["reason"] = geo::occtSourceNameReasonName(value.reason);
    result["mappings_visited"] = value.mappingsVisited;
    result["named_mappings"] = value.namedMappings;
    result["unmatched_named_mappings"] = value.unmatchedNamedMappings;
    result["unsupported_named_mappings"] = value.unsupportedNamedMappings;
    result["candidate_names"] = value.candidateNames;
    result["candidate_bytes"] = value.candidateBytes;
    result["named_entities"] = value.namedEntities;
    result["ambiguous_entities"] = value.ambiguousEntities;
    const std::array labels{"vertices", "edges", "faces", "bodies"};
    for (std::size_t dimension = 0; dimension < labels.size(); ++dimension)
        unmatched[labels[dimension]] = value.unmatchedByDimension[dimension];
    result["unmatched_by_dimension"] = unmatched;
    return result;
}
py::dict faceObject(const FaceDto& value) {
    auto result = entityObject(value.entity);
    const auto& properties = value.properties;
    result["surface"] = geo::occtSurfaceTypeName(properties.surfaceType);
    result["area"] = scalarObject(properties.area);
    result["area_centroid"] = pointObject(properties.areaCentroid);
    result["area_relative_error_estimate"] = scalarObject(properties.areaRelativeErrorEstimate);
    result["bounds"] = boundsObject(properties.bounds);
    result["metrics_reason"] = geo::occtPropertyReasonName(properties.reason);
    result["normal"] = pointObject(properties.planeReferenceNormal);
    result["normal_kind"] = properties.planeReferenceNormal ? "plane_reference" : "unavailable";
    py::list owners;
    for (const auto& valueOwner : value.owners) {
        py::dict owner;
        owner["body_tag"] = valueOwner.solidTag;
        owner["sense"] = valueOwner.sense;
        owner["outward_normal"] = properties.planeReferenceNormal && valueOwner.sense != 0
            ? pointObject(*properties.planeReferenceNormal * static_cast<double>(valueOwner.sense))
            : py::object(py::none());
        owners.append(std::move(owner));
    }
    result["owners"] = owners;
    py::dict geometry;
    if (value.descriptor) {
        const auto& descriptor = *value.descriptor;
        const bool analytic = descriptor.kind == geo::SurfaceKind::Plane ||
            descriptor.kind == geo::SurfaceKind::Cylinder || descriptor.kind == geo::SurfaceKind::Cone ||
            descriptor.kind == geo::SurfaceKind::Sphere || descriptor.kind == geo::SurfaceKind::Torus;
        geometry["origin"] = analytic ? pointObject(descriptor.origin) : py::object(py::none());
        geometry["axis"] = analytic ? pointObject(descriptor.axis) : py::object(py::none());
        geometry["x_direction"] = analytic ? pointObject(descriptor.xDirection) : py::object(py::none());
        geometry["radius"] = descriptor.kind == geo::SurfaceKind::Cylinder || descriptor.kind == geo::SurfaceKind::Sphere
            ? finiteValue(descriptor.radius) : py::object(py::none());
        geometry["reference_radius"] = descriptor.kind == geo::SurfaceKind::Cone
            ? finiteValue(descriptor.radius) : py::object(py::none());
        geometry["semi_angle"] = descriptor.kind == geo::SurfaceKind::Cone
            ? finiteValue(descriptor.semiAngle) : py::object(py::none());
        geometry["minor_radius"] = descriptor.kind == geo::SurfaceKind::Torus
            ? finiteValue(descriptor.radius) : py::object(py::none());
        geometry["major_radius"] = descriptor.kind == geo::SurfaceKind::Torus
            ? finiteValue(descriptor.majorRadius) : py::object(py::none());
        result["parameter_bounds"] = py::make_tuple(
            std::array{descriptor.parameterLower.u, descriptor.parameterLower.v},
            std::array{descriptor.parameterUpper.u, descriptor.parameterUpper.v});
        result["periodic_u"] = descriptor.periodicU;
        result["periodic_v"] = descriptor.periodicV;
        result["orientation"] = orientationName(descriptor.orientation);
    }
    result["geometry"] = geometry;
    result["edge_tags"] = value.edgeTags;
    result["wire_count"] = value.wireCount;
    result["hole_count"] = value.holeCount;
    result["adjacent_faces"] = value.adjacentFaces;
    if (value.wires) {
        py::list wires;
        for (const auto& valueWire : *value.wires) {
            py::dict wire;
            wire["outer"] = valueWire.outer;
            py::list edges;
            for (const auto& valueEdge : valueWire.edges) {
                py::dict edge;
                edge["use_id"] = valueEdge.useId;
                edge["edge_tag"] = valueEdge.edgeTag;
                edge["start_vertex_tag"] = valueEdge.startVertexTag;
                edge["end_vertex_tag"] = valueEdge.endVertexTag;
                edge["orientation"] = orientationName(valueEdge.orientation);
                edge["sense"] = valueEdge.sense;
                edge["seam"] = valueEdge.seam;
                edge["degenerated"] = valueEdge.degenerated;
                edge["parameter_range"] = std::array{valueEdge.firstParameter, valueEdge.lastParameter};
                edge["tolerance"] = finiteValue(valueEdge.tolerance);
                edges.append(std::move(edge));
            }
            wire["edges"] = edges;
            wires.append(std::move(wire));
        }
        result["wires"] = wires;
    }
    return result;
}
py::dict bodyObject(const BodyDto& value) {
    auto result = entityObject(value.entity);
    result["volume"] = finiteValue(value.properties.volume);
    result["volume_centroid"] = pointObject(value.properties.volumeCentroid);
    result["bounds"] = boundsObject(value.properties.bounds);
    result["metrics_reason"] = geo::occtPropertyReasonName(value.properties.reason);
    py::list uses;
    std::set<int> tags;
    for (const auto& valueUse : value.faces) {
        tags.insert(valueUse.faceTag);
        py::dict use;
        use["face_tag"] = valueUse.faceTag;
        use["orientation"] = orientationName(valueUse.orientation);
        use["sense"] = valueUse.sense;
        uses.append(std::move(use));
    }
    result["face_tags"] = std::vector<int>(tags.begin(), tags.end());
    result["face_uses"] = uses;
    return result;
}
py::dict edgeObject(const EdgeDto& value) {
    auto result = entityObject(value.entity);
    result["curve_type"] = curveName(value.curve);
    result["length"] = finiteValue(value.length);
    result["start_vertex_tag"] = value.startVertex;
    result["end_vertex_tag"] = value.endVertex;
    result["vertex_tags"] = std::array{value.startVertex, value.endVertex};
    result["parameter_range"] = std::isfinite(value.first) && std::isfinite(value.last)
        ? py::cast(std::array{value.first, value.last}) : py::object(py::none());
    result["periodic"] = value.periodic;
    result["degenerated"] = value.length == 0;
    result["face_tags"] = value.faceTags;
    return result;
}
py::dict vertexObject(const VertexDto& value) {
    auto result = entityObject(value.entity);
    result["position"] = pointObject(value.position);
    return result;
}
} // namespace cad_inspection_detail

py::dict namesSummary(const geo::OcctModel& model) {
    geo::OcctSourceNameSummary summary;
    {
        py::gil_scoped_release release;
        std::lock_guard lock(nativeOperationMutex);
        summary = model.sourceNameSummary();
    }
    return cad_inspection_detail::summaryObject(summary);
}
py::list inspectEntities(const CadModel& holder, const std::string& kind,
                         std::vector<int> tags, bool topology, std::uint64_t maximumUses) {
    if (tags.size() > 1000)
        throw std::invalid_argument("CAD inspection accepts at most 1000 entity tags per page");
    if (maximumUses == 0 || maximumUses > 10000)
        throw std::invalid_argument("max_uses must be between 1 and 10000");
    if (!holder.imported.model)
        throw std::invalid_argument("CAD inspection requires an imported native model");
    cad_inspection_detail::Snapshot snapshot;
    {
        py::gil_scoped_release release;
        std::lock_guard lock(nativeOperationMutex);
        snapshot = cad_inspection_detail::collect(*holder.imported.model, kind, tags, topology, maximumUses);
    }
    py::list result;
    for (const auto& value : snapshot.faces) result.append(cad_inspection_detail::faceObject(value));
    for (const auto& value : snapshot.bodies) result.append(cad_inspection_detail::bodyObject(value));
    for (const auto& value : snapshot.edges) result.append(cad_inspection_detail::edgeObject(value));
    for (const auto& value : snapshot.vertices) result.append(cad_inspection_detail::vertexObject(value));
    return result;
}
