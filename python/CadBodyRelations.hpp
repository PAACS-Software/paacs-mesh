// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Included inside bindings.cpp's anonymous namespace after CadModel. The
// geometry header is included at translation-unit scope by bindings.cpp.
py::dict bodyRelations(const CadModel& holder, const py::object& pairValues,
                       const py::object& tolerance, const py::object& maxWork,
                       const py::object& maxBytes, const py::object& maxSeconds) {
    geo::OcctBodyRelationOptions options;
    options.maximumPairs = 64;
    options.maximumWork = positiveInteger(maxWork, "max_work");
    options.maximumRetainedBytes = positiveInteger(maxBytes, "max_working_bytes");
    options.distanceTolerance = finiteNumber(tolerance, "distance_tolerance");
    const auto deadlineSeconds = positiveNumber(maxSeconds, "max_seconds");
    const double clockRange = std::chrono::duration<double>(
        std::chrono::steady_clock::duration::max()).count();
    if (deadlineSeconds > clockRange / 2)
        throw py::value_error("max_seconds exceeds the native clock range");
    if (options.distanceTolerance < 0)
        throw py::value_error("distance_tolerance must be nonnegative");
    if (options.maximumWork > 10000000 || options.maximumRetainedBytes > 1024ULL * 1024 * 1024)
        throw py::value_error("Body relation caps are 10000000 work units and 1 GiB tracked storage");
    options.maximumSourceBytes = std::max<std::uint64_t>(1, std::min<std::uint64_t>(
        64ULL * 1024 * 1024, options.maximumRetainedBytes / 4));
    if (!PySequence_Check(pairValues.ptr()) || py::isinstance<py::str>(pairValues) ||
        py::isinstance<py::bytes>(pairValues))
        throw py::type_error("pairs must be a bounded sequence of body pairs");
    const auto values = py::reinterpret_borrow<py::sequence>(pairValues);
    const auto count = py::len(values);
    if (count < 1 || count > options.maximumPairs)
        throw py::value_error("pairs must contain between one and 64 pairs");
    std::vector<geo::OcctBodyPairRequest> pairs;
    pairs.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const auto item = values[index];
        if (!PySequence_Check(item.ptr()) || py::isinstance<py::str>(item) ||
            py::isinstance<py::bytes>(item))
            throw py::type_error("Each body pair must be a sequence of two body tags");
        const auto row = py::reinterpret_borrow<py::sequence>(item);
        if (py::len(row) != 2)
            throw py::value_error("Each body pair must contain exactly two body tags");
        const auto first = positiveInteger(row[0], "first body tag");
        const auto second = positiveInteger(row[1], "second body tag");
        if (first > static_cast<std::uint64_t>((std::numeric_limits<int>::max)()) ||
            second > static_cast<std::uint64_t>((std::numeric_limits<int>::max)()))
            throw py::value_error("Body tag exceeds the native integer range");
        pairs.push_back({static_cast<int>(first), static_cast<int>(second)});
    }
    geo::OcctBodyRelationReport measured;
    double seconds = 0;
    {
        py::gil_scoped_release release;
        std::lock_guard lock(nativeOperationMutex);
        const auto started = std::chrono::steady_clock::now();
        Deadline deadline(deadlineSeconds);
        options.cancellation = &deadline.cancelled;
        measured = geo::inspectOcctBodyPairs(*holder.imported.model, pairs, options);
        seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    }
    py::dict result;
    result["schema"] = "paacs-cad-body-relations-v1";
    result["source_sha256"] = holder.sourceSha256;
    result["model_identity"] = holder.imported.model->identity();
    result["kernel_version"] = geo::occtKernelVersion();
    result["units"] = "mm";
    result["completed"] = measured.complete;
    result["reason"] = geo::nativeCadMeshReasonName(measured.reason);
    result["diagnostic"] = measured.diagnostic;
    result["source_audit_available"] = measured.sourceAuditAvailable;
    result["source_model_verified_unchanged"] = measured.sourceUnmodified;
    result["source_snapshots"] = measured.sourceSnapshots;
    result["work"] = measured.work;
    result["maximum_logical_working_bytes"] = measured.maximumLogicalWorkingBytes;
    result["seconds"] = seconds;
    result["mesh_interface_certified"] = false;
    result["requested_pairs"] = pairs.size();
    py::list rows;
    for (const auto& pair : measured.pairs) {
        py::dict row;
        row["first_body_tag"] = pair.firstSolidTag;
        row["second_body_tag"] = pair.secondSolidTag;
        row["relation"] = geo::occtBodyRelationName(pair.relation);
        row["evidence"] = geo::occtBodyRelationEvidenceName(pair.evidence);
        row["issue"] = geo::occtBodyRelationIssueName(pair.issue);
        row["classification_tolerance"] = pair.classificationTolerance;
        row["minimum_distance"] = py::cast(pair.minimumDistance);
        row["native_contact_area"] = py::cast(pair.nativeContactArea);
        row["native_contact_length"] = py::cast(pair.nativeContactLength);
        row["native_overlap_volume"] = py::cast(pair.nativeOverlapVolume);
        row["material_query_complete"] = pair.materialQueryComplete;
        row["boundary_query_complete"] = pair.boundaryQueryComplete;
        row["native_copies"] = pair.nativeCopies;
        row["native_boolean_calls"] = pair.nativeBooleanCalls;
        row["native_section_calls"] = pair.nativeSectionCalls;
        row["native_distance_calls"] = pair.nativeDistanceCalls;
        row["native_shape_validations"] = pair.nativeShapeValidations;
        row["native_volume_integrations"] = pair.nativeVolumeIntegrations;
        row["native_boundary_measures"] = pair.nativeBoundaryMeasures;
        rows.append(std::move(row));
    }
    result["pairs"] = std::move(rows);
    return result;
}
