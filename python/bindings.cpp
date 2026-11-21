// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "NumpyArrays.hpp"

#include "paacs/algo/Meshers.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include "paacs/sizing/RegionSizeField.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <vector>

namespace py = pybind11;
namespace fs = std::filesystem;
namespace {
using namespace paacs;
PyObject* importError = nullptr;
PyObject* meshingError = nullptr;
PyObject* exportError = nullptr;
// The core worker dispatcher currently permits one external caller at a time.
// Serialize Python import/mesh operations across models; internal worker counts
// still control the single active native job. Export reads only owned mesh data.
std::mutex nativeOperationMutex;

std::string utf8(const fs::path& path) {
    const auto text = path.u8string();
    return {text.begin(), text.end()};
}

py::object finiteValue(double value) {
    // Failed native analyses can retain an unset inf/NaN sentinel. JSON
    // reports use null for unavailable measurements and retain the reason.
    if (std::isfinite(value)) return py::cast(value);
    return py::none();
}

fs::path checkedPath(const fs::path& path) {
    if (path.empty() || path.native().find(fs::path::value_type{}) != fs::path::string_type::npos)
        throw py::value_error("path must be nonempty and contain no NUL characters");
    return fs::absolute(path).lexically_normal();
}

double finiteNumber(const py::object& value, const char* name) {
    if (py::isinstance<py::bool_>(value))
        throw py::type_error(std::string(name) + " must be a number, not bool");
    const double number = PyFloat_AsDouble(value.ptr());
    if (PyErr_Occurred()) throw py::error_already_set();
    if (!std::isfinite(number))
        throw py::value_error(std::string(name) + " must be finite");
    return number;
}

double positiveNumber(const py::object& value, const char* name) {
    const double number = finiteNumber(value, name);
    if (number <= 0)
        throw py::value_error(std::string(name) + " must be finite and positive");
    return number;
}

std::uint64_t positiveInteger(const py::object& value, const char* name) {
    if (!py::isinstance<py::int_>(value) || py::isinstance<py::bool_>(value))
        throw py::type_error(std::string(name) + " must be an integer, not bool");
    const auto number = PyLong_AsUnsignedLongLong(value.ptr());
    if (PyErr_Occurred()) throw py::error_already_set();
    if (number == 0)
        throw py::value_error(std::string(name) + " must be positive");
    return number;
}

double growthNumber(const py::object& value, const char* name) {
    const double number = finiteNumber(value, name);
    if (number < 1)
        throw py::value_error(std::string(name) + " must be at least one");
    return number;
}

std::optional<double> optionalGrowth(const py::object& value) {
    if (value.is_none()) return std::nullopt;
    return growthNumber(value, "growth");
}

Vec3 point3(const py::object& value, const char* name) {
    if (!PySequence_Check(value.ptr()) || py::isinstance<py::str>(value) ||
        py::isinstance<py::bytes>(value))
        throw py::type_error(std::string(name) + " must be a sequence of three finite numbers");
    const auto sequence = py::reinterpret_borrow<py::sequence>(value);
    if (py::len(sequence) != 3)
        throw py::value_error(std::string(name) + " must contain exactly three coordinates");
    return {finiteNumber(sequence[0], name), finiteNumber(sequence[1], name),
            finiteNumber(sequence[2], name)};
}

// These value objects own only C++ data. Python properties are read-only, and
// each mesh request copies the controls before entering the native operation.
struct FaceSize {
    int tag;
    double size;
    std::optional<double> growth;
};
struct BoxSize {
    Vec3 lower, upper;
    double size;
    std::optional<double> growth;
};

FaceSize faceSize(const py::object& tag, const py::object& size, const py::object& growth) {
    const auto integer = positiveInteger(tag, "tag");
    if (integer > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw py::value_error("tag exceeds the native face-tag range");
    return {static_cast<int>(integer), positiveNumber(size, "size"), optionalGrowth(growth)};
}

BoxSize boxSize(const py::object& lower, const py::object& upper,
                const py::object& size, const py::object& growth) {
    const auto lo = point3(lower, "lower"), hi = point3(upper, "upper");
    if (!(lo.x < hi.x && lo.y < hi.y && lo.z < hi.z))
        throw py::value_error("lower must be strictly below upper in every coordinate");
    return {lo, hi, positiveNumber(size, "size"), optionalGrowth(growth)};
}

template<class Control>
std::vector<Control> copiedControls(const py::object& value, const char* name) {
    if (value.is_none()) return {};
    if (!PySequence_Check(value.ptr()) || py::isinstance<py::str>(value) ||
        py::isinstance<py::bytes>(value))
        throw py::type_error(std::string(name) + " must be a sequence of sizing controls");
    const auto sequence = py::reinterpret_borrow<py::sequence>(value);
    std::vector<Control> result;
    // Do not trust an arbitrary sequence length for an eager allocation.
    // Native route admission still decides whether the control count fits.
    const auto eagerLimit = static_cast<std::size_t>(sizing::RegionSizeSpecification{}.maximumControls);
    result.reserve(std::min(py::len(sequence), eagerLimit));
    for (const auto item : sequence) {
        if (!py::isinstance<Control>(item))
            throw py::type_error(std::string(name) + " contains an incorrect control type");
        result.push_back(py::cast<Control>(item));
    }
    return result;
}

struct SizingRequest {
    double minimum = 0, growth = 1.3;
    std::optional<double> requestedMinimum;
    std::optional<bool> requestedNative;
    bool native = false;
    std::vector<FaceSize> faces;
    std::vector<BoxSize> boxes;

    py::dict report(double globalSize) const {
        py::dict result;
        result["global_size"] = globalSize;
        result["min_size"] = minimum;
        result["requested_min_size"] = py::cast(requestedMinimum);
        result["growth"] = growth;
        result["native_sizing"] = native;
        result["requested_native_sizing"] = py::cast(requestedNative);
        result["mode"] = native ? "spatial" : "standard";
        py::list faceReports, boxReports;
        for (const auto& face : faces) {
            py::dict control;
            control["tag"] = face.tag;
            control["size"] = face.size;
            control["effective_size"] = std::max(minimum, face.size);
            control["growth"] = face.growth.value_or(growth);
            control["requested_growth"] = py::cast(face.growth);
            faceReports.append(std::move(control));
        }
        for (const auto& box : boxes) {
            py::dict control;
            control["kind"] = "box";
            control["lower"] = std::array{box.lower.x, box.lower.y, box.lower.z};
            control["upper"] = std::array{box.upper.x, box.upper.y, box.upper.z};
            control["size"] = box.size;
            control["effective_size"] = std::max(minimum, box.size);
            control["growth"] = box.growth.value_or(growth);
            control["requested_growth"] = py::cast(box.growth);
            boxReports.append(std::move(control));
        }
        result["face_sizes"] = faceReports;
        result["volume_controls"] = boxReports;
        return result;
    }
};

[[noreturn]] void fail(PyObject* type, const std::string& stage,
                       const std::string& reason, py::dict report) {
    py::object error = py::reinterpret_borrow<py::object>(type)(stage + ": " + reason);
    error.attr("stage") = stage;
    error.attr("reason") = reason;
    error.attr("report") = std::move(report);
    PyErr_SetObject(type, error.ptr());
    throw py::error_already_set();
}

// The deadline requests native cooperative cancellation. It does not interrupt
// an OCCT operation or impose a process memory/runtime limit.
class Deadline {
public:
    explicit Deadline(std::optional<double> seconds) {
        if (seconds) {
            timer_ = std::thread([this, seconds] {
                std::unique_lock lock(mutex_);
                if (!condition_.wait_for(lock, std::chrono::duration<double>(*seconds),
                                         [this] { return finished_; }))
                    cancelled.store(true, std::memory_order_relaxed);
            });
        }
    }
    ~Deadline() {
        { std::lock_guard lock(mutex_); finished_ = true; }
        condition_.notify_all();
        if (timer_.joinable()) timer_.join();
    }
    std::atomic<bool> cancelled{false};
private:
    std::mutex mutex_;
    std::condition_variable condition_;
    bool finished_ = false;
    std::thread timer_;
};

py::dict importReport(const geo::NativeCadImportOutcome& cad, const fs::path& source,
                      double seconds) {
    py::dict r;
    r["source_path"] = utf8(source);
    r["reason"] = geo::cadImportReasonName(cad.reason);
    r["imported"] = cad.imported();
    r["file_bytes"] = cad.fileBytes;
    r["body_count"] = cad.solids;
    r["face_count"] = cad.faces;
    r["edge_count"] = cad.edges;
    r["vertex_count"] = cad.vertices;
    r["curved_faces"] = cad.curvedFaces;
    r["curved_edges"] = cad.curvedEdges;
    r["cad_volume"] = finiteValue(cad.cadVolume);
    r["units"] = "mm";
    r["millimetres_per_model_unit"] = finiteValue(cad.millimetresPerModelUnit);
    r["source_length_units"] = cad.sourceLengthUnits;
    r["import_seconds"] = finiteValue(seconds);
    return r;
}

struct MeshResult {
    mesh::Mesh value;
    geo::NativeCadMeshOutcome outcome;
    fs::path source;
    double importSeconds = 0;
    double requestedSize = 0;
    std::string recovery;
    std::string requestedAssembly;
    SizingRequest sizing;
    unsigned workers = 1, validationWorkers = 1;

    std::vector<int> tags(int dimension) const {
        std::set<int> found;
        for (Index i = 0; i < value.cellCount(); ++i)
            if (value.cell(i).onEntity.dim == dimension)
                found.insert(value.cell(i).onEntity.tag);
        return {found.begin(), found.end()};
    }

    py::dict report() const {
        const auto& q = outcome.quality;
        const auto& a = outcome.audit;
        py::dict r, audit, quality, times;
        r["status"] = outcome.published ? "succeeded" : "refused";
        r["reason"] = geo::nativeCadMeshReasonName(outcome.reason);
        r["mesh_published"] = outcome.published;
        r["source_path"] = utf8(source);
        r["units"] = "mm";
        r["requested_size"] = requestedSize;
        r["recovery"] = recovery;
        r["requested_assembly"] = requestedAssembly;
        r["sizing"] = sizing.report(requestedSize);
        r["workers"] = workers;
        r["validation_workers"] = validationWorkers;
        // Early refusals may retain the native outcome's default policy. With
        // no published mesh, report the request rather than that unset value.
        const bool independent = outcome.published
            ? outcome.assemblyPolicy == geo::NativeCadAssemblyPolicy::IndependentParts
            : requestedAssembly == "preserve_parts";
        r["mesh_scope"] = independent ? "independent_parts" : "conforming_domain";
        r["interbody_nonoverlap_certified"] = outcome.published && !independent;
        r["validated_parts"] = outcome.validatedParts;
        r["nodes"] = value.nodeCount();
        r["cells"] = value.cellCount();
        r["tetrahedra"] = q.tet4;
        r["body_tags"] = tags(3);
        r["face_tags"] = tags(2);
        r["audit_valid"] = a.valid;
        r["quality_warning"] = outcome.qualityWarning;
        audit["complete"] = a.complete;
        audit["valid"] = a.valid;
        audit["volume_cells"] = a.volumeCells;
        audit["boundary_cells"] = a.boundaryCells;
        audit["missing_boundary"] = a.missingBoundary;
        audit["unexpected_boundary"] = a.unexpectedBoundary;
        audit["nonmanifold_faces"] = a.nonmanifoldFaces;
        audit["inconsistent_orientation"] = a.inconsistentOrientation;
        audit["invalid_connectivity"] = a.invalidConnectivity;
        audit["nonfinite_nodes"] = a.nonfiniteNodes;
        r["audit"] = audit;
        quality["complete"] = q.complete;
        quality["valid"] = q.valid;
        quality["metric"] = "subtet_mean_ratio";
        quality["minimum"] = finiteValue(q.minimumQuality);
        quality["percentile_1"] = finiteValue(q.percentile1);
        quality["percentile_5"] = finiteValue(q.percentile5);
        quality["median"] = finiteValue(q.medianQuality);
        quality["mean"] = finiteValue(q.meanQuality);
        quality["inverted_cells"] = q.invertedCells;
        quality["total_volume"] = finiteValue(q.totalVolume);
        quality["minimum_volume"] = finiteValue(q.minimumVolume);
        r["quality"] = quality;
        times["import"] = finiteValue(importSeconds);
        times["preparation"] = finiteValue(outcome.preparationSeconds);
        times["generation"] = finiteValue(outcome.generationSeconds);
        times["validation"] = finiteValue(outcome.validationSeconds);
        times["audit"] = finiteValue(outcome.auditSeconds);
        times["quality"] = finiteValue(outcome.qualitySeconds);
        times["meshing_total"] = finiteValue(outcome.totalSeconds);
        r["seconds"] = times;
        py::list bodies, interferences;
        for (const auto& body : outcome.bodies) {
            py::dict b;
            b["solid_tag"] = body.solidTag;
            b["route"] = body.route;
            b["reason"] = geo::nativeCadMeshReasonName(body.reason);
            b["surface_reason"] = body.nativeSurfaceReason;
            b["volume_reason"] = body.nativeVolumeReason;
            b["cad_volume"] = finiteValue(body.cadVolume);
            b["mesh_volume"] = finiteValue(body.meshVolume);
            b["relative_volume_error"] = finiteValue(body.relativeVolumeError);
            b["maximum_boundary_node_error"] = finiteValue(body.maximumBoundaryNodeError);
            b["maximum_boundary_sample_error"] = finiteValue(body.maximumBoundarySampleError);
            b["boundary_embedding_certified"] = body.boundaryEmbeddingCertified;
            b["centroid_outside_distance_limit"] = finiteValue(body.centroidOutsideDistanceLimit);
            b["native_centroids_checked"] = body.nativeCentroidsChecked;
            b["local_delaunay_refinement_nodes"] = body.localDelaunayRefinementNodes;
            b["surface_chord_error_limit"] = finiteValue(body.surfaceChordErrorLimit);
            bodies.append(std::move(b));
        }
        for (const auto& pair : outcome.interferences) {
            py::dict p;
            p["first_solid_tag"] = pair.firstSolidTag;
            p["second_solid_tag"] = pair.secondSolidTag;
            p["state"] = pair.state == geo::NativeCadInterferenceState::MaterialOverlap
                ? "material_overlap" : pair.state == geo::NativeCadInterferenceState::NoMaterialOverlap
                ? "no_material_overlap" : "unresolved";
            p["overlap_volume"] = finiteValue(pair.nativeOverlapVolume);
            interferences.append(std::move(p));
        }
        r["bodies"] = bodies;
        r["interferences"] = interferences;
        return r;
    }

    fs::path outputPath(const fs::path& path) const {
        const auto destination = checkedPath(path);
        std::error_code error;
        if (fs::equivalent(source, destination, error) ||
            fs::weakly_canonical(source) == fs::weakly_canonical(destination))
            throw py::value_error("CAD input and mesh output must be different files");
        return destination;
    }

    void writeVtk(const fs::path& path, bool volumeOnly) const {
        const auto destination = outputPath(path);
        std::string message;
        { py::gil_scoped_release release;
          try { io::writeVtkFile(value, destination, volumeOnly); }
          catch (const std::bad_alloc&) { throw; }
          catch (const std::exception& error) { message = error.what(); }
        }
        if (!message.empty()) {
            py::dict r; r["path"] = utf8(destination); r["detail"] = message;
            fail(exportError, "export", "io_error", std::move(r));
        }
    }

    void writeArchive(const fs::path& path) const {
        const auto destination = outputPath(path);
        sdk::MeshArchiveOutcome saved;
        { py::gil_scoped_release release;
          saved = sdk::writeMeshArchiveFile(value, destination); }
        if (!saved.completed) {
            constexpr const char* reasons[] = {"none", "invalid_mesh", "invalid_header",
                "unsupported_version", "truncated", "count_overflow", "invalid_cell",
                "io_error", "allocation_failed", "internal_error", "resource_limit"};
            const auto index = static_cast<unsigned>(saved.reason);
            py::dict r; r["path"] = utf8(destination);
            fail(exportError, "export", index < std::size(reasons) ? reasons[index] : "unknown", r);
        }
    }

    void write(const fs::path& path) const {
        auto suffix = utf8(path.extension());
        std::transform(suffix.begin(), suffix.end(), suffix.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (suffix == ".vtk") writeVtk(path, false);
        else if (suffix == ".paacsmesh") writeArchive(path);
        else throw py::value_error("output extension must be .vtk or .paacsmesh");
    }
};

struct CadModel {
    geo::NativeCadImportOutcome imported;
    fs::path source;
    double importSeconds = 0;

    std::vector<int> faceTags() const {
        std::vector<int> result;
        for (const auto* face : imported.model->faces()) result.push_back(face->tag());
        std::sort(result.begin(), result.end());
        return result;
    }

    std::vector<int> bodyTags() const {
        std::vector<int> result;
        for (const auto* body : imported.model->solids()) result.push_back(body->tag());
        std::sort(result.begin(), result.end());
        return result;
    }

    std::unique_ptr<MeshResult> mesh(const py::object& size, const std::string& recovery,
        const std::string& assembly, const py::object& workers, const py::object& validationWorkers,
        const py::object& maxNodes, const py::object& maxCells, const py::object& maxWork,
        const py::object& maxOutputBytes, const py::object& maxWorkingBytes,
        const py::object& maxAuditBytes, const py::object& maxSeconds,
        const py::object& faceSizes, const py::object& volumeControls,
        const py::object& minimumSize, const py::object& growth,
        const py::object& nativeSizing) {
        geo::NativeCadMeshOptions options;
        options.geometry.targetSize = positiveNumber(size, "size");
        SizingRequest sizing;
        sizing.growth = growthNumber(growth, "growth");
        sizing.faces = copiedControls<FaceSize>(faceSizes, "face_sizes");
        sizing.boxes = copiedControls<BoxSize>(volumeControls, "volume_controls");
        if (!minimumSize.is_none()) {
            const double minimum = finiteNumber(minimumSize, "min_size");
            if (minimum < 0 || minimum > options.geometry.targetSize)
                throw py::value_error("min_size must be between zero and size");
            sizing.minimum = minimum;
            sizing.requestedMinimum = minimum;
        }
        if (!nativeSizing.is_none()) {
            if (!py::isinstance<py::bool_>(nativeSizing))
                throw py::type_error("native_sizing must be bool or None");
            sizing.requestedNative = nativeSizing.cast<bool>();
        }
        if (sizing.requestedNative == false && !sizing.boxes.empty())
            throw py::value_error("volume_controls require native_sizing=True or None");
        sizing.native = sizing.requestedNative.value_or(!sizing.faces.empty() || !sizing.boxes.empty());
        options.useSpatialSizing = sizing.native;
        options.minimumSize = sizing.minimum;
        for (const auto& face : sizing.faces)
            options.faceSizes.push_back({face.tag, face.size, face.growth.value_or(sizing.growth)});
        for (const auto& box : sizing.boxes)
            options.sizeBoxes.push_back({box.lower, box.upper, box.size, box.growth.value_or(sizing.growth)});
        if (recovery == "standard") options.recoveryMode = geo::NativeCadRecoveryMode::Standard;
        else if (recovery == "advanced") options.recoveryMode = geo::NativeCadRecoveryMode::Advanced;
        else throw py::value_error("recovery must be 'standard' or 'advanced'");
        if (assembly == "preserve_parts") options.assemblyPolicy = geo::NativeCadAssemblyPolicy::IndependentParts;
        else if (assembly == "conforming") options.assemblyPolicy = geo::NativeCadAssemblyPolicy::ConformingDomain;
        else throw py::value_error("assembly must be 'preserve_parts' or 'conforming'");
        const auto nw = positiveInteger(workers, "workers");
        const auto nv = positiveInteger(validationWorkers, "validation_workers");
        if (nw > 1024 || nv > 8)
            throw py::value_error("workers must be 1..1024 and validation_workers 1..8");
        options.geometry.threads = static_cast<int>(nw);
        options.maximumNativeValidationThreads = static_cast<unsigned>(nv);
        if (!maxNodes.is_none()) options.geometry.maximumNodes = positiveInteger(maxNodes, "max_nodes");
        if (!maxCells.is_none()) options.geometry.maximumCells = positiveInteger(maxCells, "max_cells");
        if (!maxWork.is_none()) {
            const auto work = positiveInteger(maxWork, "max_work");
            options.geometry.maximumMeshingWork = work;
            options.maximumValidationWork = work;
        }
        if (!maxOutputBytes.is_none()) options.geometry.maximumOutputBytes = positiveInteger(maxOutputBytes, "max_output_bytes");
        if (!maxWorkingBytes.is_none()) options.maximumPrimitiveWorkingBytes = positiveInteger(maxWorkingBytes, "max_working_bytes");
        if (!maxAuditBytes.is_none()) options.geometry.maximumAuditTemporaryBytes = positiveInteger(maxAuditBytes, "max_audit_bytes");
        std::optional<double> seconds;
        if (!maxSeconds.is_none()) {
            seconds = positiveNumber(maxSeconds, "max_seconds");
            const double clockRange = std::chrono::duration<double>(
                std::chrono::steady_clock::duration::max()).count();
            if (*seconds > clockRange / 2)
                throw py::value_error("max_seconds exceeds the native clock range");
        }
        auto result = std::make_unique<MeshResult>();
        result->source = source;
        result->importSeconds = importSeconds;
        result->requestedSize = options.geometry.targetSize;
        result->recovery = recovery;
        result->requestedAssembly = assembly;
        result->sizing = std::move(sizing);
        result->workers = static_cast<unsigned>(nw);
        result->validationWorkers = static_cast<unsigned>(nv);
        { py::gil_scoped_release release;
          std::lock_guard lock(nativeOperationMutex);
          Deadline deadline(seconds);
          if (seconds) options.geometry.cancellation = &deadline.cancelled;
          result->outcome = geo::meshNativeCad(*imported.model, options, result->value);
        }
        const auto& outcome = result->outcome;
        if (!outcome.published) {
            auto report = result->report();
            if (outcome.reason == geo::NativeCadMeshReason::UnsupportedControl && options.useSpatialSizing) {
                if (options.assemblyPolicy == geo::NativeCadAssemblyPolicy::IndependentParts)
                    report["hint"] = "Native spatial sizing requires assembly='conforming' and a single solid. "
                                     "The requested assembly policy was preserved.";
                else if (imported.solids != 1)
                    report["hint"] = "Native spatial sizing currently supports one Tet4 solid. "
                                     "Use a single-solid STEP input.";
            }
            else if (outcome.reason == geo::NativeCadMeshReason::UnsupportedControl &&
                     options.minimumSize > 0 && options.faceSizes.empty())
                report["hint"] = "A nonzero min_size requires a sizing control, or native_sizing=True "
                                 "with assembly='conforming' and a single solid.";
            fail(meshingError, "mesh", geo::nativeCadMeshReasonName(outcome.reason), std::move(report));
        }
        if (outcome.reason != geo::NativeCadMeshReason::None || !outcome.audit.complete ||
            !outcome.audit.valid || !outcome.quality.complete || !outcome.quality.valid)
            fail(meshingError, "mesh", "publication_checks_failed", result->report());
        return result;
    }
};

std::unique_ptr<CadModel> readStep(const fs::path& path, const py::object& maxFileBytes) {
    auto result = std::make_unique<CadModel>();
    result->source = checkedPath(path);
    geo::NativeCadImportOptions options;
    if (!maxFileBytes.is_none()) options.maximumFileBytes = positiveInteger(maxFileBytes, "max_file_bytes");
    { py::gil_scoped_release release;
      std::lock_guard lock(nativeOperationMutex);
      const auto start = std::chrono::steady_clock::now();
      result->imported = geo::importOcctCad(result->source, options);
      result->importSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    if (!result->imported.imported())
        fail(importError, "import", geo::cadImportReasonName(result->imported.reason),
             importReport(result->imported, result->source, result->importSeconds));
    return result;
}
} // namespace

PYBIND11_MODULE(_native, m) {
    m.doc() = "Owned Python interface to the native PAACS STEP/Tet4 engine.";
    m.attr("engine_version") = paacs::sdk::kApiVersion;
    importError = PyErr_NewException("paacs_mesh.CadImportError", PyExc_RuntimeError, nullptr);
    meshingError = PyErr_NewException("paacs_mesh.MeshingError", PyExc_RuntimeError, nullptr);
    exportError = PyErr_NewException("paacs_mesh.ExportError", PyExc_OSError, nullptr);
    m.add_object("CadImportError", py::reinterpret_steal<py::object>(importError));
    m.add_object("MeshingError", py::reinterpret_steal<py::object>(meshingError));
    m.add_object("ExportError", py::reinterpret_steal<py::object>(exportError));
    py::class_<FaceSize>(m, "FaceSize", "Read-only CAD face sizing control, in millimetres.")
        .def(py::init(&faceSize), py::arg("tag"), py::arg("size"), py::kw_only(),
             py::arg("growth") = py::none())
        .def_readonly("tag", &FaceSize::tag)
        .def_readonly("size", &FaceSize::size)
        .def_readonly("growth", &FaceSize::growth);
    py::class_<BoxSize>(m, "BoxSize", "Read-only volume refinement box, in millimetres.")
        .def(py::init(&boxSize), py::arg("lower"), py::arg("upper"), py::arg("size"),
             py::kw_only(), py::arg("growth") = py::none())
        .def_property_readonly("lower", [](const BoxSize& box) {
            return py::make_tuple(box.lower.x, box.lower.y, box.lower.z);
        })
        .def_property_readonly("upper", [](const BoxSize& box) {
            return py::make_tuple(box.upper.x, box.upper.y, box.upper.z);
        })
        .def_readonly("size", &BoxSize::size)
        .def_readonly("growth", &BoxSize::growth);
    py::class_<MeshResult>(m, "Mesh")
        .def_property_readonly("node_count", [](const MeshResult& mesh) { return mesh.value.nodeCount(); })
        .def_property_readonly("cell_count", [](const MeshResult& mesh) { return mesh.value.cellCount(); })
        .def_property_readonly("tetrahedron_count", [](const MeshResult& mesh) { return mesh.outcome.quality.tet4; })
        .def_property_readonly("face_tags", [](const MeshResult& mesh) { return mesh.tags(2); })
        .def_property_readonly("body_tags", [](const MeshResult& mesh) { return mesh.tags(3); })
        .def_property_readonly("report", &MeshResult::report)
        .def("to_numpy", [](const MeshResult& mesh, const py::object& maximumBytes) {
            return paacs::python::meshArrays(mesh.value, positiveInteger(maximumBytes, "max_bytes"));
        }, py::kw_only(), py::arg("max_bytes") = py::int_(536870912),
            "Copy owned, initially read-only NumPy arrays; max_bytes bounds their numerical payload.")
        .def("write", &MeshResult::write, py::arg("path"))
        .def("write_vtk", &MeshResult::writeVtk, py::arg("path"), py::kw_only(), py::arg("volume_only") = false)
        .def("write_archive", &MeshResult::writeArchive, py::arg("path"));
    py::class_<CadModel>(m, "CadModel")
        .def_property_readonly("info", [](const CadModel& model) {
            return importReport(model.imported, model.source, model.importSeconds);
        })
        .def_property_readonly("face_tags", &CadModel::faceTags)
        .def_property_readonly("body_tags", &CadModel::bodyTags)
        .def("mesh", &CadModel::mesh, py::arg("size"), py::kw_only(),
             py::arg("recovery") = "standard", py::arg("assembly") = "preserve_parts",
             py::arg("workers") = py::int_(1), py::arg("validation_workers") = py::int_(1),
             py::arg("max_nodes") = py::none(), py::arg("max_cells") = py::none(),
             py::arg("max_work") = py::none(), py::arg("max_output_bytes") = py::none(),
             py::arg("max_working_bytes") = py::none(), py::arg("max_audit_bytes") = py::none(),
             py::arg("max_seconds") = py::none(),
             py::arg("face_sizes") = py::none(), py::arg("volume_controls") = py::none(),
             py::arg("min_size") = py::none(), py::arg("growth") = py::float_(1.3),
             py::arg("native_sizing") = py::none());
    m.def("read_step", &readStep, py::arg("path"), py::kw_only(), py::arg("max_file_bytes") = py::none());
}
