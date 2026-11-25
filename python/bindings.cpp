// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "NumpyArrays.hpp"

#include "paacs/algo/Meshers.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "paacs/geometry/OcctBodyRelations.hpp"
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/geometry/OcctPreview.hpp"
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
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
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

double qualityThreshold(const py::object& value, const char* name) {
    const double number = finiteNumber(value, name);
    if (number < 0 || number > 1)
        throw py::value_error(std::string(name) + " must be between zero and one");
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

py::dict materialClassificationReport(const geo::NativeCadMaterialClassificationReport& source) {
    py::dict r;r["diagnostic"]=source.diagnostic;r["complete"]=source.completed;
    r["source_model_verified_unchanged"]=source.sourceModelUnmodified;
    r["attempts"]=source.attempts;r["requested_points"]=source.requestedPoints;
    r["attempted_points"]=source.attemptedPoints;r["native_queries"]=source.nativeQueries;
    r["strict_reuses"]=source.strictReuses;r["ray_admissions"]=source.rayAdmissions;
    r["work"]=source.work;r["seconds"]=finiteValue(source.seconds);
    r["profiled_batches"]=source.profiledBatches;r["ray_witness_attempts"]=source.rayWitnessAttempts;
    r["classifier_setup_seconds"]=finiteValue(source.classifierSetupSeconds);
    r["native_query_seconds"]=finiteValue(source.nativeQuerySeconds);
    r["ray_witness_seconds"]=finiteValue(source.rayWitnessSeconds);
    r["scope"]="cumulative material classification attempts; last batch status; not final publication";
    return r;
}

py::dict sourceShellReport(const geo::NativeCadSourceShellReport& source) {
    py::dict r;py::list attempts;
    const auto count=std::min<std::size_t>(source.count,source.attempts.size());
    r["available"]=count>0;r["accepted"]=count?py::cast(source.accepted):py::none();
    r["diagnostic"]=source.diagnostic;
    r["native_fallback_requested"]=source.nativeFallbackRequested;
    r["native_fallback_attempted"]=source.nativeFallbackAttempted;
    r["native_fallback_diagnostic"]=source.nativeFallbackDiagnostic;
    r["requested_global_size"]=source.requestedGlobalSize>0?finiteValue(source.requestedGlobalSize):py::none();
    r["requested_chord_error"]=source.requestedChordError>0?finiteValue(source.requestedChordError):py::none();
    r["maximum_relative_volume_error"]=source.maximumRelativeVolumeError>0?finiteValue(source.maximumRelativeVolumeError):py::none();
    for(std::size_t i=0;i<count;++i){
        const auto& a=source.attempts[i];py::dict row;
        row["attempt"]=i+1;row["chord_error_limit"]=finiteValue(a.chordErrorLimit);
        row["work"]=a.work;row["seconds"]=finiteValue(a.seconds);row["nodes"]=a.nodes;row["triangles"]=a.triangles;
        row["surface_published"]=a.published;row["surface_diagnostic"]=a.surfaceDiagnostic;
        row["source_model_verified_unchanged"]=a.sourceModelUnmodified;
        row["volume_criterion_met"]=a.published?py::cast(a.volumeCriterionMet):py::none();
        row["signed_volume_mm3"]=a.published?finiteValue(a.signedVolume):py::none();
        row["cad_volume_mm3"]=a.published?finiteValue(a.cadVolume):py::none();
        row["relative_volume_error"]=a.published?finiteValue(a.relativeVolumeError):py::none();
        row["construction_budget_limited"]=a.constructionBudgetLimited;
        row["construction_stalled"]=a.constructionStalled;
        row["provider"]=a.provider;
        row["global_sizing_applied"]=a.globalSizingApplied;
        row["curvature_points_per_turn"]=a.curvaturePointsPerTurn;
        row["angular_deflection"]=a.angularDeflection>0?finiteValue(a.angularDeflection):py::none();
        py::dict witness;
        witness["stage"]=a.witness.stage;witness["detail"]=a.witness.detail;
        witness["edge_tag"]=a.witness.edgeTag>=0?py::cast(a.witness.edgeTag):py::none();
        witness["face_tag"]=a.witness.faceTag>=0?py::cast(a.witness.faceTag):py::none();
        witness["use_id"]=py::cast(a.witness.useId);
        witness["geometric_refusal"]=a.witness.geometricRefusal;
        row["witness"]=std::move(witness);
        attempts.append(std::move(row));
    }
    r["attempts"]=std::move(attempts);
    r["scope"]="latest_source_strategy; geometry acceptance is not a sizing certificate";
    return r;
}

py::dict routeResolutionReport(const geo::NativeCadRouteResolutionReport& source) {
    py::dict r;r["reason"]=source.reason;r["candidate"]=source.candidate;
    r["preflight_attempted"]=source.preflightAttempted;r["whole_body_verified"]=source.wholeBodyVerified;
    r["preflight_work"]=source.preflightWork;r["preflight_seconds"]=finiteValue(source.preflightSeconds);
    r["fallback_chord_limit"]=source.fallbackChordLimit?finiteValue(*source.fallbackChordLimit):py::none();
    py::list attempts;
    for(std::size_t i=0;i<std::min<std::size_t>(source.count,source.attempts.size());++i){
        const auto& a=source.attempts[i];py::dict row;
        row["method"]=a.method;row["reason"]=geo::nativeCadMeshReasonName(a.reason);
        row["published"]=a.published;row["work"]=a.work;row["seconds"]=finiteValue(a.seconds);
        attempts.append(std::move(row));
    }
    r["attempts"]=std::move(attempts);
    r["scope"]="method choice and charged attempts; not an element-size or solver certificate";
    return r;
}

py::dict constraintRecoveryReport(const geo::NativeCadConstraintRecoveryReport& source) {
    py::dict r;
    r["diagnostic"]=source.diagnostic;
    r["candidate_cells"]=source.candidateCells;
    r["missing_edges_before"]=source.missingEdgesBefore;
    r["missing_edges_after"]=source.missingEdgesAfter;
    r["missing_faces_before"]=source.missingFacesBefore;
    r["missing_faces_after"]=source.missingFacesAfter;
    r["patches"]=source.patches;r["seed_cells"]=source.seedCells;
    r["work"]=source.work;r["seconds"]=finiteValue(source.seconds);
    r["edge_walks"]=source.edgeWalks;r["facet_walks"]=source.facetWalks;
    r["walk_fallbacks"]=source.walkFallbacks;
    r["steiner_attempts"]=source.steinerAttempts;r["steiner_work"]=source.steinerWork;
    r["steiner_nodes"]=source.steinerNodes;r["retained_steiner_nodes"]=source.retainedSteinerNodes;
    r["piercing_edge_patches"]=source.piercingPatches;
    r["initial_topology_used"]=source.initialTopologyUsed;
    r["initial_topology_diagnostic"]=source.initialTopologyDiagnostic;
    r["initial_topology_cells"]=source.initialTopologyCells;
    r["initial_topology_import_work"]=source.initialTopologyImportWork;
    r["initial_topology_peak_modeled_bytes"]=source.initialTopologyPeakModeledBytes;
    return r;
}

py::dict initialTopologyQualityReport(const geo::NativeCadInitialTopologyQualityReport& source) {
    py::dict r;
    r["attempted"]=source.attempted;r["complete"]=source.completed;r["changed"]=source.changed;
    r["cancelled"]=source.cancelled;r["budget_limited"]=source.budgetLimited;
    r["summary_available"]=source.summaryAvailable;r["candidate_list_truncated"]=source.candidateListTruncated;
    r["diagnostic"]=source.diagnostic;r["work"]=source.work;r["work_limit"]=source.workLimit;
    r["seconds"]=finiteValue(source.seconds);r["peak_modeled_bytes"]=source.maximumModeledBytes;
    r["search_below_mean_ratio"]=source.searchBelowMeanRatio;
    r["candidates"]=source.candidates;r["proposals"]=source.proposals;r["certifications"]=source.certifications;
    r["staged_edits"]=source.stagedEdits;r["edits"]=source.edits;
    r["flips23"]=source.flips23;r["flips32"]=source.flips32;r["flips44"]=source.flips44;
    r["live_before"]=source.liveBefore;r["live_after"]=source.liveAfter;
    r["severe_before"]=source.severeBefore;r["severe_after"]=source.severeAfter;
    r["protected_source_faces"]=source.protectedSourceFaces;r["protected_source_edges"]=source.protectedSourceEdges;
    r["minimum_mean_ratio_before"]=source.summaryAvailable?finiteValue(source.minimumMeanRatioBefore):py::none();
    r["minimum_mean_ratio_after"]=source.summaryAvailable?finiteValue(source.minimumMeanRatioAfter):py::none();
    r["minimum_scaled_jacobian_before"]=source.summaryAvailable?finiteValue(source.minimumScaledJacobianBefore):py::none();
    r["minimum_scaled_jacobian_after"]=source.summaryAvailable?finiteValue(source.minimumScaledJacobianAfter):py::none();
    r["delaunay_certified"]=false;
    r["scope"]="nonguard numerical cells including exterior; final Native audits remain required";
    return r;
}
py::dict initialFallbackReport(const geo::NativeCadInitialFrontFallbackReport& fallback) {
    py::dict r;
    r["attempted"] = fallback.attempted;
    r["reason"] = geo::nativeCadMeshReasonName(fallback.reason);
    r["work"] = fallback.work;
    r["seconds"] = finiteValue(fallback.seconds);
    r["surface_diagnostic"] = fallback.surfaceDiagnostic;
    r["volume_diagnostic"] = fallback.volumeDiagnostic;
    r["source_shell"] = sourceShellReport(fallback.sourceShell);
    return r;
}

py::dict planarBodyReport(const algo::GeometryBodyReport& body) {
    py::dict r, options, estimates, engine, output, conversion, boundary, surface;
    r["solid_tag"]=body.solidTag;
    r["strategy"]=algo::geometryStrategyName(body.strategy);
    r["reason"]=algo::geometryReasonName(body.reason);
    options["elements"]=body.elements==algo::GeometryElements::Tet4?"tet4":"native";
    options["target_size"]=finiteValue(body.targetSize);
    options["source_face_tag"]=body.sourceFaceTag;options["target_face_tag"]=body.targetFaceTag;
    options["layers"]=body.layers;options["end_to_start_layer_width_ratio"]=finiteValue(body.endToStartLayerWidthRatio);
    options["used_spatial_field"]=body.usedSpatialField;options["used_size_profile"]=body.usedSizeProfile;
    options["spatial_minimum_bound"]=finiteValue(body.spatialMinimumBound);
    r["resolved_options"]=options;
    estimates["nodes"]=body.estimatedNodes;estimates["total_cells"]=body.estimatedCells;
    estimates["output_bytes"]=body.estimatedOutputBytes;estimates["exact"]=body.exactCountEstimate;
    r["estimates"]=estimates;
    engine["status_code"]=body.engineStatus;engine["reason_code"]=body.engineReason;
    engine["phase_code"]=body.enginePhase;r["engine"]=engine;
    output["nodes"]=body.outputNodes;output["total_cells"]=body.outputCells;output["bytes"]=body.outputBytes;
    r["output"]=output;r["generation_seconds"]=finiteValue(body.generationSeconds);
    const auto& c=body.conversion;
    conversion["published"]=c.published;conversion["reason"]=algo::convexToTetReasonName(c.reason);
    conversion["tetrahedra"]=c.tetrahedra;conversion["boundary_triangles"]=c.boundaryTriangles;
    conversion["output_bytes"]=c.outputBytes;conversion["working_bytes"]=c.workingBytes;
    conversion["minimum_quality"]=finiteValue(c.minimumQuality);conversion["low_quality_cells"]=c.lowQualityCells;
    r["conversion"]=conversion;
    boundary["available"]=body.hasBoundaryDiagnostic;boundary["status_code"]=body.boundaryStatus;
    boundary["reason_code"]=body.boundaryReason;boundary["phase_code"]=body.boundaryPhase;
    r["boundary"]=boundary;
    surface["available"]=body.hasSurfaceDiagnostic;surface["status_code"]=body.surfaceStatus;
    surface["reason_code"]=body.surfaceReason;surface["failure_face_id"]=body.surfaceFailureFacePersistentId;
    surface["processed_faces"]=body.surfaceProcessedFaces;surface["unattempted_faces"]=body.surfaceUnattemptedFaces;
    r["surface"]=surface;
    return r;
}

py::dict planarReport(const geo::NativeCadBodyMeshReport& body) {
    py::dict r;
    r["attempted"]=body.planarAttempted;
    r["scope"]="retained geometry-core receipt; default substage codes do not prove stage entry; separate from Native CAD certification";
    if(!body.planarAttempted){r["published"]=false;r["reason"]=py::none();return r;}
    const auto& p=body.planar;
    r["published"]=p.published;r["reason"]=algo::geometryReasonName(p.reason);
    py::dict analysis,times;
    analysis["complete"]=p.analysis.complete;analysis["reason"]=algo::geometryReasonName(p.analysis.reason);
    analysis["work"]=p.analysis.work;r["analysis"]=analysis;
    times["analysis"]=finiteValue(p.analysisSeconds);times["generation"]=finiteValue(p.generationSeconds);
    times["audit"]=finiteValue(p.auditSeconds);times["quality"]=finiteValue(p.qualitySeconds);
    r["seconds"]=times;
    constexpr std::size_t maximumBodyRows=16;
    r["bodies_total"]=p.analysis.bodies.size();r["bodies_truncated"]=p.analysis.bodies.size()>maximumBodyRows;
    py::list bodies;
    for(std::size_t i=0;i<std::min(p.analysis.bodies.size(),maximumBodyRows);++i)
        bodies.append(planarBodyReport(p.analysis.bodies[i]));
    r["bodies"]=bodies;
    return r;
}

struct MeshResult {
    mesh::Mesh value;
    geo::NativeCadMeshOutcome outcome;
    fs::path source;
    std::string sourceSha256;
    double importSeconds = 0;
    double requestedSize = 0;
    std::string recovery;
    std::string requestedAssembly;
    SizingRequest sizing;
    std::optional<double> requestedWarningQuality, requestedMinimumQuality;
    std::optional<double> requestedChordError, requestedRelativeVolumeError;
    double maximumRelativeVolumeError=.01;
    std::optional<unsigned> requestedCurvaturePoints;
    unsigned curvaturePoints=0;
    std::optional<bool> requestedInitialTopologyQuality;
    bool requestedNativeSourceFallback=false;
    bool requestedBoundaryFirst=false;
    std::optional<bool> requestedAdvancingInteriorFront;
    double warningQuality = 0, minimumQuality = 0;
    unsigned workers = 1, validationWorkers = 1;
    geo::CadConstructionBudget construction;
    std::optional<std::uint64_t> requestedWork;
    std::uint64_t effectiveWork=0,workingBytes=0;

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
        r["source_sha256"] = sourceSha256;
        r["kernel_version"] = geo::occtKernelVersion();
        r["units"] = "mm";
        r["requested_size"] = requestedSize;
        r["recovery"] = recovery;
        r["requested_assembly"] = requestedAssembly;
        r["requested_advancing_interior_front"] = py::cast(requestedAdvancingInteriorFront);
        r["initial_front_fallback"] = initialFallbackReport(outcome.initialFrontFallback);
        r["sizing"] = sizing.report(requestedSize);
        py::dict constructionReport;
        constructionReport["effort"] = construction.effort==geo::CadConstructionEffort::Legacy?"legacy":
            construction.effort==geo::CadConstructionEffort::Automatic?"auto":"thorough";
        const auto requestedCap=[](std::uint64_t value)->py::object {
            return value?py::object(py::int_(value)):py::object(py::none());
        };
        constructionReport["max_surface_points"] = requestedCap(construction.surfacePointsPerFace);
        constructionReport["max_surface_chord_points"] = requestedCap(construction.surfaceChordPointsPerFace);
        constructionReport["max_interior_candidates"] = requestedCap(construction.interiorCandidates);
        constructionReport["max_interior_sites"] = requestedCap(construction.interiorSites);
        constructionReport["scope"] = "stage ceilings subordinate to overall resource limits; not a sizing certificate";
        r["construction"] = std::move(constructionReport);
        py::dict resources;
        resources["policy"] = construction.resourcePolicy==geo::CadResourcePolicy::Normal?"normal":"bounded";
        resources["requested_max_work"] = requestedWork?py::cast(*requestedWork):py::none();
        resources["effective_max_work"] = effectiveWork;
        resources["logical_work_cap_disabled"] = construction.resourcePolicy==geo::CadResourcePolicy::Normal&&!requestedWork;
        resources["max_working_bytes"] = workingBytes;
        resources["scope"] = "logical native admission and modeled memory; not CPU instructions or a process RSS limit";
        r["resources"] = std::move(resources);
        py::dict fidelity;
        fidelity["requested_max_chord_error"]=py::cast(requestedChordError);
        fidelity["requested_max_relative_volume_error"]=py::cast(requestedRelativeVolumeError);
        fidelity["effective_max_relative_volume_error"]=maximumRelativeVolumeError;
        fidelity["requested_curvature_points_per_turn"]=py::cast(requestedCurvaturePoints);
        fidelity["effective_curvature_points_per_turn"]=curvaturePoints;
        fidelity["chord_default"]=requestedChordError?"explicit":"route_specific";
        fidelity["scope"]="native geometry acceptance limits; measured errors and achieved sizes are separate";
        r["fidelity"]=std::move(fidelity);
        r["route_resolution"]=routeResolutionReport(outcome.routeResolution);
        r["requested_initial_topology_quality"]=py::cast(requestedInitialTopologyQuality);
        r["requested_native_source_fallback"]=requestedNativeSourceFallback;
        r["requested_boundary_first"]=requestedBoundaryFirst;
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
        quality["warning_threshold"] = warningQuality;
        quality["warning_enabled"] = warningQuality > 0;
        quality["minimum_required"] = minimumQuality;
        quality["requested_warning_threshold"] = py::cast(requestedWarningQuality);
        quality["requested_minimum_required"] = py::cast(requestedMinimumQuality);
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
        py::list bodies, interferences, bodyRoutes;
        for (const auto& body : outcome.bodies) {
            py::dict b;
            b["planar"]=planarReport(body);
            b["solid_tag"] = body.solidTag;
            b["route"] = body.route;
            b["reason"] = geo::nativeCadMeshReasonName(body.reason);
            b["surface_reason"] = body.nativeSurfaceReason;
            b["volume_reason"] = body.nativeVolumeReason;
            b["surface_triangles"] = body.surfaceTriangles;
            b["surface_attempts"] = body.surfaceAttempts;
            b["recovery_source_facet_nodes"] = body.recoverySourceFacetNodes;
            b["native_surface_work"] = body.nativeSurfaceWork;
            b["native_volume_work"] = body.nativeVolumeWork;
            b["native_surface_seconds"] = finiteValue(body.nativeSurfaceSeconds);
            b["native_volume_seconds"] = finiteValue(body.nativeVolumeSeconds);
            b["native_partition_seconds"] = finiteValue(body.nativePartitionSeconds);
            b["native_assembly_seconds"] = finiteValue(body.nativeAssemblySeconds);
            b["partition_reason"]=body.partitionReason;
            b["partition_policy"]=body.partitionPolicy;
            b["partition_attempts"]=body.partitionAttempts;
            b["partition_pieces"]=body.partitionPieces;
            b["partition_interfaces"]=body.partitionInterfaces;
            b["partition_work"]=body.partitionWork;
            b["partition_volume_change"]=finiteValue(body.partitionVolumeChange);
            b["partition_volume_admission"]=finiteValue(body.partitionVolumeAdmission);
            b["sampled_partition_trim_witness"]=body.sampledPartitionTrimWitness;
            b["partition_trim_samples"]=body.partitionTrimSamples;
            b["maximum_partition_trim_drift"]=finiteValue(body.maximumPartitionTrimDrift);
            b["maximum_partition_trim_tolerance_ratio"]=finiteValue(body.maximumPartitionTrimToleranceRatio);
            b["surface_diagonal_repairs"]=body.surfaceDiagonalRepairs;
            b["exact_surface_pairs"]=body.exactSurfacePairs;
            b["kernel_trials"]=body.kernelTrials;
            py::list partitionAttempts;
            for(const auto& attempt:body.partitionAttemptReports){
                py::dict row;row["policy"]=attempt.policy;row["diagnostic"]=attempt.diagnostic;
                row["failed_piece_tag"]=attempt.failedPieceTag;row["pieces"]=attempt.pieces;
                row["interfaces"]=attempt.interfaces;row["work"]=attempt.work;
                row["seconds"]=finiteValue(attempt.seconds);row["surface_deflection"]=finiteValue(attempt.surfaceDeflection);
                partitionAttempts.append(std::move(row));
            }
            b["partition_attempt_reports"]=std::move(partitionAttempts);
            b["cad_volume"] = finiteValue(body.cadVolume);
            b["mesh_volume"] = finiteValue(body.meshVolume);
            b["relative_volume_error"] = finiteValue(body.relativeVolumeError);
            b["maximum_boundary_node_error"] = finiteValue(body.maximumBoundaryNodeError);
            b["maximum_boundary_sample_error"] = finiteValue(body.maximumBoundarySampleError);
            b["boundary_embedding_certified"] = body.boundaryEmbeddingCertified;
            b["centroid_outside_distance_limit"] = finiteValue(body.centroidOutsideDistanceLimit);
            b["native_centroids_checked"] = body.nativeCentroidsChecked;
            b["native_candidate_queries"] = body.nativeCandidateQueries;
            b["native_candidate_ray_admissions"] = body.nativeCandidateRayAdmissions;
            b["local_delaunay_refinement_nodes"] = body.localDelaunayRefinementNodes;
            b["surface_chord_error_limit"] = finiteValue(body.surfaceChordErrorLimit);
            b["source_shell"] = sourceShellReport(body.sourceShell);
            b["initial_topology_quality"]=initialTopologyQualityReport(body.initialTopologyQuality);
            b["route_resolution"]=routeResolutionReport(body.routeResolution);
            py::dict bodyFidelity;
            bodyFidelity["effective_max_chord_error"]=body.effectiveMaximumChordError?
                finiteValue(*body.effectiveMaximumChordError):py::none();
            bodyFidelity["chord_policy"]=body.chordErrorPolicy;
            bodyFidelity["scope"]="selected construction limit; not a complete CAD distance certificate";
            b["fidelity"]=std::move(bodyFidelity);
            b["constraint_recovery"] = constraintRecoveryReport(body.constraintRecovery);
            b["material_classification"] = materialClassificationReport(body.materialClassification);
            py::dict sizingApplication;
            sizingApplication["applied"]=py::cast(body.volumeSizingApplied);
            sizingApplication["diagnostic"]=body.volumeSizingDiagnostic;
            sizingApplication["requested_global_size"]=finiteValue(requestedSize);
            sizingApplication["sizing_completion_certified"]=false;
            sizingApplication["basis"]="reported_native_construction";
            b["volume_sizing"]=std::move(sizingApplication);
            const auto& initial = body.initialFront;
            py::dict placement;
            placement["diagnostic"] = initial.diagnostic;
            placement["nodes"] = initial.nodes;
            placement["wall_nodes"] = initial.wallNodes;
            placement["void_nodes"] = initial.voidNodes;
            placement["void_candidates"] = initial.voidCandidates;
            placement["source_sweep_complete"] = initial.sourceSweepComplete;
            placement["void_queue_exhausted"] = initial.voidQueueExhausted;
            placement["maximum_unresolved_void_ratio"] = finiteValue(initial.maximumUnresolvedVoidRatio);
            placement["candidates"] = initial.candidates;
            placement["work"] = initial.work;
            placement["queries"] = initial.queries;
            placement["delaunay_commits"] = initial.delaunayCommits;
            placement["full_certifications"] = initial.fullCertifications;
            placement["seconds"] = finiteValue(initial.seconds);
            placement["budget_reason"] = initial.budgetReason;
            placement["candidate_limit"] = initial.candidateLimit;
            placement["site_limit"] = initial.siteLimit;
            placement["query_limit"] = initial.queryLimit;
            placement["work_limit"] = initial.workLimit;
            placement["reserved_bulk_candidates"] = initial.reservedBulkCandidates;
            placement["reserved_bulk_sites"] = initial.reservedBulkSites;
            placement["source_candidates"] = initial.sourceCandidates;
            placement["source_deferred_faces"] = initial.sourceDeferredFaces;
            placement["source_phase_deferred"] = initial.sourcePhaseDeferred;
            placement["cones"] = initial.cones;
            placement["propagated_sites"] = initial.propagatedSites;
            placement["remaining_faces"] = initial.remainingFaces;
            placement["storage_rejections"] = initial.storageRejections;
            placement["compactions"] = initial.compactions;
            placement["retired_cells_reclaimed"] = initial.retiredCellsReclaimed;
            placement["compaction_work"] = initial.compactionWork;
            placement["auxiliary_stored_cells"] = initial.auxiliaryStoredCells;
            placement["auxiliary_live_cells"] = initial.auxiliaryLiveCells;
            placement["auxiliary_stored_peak"] = initial.auxiliaryStoredPeak;
            placement["auxiliary_storage_capacity"] = initial.auxiliaryStorageCapacity;
            placement["compaction_peak_modeled_bytes"] = initial.compactionPeakModeledBytes;
            placement["auxiliary_working_bytes"] = initial.auxiliaryWorkingBytes;
            placement["auxiliary_imported_cells"] = initial.auxiliaryImportedCells;
            placement["auxiliary_import_work"] = initial.auxiliaryImportWork;
            placement["auxiliary_import_modeled_bytes"] = initial.auxiliaryImportModeledBytes;
            placement["auxiliary_initial_capacity"] = initial.auxiliaryInitialCapacity;
            placement["auxiliary_hard_capacity"] = initial.auxiliaryHardCapacity;
            placement["auxiliary_growth_attempts"] = initial.auxiliaryGrowthAttempts;
            placement["auxiliary_growths"] = initial.auxiliaryGrowths;
            placement["auxiliary_growth_work"] = initial.auxiliaryGrowthWork;
            placement["maximum_combined_modeled_bytes"] = initial.maximumCombinedModeledBytes;
            placement["maximum_heap_capacity"] = initial.maximumHeapCapacity;
            placement["resource_retries"] = initial.resourceRetries;
            placement["geometry_refusal_cache_hits"] = initial.geometryRefusalCacheHits;
            placement["cached_geometry_refusals"] = initial.cachedGeometryRefusals;
            placement["auxiliary_resource_diagnostic"] = initial.auxiliaryResourceDiagnostic;
            b["initial_front"] = std::move(placement);
            py::list surfaceSpacing;
            for(const auto& face:body.surfaceSpacing){
                py::dict row;
                row["face_tag"]=face.faceTag;
                row["point_limit"]=face.pointLimit;
                row["chord_point_limit"]=face.chordPointLimit;
                row["chart_node_limit"]=face.chartNodeLimit;
                row["added_points"]=face.addedPoints;
                row["chord_added_points"]=face.chordAddedPoints;
                row["pending_entries"]=face.pendingEntries;
                row["blocked_candidates"]=face.blockedCandidates;
                row["termination"]=face.termination;
                row["queue_exhausted"]=face.queueExhausted;
                row["budget_limited"]=face.budgetLimited;
                surfaceSpacing.append(std::move(row));
            }
            b["surface_spacing"]=std::move(surfaceSpacing);
            b["initial_front_fallback"] = initialFallbackReport(body.initialFrontFallback);
            // This is the native adaptive stage after the initial point set.
            // Its diagnostic distinguishes completion from a bounded stop;
            // a route name alone does not certify achieved sizing.
            py::dict interior;
            interior["diagnostic"] = body.interiorFrontDiagnostic;
            interior["nodes"] = body.interiorFrontNodes;
            interior["candidates"] = body.interiorFrontCandidates;
            interior["work"] = body.interiorFrontWork;
            interior["queries"] = body.interiorFrontQueries;
            interior["remaining_oversized"] = body.interiorFrontRemainingOversized;
            interior["front_activations"] = body.interiorFrontActivations;
            interior["void_activations"] = body.interiorVoidActivations;
            interior["front_advances"] = body.interiorFrontAdvances;
            interior["seconds"] = finiteValue(body.interiorFrontSeconds);
            b["interior_front"] = std::move(interior);
            py::dict route;
            route["solid_tag"] = body.solidTag;
            route["route"] = body.route;
            route["reason"] = geo::nativeCadMeshReasonName(body.reason);
            bodyRoutes.append(std::move(route));
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
        // Preserve the legacy sizing report, which describes wrapper policy.
        // AUTO may select a different construction route in the native core.
        // Retain requests and actual body reports separately, including refused
        // attempts; none of these fields asserts that size targets were met.
        py::dict policy;
        policy["requested_native_sizing"] = py::cast(sizing.requestedNative);
        policy["requested_advancing_interior_front"] = py::cast(requestedAdvancingInteriorFront);
        policy["requested_native_source_fallback"] = requestedNativeSourceFallback;
        policy["requested_boundary_first"] = requestedBoundaryFirst;
        policy["wrapper_native_sizing"] = sizing.native;
        policy["wrapper_sizing_mode"] = sizing.native ? "spatial" : "standard";
        policy["body_routes"] = std::move(bodyRoutes);
        r["effective_policy"] = std::move(policy);
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
    std::shared_ptr<const std::string> snapshot;
    geo::NativeCadImportOutcome imported;
    fs::path source;
    std::string sourceSha256;
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
        const py::object& nativeSizing, const py::object& warningQuality,
        const py::object& minimumQuality, const py::object& advancingInteriorFront,
        const std::string& constructionEffort, const py::object& maxSurfacePoints,
        const py::object& maxSurfaceChordPoints, const py::object& maxInteriorCandidates,
        const py::object& maxInteriorSites, const std::string& resourcePolicy,
        const py::object& maxChordError, const py::object& maxRelativeVolumeError,
        const py::object& curvaturePointsPerTurn, const py::object& initialTopologyQuality,
        const py::object& nativeSourceFallback, const py::object& boundaryFirst) {
        if(resourcePolicy!="normal"&&resourcePolicy!="bounded")
            throw py::value_error("resource_policy must be 'normal' or 'bounded'");
        auto options=resourcePolicy=="normal"?geo::normalNativeCadMeshOptions():geo::NativeCadMeshOptions{};
        options.geometry.targetSize = positiveNumber(size, "size");
        if(!py::isinstance<py::bool_>(nativeSourceFallback))
            throw py::type_error("native_source_fallback must be bool");
        options.allowNativeSourceFallback=nativeSourceFallback.cast<bool>();
        if(!py::isinstance<py::bool_>(boundaryFirst))
            throw py::type_error("boundary_first must be bool");
        options.recoverBoundaryBeforeInteriorSites=boundaryFirst.cast<bool>();
        if(!initialTopologyQuality.is_none()){
            if(!py::isinstance<py::bool_>(initialTopologyQuality))
                throw py::type_error("initial_topology_quality must be bool or None");
            options.initialTopologyQuality=initialTopologyQuality.cast<bool>();
        }
        std::optional<unsigned> requestedCurvaturePoints;
        if(!curvaturePointsPerTurn.is_none()){
            if(py::isinstance<py::bool_>(curvaturePointsPerTurn)||!PyIndex_Check(curvaturePointsPerTurn.ptr()))
                throw py::type_error("curvature_points_per_turn must be an integer");
            auto index=py::reinterpret_steal<py::object>(PyNumber_Index(curvaturePointsPerTurn.ptr()));
            if(!index)throw py::error_already_set();
            const auto value=PyLong_AsLong(index.ptr());
            if(PyErr_Occurred())throw py::error_already_set();
            if(value<0||value>64)throw py::value_error("curvature_points_per_turn must be between zero and 64");
            requestedCurvaturePoints=static_cast<unsigned>(value);
            options.curvaturePointsPerTurn=*requestedCurvaturePoints;
        }
        std::optional<double> requestedChordError, requestedRelativeVolumeError;
        if(!maxChordError.is_none()){
            requestedChordError=positiveNumber(maxChordError,"max_chord_error");
            options.maximumChordError=*requestedChordError;
        }
        if(!maxRelativeVolumeError.is_none()){
            requestedRelativeVolumeError=positiveNumber(maxRelativeVolumeError,"max_relative_volume_error");
            if(*requestedRelativeVolumeError>=1)throw py::value_error("max_relative_volume_error must be less than one");
            options.maximumRelativeVolumeError=*requestedRelativeVolumeError;
        }
        if(constructionEffort=="legacy")options.construction.effort=geo::CadConstructionEffort::Legacy;
        else if(constructionEffort=="auto")options.construction.effort=geo::CadConstructionEffort::Automatic;
        else if(constructionEffort=="thorough")options.construction.effort=geo::CadConstructionEffort::Thorough;
        else throw py::value_error("construction_effort must be 'legacy', 'auto' or 'thorough'");
        if(!maxSurfacePoints.is_none())options.construction.surfacePointsPerFace=positiveInteger(maxSurfacePoints,"max_surface_points");
        if(!maxSurfaceChordPoints.is_none())options.construction.surfaceChordPointsPerFace=positiveInteger(maxSurfaceChordPoints,"max_surface_chord_points");
        if(!maxInteriorCandidates.is_none())options.construction.interiorCandidates=positiveInteger(maxInteriorCandidates,"max_interior_candidates");
        if(!maxInteriorSites.is_none())options.construction.interiorSites=positiveInteger(maxInteriorSites,"max_interior_sites");
        if (!advancingInteriorFront.is_none()) {
            if (!py::isinstance<py::bool_>(advancingInteriorFront))
                throw py::type_error("advancing_interior_front must be bool or None");
            options.useAdvancingInteriorFront = advancingInteriorFront.cast<bool>();
        }
        const auto requestedPlacement = options.useAdvancingInteriorFront;
        std::optional<double> requestedWarningQuality, requestedMinimumQuality;
        if (!warningQuality.is_none()) {
            requestedWarningQuality = qualityThreshold(warningQuality, "warning_quality");
            options.geometry.warningElementQuality = *requestedWarningQuality;
        }
        if (!minimumQuality.is_none()) {
            requestedMinimumQuality = qualityThreshold(minimumQuality, "minimum_quality");
            options.geometry.minimumElementQuality = *requestedMinimumQuality;
        }
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
        if (sizing.requestedNative == false) {
            if (requestedPlacement == true)
                throw py::value_error("advancing_interior_front=True conflicts with native_sizing=False");
            options.useAdvancingInteriorFront = false;
        }
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
        result->sourceSha256 = sourceSha256;
        result->importSeconds = importSeconds;
        result->requestedSize = options.geometry.targetSize;
        result->recovery = recovery;
        result->requestedAssembly = assembly;
        result->sizing = std::move(sizing);
        result->construction = options.construction;
        if(!maxWork.is_none())result->requestedWork=positiveInteger(maxWork,"max_work");
        result->effectiveWork=std::min(options.maximumValidationWork,options.geometry.maximumMeshingWork);
        result->workingBytes=options.maximumPrimitiveWorkingBytes;
        result->requestedWarningQuality = requestedWarningQuality;
        result->requestedMinimumQuality = requestedMinimumQuality;
        result->requestedChordError=requestedChordError;
        result->requestedRelativeVolumeError=requestedRelativeVolumeError;
        result->maximumRelativeVolumeError=options.maximumRelativeVolumeError;
        result->requestedCurvaturePoints=requestedCurvaturePoints;
        result->curvaturePoints=options.curvaturePointsPerTurn;
        result->requestedInitialTopologyQuality=options.initialTopologyQuality;
        result->requestedNativeSourceFallback=options.allowNativeSourceFallback;
        result->requestedBoundaryFirst=options.recoverBoundaryBeforeInteriorSites;
        result->requestedAdvancingInteriorFront = requestedPlacement;
        result->warningQuality = options.geometry.warningElementQuality;
        result->minimumQuality = options.geometry.minimumElementQuality;
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


// The source identity covers the exact immutable bytes passed to the parser.
// Read in bounded chunks; a caller-supplied large limit does not eagerly
// allocate that amount. Hashing uses Python's standard SHA-256 implementation.
std::unique_ptr<CadModel> readStep(const fs::path& path, const py::object& maxFileBytes,
                                  bool names, bool retainSource) {
    auto result = std::make_unique<CadModel>();
    result->source = checkedPath(path);
    geo::NativeCadImportOptions options;
    if (!maxFileBytes.is_none()) options.maximumFileBytes = positiveInteger(maxFileBytes, "max_file_bytes");
    options.retainRepresentationNames = names;
    auto snapshot = std::make_shared<std::string>();
    const auto start = std::chrono::steady_clock::now();
    { py::gil_scoped_release release;
      result->imported.reason = geo::CadImportReason::FileUnavailable;
      std::error_code fileError;
      const auto advertisedBytes = fs::file_size(result->source, fileError);
      if (!fileError) result->imported.fileBytes = advertisedBytes;
      if (!fileError && advertisedBytes > options.maximumFileBytes)
          result->imported.reason = geo::CadImportReason::FileBudget;
      std::ifstream input;
      if (!fileError && advertisedBytes <= options.maximumFileBytes)
          input.open(result->source, std::ios::binary);
      if (input.is_open()) {
          std::array<char, 65536> buffer;
          bool budget = false;
          while (input) {
              input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
              const auto count = static_cast<std::size_t>(input.gcount());
              if (count > options.maximumFileBytes - snapshot->size()) {
                  budget = true;
                  break;
              }
              snapshot->append(buffer.data(), count);
          }
          result->imported.fileBytes = snapshot->size();
          if (budget) result->imported.reason = geo::CadImportReason::FileBudget;
          else if (!input.bad()) result->imported.reason = geo::CadImportReason::None;
      }
    }
    if (result->imported.reason == geo::CadImportReason::None) {
        auto hash = py::module_::import("hashlib").attr("sha256")();
        for (std::size_t offset = 0; offset < snapshot->size(); offset += 65536)
            hash.attr("update")(py::bytes(snapshot->data() + offset,
                std::min<std::size_t>(65536, snapshot->size() - offset)));
        result->sourceSha256 = hash.attr("hexdigest")().cast<std::string>();
        { py::gil_scoped_release release;
          std::lock_guard lock(nativeOperationMutex);
          result->imported = geo::importOcctStepData(
              std::span<const char>(snapshot->data(), snapshot->size()), utf8(result->source), options);
        }
    }
    result->importSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!result->imported.imported()) {
        auto report = importReport(result->imported, result->source, result->importSeconds);
        if (!result->sourceSha256.empty()) report["source_sha256"] = result->sourceSha256;
        fail(importError, "import", geo::cadImportReasonName(result->imported.reason), std::move(report));
    }
    if (retainSource) result->snapshot = std::move(snapshot);
    return result;
}

#include "CadInspection.hpp"
#include "CadPreview.hpp"
#include "CadBodyRelations.hpp"

// Private diagnostic query: distances are to trimmed faces, including their
// edges/vertices, not merely projections onto the underlying infinite surface.
// PAACS bounds samples/output and checks time between OCCT calls. One OCCT
// distance solve is not preemptible and its internal allocations are not capped.
py::dict faceDistances(const CadModel& holder, const py::object& pointValues,
                      const py::object& faceValues, const py::object& maxWork,
                      const py::object& maxSeconds, bool includeFaces,
                      const py::object& maxResultBytes) {
    const auto maximumWork = positiveInteger(maxWork, "max_work");
    const auto maximumBytes = positiveInteger(maxResultBytes, "max_result_bytes");
    const auto seconds = positiveNumber(maxSeconds, "max_seconds");
    if (maximumWork > 100000 || maximumBytes > 8ULL * 1024 * 1024)
        throw std::invalid_argument("CAD proximity caps are 100000 evaluations and 8 MiB output scratch");
    auto sequence = [](const py::object& value, const char* name, py::ssize_t maximum) {
        if (!PySequence_Check(value.ptr()) || py::isinstance<py::str>(value) || py::isinstance<py::bytes>(value))
            throw std::invalid_argument(std::string(name) + " must be a bounded sequence");
        auto result = py::reinterpret_borrow<py::sequence>(value);
        if (py::len(result) < 1 || py::len(result) > maximum)
            throw std::invalid_argument(std::string(name) + " count exceeds its diagnostic cap");
        return result;
    };
    const auto pointSequence = sequence(pointValues, "points", 256);
    const auto faceSequence = sequence(faceValues, "face_tags", 1000);
    const auto evaluations = static_cast<std::uint64_t>(py::len(pointSequence)) * py::len(faceSequence);
    if (evaluations > maximumWork)
        throw std::length_error("CAD proximity exceeds max_work; restrict cells or face tags");
    const auto scratch = evaluations * (includeFaces ? 1024 : 160) + py::len(pointSequence) * 4096 + 65536;
    if (scratch > maximumBytes)
        throw std::length_error("CAD proximity output exceeds max_result_bytes; request a smaller chunk");
    std::vector<Vec3> points;
    std::vector<int> tags;
    points.reserve(py::len(pointSequence));
    tags.reserve(py::len(faceSequence));
    for (const auto point : pointSequence) points.push_back(point3(py::reinterpret_borrow<py::object>(point), "point"));
    for (const auto value : faceSequence) {
        const auto tag = positiveInteger(py::reinterpret_borrow<py::object>(value), "face tag");
        if (tag > static_cast<std::uint64_t>((std::numeric_limits<int>::max)()))
            throw std::invalid_argument("Face tag exceeds the native integer range");
        tags.push_back(static_cast<int>(tag));
    }
    std::sort(tags.begin(), tags.end());
    if (std::adjacent_find(tags.begin(), tags.end()) != tags.end())
        throw std::invalid_argument("face_tags must be unique");
    struct Distance { Vec3 closest; double distance, tolerance; int tag; };
    std::vector<Distance> measured;
    measured.reserve(evaluations);
    const auto started = std::chrono::steady_clock::now();
    const auto checkTime = [&] {
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() > seconds)
            throw std::runtime_error("CAD proximity exceeded cooperative max_seconds; no partial result published");
    };
    {
        py::gil_scoped_release release;
        std::unique_lock lock(nativeOperationMutex, std::defer_lock);
        while (!lock.try_lock()) { checkTime(); std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
        if (!holder.imported.model) throw std::invalid_argument("CAD proximity requires an imported model");
        for (const auto tag : tags)
            if (!holder.imported.model->face(tag)) throw std::invalid_argument("CAD proximity has an unknown face tag");
        for (const auto point : points) for (const auto tag : tags) {
            checkTime();
            Vec3 closest;
            if (!geo::closestPointOnOcctFace(*holder.imported.model, tag, point, closest))
                throw std::runtime_error("Trimmed CAD face distance failed; no partial result published");
            const auto distance = std::hypot(point.x - closest.x, point.y - closest.y, point.z - closest.z);
            const auto tolerance = holder.imported.model->tolerance(2, tag);
            if (!std::isfinite(distance) || !std::isfinite(tolerance) || tolerance < 0)
                throw std::runtime_error("Trimmed CAD distance is outside the supported numerical range");
            measured.push_back({closest, distance, tolerance, tag});
            checkTime();
        }
    }
    py::list rows;
    for (std::size_t index = 0; index < points.size(); ++index) {
        const auto begin = measured.begin() + index * tags.size(), end = begin + tags.size();
        const auto best = std::min_element(begin, end, [](const auto& a, const auto& b) { return a.distance < b.distance; });
        py::dict row;
        row["point"] = cad_inspection_detail::pointObject(points[index]);
        row["minimum_distance"] = best->distance;
        row["closest_face_tag"] = best->tag;
        row["closest_point"] = cad_inspection_detail::pointObject(best->closest);
        row["face_tolerance"] = best->tolerance;
        row["measured_face_count"] = tags.size();
        py::list nearest, faces;
        const auto roundoff = 32 * std::numeric_limits<double>::epsilon() * std::max(1., best->distance);
        for (auto current = begin; current != end; ++current) {
            if (current->distance - best->distance <= std::max({roundoff, current->tolerance, best->tolerance}))
                nearest.append(current->tag);
            if (includeFaces) {
                py::dict face;
                face["face_tag"] = current->tag; face["distance"] = current->distance;
                face["closest_point"] = cad_inspection_detail::pointObject(current->closest);
                face["tolerance"] = current->tolerance;
                faces.append(std::move(face));
            }
        }
        row["nearest_face_tags"] = std::move(nearest);
        if (includeFaces) row["faces"] = std::move(faces);
        rows.append(std::move(row));
        checkTime();
    }
    py::dict result;
    result["schema"] = "paacs-cad-face-distances-v1";
    result["source_sha256"] = holder.sourceSha256;
    result["kernel_version"] = geo::occtKernelVersion();
    result["evaluations"] = evaluations;
    result["rows"] = std::move(rows);
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
            auto report = importReport(model.imported, model.source, model.importSeconds);
            report["source_sha256"] = model.sourceSha256;
            report["kernel_version"] = geo::occtKernelVersion();
            report["representation_names"] = namesSummary(*model.imported.model);
            return report;
        })
        .def_property_readonly("face_tags", &CadModel::faceTags)
        .def_property_readonly("body_tags", &CadModel::bodyTags)
        .def_property_readonly("edge_tags", [](const CadModel& model) {
            std::vector<int> tags;
            for (const auto* entity : model.imported.model->edges()) tags.push_back(entity->tag());
            std::sort(tags.begin(), tags.end());
            return tags;
        })
        .def_property_readonly("vertex_tags", [](const CadModel& model) {
            std::vector<int> tags;
            for (const auto* entity : model.imported.model->vertices()) tags.push_back(entity->tag());
            std::sort(tags.begin(), tags.end());
            return tags;
        })
        .def_property_readonly("model_identity", [](const CadModel& model) { return model.imported.model->identity(); })
        .def_property_readonly("source_sha256", [](const CadModel& model) { return model.sourceSha256; })
        .def_property_readonly("kernel_version", [](const CadModel&) { return geo::occtKernelVersion(); })
        .def_property_readonly("_source_bytes", [](const CadModel& model) -> py::object {
            return model.snapshot ? py::object(py::bytes(*model.snapshot)) : py::object(py::none());
        })
        .def("_inspect_entities", &inspectEntities, py::arg("kind"), py::arg("tags"), py::kw_only(),
             py::arg("topology") = false, py::arg("max_uses") = 10000)
        .def("_preview", &previewCad, py::kw_only(), py::arg("linear_deflection"),
             py::arg("max_points") = py::int_(50000), py::arg("max_triangles") = py::int_(100000))
        .def("_body_relations", &bodyRelations, py::arg("pairs"), py::kw_only(),
             py::arg("distance_tolerance") = py::float_(1e-7),
             py::arg("max_work") = py::int_(100000),
             py::arg("max_working_bytes") = py::int_(64 * 1024 * 1024),
             py::arg("max_seconds") = py::float_(30.0),
             "Query selected native body relations without constructing a mesh interface.")
        .def("_face_distances", &faceDistances, py::arg("points"), py::arg("face_tags"), py::kw_only(),
             py::arg("max_work") = py::int_(100000), py::arg("max_seconds") = py::float_(10.0),
             py::arg("include_faces") = false, py::arg("max_result_bytes") = py::int_(8 * 1024 * 1024),
             "Private centroid diagnostics use true trimmed-face distances; one OCCT call is not preemptible.")
        .def("mesh", &CadModel::mesh, py::arg("size"), py::kw_only(),
             py::arg("recovery") = "standard", py::arg("assembly") = "preserve_parts",
             py::arg("workers") = py::int_(1), py::arg("validation_workers") = py::int_(1),
             py::arg("max_nodes") = py::none(), py::arg("max_cells") = py::none(),
             py::arg("max_work") = py::none(), py::arg("max_output_bytes") = py::none(),
             py::arg("max_working_bytes") = py::none(), py::arg("max_audit_bytes") = py::none(),
             py::arg("max_seconds") = py::none(),
             py::arg("face_sizes") = py::none(), py::arg("volume_controls") = py::none(),
             py::arg("min_size") = py::none(), py::arg("growth") = py::float_(1.3),
             py::arg("native_sizing") = py::none(),
             py::arg("warning_quality") = py::none(), py::arg("minimum_quality") = py::none(),
             py::arg("advancing_interior_front") = py::none(),
             py::arg("construction_effort") = "legacy",
             py::arg("max_surface_points") = py::none(), py::arg("max_surface_chord_points") = py::none(),
             py::arg("max_interior_candidates") = py::none(), py::arg("max_interior_sites") = py::none(),
             py::arg("resource_policy") = "normal",
             py::arg("max_chord_error") = py::none(), py::arg("max_relative_volume_error") = py::none(),
             py::arg("curvature_points_per_turn") = py::none(), py::arg("initial_topology_quality") = py::none(),
             py::arg("native_source_fallback") = py::bool_(false), py::arg("boundary_first") = py::bool_(false));
    m.def("read_step", &readStep, py::arg("path"), py::kw_only(), py::arg("max_file_bytes") = py::none(), py::arg("names") = true, py::arg("retain_source") = false);
}
