// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/geometry/OcctImport.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include "NumberArgument.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <new>
#include <stdexcept>
#include <string_view>

namespace {
template <class T> T number(const char *text) {
    const std::string_view s(text);
    T value{};
    if (!paacs::cli::parseNumber(s, value))
        throw std::invalid_argument("invalid number");
    return value;
}
void jsonString(std::string_view text) {
    std::cout << '"';
    constexpr char hex[] = "0123456789abcdef";
    for (const unsigned char c : text) {
        if (c == '"' || c == '\\')
            std::cout << '\\' << static_cast<char>(c);
        else if (c < 32)
            std::cout << "\\u00" << hex[c >> 4] << hex[c & 15];
        else
            std::cout << static_cast<char>(c);
    }
    std::cout << '"';
}
const char* surfaceReasonName(std::uint32_t reason) noexcept {
    // SurfaceReason order is the numeric receipt retained by GeometryBodyReport.
    constexpr const char* names[] = {
        "none", "invalid_loop_topology", "invalid_constraint_embedding",
        "invalid_segment_identity", "invalid_physical_policy", "stale_geometry_provider",
        "stale_candidate_provider", "stale_face_artifact", "stale_curve_artifact",
        "duplicate_face_identity", "duplicate_shared_binding", "missing_shared_binding",
        "unused_shared_binding", "unknown_shared_global_index", "shared_classification_conflict",
        "shared_identity_conflict", "shared_mapped_position_mismatch", "invalid_classification",
        "nonfinite_parameter", "nonfinite_mapped_position", "nonfinite_mapped_normal",
        "provider_domain_rejected", "unsupported_periodic_chart", "unsupported_seam_chart",
        "unsupported_singular_chart", "boundary_steiner_prohibited",
        "curve_artifact_subdivision_required", "minimum_feature_reached", "split_budget_exceeded",
        "geometry_provider_error", "candidate_provider_error", "candidate_validation_failed",
        "recovery_stopped", "domain_extraction_failed", "physical_orientation_failed",
        "local_audit_failed", "publication_validation_failed", "allocation_failed", "internal_error",
    };
    return reason < std::size(names) ? names[reason] : "unknown";
}
const char* recoveryReasonName(std::uint32_t kind, std::uint32_t detail) noexcept {
    constexpr const char* rejects[] = {
        "none", "stale_geometry_provider", "nonfinite_mapped_position", "provider_domain_rejected",
        "non_injective_staged_id", "unknown_staged_reference", "stale_candidate_stamp",
        "chart_mismatch", "invalid_split_policy", "invalid_classification", "pslg_rejected",
    };
    constexpr const char* refusals[] = {
        "none", "steiner_disabled", "split_budget_exceeded", "minimum_feature_reached",
        "unsupported_chart", "unsupported_seam", "crossed_required_edge", "partial_chain",
        "strip_not_recoverable",
    };
    constexpr const char* failures[] = {
        "none", "invalid_steiner_geometry", "cavity_rejected", "candidate_allocation_failed",
        "geometry_provider_error", "staged_node_allocation_failed", "cavity_threw", "index_rebuild_failed",
    };
    if (kind == 0) return "none";
    if (kind == 1 && detail < std::size(rejects)) return rejects[detail];
    if (kind == 2 && detail < std::size(refusals)) return refusals[detail];
    if (kind == 3 && detail < std::size(failures)) return failures[detail];
    return "unknown";
}
} // namespace
int main(int argc, char **argv) {
    for(int i=1;i<argc;++i) if(std::string_view(argv[i])=="--native-cad") {
        extern int runNativeCadCli(int, char**);
        return runNativeCadCli(argc,argv);
    }
    const auto started = std::chrono::steady_clock::now();
    bool reportJson = false;
    for (int i = 1; i < argc; ++i)
        reportJson |= std::string_view(argv[i]) == "--report-json";
    paacs::geo::CadImportOptions import;
    paacs::geo::CadImportOutcome cad;
    paacs::algo::GeometryMeshingOptions options;
    paacs::algo::GeometryMeshingOutcome result;
    paacs::algo::MixedMeshQualityReport quality;
    paacs::mesh::Mesh mesh;
    double importSeconds = 0, archiveSeconds = 0, sizeFactor = 0, diagonal = 0;
    bool archivePublished = false;
    auto archiveReason = paacs::sdk::MeshArchiveReason::InternalError;
    const char *activeStage = "arguments";
    const auto seconds = [](const auto &begin) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    };
    const auto report = [&](const char *status, const char *stage, const char *reason) {
        if (!reportJson)
            return;
        const auto value = [](double v) {
            if (std::isfinite(v))
                std::cout << v;
            else
                std::cout << "null";
        };
        std::cout << std::setprecision(17)
                  << "{\"schema\":\"paacs-cad-report-v1\",\"status\":\"" << status
                  << "\",\"stage\":\"" << stage << "\",\"reason\":\"" << reason
                  << "\",\"mesh_published\":" << (result.published ? "true" : "false")
                  << ",\"archive_published\":" << (archivePublished ? "true" : "false")
                  << ",\"cad_import_reason\":\"" << paacs::geo::cadImportReasonName(cad.reason)
                  << "\",\"polyhedral_reason\":\""
                  << paacs::geo::polyhedralReasonName(cad.polyhedralReason)
                  << "\",\"file_bytes\":" << cad.fileBytes << ",\"solids\":" << cad.solids
                  << ",\"faces\":" << cad.faces << ",\"edges\":" << cad.edges
                  << ",\"vertices\":" << cad.vertices << ",\"curved_faces\":" << cad.curvedFaces
                  << ",\"curved_edges\":" << cad.curvedEdges
                  << ",\"working_length_unit\":";
        if (cad.millimetresPerModelUnit == 1)
            jsonString("millimetre");
        else
            std::cout << "null";
        std::cout << ",\"millimetres_per_model_unit\":";
        value(cad.millimetresPerModelUnit);
        std::cout << ",\"source_length_units\":[";
        for (std::size_t i = 0; i < cad.sourceLengthUnits.size(); ++i) {
            if (i)
                std::cout << ',';
            jsonString(cad.sourceLengthUnits[i]);
        }
        std::cout << "],\"unit_policy\":\""
                  << (import.format == paacs::geo::CadFileFormat::Step
                          ? "STEP transferred explicitly to millimetres"
                          : "BREP coordinate units; physical length unit not specified")
                  << "\",\"target_size\":";
        if (sizeFactor > 0 && !cad.imported())
            std::cout << "null";
        else
            value(options.targetSize);
        std::cout << ",\"size_factor\":";
        value(sizeFactor);
        std::cout << ",\"bounding_box_diagonal\":";
        value(diagonal);
        std::cout << ",\"nodes\":" << mesh.nodeCount() << ",\"cells\":" << mesh.cellCount()
                  << ",\"cad_volume\":";
        value(cad.cadVolume);
        std::cout << ",\"mesh_volume\":";
        value(quality.totalVolume);
        std::cout << ",\"minimum_quality\":";
        value(quality.minimumQuality);
        std::cout << ",\"output_bytes\":" << result.outputBytes
                  << ",\"import_seconds\":" << importSeconds
                  << ",\"preparation_seconds\":" << result.analysisSeconds
                  << ",\"generation_seconds\":" << result.generationSeconds
                  << ",\"audit_seconds\":" << result.auditSeconds
                  << ",\"quality_seconds\":" << result.qualitySeconds
                  << ",\"archive_seconds\":" << archiveSeconds
                  << ",\"total_seconds\":" << seconds(started)
                  << ",\"audit_temporary_bytes\":" << result.audit.temporaryBytes
                  << ",\"audit_valid\":" << (result.audit.valid ? "true" : "false")
                  << ",\"quality_valid\":" << (quality.valid ? "true" : "false")
                  << ",\"archive_reason\":" << static_cast<unsigned>(archiveReason)
                  << ",\"work_budget_scope\":\"Per-body structured node+cell construction floor, planar source-recovery and optional tet refinement/Tet10 phase counters; baseline constrained Tet4 construction has no work counter\""
                  << ",\"limits\":{\"file_bytes\":" << import.maximumFileBytes
                  << ",\"meshing_work\":" << options.maximumMeshingWork
                  << ",\"nodes\":" << options.maximumNodes
                  << ",\"cells\":" << options.maximumCells
                  << ",\"logical_output_bytes\":" << options.maximumOutputBytes
                  << ",\"audit_temporary_bytes\":" << options.maximumAuditTemporaryBytes << '}';
        if (!result.analysis.bodies.empty()) {
            const auto &body = result.analysis.bodies.front();
            std::cout << ",\"strategy\":\"" << paacs::algo::geometryStrategyName(body.strategy)
                      << "\",\"source_face\":" << body.sourceFaceTag
                      << ",\"target_face\":" << body.targetFaceTag
                      << ",\"layers\":" << body.layers << ",\"layer_ratio\":";
            value(body.endToStartLayerWidthRatio);
            std::cout << ",\"engine_status\":" << body.engineStatus
                      << ",\"engine_reason\":" << body.engineReason
                      << ",\"engine_phase\":" << body.enginePhase
                      << ",\"boundary_reason\":" << body.boundaryReason
                      << ",\"boundary_phase\":" << body.boundaryPhase
                      << ",\"recovery_kind\":" << body.recoveryKind
                      << ",\"recovery_detail\":" << body.recoveryDetail
                      << ",\"recovery_reason\":\""
                      << recoveryReasonName(body.recoveryKind, body.recoveryDetail) << '"';
            if (body.hasSurfaceDiagnostic)
                std::cout << ",\"surface_status\":" << body.surfaceStatus
                          << ",\"surface_reason\":" << body.surfaceReason
                          << ",\"surface_reason_name\":\"" << surfaceReasonName(body.surfaceReason)
                          << "\",\"surface_failure_face_persistent_id\":"
                          << body.surfaceFailureFacePersistentId
                          << ",\"surface_processed_faces\":" << body.surfaceProcessedFaces
                          << ",\"surface_unattempted_faces\":" << body.surfaceUnattemptedFaces;
        }
        std::cout << "}\n";
    };
    if (argc < 4) {
        std::cerr << "PAACS Mesh " << paacs::sdk::kApiVersion << "\n"
                  << "usage: paacs-cad-mesh INPUT.step SIZE OUTPUT.paacsmesh\n"
                  << "  [--brep] [--threads N] [--strategy auto|hex|sweep|tet]\n"
                  << "  [--max-nodes N] [--max-cells N] [--max-bytes N] [--min-quality Q]\n"
                  << "  [--max-audit-bytes N] [--layers N] [--layer-ratio R]\n"
                  << "  [--source-face TAG] [--max-layer-growth R]\n"
                  << "  [--report-json] [--max-file-bytes N] [--max-work N]\n"
                  << "  [--size-factor F] (SIZE replaced by F * solid bounding-box diagonal)\n";
        report("rejected", "arguments", "invalid-request");
        return 2;
    }
    try {
        const std::filesystem::path path(argv[1]);
        if (path.extension() == ".brep")
            import.format = paacs::geo::CadFileFormat::Brep;
        options.targetSize = number<double>(argv[2]);
        paacs::algo::GeometryBodyControl control;
        bool controlled = false;
        for (int i = 4; i < argc; ++i) {
            const std::string_view a(argv[i]);
            if (a == "--report-json")
                reportJson = true;
            else if (a == "--brep")
                import.format = paacs::geo::CadFileFormat::Brep;
            else if (a == "--max-file-bytes" && i + 1 < argc)
                import.maximumFileBytes = number<std::uint64_t>(argv[++i]);
            else if (a == "--max-work" && i + 1 < argc)
                options.maximumMeshingWork = number<std::uint64_t>(argv[++i]);
            else if (a == "--size-factor" && i + 1 < argc) {
                sizeFactor = number<double>(argv[++i]);
                if (!(sizeFactor > 0) || !std::isfinite(sizeFactor))
                    throw std::invalid_argument("size factor must be finite and positive");
            }
            else if (a == "--threads" && i + 1 < argc) {
                options.threads = number<int>(argv[++i]);
                if (options.threads < 1)
                    throw std::invalid_argument("threads must be positive");
            } else if (a == "--max-nodes" && i + 1 < argc)
                options.maximumNodes = number<std::uint64_t>(argv[++i]);
            else if (a == "--max-cells" && i + 1 < argc)
                options.maximumCells = number<std::uint64_t>(argv[++i]);
            else if (a == "--max-bytes" && i + 1 < argc)
                options.maximumOutputBytes = number<std::uint64_t>(argv[++i]);
            else if (a == "--max-audit-bytes" && i + 1 < argc)
                options.maximumAuditTemporaryBytes = number<std::uint64_t>(argv[++i]);
            else if (a == "--layers" && i + 1 < argc) {
                control.layers = number<std::uint32_t>(argv[++i]);
                controlled = true;
                if (!control.layers)
                    throw std::invalid_argument("layers must be positive");
            } else if (a == "--layer-ratio" && i + 1 < argc) {
                control.endToStartLayerWidthRatio = number<double>(argv[++i]);
                controlled = true;
            } else if (a == "--max-layer-growth" && i + 1 < argc) {
                control.maximumAdjacentLayerRatio = number<double>(argv[++i]);
                controlled = true;
            } else if (a == "--source-face" && i + 1 < argc) {
                control.sourceFaceTag = number<int>(argv[++i]);
                controlled = true;
                if (control.sourceFaceTag < 0)
                    throw std::invalid_argument("source face must be nonnegative");
            } else if (a == "--min-quality" && i + 1 < argc)
                options.minimumElementQuality = number<double>(argv[++i]);
            else if (a == "--strategy" && i + 1 < argc) {
                const std::string_view s(argv[++i]);
                using S = paacs::algo::GeometryStrategy;
                if (s == "auto")
                    options.strategy = S::Auto;
                else if (s == "hex")
                    options.strategy = S::MappedHex8;
                else if (s == "sweep")
                    options.strategy = S::LinearPrismSweep;
                else if (s == "tet")
                    options.strategy = S::BoundaryConformingTet4;
                else
                    throw std::invalid_argument("unknown strategy");
            } else
                throw std::invalid_argument("unknown option");
        }
        if (!(options.targetSize > 0) || !std::isfinite(options.targetSize))
            throw std::invalid_argument("size must be finite and positive");
        const auto importStart = std::chrono::steady_clock::now();
        activeStage = "import";
        cad = paacs::geo::importPlanarCad(path, import);
        importSeconds = seconds(importStart);
        if (!cad.imported()) {
            std::cerr << "CAD import: " << paacs::geo::cadImportReasonName(cad.reason) << '\n';
            report("refused", "import", paacs::geo::cadImportReasonName(cad.reason));
            return 3;
        }
        activeStage = "preparation";
        paacs::Vec3 lo, hi;
        cad.model->solids().at(0)->bbox(lo, hi);
        const auto extent = hi - lo;
        diagonal = std::hypot(extent.x, extent.y, extent.z);
        if (sizeFactor > 0) {
            options.targetSize = diagonal * sizeFactor;
            if (!(options.targetSize > 0) || !std::isfinite(options.targetSize))
                throw std::invalid_argument("relative size produced an invalid target");
        }
        if (controlled) {
            control.solidTag = cad.model->solids().at(0)->tag();
            options.bodyControls.push_back(control);
        }
        activeStage = "meshing";
        result = paacs::algo::meshGeometryModel(*cad.model, options, mesh);
        if (!result.published) {
            std::cerr << "meshing: " << paacs::algo::geometryReasonName(result.reason);
            for (const auto& body : result.analysis.bodies)
                if (body.hasSurfaceDiagnostic && body.surfaceReason != 0) {
                    std::cerr << "; surface: " << surfaceReasonName(body.surfaceReason)
                              << "; face persistent id: " << body.surfaceFailureFacePersistentId;
                    if (body.recoveryKind != 0)
                        std::cerr << "; recovery: "
                                  << recoveryReasonName(body.recoveryKind, body.recoveryDetail);
                    break;
                }
            std::cerr << '\n';
            report("refused", "meshing", paacs::algo::geometryReasonName(result.reason));
            return 4;
        }
        activeStage = "quality";
        const auto qualityStart = std::chrono::steady_clock::now();
        quality =
            result.qualityAssessed ? result.quality : paacs::algo::analyzeMixedMeshQuality(mesh);
        if (!result.qualityAssessed)
            result.qualitySeconds = seconds(qualityStart);
        if (!quality.valid) {
            std::cerr << "quality audit failed\n";
            report("failed", "quality", "quality-audit-failed");
            return 5;
        }
        activeStage = "archive";
        const auto archiveStart = std::chrono::steady_clock::now();
        const auto archive = paacs::sdk::writeMeshArchive(mesh, argv[3]);
        archiveReason = archive.reason;
        archiveSeconds = seconds(archiveStart);
        if (!archive.completed) {
            std::cerr << "archive write failed\n";
            report("failed", "archive",
                   archive.reason == paacs::sdk::MeshArchiveReason::AllocationFailed
                       ? "allocation-failed"
                       : "archive-write-failed");
            return 6;
        }
        archivePublished = true;
        if (reportJson) {
            report("complete", "archive", "none");
            return 0;
        }
        const auto &body = result.analysis.bodies.front();
        std::cout << std::setprecision(17) << "{\"status\":\"complete\",\"strategy\":\""
                  << paacs::algo::geometryStrategyName(result.analysis.bodies.front().strategy)
                  << "\",\"nodes\":" << mesh.nodeCount() << ",\"cells\":" << mesh.cellCount()
                  << ",\"cad_volume\":" << cad.cadVolume
                  << ",\"mesh_volume\":" << quality.totalVolume
                  << ",\"minimum_quality\":" << quality.minimumQuality
                  << ",\"output_bytes\":" << result.outputBytes
                  << ",\"generation_seconds\":" << result.generationSeconds
                  << ",\"audit_seconds\":" << result.auditSeconds
                  << ",\"audit_temporary_bytes\":" << result.audit.temporaryBytes
                  << ",\"source_face\":" << body.sourceFaceTag
                  << ",\"target_face\":" << body.targetFaceTag << ",\"layers\":" << body.layers
                  << ",\"layer_ratio\":" << body.endToStartLayerWidthRatio << "}\n";
        return 0;
    } catch (const std::invalid_argument &e) {
        std::cerr << e.what() << '\n';
        report("rejected", "arguments", "invalid-request");
        return 2;
    } catch (const std::bad_alloc &e) {
        std::cerr << e.what() << '\n';
        report("failed", activeStage, "allocation-failed");
        return 2;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        report("failed", activeStage, "kernel-failure");
        return 2;
    }
}
