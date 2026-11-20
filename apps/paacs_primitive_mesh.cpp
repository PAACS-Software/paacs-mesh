// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/PrimitiveMesher.hpp"
#include "paacs/algo/SurfaceFeatures.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include "NumberArgument.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>
namespace {
template <class T> T number(const char *text) {
    const std::string_view s(text);
    T v{};
    if (!paacs::cli::parseNumber(s, v))
        throw std::invalid_argument("invalid number");
    return v;
}
} // namespace
int main(int argc, char **argv) {
    using namespace paacs;
    if (argc < 4) {
        std::cerr << "PAACS " << sdk::kApiVersion
                  << "\nusage: paacs-primitive-mesh SHAPE SIZE OUTPUT.paacsmesh\n"
                  << "  SHAPE: cylinder | cone | frustum | tube | torus | elbow\n"
                  << "  [--radius R] [--height H] [--end-radius R] [--inner-radius R]\n"
                  << "  [--major-radius R] [--angle-degrees A] [--threads N]\n"
                  << "  [--elements hex-dominant|tet4|tet10] [--cone-strategy auto|swept|graded]\n"
                  << "  [--cylinder-strategy auto|swept|uniform]\n"
                  << "  [--circular-strategy auto|swept|uniform]\n"
                  << "  [--detect-features on|off] [--feature-angle DEGREES]\n"
                  << "  [--feature-angle-tolerance DEGREES] [--feature-min-length LENGTH]\n"
                  << "  [--feature-tag-boundaries on|off] [--feature-max-edges N]\n"
                  << "  [--warn-quality Q] [--quality-policy warn|reject]\n"
                  << "  [--chord-error E] [--volume-error E] [--min-quality Q]\n"
                  << "  [--max-nodes N] [--max-cells N] [--max-bytes N] [--max-working-bytes N]\n";
        return 2;
    }
    try {
        algo::PrimitiveShape s;
        bool found = false;
        for (auto k : {algo::PrimitiveKind::Cylinder, algo::PrimitiveKind::Cone,
                       algo::PrimitiveKind::Frustum, algo::PrimitiveKind::Tube,
                       algo::PrimitiveKind::Torus, algo::PrimitiveKind::Elbow})
            if (std::string_view(argv[1]) == algo::primitiveKindName(k)) {
                s.kind = k;
                found = true;
            }
        if (!found)
            throw std::invalid_argument("unknown shape");
        if (s.kind == algo::PrimitiveKind::Tube)
            s.innerRadius = .8;
        algo::PrimitiveMeshingOptions o;
        algo::SurfaceFeatureOptions features;
        o.targetSize = number<double>(argv[2]);
        for (int i = 4; i < argc; ++i) {
            const std::string_view a(argv[i]);
            if (i + 1 >= argc)
                throw std::invalid_argument("missing option value");
            const char *v = argv[++i];
            if (a == "--radius")
                s.radius = number<double>(v);
            else if (a == "--height")
                s.height = number<double>(v);
            else if (a == "--end-radius")
                s.endRadius = number<double>(v);
            else if (a == "--inner-radius")
                s.innerRadius = number<double>(v);
            else if (a == "--major-radius")
                s.majorRadius = number<double>(v);
            else if (a == "--angle-degrees")
                s.sweepAngle = number<double>(v) * std::numbers::pi / 180;
            else if (a == "--threads")
                o.threads = number<int>(v);
            else if (a == "--chord-error")
                o.maximumChordError = number<double>(v);
            else if (a == "--volume-error")
                o.maximumRelativeVolumeError = number<double>(v);
            else if (a == "--min-quality")
                o.minimumQuality = number<double>(v);
            else if (a == "--elements") {
                bool selected = false;
                for (auto e : {algo::PrimitiveElements::HexDominant, algo::PrimitiveElements::Tet4,
                               algo::PrimitiveElements::Tet10})
                    if (std::string_view(v) == algo::primitiveElementsName(e)) {
                        o.elements = e;
                        selected = true;
                    }
                if (!selected)
                    throw std::invalid_argument("unknown element family");
            } else if (a == "--cone-strategy") {
                const std::string_view value(v);
                if (value == "auto")
                    o.coneStrategy = algo::PrimitiveConeStrategy::Automatic;
                else if (value == "swept")
                    o.coneStrategy = algo::PrimitiveConeStrategy::Swept;
                else if (value == "graded")
                    o.coneStrategy = algo::PrimitiveConeStrategy::Graded;
                else
                    throw std::invalid_argument("unknown cone strategy");
            } else if (a == "--cylinder-strategy") {
                const std::string_view value(v);
                if (value == "auto")
                    o.cylinderStrategy = algo::PrimitiveCylinderStrategy::Automatic;
                else if (value == "swept")
                    o.cylinderStrategy = algo::PrimitiveCylinderStrategy::Swept;
                else if (value == "uniform")
                    o.cylinderStrategy = algo::PrimitiveCylinderStrategy::Uniform;
                else
                    throw std::invalid_argument("unknown cylinder strategy");
            } else if (a == "--circular-strategy") {
                const std::string_view value(v);
                if (value == "auto")
                    o.circularStrategy = algo::PrimitiveCircularStrategy::Automatic;
                else if (value == "swept")
                    o.circularStrategy = algo::PrimitiveCircularStrategy::Swept;
                else if (value == "uniform")
                    o.circularStrategy = algo::PrimitiveCircularStrategy::Uniform;
                else
                    throw std::invalid_argument("unknown circular strategy");
            } else if (a == "--detect-features" || a == "--feature-tag-boundaries") {
                const std::string_view value(v);
                if (value != "on" && value != "off")
                    throw std::invalid_argument("expected on or off");
                (a == "--detect-features" ? features.enabled : features.includeTagBoundaries) =
                    value == "on";
            } else if (a == "--feature-angle")
                features.angleDegrees = number<double>(v);
            else if (a == "--feature-angle-tolerance")
                features.angularToleranceDegrees = number<double>(v);
            else if (a == "--feature-min-length")
                features.minimumLength = number<double>(v);
            else if (a == "--feature-max-edges")
                features.maximumEdges = number<std::uint64_t>(v);
            else if (a == "--warn-quality")
                o.warningQuality = number<double>(v);
            else if (a == "--quality-policy") {
                const std::string_view value(v);
                if (value == "warn")
                    o.qualityPolicy = algo::PrimitiveQualityPolicy::Warn;
                else if (value == "reject")
                    o.qualityPolicy = algo::PrimitiveQualityPolicy::Reject;
                else
                    throw std::invalid_argument("unknown quality policy");
            } else if (a == "--max-nodes")
                o.maximumNodes = number<std::uint64_t>(v);
            else if (a == "--max-cells")
                o.maximumCells = number<std::uint64_t>(v);
            else if (a == "--max-bytes")
                o.maximumOutputBytes = number<std::uint64_t>(v);
            else if (a == "--max-working-bytes")
                o.maximumWorkingBytes = number<std::uint64_t>(v);
            else
                throw std::invalid_argument("unknown option");
        }
        mesh::Mesh m;
        const auto r = algo::meshPrimitive(s, o, m);
        if (!r.published) {
            std::cerr << "meshing: " << algo::primitiveReasonName(r.reason) << '\n';
            if (r.reason == algo::PrimitiveReason::QualityRequirement)
                std::cerr << "minimum mean ratio " << r.minimumQuality << "; worst cell "
                          << r.quality.worstCell << "; " << r.quality.lowQualityCells
                          << " cells below warning threshold\n";
            return 1;
        }
        if (r.qualityWarning)
            std::cerr << "quality warning: " << r.quality.lowQualityCells << " volume cells below "
                      << o.warningQuality << "; worst cell " << r.quality.worstCell
                      << " has mean ratio " << r.minimumQuality << '\n';
        // Optional analysis precedes archive publication, so a failed feature
        // request cannot replace an existing user's output file.
        features.maximumWorkingBytes = o.maximumWorkingBytes;
        features.maximumOutputBytes = std::min(features.maximumOutputBytes, o.maximumOutputBytes);
        const auto featureBegin = std::chrono::steady_clock::now();
        const auto featureResult = algo::detectSurfaceFeatures(m, features);
        const double featureSeconds =
            features.enabled
                ? std::chrono::duration<double>(std::chrono::steady_clock::now() - featureBegin)
                      .count()
                : 0;
        if (!featureResult.completed) {
            std::cerr << "feature analysis: "
                      << algo::surfaceFeatureReasonName(featureResult.reason) << '\n';
            return 1;
        }
        if (!sdk::writeMeshArchive(m, argv[3]).completed) {
            std::cerr << "cannot write archive\n";
            return 1;
        }
        std::cout
            << std::setprecision(17) << "{\"shape\":\"" << algo::primitiveKindName(s.kind)
            << "\",\"nodes\":" << r.nodes << ",\"hexes\":" << r.hexes << ",\"element_family\":\""
            << algo::primitiveElementsName(o.elements) << "\""
            << ",\"tet4s\":" << r.tet4s << ",\"tet10s\":" << r.tet10s
            << ",\"certified_tet10s\":" << r.certifiedTet10s
            << ",\"quadratic_connectivity_verified\":"
            << (r.quadraticConnectivityVerified ? "true" : "false")
            << ",\"graded_cone\":" << (r.gradedCone ? "true" : "false")
            << ",\"uniform_cylinder\":" << (r.uniformCylinder ? "true" : "false")
            << ",\"cylinder_geometry_limited\":" << (r.cylinderGeometryLimited ? "true" : "false")
            << ",\"cylinder_radial_spacing\":" << r.cylinderRadialSpacing
            << ",\"cylinder_axial_spacing\":" << r.cylinderAxialSpacing
            << ",\"uniform_circular\":" << (r.uniformCircular ? "true" : "false")
            << ",\"circular_geometry_limited\":" << (r.circularGeometryLimited ? "true" : "false")
            << ",\"circular_radial_spacing\":" << r.circularRadialSpacing
            << ",\"circular_inner_spacing\":" << r.circularInnerSpacing
            << ",\"circular_outer_spacing\":" << r.circularOuterSpacing
            << ",\"feature_analysis_enabled\":" << (features.enabled ? "true" : "false")
            << ",\"feature_edges\":" << featureResult.metrics.selectedEdges
            << ",\"sharp_edges\":" << featureResult.metrics.sharpEdges
            << ",\"feature_filtered_short_edges\":" << featureResult.metrics.filteredShortEdges
            << ",\"feature_working_bytes\":" << featureResult.metrics.workingBytes
            << ",\"feature_output_bytes\":" << featureResult.metrics.outputBytes
            << ",\"feature_analysis_seconds\":" << featureSeconds
            << ",\"cone_refinement_levels\":" << r.coneRefinementLevels
            << ",\"quality_warning\":" << (r.qualityWarning ? "true" : "false")
            << ",\"low_quality_cells\":" << r.quality.lowQualityCells
            << ",\"sharp_feature_low_quality_cells\":" << r.quality.sharpFeatureLowQualityCells
            << ",\"worst_cell\":" << r.quality.worstCell << ",\"pyramids\":" << r.pyramids
            << ",\"boundary_quads\":" << r.boundaryQuads
            << ",\"boundary_triangles\":" << r.boundaryTriangles
            << ",\"output_bytes\":" << r.outputBytes
            << ",\"planned_working_bytes\":" << r.plannedWorkingBytes
            << ",\"mesh_volume\":" << r.meshVolume << ",\"analytic_volume\":" << r.analyticVolume
            << ",\"relative_volume_error\":" << r.relativeVolumeError
            << ",\"chord_error_bound\":" << r.chordErrorBound
            << ",\"minimum_quality\":" << r.minimumQuality
            << ",\"generation_seconds\":" << r.generationSeconds
            << ",\"validation_seconds\":" << r.validationSeconds
            << ",\"audit_seconds\":" << r.auditSeconds
            << ",\"audit_passed\":true,\"quality_examples\":[";
        for (unsigned i = 0; i < r.quality.exampleCount; ++i) {
            if (i)
                std::cout << ',';
            const auto &q = r.quality.examples[i];
            std::cout << "{\"cell\":" << q.cell << ",\"quality\":" << q.quality
                      << ",\"sharp_feature\":" << (q.touchesSharpFeature ? "true" : "false") << '}';
        }
        std::cout << "],\"quality_histogram\":[";
        for (unsigned i = 0; i < r.quality.histogram.size(); ++i) {
            if (i)
                std::cout << ',';
            std::cout << r.quality.histogram[i];
        }
        std::cout << "],\"quality_examples_truncated\":"
                  << (r.quality.examplesTruncated ? "true" : "false") << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
