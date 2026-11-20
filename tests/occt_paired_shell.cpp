// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctModelInternal.hpp"
#include "OcctPairedShell.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <cmath>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace paacs;
using namespace paacs::geo::detail;
void need(bool ok, const char *why) {
    if (!ok)
        throw std::runtime_error(why);
}
std::string snapshot(const geo::OcctModel &m) {
    std::ostringstream out;
    BRepTools::Write(OcctModelAccess::shape(m), out);
    return out.str();
}
mesh::Mesh sentinel() {
    mesh::Mesh m;
    m.addNode({17, 19, 23}, {0, 777});
    return m;
}
void unchanged(const mesh::Mesh &m) {
    need(m.nodeCount() == 1 && m.cellCount() == 0 && m.node(0).x == 17 && m.node(0).y == 19 &&
             m.node(0).z == 23 && m.nodeClass(0) == mesh::Classification{0, 777},
         "paired-shell rollback");
}
void validate(const geo::OcctModel &m, const OcctPairedShellOutcome &o, const mesh::Mesh &mesh,
              double chord) {
    need(o.published && !o.destinationUnchanged && o.reason == OcctPairedShellReason::None,
         o.diagnostic);
    need(o.sourceImmutabilityVerified && o.audit.complete && o.audit.valid &&
             o.exteriorCertificate.certified && o.quality.complete && o.quality.valid,
         "paired-shell full final gates");
    need(o.workUsed == o.recognitionWork + o.surfaceWork + o.sweepWork + o.validationWork,
         "paired-shell shared work receipt");
    need(o.maximumBoundarySampleError <= chord && o.maximumBoundaryNodeError <= 1e-6 &&
             o.relativeVolumeError <= .01,
         "paired-shell original fidelity");
    need(o.nativeInteriorQueries >= mesh.nodeCount() + mesh.countCells(mesh::CellType::Tet4),
         "all native nodes and tet centroids checked");
    for (Index i = 0; i < mesh.nodeCount(); ++i) {
        const auto c = mesh.nodeClass(i);
        need(c.dim == 0   ? m.vertex(c.tag) != nullptr
             : c.dim == 1 ? m.edge(c.tag) != nullptr
             : c.dim == 2 ? m.face(c.tag) != nullptr
                          : c.dim == 3 && m.solid(c.tag) != nullptr,
             "paired-shell original classification");
    }
}
int main(int argc, char **argv) {
    try {
        OcctPairedShellOptions options;
        options.boundary.linearDeflection = .04;
        options.boundary.maximumChordError = .1;
        options.maximumWork = 20000000;
        auto cube = OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2, 3, 4).Shape());
        need(cube.imported(), "native cube");
        auto dst = sentinel();
        auto unsupported = buildOcctPairedShellVolume(*cube.model, 1, options, dst);
        need(!unsupported.published &&
                 unsupported.reason == OcctPairedShellReason::UnsupportedGeometry,
             "plain solid unsupported");
        unchanged(dst);
        const double pi = std::acos(-1.);
        auto outer = BRepPrimAPI_MakeSphere(10, 0., pi / 2, 2 * pi).Shape();
        auto inner = BRepPrimAPI_MakeSphere(8, 0., pi / 2, 2 * pi).Shape();
        BRepAlgoAPI_Cut cut(outer, inner);
        need(cut.IsDone(), "native spherical shell construction");
        auto source = OcctModelAccess::fromShape(cut.Shape());
        need(source.imported(), "native spherical shell");
        const auto before = snapshot(*source.model);
        mesh::Mesh mesh;
        auto result = buildOcctPairedShellVolume(*source.model, 1, options, mesh);
        validate(*source.model, result, mesh, options.boundary.maximumChordError);
        need(result.nativeValidationThreads == 1, "paired shell default serial native cap");
        need(snapshot(*source.model) == before, "paired shell source unchanged");
        need(std::abs(result.cadVolume - 2 * pi * (1000 - 512) / 3) < 1e-6,
             "spherical shell accurate mass");
        auto reject = [&](OcctPairedShellOptions limits, OcctPairedShellReason reason) {
            auto dst = sentinel();
            auto out = buildOcctPairedShellVolume(*source.model, 1, limits, dst);
            need(!out.published && out.reason == reason, out.diagnostic);
            unchanged(dst);
            need(snapshot(*source.model) == before, "refusal source unchanged");
        };
        auto limited = options;
        limited.maximumWork = 0;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumWork = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumWorkingBytes = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumSnapshotBytes = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumNodes = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumCells = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.maximumOutputBytes = 1;
        reject(limited, OcctPairedShellReason::BudgetExceeded);
        limited = options;
        limited.boundary.linearDeflection = std::numeric_limits<double>::quiet_NaN();
        reject(limited, OcctPairedShellReason::InvalidRequest);
        limited = options;
        limited.maximumNativeValidationThreads = 0;
        reject(limited, OcctPairedShellReason::InvalidRequest);
        std::atomic_bool cancelled{true};
        limited = options;
        limited.cancellation = &cancelled;
        reject(limited, OcctPairedShellReason::Cancelled);
        limited = options;
        limited.minimumQuality = 1;
        reject(limited, OcctPairedShellReason::QualityRequirement);
        limited = options;
        limited.maximumLayerThickness = .9;
        limited.maximumNativeValidationThreads = 2;
        mesh::Mesh layered;
        auto layeredResult = buildOcctPairedShellVolume(*source.model, 1, limited, layered);
        validate(*source.model, layeredResult, layered, limited.boundary.maximumChordError);
        need(layeredResult.layers >= 3, "automatic shell thickness layers");
        need(layeredResult.nativeValidationThreads == 2, "paired shell explicit native cap");
        gp_Trsf transform;
        transform.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), .43);
        transform.SetTranslationPart(gp_Vec(17, -29, 43));
        BRepBuilderAPI_Transform moved(cut.Shape(), transform, true);
        auto rotated = OcctModelAccess::fromShape(moved.Shape());
        need(rotated.imported(), "rotated shell");
        mesh::Mesh rotatedMesh;
        validate(*rotated.model,
                 buildOcctPairedShellVolume(*rotated.model, 1, options, rotatedMesh), rotatedMesh,
                 options.boundary.maximumChordError);
        if (argc >= 2) {
            auto actual = geo::importOcctCad(argv[1]);
            need(actual.imported(), "actual paired native shell import");
            Vec3 lo, hi;
            actual.model->solid(1)->bbox(lo, hi);
            auto actualOptions = options;
            actualOptions.maximumWork = 100000000;
            actualOptions.maximumNativeValidationThreads = 4;
            actualOptions.maximumWorkingBytes = 1024ULL * 1024 * 1024;
            actualOptions.maximumSnapshotBytes = 64ULL * 1024 * 1024;
            actualOptions.boundary.maximumWork = 100000000;
            actualOptions.boundary.maximumChordError = (hi - lo).norm() * .0025;
            actualOptions.boundary.linearDeflection = actualOptions.boundary.maximumChordError * .1;
            mesh::Mesh actualMesh;
            auto out = buildOcctPairedShellVolume(*actual.model, 1, actualOptions, actualMesh);
            std::cout << std::setprecision(17)
                      << "{\"published\":" << (out.published ? "true" : "false")
                      << ",\"reason\":" << static_cast<unsigned>(out.reason) << ",\"diagnostic\":\""
                      << out.diagnostic << "\",\"source_unmodified\":"
                      << (out.sourceImmutabilityVerified ? "true" : "false")
                      << ",\"nodes\":" << out.nodes << ",\"cells\":" << out.cells
                      << ",\"work\":" << out.workUsed
                      << ",\"native_queries\":" << out.nativeInteriorQueries
                      << ",\"native_threads\":" << out.nativeValidationThreads
                      << ",\"cad_volume\":" << out.cadVolume
                      << ",\"mesh_volume\":" << out.meshVolume
                      << ",\"relative_volume_error\":" << out.relativeVolumeError
                      << ",\"node_error\":" << out.maximumBoundaryNodeError
                      << ",\"sample_error\":" << out.maximumBoundarySampleError
                      << ",\"chord\":" << actualOptions.boundary.maximumChordError
                      << ",\"quality\":" << out.quality.minimumQuality
                      << ",\"recognition_seconds\":" << out.recognitionSeconds
                      << ",\"surface_seconds\":" << out.surfaceSeconds
                      << ",\"construction_seconds\":" << out.sweepSeconds
                      << ",\"validation_seconds\":" << out.validationSeconds
                      << ",\"total_seconds\":" << out.totalSeconds << "}\n"
                      << std::flush;
            validate(*actual.model, out, actualMesh, actualOptions.boundary.maximumChordError);
            if (argc >= 3) {
                need(sdk::writeMeshArchive(actualMesh, argv[2]).completed, "actual paired archive");
                mesh::Mesh reopened;
                need(sdk::readMeshArchive(argv[2], reopened).completed,
                     "actual paired archive read");
                std::vector<std::uint8_t> originalBytes, reopenedBytes;
                need(sdk::serializeMeshArchive(actualMesh, originalBytes).completed &&
                         sdk::serializeMeshArchive(reopened, reopenedBytes).completed &&
                         originalBytes == reopenedBytes,
                     "actual paired archive exact roundtrip");
            }
            std::cout << "actual paired shell nodes=" << out.nodes << " cells=" << out.cells
                      << " work=" << out.workUsed
                      << " chord=" << actualOptions.boundary.maximumChordError
                      << " error=" << out.maximumBoundarySampleError
                      << " quality=" << out.quality.minimumQuality << '\n';
        }
        std::cout << "Native paired shell topology, fidelity, ownership, work, and rollback checks "
                     "passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    } catch (...) {
        return 1;
    }
}
