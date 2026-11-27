// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/geometry/ExperimentalCadVolume.hpp"
#include "paacs/geometry/OcctBoundary.hpp"
#include "../src/geometry/OcctModelInternal.hpp"
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <gp_Pnt.hxx>
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace paacs;
namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
auto key(Vec3 p) {
    return std::array{std::bit_cast<std::uint64_t>(p.x),
        std::bit_cast<std::uint64_t>(p.y), std::bit_cast<std::uint64_t>(p.z)};
}
void sameMesh(const mesh::Mesh& a, const mesh::Mesh& b) {
    require(a.nodeCount() == b.nodeCount() && a.cellCount() == b.cellCount(), "mesh counts changed");
    for (Index i = 0; i < a.nodeCount(); ++i)
        require(key(a.node(i)) == key(b.node(i)) && a.nodeClass(i) == b.nodeClass(i),
            "original point or classification changed");
    for (Index i = 0; i < a.cellCount(); ++i) {
        require(a.cell(i).type == b.cell(i).type && a.cell(i).onEntity == b.cell(i).onEntity,
            "original cell classification changed");
        require(std::equal(a.cellNodes(i), a.cellNodes(i) + mesh::nodesPerCell(a.cell(i).type),
            b.cellNodes(i)), "original connectivity changed");
    }
}
void qualified(const geo::ExperimentalCadVolumeOutcome& out) {
    if (!out.native.published)
        throw std::runtime_error(std::string("volume refused: ") + out.experimental.diagnostic +
            " / " + geo::nativeCadMeshReasonName(out.native.reason));
    require(out.native.audit.complete && out.native.audit.valid, "Native audit not qualified");
    require(out.native.quality.complete && out.native.quality.valid, "Native quality not qualified");
    require(out.experimental.boundaryRevalidated && out.experimental.boundaryUnchanged,
        "boundary not revalidated and preserved");
    require(std::abs(out.native.quality.totalVolume - 1.) < 1e-8, "independent box volume mismatch");
}
void placed(const geo::ExperimentalCadInitialPointsOutcome& out) {
    require(out.generated && out.reason == geo::NativeCadMeshReason::None,
        "initial placement checkpoint refused");
    require(out.stages.initialPointsAvailable && out.stages.boundaryRevalidated && out.stages.boundaryUnchanged,
        "initial placement did not revalidate and preserve its boundary");
    require(out.stages.recoveryWork == 0 && out.stages.recoverySeconds == 0 &&
        out.stages.fillSeconds == 0 && out.stages.validationSeconds == 0,
        "point generation invoked downstream volume stages");
    require(out.work >= out.stages.sourceWork + out.stages.placementWork &&
        out.seconds >= out.stages.placementSeconds && out.stages.totalSeconds == out.seconds,
        "point generation cost scopes are inconsistent");
    for (const auto point : out.stages.initialPoints)
        require(std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) &&
            point.x > 0 && point.x < 1 && point.y > 0 && point.y < 1 && point.z > 0 && point.z < 1,
            "generated box point is nonfinite or not strictly inside");
}
}
int main() {
    try {
        auto imported = geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1., 1., 1.).Shape());
        require(imported.imported(), "owned box import failed");
        auto options = geo::normalNativeCadMeshOptions();
        options.geometry.targetSize = .4;
        options.geometry.threads = 1;
        options.maximumNativeValidationThreads = 1;
        options.recoveryMode = geo::NativeCadRecoveryMode::Advanced;
        options.useAdvancingInteriorFront = true;
        options.useSpatialSizing = true;
        options.construction.interiorSites = 16;
        options.initialTopologyQuality = false;
        geo::ExperimentalCadBoundary boundary;
        auto prepared = geo::prepareExperimentalCadBoundary(*imported.model, options, boundary);
        if (!prepared.prepared) throw std::runtime_error(prepared.diagnostic);
        require(prepared.nodes > 0 && prepared.triangles > 0, "empty prepared shell");
        const auto frozen = boundary.boundary;

        const auto generated = geo::generateExperimentalCadInitialPoints(*imported.model, boundary, options);
        placed(generated);
        require(!generated.stages.initialPoints.empty(), "ordinary owned placement unexpectedly empty");
        sameMesh(frozen, boundary.boundary);

        mesh::Mesh classical;
        auto reference = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Classical, {}, classical);
        qualified(reference);
        require(reference.experimental.initialPointsAvailable && reference.experimental.initialPoints.size() ==
            generated.stages.initialPoints.size(), "early checkpoint differs from classical cloud count");
        for (std::size_t i = 0; i < generated.stages.initialPoints.size(); ++i)
            require(key(generated.stages.initialPoints[i]) == key(reference.experimental.initialPoints[i]),
                "early checkpoint changed classical point bits or ordering");
        sameMesh(frozen, boundary.boundary);

        mesh::Mesh generatedReplay;
        const auto reconstructed = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Supplied, generated.stages.initialPoints, generatedReplay);
        qualified(reconstructed);
        require(reconstructed.experimental.allSuppliedPointsRetained && reconstructed.experimental.retainedPoints ==
            generated.stages.initialPoints.size(), "generated cloud cannot replay with all points retained");
        sameMesh(frozen, boundary.boundary);

        auto boundaryOnly = options;
        boundaryOnly.geometry.maximumNodes = prepared.nodes;
        const auto empty = geo::generateExperimentalCadInitialPoints(*imported.model, boundary, boundaryOnly);
        placed(empty);
        require(empty.stages.initialPoints.empty(), "zero available interior capacity did not return an explicit empty cloud");
        sameMesh(frozen, boundary.boundary);

        const std::array<Vec3, 1> center{{{.5, .5, .5}}};
        mesh::Mesh supplied;
        auto seeded = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Supplied, center, supplied);
        qualified(seeded);
        require(seeded.experimental.allSuppliedPointsRetained && seeded.experimental.retainedPoints == 1,
            "explicit center was dropped");
        sameMesh(frozen, boundary.boundary);

        auto second = geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(1., 1., 1.).Shape());
        require(second.imported() && second.model->identity() != imported.model->identity(),
            "fresh import did not get a distinct identity");
        geo::ExperimentalCadBoundary rebound;
        auto reloaded = geo::rebindExperimentalCadBoundary(*second.model, frozen, options, rebound);
        require(reloaded.prepared && reloaded.rebound, "valid shell cannot rebind to fresh identical CAD");
        mesh::Mesh replayed;
        qualified(geo::replayExperimentalCadVolume(*second.model, rebound, options,
            geo::ExperimentalInteriorMode::Supplied, center, replayed));

        const auto previous = supplied;
        const std::array<Vec3, 1> outside{{{1.5, .5, .5}}};
        auto failed = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Supplied, outside, supplied);
        require(!failed.native.published, "outside point accepted");
        sameMesh(previous, supplied);
        const std::array<Vec3, 2> duplicate{{center[0], center[0]}};
        failed = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Supplied, duplicate, supplied);
        require(!failed.native.published, "duplicate points accepted");
        sameMesh(previous, supplied);
        const std::array<Vec3, 1> nonfinite{{{std::numeric_limits<double>::quiet_NaN(), .5, .5}}};
        failed = geo::replayExperimentalCadVolume(*imported.model, boundary, options,
            geo::ExperimentalInteriorMode::Supplied, nonfinite, supplied);
        require(!failed.native.published, "nonfinite point accepted");
        sameMesh(previous, supplied);

        auto altered = geo::detail::OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(.9, 1., 1.).Shape());
        require(altered.imported(), "altered owned CAD import failed");
        auto rejected = geo::rebindExperimentalCadBoundary(*altered.model, frozen, options, rebound);
        require(!rejected.prepared, "shell rebound to geometrically incompatible CAD");
        sameMesh(frozen, rebound.boundary);

        auto tiny = options;
        tiny.geometry.maximumMeshingWork = 1;
        tiny.maximumValidationWork = 1;
        rejected = geo::prepareExperimentalCadBoundary(*imported.model, tiny, boundary);
        require(!rejected.prepared, "exhausted preparation budget accepted");
        sameMesh(frozen, boundary.boundary);

        auto declined = geo::generateExperimentalCadInitialPoints(*imported.model, boundary, tiny);
        require(!declined.generated && !declined.stages.initialPointsAvailable && declined.stages.initialPoints.empty() &&
            declined.reason == geo::NativeCadMeshReason::BudgetExceeded, "exhausted point stage published a cloud");
        sameMesh(frozen, boundary.boundary);
        auto mismatched = options;
        mismatched.geometry.targetSize = .41;
        declined = geo::generateExperimentalCadInitialPoints(*imported.model, boundary, mismatched);
        require(!declined.generated && !declined.stages.initialPointsAvailable && declined.stages.initialPoints.empty(),
            "point generator changed the frozen sizing contract");
        std::atomic_bool cancelled{true};
        auto interrupted = options;
        interrupted.geometry.cancellation = &cancelled;
        declined = geo::generateExperimentalCadInitialPoints(*imported.model, boundary, interrupted);
        require(!declined.generated && declined.reason == geo::NativeCadMeshReason::Cancelled &&
            !declined.stages.initialPointsAvailable && declined.stages.initialPoints.empty(),
            "cancelled point generator published a partial cloud");
        declined = geo::generateExperimentalCadInitialPoints(*altered.model, boundary, options);
        require(!declined.generated && !declined.stages.initialPointsAvailable,
            "point generator admitted a shell from incompatible CAD");
        sameMesh(frozen, boundary.boundary);

        // Reversing a tiny cavity changes volume by only .2%, below the 1%
        // acceptance bound. CAD face-use orientation must reject it directly.
        const auto cavityShape = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(1., 1., 1.).Shape(),
            BRepPrimAPI_MakeBox(gp_Pnt(.45, .45, .45), .1, .1, .1).Shape()).Shape();
        auto cavity = geo::detail::OcctModelAccess::fromShape(cavityShape);
        require(cavity.imported(), "owned enclosed-cavity import failed");
        geo::OcctBoundaryOptions shellOptions;
        mesh::Mesh cavityShell;
        const auto shell = geo::buildOcctBoundary(*cavity.model, cavity.model->solids().front()->tag(),
            shellOptions, cavityShell);
        require(shell.published, "owned enclosed-cavity boundary failed");
        geo::ExperimentalCadBoundary cavitySnapshot;
        require(geo::rebindExperimentalCadBoundary(*cavity.model, cavityShell, options, cavitySnapshot).prepared,
            "correctly oriented enclosed-cavity shell failed rebind");
        mesh::Mesh reversedCavity;
        for (Index i = 0; i < cavityShell.nodeCount(); ++i)
            reversedCavity.addNode(cavityShell.node(i), cavityShell.nodeClass(i));
        unsigned reversed = 0;
        for (Index i = 0; i < cavityShell.cellCount(); ++i) {
            const auto* ids = cavityShell.cellNodes(i);
            Vec3 centroid{};
            for (unsigned j = 0; j < 3; ++j) centroid = centroid + cavityShell.node(ids[j]) / 3.;
            const bool interior = centroid.x > .4 && centroid.x < .6 && centroid.y > .4 &&
                centroid.y < .6 && centroid.z > .4 && centroid.z < .6;
            if (interior) ++reversed;
            reversedCavity.addCell(mesh::CellType::Tri3,
                {ids[0], ids[interior ? 2 : 1], ids[interior ? 1 : 2]}, cavityShell.cell(i).onEntity);
        }
        require(reversed > 0, "cavity reversal exercised no triangles");
        const auto cavityBefore = cavitySnapshot.boundary;
        rejected = geo::rebindExperimentalCadBoundary(*cavity.model, reversedCavity, options, cavitySnapshot);
        require(!rejected.prepared, "reversed cavity hidden by small relative volume was accepted");
        sameMesh(cavityBefore, cavitySnapshot.boundary);
        std::cout << "Experimental boundary capture, early ordered point checkpoint, fresh rebind, classical/supplied volume, "
                     "original CAD volume, fixed boundary and transactional refusals passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
