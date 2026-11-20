// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/algo/LayeredPatchMesher.hpp>
#include <paacs/sdk/MeshArchive.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
constexpr double side = 10.0;
constexpr double length = 40.0;
constexpr double growth = 1.22;

std::vector<double> symmetricStations(double firstHeight, unsigned wallRows, double coreCap) {
    paacs::algo::LayerScheduleRequest request;
    request.firstLayerThickness = firstHeight;
    request.wallLayers = wallRows;
    request.growthRatio = growth;
    request.maximumAllowedGrowth = growth + 1e-12;
    request.totalDepth = side / 2;
    request.maximumCoreThickness = coreCap;
    request.maximumLayers = 48;
    const auto half = paacs::algo::planLayerSchedule(request);
    if (!half.planned) throw std::runtime_error("Opposing layer schedule was refused");
    auto stations = half.offsets;
    for (std::size_t i = half.offsets.size() - 1; i > 0; --i)
        stations.push_back(side - half.offsets[i - 1]);
    return stations;
}

struct WallPatch {
    paacs::mesh::Mesh surface;
    std::vector<paacs::algo::SweepSideClassification> boundary;
};

WallPatch makeWallPatch(const std::vector<double>& stations, unsigned axialIntervals) {
    using paacs::Index;
    using paacs::mesh::CellType;
    constexpr paacs::mesh::Classification lowerWall{2, 3};
    WallPatch wall;
    const auto lastZ = stations.size() - 1;
    for (double z : stations) {
        for (unsigned x = 0; x <= axialIntervals; ++x)
            wall.surface.addNode({length * x / axialIntervals, -side / 2, z - side / 2}, lowerWall);
    }
    const auto index = [axialIntervals](unsigned x, std::size_t z) {
        return static_cast<Index>(z * (axialIntervals + 1) + x);
    };
    for (std::size_t z = 0; z < lastZ; ++z) {
        for (unsigned x = 0; x < axialIntervals; ++x) {
            const auto a = index(x, z), b = index(x + 1, z);
            const auto c = index(x + 1, z + 1), d = index(x, z + 1);
            // Orient the source triangles into the fluid along +y.
            const std::array<Index, 3> first{a, c, b}, second{a, d, c};
            wall.surface.addCell(CellType::Tri3, first, lowerWall);
            wall.surface.addCell(CellType::Tri3, second, lowerWall);
        }
        wall.boundary.push_back({index(0, z), index(0, z + 1), {2, 1}});
        wall.boundary.push_back({index(axialIntervals, z), index(axialIntervals, z + 1), {2, 2}});
    }
    for (unsigned x = 0; x < axialIntervals; ++x) {
        wall.boundary.push_back({index(x, 0), index(x + 1, 0), {2, 5}});
        wall.boundary.push_back({index(x, lastZ), index(x + 1, lastZ), {2, 6}});
    }
    return wall;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: duct-opposed-wall-layers output.paacsmesh first-mm rows core-cap-mm axial-layers\n";
        return 2;
    }
    try {
        const auto firstHeight = std::stod(argv[2]), coreCap = std::stod(argv[4]);
        const auto wallRows = std::stoul(argv[3]), axialIntervals = std::stoul(argv[5]);
        if (!std::isfinite(firstHeight) || firstHeight <= 0 || !std::isfinite(coreCap) || coreCap <= 0 ||
            wallRows < 1 || wallRows > 16 || axialIntervals < 1 || axialIntervals > 64)
            throw std::runtime_error("Invalid opposing wall layer settings");
        const auto stations = symmetricStations(firstHeight, static_cast<unsigned>(wallRows), coreCap);
        auto wall = makeWallPatch(stations, static_cast<unsigned>(axialIntervals));
        paacs::algo::PreparedMixedSweepRequest request;
        request.layers = static_cast<std::uint32_t>(stations.size() - 1);
        request.direction = {0., 1., 0.};
        request.totalThickness = side;
        request.maximumAllowedGrowth = growth + 1e-12;
        request.targetClassification = {2, 4};
        request.sideEdgeClassifications = std::move(wall.boundary);
        for (const double station : stations) request.layerFractions.push_back(station / side);
        request.maximumNodes = 100'000;
        request.maximumCells = 100'000;
        request.executionThreads = 1;
        paacs::mesh::Mesh mesh;
        const auto result = paacs::algo::generatePreparedMixedSweep(wall.surface, request, mesh);
        const auto& metrics = result.metrics;
        if (!result.published() || !metrics.usedFullGeometricAudit || !metrics.audit.valid ||
            !paacs::mesh::auditLinearVolumeBoundary(mesh).valid ||
            metrics.wedges != wall.surface.cellCount() * request.layers || metrics.hexes ||
            std::abs(metrics.totalVolume - length * side * side) > 1e-9)
            throw std::runtime_error("Opposing wall-normal prism mesh failed validation");
        if (!paacs::sdk::writeMeshArchive(mesh, argv[1]).completed)
            throw std::runtime_error("Could not write the wall-normal mesh archive");
        std::cout << std::setprecision(16)
                  << "{\"prism6\":" << metrics.wedges << ",\"wall_normal_intervals\":" << request.layers
                  << ",\"prescribed_rows_per_wall\":" << wallRows
                  << ",\"first_height_mm\":" << firstHeight << ",\"maximum_growth\":" << metrics.maximumAdjacentGrowth
                  << ",\"volume_mm3\":" << metrics.totalVolume << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
