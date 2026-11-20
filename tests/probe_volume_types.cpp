// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/StructuredPolyhedralMesher.hpp"
#include "paacs/mesh/PolyhedralMesh.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

[[noreturn]] void fail(const char* message) {
    std::fprintf(stderr, "probe_volume_types: %s\n", message);
    std::exit(1);
}

void require(bool value, const char* message) {
    if (!value) fail(message);
}

paacs::algo::StructuredBlockGeometry geometry(std::uint32_t n) {
    paacs::algo::StructuredBlockGeometry out;
    out.origin = {-2.0, 3.0, 5.0};
    out.axisU = {2.0, 0.0, 0.0};
    out.axisV = {0.2, 3.0, 0.0};
    out.axisW = {0.1, 0.15, 4.0};
    out.cellsU = n;
    out.cellsV = n;
    out.cellsW = n;
    out.volumeClassification = {3, 77};
    for (std::size_t i = 0; i < out.boundaryClassifications.size(); ++i)
        out.boundaryClassifications[i] =
            {2, static_cast<int>(301 + i)};
    return out;
}

}  // namespace

int main() {
    paacs::algo::MappedHexBlockRequest hexRequest;
    hexRequest.geometry = geometry(20);
    hexRequest.validation =
        paacs::algo::MappedHexValidationMode::FullTopologyAudit;
    paacs::mesh::Mesh hex;
    const auto hexOut = paacs::algo::generateMappedHexBlock(
        hexRequest, hex);
    require(hexOut.published(), "full-audit Hex8 generation failed");
    require(hexOut.metrics.audit.valid, "Hex8 audit failed");
    require(hexOut.metrics.hexes == 8000, "wrong Hex8 count");
    require(hexOut.metrics.boundaryQuads == 2400,
            "wrong Quad4 boundary count");
    require(std::abs(hexOut.metrics.audit.totalVolume - 24.0) < 1.0e-9,
            "wrong Hex8 volume");

    paacs::mesh::PolyhedralMesh converted;
    const auto convertedOut = paacs::algo::convertHex8ToPolyhedral(
        hex, converted);
    require(convertedOut.published(), "Hex8 conversion failed");
    require(convertedOut.metrics.audit.valid,
            "converted polyhedral audit failed");
    require(converted.cellCount() == 8000,
            "wrong converted polyhedral count");

    paacs::mesh::PolyhedralMesh direct;
    const auto directOut = paacs::algo::generateStructuredPolyhedralBlock(
        geometry(16), direct,
        paacs::algo::StructuredPolyValidationMode::FullAudit);
    require(directOut.published(), "full-audit poly generation failed");
    require(directOut.metrics.audit.valid, "polyhedral audit failed");
    require(direct.cellCount() == 4096, "wrong polyhedral count");
    require(std::abs(directOut.metrics.audit.totalVolume - 24.0) < 1.0e-9,
            "wrong polyhedral volume");

    paacs::mesh::Mesh sentinel;
    (void)sentinel.addNode({9.0, 8.0, 7.0}, {0, 99});
    hexRequest.geometry.axisU = {0.0, 1.0, 0.0};
    hexRequest.geometry.axisV = {1.0, 0.0, 0.0};
    const auto rejected = paacs::algo::generateMappedHexBlock(
        hexRequest, sentinel);
    require(!rejected.published(), "negative frame was accepted");
    require(rejected.destinationUnchanged,
            "rejected Hex8 request mutated destination");
    require(sentinel.nodeCount() == 1 && sentinel.cellCount() == 0,
            "sentinel destination changed");

    std::puts("probe_volume_types: PASS");
    return 0;
}
