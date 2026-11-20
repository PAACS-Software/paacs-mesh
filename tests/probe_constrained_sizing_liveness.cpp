// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <cstdio>
#include <cstdint>
#include <memory>
#include <vector>

#include "paacs/algo/ProductionModelMesher.hpp"

namespace {
using namespace paacs;

struct Snapshot {
    std::vector<double> coordinates;
    std::vector<std::int64_t> topology;
    friend bool operator==(const Snapshot&, const Snapshot&) = default;
};

Snapshot snapshot(const mesh::Mesh& value) {
    Snapshot out;
    out.topology.push_back(value.nodeCount());
    out.topology.push_back(value.cellCount());
    out.topology.push_back(
        static_cast<std::int64_t>(value.dedupTolerance() * 1.0e15));
    for (Index i = 0; i < value.nodeCount(); ++i) {
        const Vec3& p = value.node(i);
        out.coordinates.insert(out.coordinates.end(), {p.x, p.y, p.z});
        out.topology.push_back(value.nodeClass(i).dim);
        out.topology.push_back(value.nodeClass(i).tag);
    }
    for (Index i = 0; i < value.cellCount(); ++i) {
        const auto& cell = value.cell(i);
        out.topology.push_back(static_cast<std::int64_t>(cell.type));
        out.topology.push_back(cell.onEntity.dim);
        out.topology.push_back(cell.onEntity.tag);
        for (int n = 0; n < mesh::nodesPerCell(cell.type); ++n)
            out.topology.push_back(value.cellNodes(i)[n]);
    }
    return out;
}

class CenterZeroField final : public sizing::SizeField {
public:
    double eval(const Vec3& p) const override {
        return p.x == 0.5 && p.y == 0.5 && p.z == 0.5 ? 0.0 : 0.3;
    }
};

void seed(mesh::Mesh& value) {
    value.setDedupTolerance(0.125);
    const Index a = value.addNode({-3, -2, -1}, {0, 91});
    const Index b = value.addNode({-2, -2, -1}, {0, 92});
    value.addCell(mesh::CellType::Line2, {a, b}, {1, 93});
}

}  // namespace

int main() {
    using namespace paacs;
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0, 0, 0}, {1, 1, 1});
    sizing::MeshSpec invalid;
    invalid.field = std::make_shared<CenterZeroField>();
    invalid.threads = 1;
    mesh::Mesh destination;
    seed(destination);
    const Snapshot before = snapshot(destination);
    const auto failed =
        algo::meshModelConstrained(*model, invalid, destination);
    if (failed.published() ||
        failed.reason != algo::ProductionModelReason::SizingEvaluationFailed ||
        failed.phase != algo::ProductionModelPhase::VolumeConstruction ||
        !failed.destinationUnchanged || !(snapshot(destination) == before)) {
        std::printf("probe_constrained_sizing_liveness FAILED: bad outcome\n");
        return 1;
    }

    sizing::MeshSpec valid;
    valid.field = std::make_shared<sizing::ConstantField>(0.3);
    valid.threads = 1;
    const auto retry = algo::meshModelConstrained(*model, valid, destination);
    if (!retry.published()) {
        std::printf("probe_constrained_sizing_liveness FAILED: retry\n");
        return 1;
    }
    std::printf("probe_constrained_sizing_liveness PASSED: center-zero "
                "contained before lattice loop; exact rollback and retry\n");
    return 0;
}
