// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/accel/Backend.hpp"
#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/Quality.hpp"

namespace paacs::accel {

namespace {

class CpuBackend final : public Backend {
public:
    std::string name() const override { return "cpu-threads"; }

    void evalField(const sizing::BakedGridField& field,
                   const std::vector<Vec3>& pts,
                   std::vector<float>& out) override {
        out.resize(pts.size());
        core::parallelFor(0, pts.size(), [&](size_t i) {
            out[i] = static_cast<float>(field.eval(pts[i]));
        });
    }

    void tetQuality(const std::vector<Vec3>& coords,
                    const std::vector<Index>& tets,
                    std::vector<float>& out) override {
        const size_t n = tets.size() / 4;
        out.resize(n);
        core::parallelFor(0, n, [&](size_t i) {
            out[i] = static_cast<float>(mesh::tetQuality(
                coords[static_cast<size_t>(tets[4 * i + 0])],
                coords[static_cast<size_t>(tets[4 * i + 1])],
                coords[static_cast<size_t>(tets[4 * i + 2])],
                coords[static_cast<size_t>(tets[4 * i + 3])]));
        });
    }
};

} // namespace

std::unique_ptr<Backend> makeCpuBackend() {
    return std::make_unique<CpuBackend>();
}

#ifndef PAACS_MESH_HAS_OPENCL
std::unique_ptr<Backend> makeOpenCLBackend() { return nullptr; }
#endif

std::unique_ptr<Backend> makeDefaultBackend() {
    if (auto cl = makeOpenCLBackend()) return cl;
    return makeCpuBackend();
}

} // namespace paacs::accel
