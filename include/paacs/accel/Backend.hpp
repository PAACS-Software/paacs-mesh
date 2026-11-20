// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// Accelerator abstraction for batch-parallel meshing kernels.
//
// Strategy (why this shape): Delaunay insertion is control-flow heavy and
// hard to accelerate; the profitable batch stages are (a) size-field
// evaluation over millions of candidate points, (b) inside/outside and
// clearance filtering, (c) per-element quality metrics. Those are pure
// data-parallel maps over flat arrays, so one narrow interface covers CPU
// threads today and GPU backends tomorrow. Fields cross the boundary as
// BakedGridField (flat float grid + trilinear lookup) because polymorphic
// C++ field objects cannot execute on a GPU.
//
// Backends:
//   CpuBackend    - std::thread pool; always available; unit-tested.
//   OpenCLBackend - vendor-agnostic GPU path (NVIDIA/AMD/Intel/Apple via
//                   OpenCL 1.2). Built with -DPAACS_MESH_WITH_OPENCL=ON.
//                   NOTE: scaffolded and compile-guarded; this environment
//                   has no GPU, so it is NOT yet validated - treat as a
//                   starting point, not shipped functionality.
//
#include <memory>
#include <string>
#include <vector>

#include "paacs/core/Types.hpp"
#include "paacs/sizing/SizeField.hpp"

namespace paacs::accel {

class Backend {
public:
    virtual ~Backend() = default;
    virtual std::string name() const = 0;

    // out[i] = trilinear lookup of the baked field at pts[i].
    virtual void evalField(const sizing::BakedGridField& field,
                           const std::vector<Vec3>& pts,
                           std::vector<float>& out) = 0;

    // out[i] = per-tet radius-ratio quality; tets as flat index quadruples
    // into coords.
    virtual void tetQuality(const std::vector<Vec3>& coords,
                            const std::vector<Index>& tets,
                            std::vector<float>& out) = 0;
};

// Always available; uses core::parallelFor.
std::unique_ptr<Backend> makeCpuBackend();

// Returns nullptr when compiled without OpenCL or when no device is found.
std::unique_ptr<Backend> makeOpenCLBackend();

// Best available backend: OpenCL device if compiled in and present,
// otherwise CPU.
std::unique_ptr<Backend> makeDefaultBackend();

} // namespace paacs::accel
