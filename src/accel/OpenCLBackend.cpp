// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// OpenCL backend: vendor-agnostic GPU path (built only with
// -DPAACS_MESH_WITH_OPENCL=ON, which defines PAACS_MESH_HAS_OPENCL).
//
// STATUS: scaffold. Written against OpenCL 1.2 C API for maximum device
// coverage (NVIDIA, AMD, Intel, Mali). This build environment has no GPU
// or OpenCL runtime, so the code compiles in principle but has NOT been
// executed - validate on real hardware before relying on it. The kernels
// mirror CpuBackend semantics exactly so the unit tests can compare
// backends bit-tolerantly once a device is available.
#ifdef PAACS_MESH_HAS_OPENCL

#include "paacs/accel/Backend.hpp"

#include <CL/cl.h>
#include <cstring>
#include <stdexcept>

namespace paacs::accel {

namespace {

const char* kSource = R"CLC(
__kernel void eval_field(__global const float* grid,
                         const float3 lo, const float3 hi,
                         const int nx, const int ny, const int nz,
                         __global const float* px, __global const float* py,
                         __global const float* pz, __global float* out,
                         const int n) {
    int gid = get_global_id(0);
    if (gid >= n) return;
    float3 p = (float3)(px[gid], py[gid], pz[gid]);
    float3 t = (p - lo) / (hi - lo);
    float tx = clamp(t.x * (nx - 1), 0.0f, (float)(nx - 1));
    float ty = clamp(t.y * (ny - 1), 0.0f, (float)(ny - 1));
    float tz = clamp(t.z * (nz - 1), 0.0f, (float)(nz - 1));
    int i = min(nx - 2, (int)tx), j = min(ny - 2, (int)ty),
        k = min(nz - 2, (int)tz);
    float fx = tx - i, fy = ty - j, fz = tz - k;
    #define AT(a,b,c) grid[(a) + nx * ((b) + ny * (c))]
    float c00 = AT(i,j,k)*(1-fx)   + AT(i+1,j,k)*fx;
    float c10 = AT(i,j+1,k)*(1-fx) + AT(i+1,j+1,k)*fx;
    float c01 = AT(i,j,k+1)*(1-fx)   + AT(i+1,j,k+1)*fx;
    float c11 = AT(i,j+1,k+1)*(1-fx) + AT(i+1,j+1,k+1)*fx;
    out[gid] = mix(mix(c00, c10, fy), mix(c01, c11, fy), fz);
    #undef AT
}

__kernel void tet_quality(__global const float* cx, __global const float* cy,
                          __global const float* cz, __global const int* tets,
                          __global float* out, const int n) {
    int gid = get_global_id(0);
    if (gid >= n) return;
    int a = tets[4*gid], b = tets[4*gid+1], c = tets[4*gid+2],
        d = tets[4*gid+3];
    float3 A = (float3)(cx[a], cy[a], cz[a]);
    float3 B = (float3)(cx[b], cy[b], cz[b]);
    float3 C = (float3)(cx[c], cy[c], cz[c]);
    float3 D = (float3)(cx[d], cy[d], cz[d]);
    float3 ba = B - A, ca = C - A, da = D - A;
    float V = fabs(dot(cross(ba, ca), da)) / 6.0f;
    if (V < 1e-30f) { out[gid] = 0.0f; return; }
    float area = 0.5f * (length(cross(B-A, C-A)) + length(cross(B-A, D-A)) +
                         length(cross(C-A, D-A)) + length(cross(C-B, D-B)));
    float rin = 3.0f * V / area;
    float3 num = dot(ba,ba)*cross(ca,da) + dot(ca,ca)*cross(da,ba) +
                 dot(da,da)*cross(ba,ca);
    float R = length(num) / (12.0f * V);
    out[gid] = clamp(3.0f * rin / R, 0.0f, 1.0f);
}
)CLC";

class OpenCLBackend final : public Backend {
public:
    OpenCLBackend() {
        cl_uint nPlat = 0;
        if (clGetPlatformIDs(0, nullptr, &nPlat) != CL_SUCCESS || !nPlat)
            throw std::runtime_error("no OpenCL platform");
        std::vector<cl_platform_id> plats(nPlat);
        clGetPlatformIDs(nPlat, plats.data(), nullptr);
        for (auto p : plats)
            if (clGetDeviceIDs(p, CL_DEVICE_TYPE_GPU, 1, &dev_, nullptr) ==
                CL_SUCCESS)
                break;
        if (!dev_ &&
            clGetDeviceIDs(plats[0], CL_DEVICE_TYPE_DEFAULT, 1, &dev_,
                           nullptr) != CL_SUCCESS)
            throw std::runtime_error("no OpenCL device");

        cl_int err = 0;
        ctx_ = clCreateContext(nullptr, 1, &dev_, nullptr, nullptr, &err);
        queue_ = clCreateCommandQueue(ctx_, dev_, 0, &err);
        const char* src = kSource;
        prog_ = clCreateProgramWithSource(ctx_, 1, &src, nullptr, &err);
        if (clBuildProgram(prog_, 1, &dev_, "", nullptr, nullptr) !=
            CL_SUCCESS)
            throw std::runtime_error("OpenCL build failed");
        kField_ = clCreateKernel(prog_, "eval_field", &err);
        kQual_  = clCreateKernel(prog_, "tet_quality", &err);
    }
    ~OpenCLBackend() override {
        clReleaseKernel(kField_);
        clReleaseKernel(kQual_);
        clReleaseProgram(prog_);
        clReleaseCommandQueue(queue_);
        clReleaseContext(ctx_);
    }

    std::string name() const override {
        char buf[256] = {};
        clGetDeviceInfo(dev_, CL_DEVICE_NAME, sizeof buf, buf, nullptr);
        return std::string("opencl:") + buf;
    }

    void evalField(const sizing::BakedGridField& f,
                   const std::vector<Vec3>& pts,
                   std::vector<float>& out) override {
        const int n = static_cast<int>(pts.size());
        out.resize(pts.size());
        std::vector<float> px(pts.size()), py(pts.size()), pz(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            px[i] = float(pts[i].x); py[i] = float(pts[i].y);
            pz[i] = float(pts[i].z);
        }
        auto buf = [&](const void* d, size_t bytes) {
            cl_int err = 0;
            return clCreateBuffer(ctx_,
                                  CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  bytes, const_cast<void*>(d), &err);
        };
        cl_mem bg = buf(f.data().data(), f.data().size() * 4);
        cl_mem bx = buf(px.data(), px.size() * 4);
        cl_mem by = buf(py.data(), py.size() * 4);
        cl_mem bz = buf(pz.data(), pz.size() * 4);
        cl_int err = 0;
        cl_mem bo = clCreateBuffer(ctx_, CL_MEM_WRITE_ONLY, out.size() * 4,
                                   nullptr, &err);
        cl_float3 lo{{float(f.gridLo().x), float(f.gridLo().y),
                      float(f.gridLo().z)}};
        cl_float3 hi{{float(f.gridHi().x), float(f.gridHi().y),
                      float(f.gridHi().z)}};
        int nx = f.nx(), ny = f.ny(), nz = f.nz();
        int arg = 0;
        clSetKernelArg(kField_, arg++, sizeof(cl_mem), &bg);
        clSetKernelArg(kField_, arg++, sizeof(cl_float3), &lo);
        clSetKernelArg(kField_, arg++, sizeof(cl_float3), &hi);
        clSetKernelArg(kField_, arg++, sizeof(int), &nx);
        clSetKernelArg(kField_, arg++, sizeof(int), &ny);
        clSetKernelArg(kField_, arg++, sizeof(int), &nz);
        clSetKernelArg(kField_, arg++, sizeof(cl_mem), &bx);
        clSetKernelArg(kField_, arg++, sizeof(cl_mem), &by);
        clSetKernelArg(kField_, arg++, sizeof(cl_mem), &bz);
        clSetKernelArg(kField_, arg++, sizeof(cl_mem), &bo);
        clSetKernelArg(kField_, arg++, sizeof(int), &n);
        size_t gsz = (pts.size() + 63) / 64 * 64;
        clEnqueueNDRangeKernel(queue_, kField_, 1, nullptr, &gsz, nullptr, 0,
                               nullptr, nullptr);
        clEnqueueReadBuffer(queue_, bo, CL_TRUE, 0, out.size() * 4,
                            out.data(), 0, nullptr, nullptr);
        for (cl_mem b : {bg, bx, by, bz, bo}) clReleaseMemObject(b);
    }

    void tetQuality(const std::vector<Vec3>& coords,
                    const std::vector<Index>& tets,
                    std::vector<float>& out) override {
        // Same buffer pattern as evalField with the tet_quality kernel;
        // omitted from the scaffold until it can be exercised on hardware.
        // Falls back to CPU semantics via the default backend chooser.
        (void)coords; (void)tets; (void)out;
        throw std::runtime_error("OpenCL tetQuality: validate on hardware");
    }

private:
    cl_device_id dev_ = nullptr;
    cl_context ctx_ = nullptr;
    cl_command_queue queue_ = nullptr;
    cl_program prog_ = nullptr;
    cl_kernel kField_ = nullptr, kQual_ = nullptr;
};

} // namespace

std::unique_ptr<Backend> makeOpenCLBackend() {
    try {
        return std::make_unique<OpenCLBackend>();
    } catch (...) {
        return nullptr;   // no platform/device: caller falls back to CPU
    }
}

} // namespace paacs::accel

#endif // PAACS_MESH_HAS_OPENCL
