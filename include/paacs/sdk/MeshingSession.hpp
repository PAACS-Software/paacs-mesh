// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/HybridMeshingExecutor.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/PreOcctMeshingCore.hpp"
#include "paacs/sdk/Version.hpp"

#include <atomic>
#include <cstdint>
#include <functional>
#include <string>

namespace paacs::sdk {

class CancellationToken {
public:
    void cancel() noexcept { cancelled_.store(true,std::memory_order_release); }
    bool cancelled() const noexcept { return cancelled_.load(std::memory_order_acquire); }
private: std::atomic<bool> cancelled_{false};
};
enum class SessionStage : std::uint8_t { Preflight, Planning, Generation, Audit, Publication };
using ProgressCallback=std::function<void(SessionStage,std::uint64_t,std::uint64_t)>;
struct SessionOptions { int threads=0; ProgressCallback progress; CancellationToken* cancellation=nullptr; };
struct SessionManifest { std::string apiVersion=kApiVersion;std::uint64_t planHash=0,meshHash=0;bool deterministic=true;int requestedThreads=0; };
struct SessionResult { bool completed=false,cancelled=false;std::string diagnostic;SessionManifest manifest{}; };
class MeshingSession {
public:
    explicit MeshingSession(SessionOptions options={}) : options_(std::move(options)) {}
    SessionResult executeHybrid(const algo::HybridExecutionRequest& request,
                                algo::HybridMesh& destination) const noexcept;
    SessionResult generateMappedHex(algo::MappedHexBlockRequest request,
                                    mesh::Mesh& destination) const noexcept;
    // Repeated-generation overload. The caller owns the workspace so one
    // session can drive independent remeshing streams without hidden global
    // state or allocator churn. The explicit Mesh result is unchanged.
    SessionResult generateMappedHex(algo::MappedHexBlockRequest request,
                                    algo::MappedHexWorkspace& workspace,
                                    mesh::Mesh& destination) const noexcept;
    // Complete geometry-independent Tet production route.  Geometry and
    // topology arrive through the PAACS Model/provider contracts; the compact
    // construction arena remains private to the Delaunay implementation.
    SessionResult meshPreparedTet(
        const geo::Model& model,
        const sizing::MeshSpec& spec,
        const algo::PreOcctTetRequest& request,
        mesh::Mesh& destination,
        const algo::VolumeBoundaryOwnershipProvider& ownershipProvider,
        const algo::Tet10PlacementProvider* placementProvider = nullptr,
        algo::PreOcctMeshingOutcome* detailedOutcome = nullptr) const noexcept;
private:
    bool cancelled() const noexcept { return options_.cancellation&&options_.cancellation->cancelled(); }
    void report(SessionStage stage,std::uint64_t done,std::uint64_t total) const;
    SessionOptions options_;
};

}  // namespace paacs::sdk
