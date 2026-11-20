// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "paacs/algo/AdaptiveBoxMesher.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/algo/LayeredPatchMesher.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include "paacs/algo/PrimitiveMesher.hpp"
#include <array>
#include <chrono>
#include <functional>
#include <memory>
#include <variant>

namespace paacs::sdk {
// Poll snapshots on the UI thread; no callbacks are invoked on a UI object.
enum class JobState : std::uint8_t {
    Queued,
    WaitingForEngine,
    Running,
    Succeeded,
    Failed,
    Cancelled
};
enum class JobStage : std::uint8_t { Queued, Preparing, Generating, Finalizing, Complete };
enum class JobDomain : std::uint8_t {
    Service,
    Primitive,
    AdaptiveBox,
    Geometry,
    LayeredPatch,
    MappedHex,
    Custom
};
enum class JobError : std::uint32_t {
    None,
    EngineRejected,
    Cancelled,
    ResultBudget,
    AllocationFailed,
    Exception
};
enum class SubmitReason : std::uint8_t { None, EmptyTask, JobLimit, Stopped, AllocationFailed };
bool terminal(JobState) noexcept;
const char *jobStateName(JobState) noexcept;
const char *jobStageName(JobStage) noexcept;

struct JobSnapshot {
    std::uint64_t id = 0;
    JobState state = JobState::Queued;
    JobStage stage = JobStage::Queued;
    bool cancellationRequested = false;
    // total == 0 means indeterminate. Never a fabricated percentage.
    std::uint64_t completed = 0, total = 0;
    double elapsedSeconds = 0;
    JobError error = JobError::None;
    JobDomain domain = JobDomain::Service;
    std::uint32_t engineReason = 0;
    std::array<char, 192> diagnostic{};
};
using EngineReport =
    std::variant<std::monostate, algo::PrimitiveMeshingResult, algo::AdaptiveBoxResult,
                 algo::GeometryMeshingOutcome, algo::LayeredPatchResult, algo::MappedHexOutcome>;
struct TaskOutcome {
    bool published = false;
    JobDomain domain = JobDomain::Custom;
    std::uint32_t reason = 0;
    std::array<char, 192> diagnostic{};
    EngineReport report;
    static TaskOutcome status(bool published, JobDomain, std::uint32_t reason,
                              const char *diagnostic) noexcept;
};
struct MeshSummary {
    std::uint64_t nodes = 0, cells = 0, connectivity = 0, logicalBytes = 0, hash = 0;
    bool hashAvailable = false;
    std::array<std::uint64_t, 10> cellTypes{};
    Vec3 lower{}, upper{};
};
namespace detail {
struct JobRecord;
struct ResultBudget;
} // namespace detail
class JobResult {
  public:
    mesh::Mesh mesh;
    MeshSummary summary;
    EngineReport report;
    ~JobResult();
    JobResult() = default;
    JobResult(const JobResult &) = delete;
    JobResult &operator=(const JobResult &) = delete;

  private:
    std::shared_ptr<detail::ResultBudget> budget_;
    std::uint64_t chargedBytes_ = 0;
    friend class JobService;
};
class JobContext {
  public:
    bool cancelled() const noexcept;
    const std::atomic<bool> &cancellation() const noexcept;
    std::uint64_t maximumResultBytes() const noexcept { return maximumBytes_; }
    // Custom native jobs can report their measured work. No callback/re-entry.
    void progress(JobStage, std::uint64_t completed = 0, std::uint64_t total = 0) noexcept;

  private:
    JobContext(detail::JobRecord &, std::uint64_t) noexcept;
    detail::JobRecord &record_;
    std::uint64_t maximumBytes_;
    friend class JobService;
};
// A native extension can invoke any existing PAACS operation here. Capture
// input ownership, use the cancellation flag, honor the output budget and only
// return published=true for a validated result. Never wait on another job or
// destroy/shut down its own service from this function.
using MeshTask = std::function<TaskOutcome(JobContext &, mesh::Mesh &)>;
class JobHandle {
  public:
    JobHandle() = default;
    explicit operator bool() const noexcept { return bool(record_); }
    JobSnapshot snapshot() const;
    // Linearized against publication: true means no result will be published.
    // A late cancellation of a terminal job returns false and changes nothing.
    bool cancel() const;
    bool waitFor(std::chrono::milliseconds timeout) const;
    std::shared_ptr<const JobResult> result() const;
    // Includes structured refusal details when result() is null.
    EngineReport report() const;

  private:
    std::shared_ptr<detail::JobRecord> record_;
    explicit JobHandle(std::shared_ptr<detail::JobRecord> record) : record_(std::move(record)) {}
    friend class JobService;
};
struct JobServiceLimits {
    // Includes terminal handles retained by the application. Release handles
    // when done; result references independently retain their byte charge.
    std::size_t maximumJobs = 8;
    std::uint64_t maximumResultBytes = 512ULL * 1024 * 1024;
    std::uint64_t maximumRetainedResultBytes = 512ULL * 1024 * 1024;
    // Strict hashing scans every coordinate/classification/connectivity byte.
    // Opt in for reproducibility; a viewport does not need this extra pass.
    bool computeResultHash = false;
};
struct JobSubmission {
    JobHandle job;
    SubmitReason reason = SubmitReason::None;
};
class JobService {
  public:
    explicit JobService(JobServiceLimits = {});
    ~JobService();
    JobService(const JobService &) = delete;
    JobService &operator=(const JobService &) = delete;
    JobSubmission submit(MeshTask) noexcept;
    // Cooperatively cancel and join. Call from the application, never a task.
    // Mapped Hex8 has no internal cancellation; its result is discarded on cancel.
    void shutdown() noexcept;
    std::uint64_t retainedResultBytes() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    void run() noexcept;
};

// Owning adapters preserve every scalar control and the complete engine report.
// Borrowed profile/position/audit spans are copied before these factories return.
// All caller cancellation pointers are replaced by the job's owned flag.
MeshTask primitiveTask(algo::PrimitiveShape, const algo::PrimitiveMeshingOptions & = {});
MeshTask adaptiveBoxTask(algo::AdaptiveBoxRequest,
                         std::shared_ptr<const sizing::BoundedSizeField> field = {});
MeshTask layeredPatchTask(std::shared_ptr<const mesh::Mesh> source,
                          const algo::LayeredPatchRequest & = {});
MeshTask mappedHexTask(algo::MappedHexBlockRequest);
struct GeometryJobInput {
    std::shared_ptr<const geo::Model> model;
    algo::GeometryMeshingOptions options;
    std::shared_ptr<const sizing::BoundedSizeField> field;
    // Empty, or one entry per options.bodyControls; null inherits global field.
    std::vector<std::shared_ptr<const sizing::BoundedSizeField>> bodyFields;
};
// Shared models/fields/source patches must remain immutable, including aliases.
// Raw field pointers must be null or match the supplied owner; mismatch throws.
MeshTask geometryTask(GeometryJobInput);
} // namespace paacs::sdk
