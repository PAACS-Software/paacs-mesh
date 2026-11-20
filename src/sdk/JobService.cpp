// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/JobService.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include <condition_variable>
#include <deque>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>

namespace paacs::sdk {
namespace detail {
using Clock = std::chrono::steady_clock;
struct ResultBudget {
    mutable std::mutex mutex;
    std::uint64_t used = 0, maximum = 0;
};
struct JobRecord {
    mutable std::mutex mutex;
    std::condition_variable changed;
    std::atomic<bool> cancelled{false};
    JobSnapshot snapshot;
    Clock::time_point submitted = Clock::now(), finished{};
    MeshTask task;
    std::shared_ptr<const JobResult> result;
    EngineReport failureReport;
};
} // namespace detail
namespace {
std::timed_mutex &engineGate() {
    static std::timed_mutex gate;
    return gate;
}
void message(std::array<char, 192> &out, const char *text) noexcept {
    out.fill(0);
    if (text)
        for (std::size_t i = 0; i + 1 < out.size() && text[i]; ++i)
            out[i] = text[i];
}
void finish(detail::JobRecord &r, JobState state, JobError error, const char *text) noexcept {
    r.snapshot.state = state;
    r.snapshot.stage = JobStage::Complete;
    r.snapshot.error = error;
    r.snapshot.completed = state == JobState::Succeeded ? 1 : 0;
    r.snapshot.total = 1;
    message(r.snapshot.diagnostic, text);
    r.finished = detail::Clock::now();
    r.changed.notify_all();
}
void summarize(const mesh::Mesh &m, MeshSummary &s, bool computeHash) {
    s.nodes = static_cast<std::uint64_t>(m.nodeCount());
    s.cells = static_cast<std::uint64_t>(m.cellCount());
    s.connectivity = static_cast<std::uint64_t>(m.connectivityCount());
    s.logicalBytes = 0;
    auto add = [&](std::uint64_t count, std::uint64_t width) {
        if (count > (UINT64_MAX - s.logicalBytes) / width)
            throw std::length_error("result size overflow");
        s.logicalBytes += count * width;
    };
    add(s.nodes, sizeof(Vec3) + sizeof(mesh::Classification));
    add(s.cells, sizeof(mesh::Cell));
    add(s.connectivity, sizeof(Index));
    if (m.nodeCount())
        s.lower = s.upper = m.node(0);
    for (const auto &p : m.nodes()) {
        s.lower = {std::min(s.lower.x, p.x), std::min(s.lower.y, p.y), std::min(s.lower.z, p.z)};
        s.upper = {std::max(s.upper.x, p.x), std::max(s.upper.y, p.y), std::max(s.upper.z, p.z)};
    }
    for (Index i = 0; i < m.cellCount(); ++i) {
        auto type = static_cast<std::size_t>(m.cell(i).type);
        if (type >= s.cellTypes.size())
            throw std::invalid_argument("invalid result cell type");
        ++s.cellTypes[type];
    }
    if (computeHash) {
        s.hash = mesh::meshHash(m);
        s.hashAvailable = true;
    }
}
template <class T> std::optional<T> copyOptional(const T *p) {
    return p ? std::optional<T>(*p) : std::nullopt;
}
void checkField(const sizing::BoundedSizeField *raw,
                const std::shared_ptr<const sizing::BoundedSizeField> &owner) {
    if (raw && raw != owner.get())
        throw std::invalid_argument("field needs matching shared ownership");
}
} // namespace
bool terminal(JobState s) noexcept {
    return s == JobState::Succeeded || s == JobState::Failed || s == JobState::Cancelled;
}
const char *jobStateName(JobState s) noexcept {
    switch (s) {
    case JobState::Queued:
        return "queued";
    case JobState::WaitingForEngine:
        return "waiting_for_engine";
    case JobState::Running:
        return "running";
    case JobState::Succeeded:
        return "succeeded";
    case JobState::Failed:
        return "failed";
    case JobState::Cancelled:
        return "cancelled";
    }
    return "unknown";
}
const char *jobStageName(JobStage s) noexcept {
    switch (s) {
    case JobStage::Queued:
        return "queued";
    case JobStage::Preparing:
        return "preparing";
    case JobStage::Generating:
        return "generating";
    case JobStage::Finalizing:
        return "finalizing";
    case JobStage::Complete:
        return "complete";
    }
    return "unknown";
}
TaskOutcome TaskOutcome::status(bool published, JobDomain domain, std::uint32_t reason,
                                const char *text) noexcept {
    TaskOutcome out;
    out.published = published;
    out.domain = domain;
    out.reason = reason;
    message(out.diagnostic, text);
    return out;
}
JobResult::~JobResult() {
    if (budget_) {
        std::lock_guard lock(budget_->mutex);
        budget_->used -= chargedBytes_;
    }
}
JobContext::JobContext(detail::JobRecord &r, std::uint64_t maximum) noexcept
    : record_(r), maximumBytes_(maximum) {}
bool JobContext::cancelled() const noexcept {
    return record_.cancelled.load(std::memory_order_acquire);
}
const std::atomic<bool> &JobContext::cancellation() const noexcept {
    return record_.cancelled;
}
void JobContext::progress(JobStage stage, std::uint64_t done, std::uint64_t total) noexcept {
    std::lock_guard lock(record_.mutex);
    if (record_.snapshot.state != JobState::Running || stage == JobStage::Complete ||
        stage == JobStage::Queued)
        return;
    record_.snapshot.stage = stage;
    record_.snapshot.total = total;
    record_.snapshot.completed = total ? std::min(done, total) : done;
}
JobSnapshot JobHandle::snapshot() const {
    if (!record_)
        throw std::logic_error("empty job handle");
    std::lock_guard lock(record_->mutex);
    auto s = record_->snapshot;
    s.cancellationRequested = record_->cancelled.load(std::memory_order_acquire);
    const auto end = terminal(s.state) ? record_->finished : detail::Clock::now();
    s.elapsedSeconds = std::chrono::duration<double>(end - record_->submitted).count();
    return s;
}
bool JobHandle::cancel() const {
    if (!record_)
        return false;
    MeshTask discarded;
    {
        std::lock_guard lock(record_->mutex);
        if (terminal(record_->snapshot.state))
            return false;
        record_->cancelled.store(true, std::memory_order_release);
        if (record_->snapshot.state != JobState::Running) {
            discarded = std::move(record_->task);
            finish(*record_, JobState::Cancelled, JobError::Cancelled,
                   "cancelled before execution");
        }
    }
    return true;
}
bool JobHandle::waitFor(std::chrono::milliseconds timeout) const {
    if (!record_)
        return false;
    std::unique_lock lock(record_->mutex);
    return record_->changed.wait_for(lock, timeout,
                                     [&] { return terminal(record_->snapshot.state); });
}
std::shared_ptr<const JobResult> JobHandle::result() const {
    if (!record_)
        return {};
    std::lock_guard lock(record_->mutex);
    return record_->result;
}
EngineReport JobHandle::report() const {
    if (!record_)
        return {};
    std::lock_guard lock(record_->mutex);
    return record_->result ? record_->result->report : record_->failureReport;
}
struct JobService::Impl {
    JobServiceLimits limits;
    std::shared_ptr<detail::ResultBudget> budget = std::make_shared<detail::ResultBudget>();
    std::mutex mutex, shutdownMutex;
    std::condition_variable changed;
    std::deque<std::shared_ptr<detail::JobRecord>> queue;
    std::vector<std::weak_ptr<detail::JobRecord>> records;
    bool stopped = false;
    std::uint64_t nextId = 1;
    std::thread worker;
    explicit Impl(JobServiceLimits l) : limits(l) {
        budget->maximum = l.maximumRetainedResultBytes;
    }
};
JobService::JobService(JobServiceLimits limits) : impl_(std::make_unique<Impl>(limits)) {
    if (!limits.maximumJobs || !limits.maximumResultBytes || !limits.maximumRetainedResultBytes)
        throw std::invalid_argument("job limits must be positive");
    (void)engineGate(); // Initialized before the worker; safe static destruction order.
    impl_->worker = std::thread([this] { run(); });
}
JobService::~JobService() {
    shutdown();
}
JobSubmission JobService::submit(MeshTask task) noexcept {
    if (!task)
        return {{}, SubmitReason::EmptyTask};
    try {
        std::lock_guard lock(impl_->mutex);
        if (impl_->stopped)
            return {{}, SubmitReason::Stopped};
        std::erase_if(impl_->records, [](const auto &p) { return p.expired(); });
        // Cancelled queued entries must not retain input storage/admission slots.
        std::erase_if(impl_->queue, [](const auto &p) {
            std::lock_guard guard(p->mutex);
            return terminal(p->snapshot.state);
        });
        std::erase_if(impl_->records, [](const auto &p) { return p.expired(); });
        if (impl_->records.size() >= impl_->limits.maximumJobs || impl_->nextId == UINT64_MAX)
            return {{}, SubmitReason::JobLimit};
        auto record = std::make_shared<detail::JobRecord>();
        record->task = std::move(task);
        record->snapshot.id = impl_->nextId;
        impl_->records.emplace_back(record);
        try {
            impl_->queue.push_back(record);
        } catch (...) {
            impl_->records.pop_back();
            throw;
        }
        ++impl_->nextId;
        impl_->changed.notify_one();
        return {JobHandle(std::move(record)), SubmitReason::None};
    } catch (...) {
        return {{}, SubmitReason::AllocationFailed};
    }
}
void JobService::shutdown() noexcept {
    std::lock_guard shutdownLock(impl_->shutdownMutex);
    std::deque<std::shared_ptr<detail::JobRecord>> discarded;
    {
        std::lock_guard lock(impl_->mutex);
        impl_->stopped = true;
        for (auto &weak : impl_->records)
            if (auto r = weak.lock()) {
                std::lock_guard guard(r->mutex);
                if (!terminal(r->snapshot.state)) {
                    r->cancelled.store(true, std::memory_order_release);
                    if (r->snapshot.state != JobState::Running)
                        finish(*r, JobState::Cancelled, JobError::Cancelled, "service shut down");
                }
            }
        discarded.swap(impl_->queue);
        impl_->changed.notify_all();
    }
    // Release captured input storage outside service locks, even if handles survive.
    for (auto &r : discarded) {
        MeshTask task;
        {
            std::lock_guard lock(r->mutex);
            task = std::move(r->task);
        }
    }
    if (impl_->worker.joinable())
        impl_->worker.join();
}
std::uint64_t JobService::retainedResultBytes() const noexcept {
    std::lock_guard lock(impl_->budget->mutex);
    return impl_->budget->used;
}
void JobService::run() noexcept {
    for (;;) {
        std::shared_ptr<detail::JobRecord> r;
        {
            std::unique_lock lock(impl_->mutex);
            impl_->changed.wait(lock, [&] { return impl_->stopped || !impl_->queue.empty(); });
            if (impl_->stopped)
                return;
            r = std::move(impl_->queue.front());
            impl_->queue.pop_front();
        }
        {
            std::lock_guard lock(r->mutex);
            if (terminal(r->snapshot.state))
                continue;
            r->snapshot.state = JobState::WaitingForEngine;
        }
        std::unique_lock gate(engineGate(), std::defer_lock);
        while (!r->cancelled.load(std::memory_order_acquire) &&
               !gate.try_lock_for(std::chrono::milliseconds(10))) {
        }
        MeshTask task;
        {
            std::lock_guard lock(r->mutex);
            task = std::move(r->task);
            if (terminal(r->snapshot.state))
                continue;
            if (r->cancelled.load(std::memory_order_acquire)) {
                finish(*r, JobState::Cancelled, JobError::Cancelled, "cancelled before execution");
                continue;
            }
            r->snapshot.state = JobState::Running;
            r->snapshot.stage = JobStage::Preparing;
        }
        try {
            std::uint64_t maximum;
            {
                std::lock_guard lock(impl_->budget->mutex);
                maximum = std::min(impl_->limits.maximumResultBytes,
                                   impl_->budget->maximum - impl_->budget->used);
            }
            if (!maximum) {
                std::lock_guard lock(r->mutex);
                if (r->cancelled.load(std::memory_order_acquire))
                    finish(*r, JobState::Cancelled, JobError::Cancelled,
                           "cancelled before execution");
                else
                    finish(*r, JobState::Failed, JobError::ResultBudget,
                           "retained result budget exhausted");
                continue;
            }
            auto result = std::make_shared<JobResult>();
            JobContext context(*r, maximum);
            auto outcome = task(context, result->mesh);
            if (outcome.published && !context.cancelled()) {
                context.progress(JobStage::Finalizing);
                summarize(result->mesh, result->summary, impl_->limits.computeResultHash);
            }
            result->report = std::move(outcome.report);
            std::lock_guard lock(r->mutex);
            r->snapshot.domain = outcome.domain;
            r->snapshot.engineReason = outcome.reason;
            if (context.cancelled()) {
                finish(*r, JobState::Cancelled, JobError::Cancelled,
                       "cancelled; no result published");
            } else if (!outcome.published) {
                r->failureReport = std::move(result->report);
                finish(*r, JobState::Failed, JobError::EngineRejected, outcome.diagnostic.data());
            } else {
                std::lock_guard budgetLock(impl_->budget->mutex);
                const auto bytes = result->summary.logicalBytes;
                if (bytes > maximum || bytes > impl_->budget->maximum - impl_->budget->used) {
                    finish(*r, JobState::Failed, JobError::ResultBudget,
                           "result byte budget exceeded");
                } else {
                    result->budget_ = impl_->budget;
                    result->chargedBytes_ = bytes;
                    impl_->budget->used += bytes;
                    r->result = std::move(result);
                    finish(*r, JobState::Succeeded, JobError::None, "complete");
                }
            }
        } catch (const std::bad_alloc &) {
            std::lock_guard lock(r->mutex);
            finish(*r, r->cancelled ? JobState::Cancelled : JobState::Failed,
                   r->cancelled ? JobError::Cancelled : JobError::AllocationFailed,
                   "allocation failed");
        } catch (...) {
            std::lock_guard lock(r->mutex);
            finish(*r, r->cancelled ? JobState::Cancelled : JobState::Failed,
                   r->cancelled ? JobError::Cancelled : JobError::Exception,
                   "native task threw an exception");
        }
    }
}

MeshTask primitiveTask(algo::PrimitiveShape shape, const algo::PrimitiveMeshingOptions &input) {
    auto profile = copyOptional(input.sizeProfile);
    std::vector<double> positions(input.axialPositions.begin(), input.axialPositions.end());
    std::vector<double> volumes(input.auditLimits.cellVolumes.begin(),
                                input.auditLimits.cellVolumes.end());
    return [shape, input, profile = std::move(profile), positions = std::move(positions),
            volumes = std::move(volumes)](JobContext &c, mesh::Mesh &m) {
        auto o = input;
        o.sizeProfile = profile ? &*profile : nullptr;
        o.axialPositions = positions;
        o.auditLimits.cellVolumes = volumes;
        o.cancellation = o.auditLimits.cancellation = o.profileLimits.cancellation =
            &c.cancellation();
        o.maximumOutputBytes = std::min(o.maximumOutputBytes, c.maximumResultBytes());
        c.progress(JobStage::Generating);
        auto r = algo::meshPrimitive(shape, o, m);
        auto out = TaskOutcome::status(r.published, JobDomain::Primitive,
                                       static_cast<std::uint32_t>(r.reason),
                                       algo::primitiveReasonName(r.reason));
        out.report = std::move(r);
        return out;
    };
}
MeshTask adaptiveBoxTask(algo::AdaptiveBoxRequest input,
                         std::shared_ptr<const sizing::BoundedSizeField> field) {
    checkField(input.field, field);
    std::vector<double> volumes(input.auditLimits.cellVolumes.begin(),
                                input.auditLimits.cellVolumes.end());
    return [input, field = std::move(field), volumes = std::move(volumes)](JobContext &c,
                                                                           mesh::Mesh &m) {
        auto o = input;
        o.field = field.get();
        o.auditLimits.cellVolumes = volumes;
        o.cancellation = o.auditLimits.cancellation = &c.cancellation();
        o.maximumOutputBytes = std::min(o.maximumOutputBytes, c.maximumResultBytes());
        c.progress(JobStage::Generating);
        auto r = algo::meshAdaptiveBox(o, m);
        auto out = TaskOutcome::status(r.published, JobDomain::AdaptiveBox,
                                       static_cast<std::uint32_t>(r.reason),
                                       algo::adaptiveBoxReasonName(r.reason));
        out.report = std::move(r);
        return out;
    };
}
MeshTask layeredPatchTask(std::shared_ptr<const mesh::Mesh> source,
                          const algo::LayeredPatchRequest &input) {
    if (!source)
        throw std::invalid_argument("source patch required");
    auto profile = copyOptional(input.sizeProfile);
    std::vector<double> volumes(input.auditLimits.cellVolumes.begin(),
                                input.auditLimits.cellVolumes.end());
    return [source = std::move(source), input, profile = std::move(profile),
            volumes = std::move(volumes)](JobContext &c, mesh::Mesh &m) {
        auto o = input;
        o.sizeProfile = profile ? &*profile : nullptr;
        o.auditLimits.cellVolumes = volumes;
        o.cancellation = o.auditLimits.cancellation = o.profileLimits.cancellation =
            &c.cancellation();
        o.maximumOutputBytes = std::min(o.maximumOutputBytes, c.maximumResultBytes());
        c.progress(JobStage::Generating);
        auto r = algo::meshLayeredPatch(*source, o, m);
        auto out = TaskOutcome::status(r.published, JobDomain::LayeredPatch,
                                       static_cast<std::uint32_t>(r.reason),
                                       "layered patch result; see report.reason");
        out.report = std::move(r);
        return out;
    };
}
MeshTask geometryTask(GeometryJobInput input) {
    if (!input.model)
        throw std::invalid_argument("model required");
    checkField(input.options.spatialField, input.field);
    const auto n = input.options.bodyControls.size();
    if (input.bodyFields.empty())
        input.bodyFields.resize(n);
    if (input.bodyFields.size() != n)
        throw std::invalid_argument("one field owner per body control required");
    std::vector<std::optional<sizing::CappedSizeProfile>> profiles;
    profiles.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        checkField(input.options.bodyControls[i].spatialField, input.bodyFields[i]);
        profiles.push_back(copyOptional(input.options.bodyControls[i].sizeProfile));
    }
    return
        [input = std::move(input), profiles = std::move(profiles)](JobContext &c, mesh::Mesh &m) {
            auto o = input.options;
            o.cancellation = &c.cancellation();
            o.spatialField = input.field.get();
            o.maximumOutputBytes = std::min(o.maximumOutputBytes, c.maximumResultBytes());
            for (std::size_t i = 0; i < o.bodyControls.size(); ++i) {
                o.bodyControls[i].sizeProfile = profiles[i] ? &*profiles[i] : nullptr;
                o.bodyControls[i].spatialField = input.bodyFields[i].get();
                o.bodyControls[i].profileLimits.cancellation = &c.cancellation();
            }
            c.progress(JobStage::Generating);
            auto r = algo::meshGeometryModel(*input.model, o, m);
            auto out = TaskOutcome::status(r.published, JobDomain::Geometry,
                                           static_cast<std::uint32_t>(r.reason),
                                           algo::geometryReasonName(r.reason));
            out.report = std::move(r);
            return out;
        };
}
MeshTask mappedHexTask(algo::MappedHexBlockRequest input) {
    std::array<std::vector<double>, 3> fractions;
    const std::array<std::span<const double>, 3> views{input.growth.u.boundaryFractions,
                                                       input.growth.v.boundaryFractions,
                                                       input.growth.w.boundaryFractions};
    for (std::size_t i = 0; i < 3; ++i)
        fractions[i].assign(views[i].begin(), views[i].end());
    return [input, fractions = std::move(fractions)](JobContext &c, mesh::Mesh &m) {
        auto o = input;
        o.growth.u.boundaryFractions = fractions[0];
        o.growth.v.boundaryFractions = fractions[1];
        o.growth.w.boundaryFractions = fractions[2];
        // Checked exact flat-output payload, before the legacy unbounded builder.
        const std::uint64_t u = o.geometry.cellsU, v = o.geometry.cellsV, w = o.geometry.cellsW;
        auto remaining = c.maximumResultBytes();
        auto charge = [&](std::uint64_t width, std::initializer_list<std::uint64_t> factors) {
            std::uint64_t n = width;
            for (auto f : factors) {
                if (f && n > remaining / f)
                    return false;
                n *= f;
            }
            if (n > remaining)
                return false;
            remaining -= n;
            return true;
        };
        const auto cellBytes = sizeof(mesh::Cell);
        if (!charge(sizeof(Vec3) + sizeof(mesh::Classification), {u + 1, v + 1, w + 1}) ||
            !charge(cellBytes + 8 * sizeof(Index), {u, v, w}) ||
            (o.emitBoundaryQuads && (!charge(2 * (cellBytes + 4 * sizeof(Index)), {u, v}) ||
                                     !charge(2 * (cellBytes + 4 * sizeof(Index)), {u, w}) ||
                                     !charge(2 * (cellBytes + 4 * sizeof(Index)), {v, w}))))
            return TaskOutcome::status(
                false, JobDomain::MappedHex,
                static_cast<std::uint32_t>(algo::MappedHexReason::BudgetExceeded),
                "mapped output budget exceeded");
        c.progress(JobStage::Generating);
        auto r = algo::generateMappedHexBlock(o, m);
        auto out = TaskOutcome::status(r.published(), JobDomain::MappedHex,
                                       static_cast<std::uint32_t>(r.reason),
                                       "mapped Hex8 result; see report.reason");
        out.report = std::move(r);
        return out;
    };
}
} // namespace paacs::sdk
