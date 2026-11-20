// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "paacs/sdk/JobService.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sizing/ControlField.hpp"
#include <atomic>
#include <chrono>
#include <future>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <thread>

using namespace paacs;
using namespace paacs::sdk;
using namespace std::chrono_literals;
namespace {
int checks = 0;
void require(bool ok, const char *message) {
    ++checks;
    if (!ok)
        throw std::runtime_error(message);
}
std::shared_ptr<const JobResult> complete(const JobHandle &job) {
    require(bool(job), "job admitted");
    require(job.waitFor(20s), "job completed within deadline");
    const auto s = job.snapshot();
    if (s.state != JobState::Succeeded)
        std::cerr << s.diagnostic.data() << " domain=" << int(s.domain)
                  << " reason=" << s.engineReason << '\n';
    require(s.state == JobState::Succeeded, "job succeeded");
    require(!job.cancel(), "late cancellation refused");
    auto r = job.result();
    require(bool(r), "immutable result available");
    return r;
}
void same(const mesh::Mesh &direct, const JobResult &job) {
    require(mesh::meshHash(direct) == job.summary.hash, "same direct-call hash");
    std::vector<std::uint8_t> a, b;
    require(serializeMeshArchive(direct, a).completed &&
                serializeMeshArchive(job.mesh, b).completed,
            "serialize comparison");
    require(a == b, "byte-identical archive");
}
struct Blocker {
    std::shared_ptr<std::promise<void>> entered = std::make_shared<std::promise<void>>();
    std::promise<void> release;
    std::shared_future<void> gate = release.get_future().share();
    MeshTask task() {
        return [entered = entered, gate = gate](JobContext &c, mesh::Mesh &) {
            c.progress(JobStage::Generating, 2, 7);
            entered->set_value();
            gate.wait();
            return TaskOutcome::status(false, JobDomain::Custom, 99, "test gate released");
        };
    }
    void wait() {
        require(entered->get_future().wait_for(10s) == std::future_status::ready, "blocker entered");
    }
};
void blockerLifetime() {
    JobService service;
    JobHandle running;
    {
        Blocker block;
        running = service.submit(block.task()).job;
        block.wait();
        const auto snapshot = running.snapshot();
        require(snapshot.completed == 2 && snapshot.total == 7,
                "blocker signals entry after progress is available");
    }
    require(running.waitFor(10s), "destroyed blocker releases its waiting task");
    const auto snapshot = running.snapshot();
    require(snapshot.state == JobState::Failed && snapshot.error == JobError::EngineRejected &&
                snapshot.domain == JobDomain::Custom && snapshot.engineReason == 99,
            "blocker task owns its synchronization state through teardown");
}
void adapters() {
    JobService service({16, 64 * 1024 * 1024, 128 * 1024 * 1024, true});
    Blocker block;
    auto waiting = service.submit(block.task()).job;
    block.wait();
    require(waiting.snapshot().completed == 2 && waiting.snapshot().total == 7,
            "measured progress exposed");
    mesh::Mesh primitive, adaptive, mapped, geometry, layers;
    MeshTask p, a, h, g, l;
    // Each factory must own all borrowed arrays/profiles after this scope exits.
    {
        algo::PrimitiveShape shape;
        shape.height = 1;
        sizing::CappedSizeSpecification spec;
        spec.globalSize = .4;
        spec.controls = {{0, 0, .2, 1.3}};
        auto profile = sizing::compileCappedSizeProfile(spec);
        require(profile.compiled, "profile compiled");
        algo::PrimitiveMeshingOptions o;
        o.elements = algo::PrimitiveElements::Tet4;
        o.targetSize = .4;
        o.maximumChordError = .04;
        o.sizeProfile = &profile.profile;
        o.threads = 4;
        p = primitiveTask(shape, o);
        // Direct calls run only once no service dispatcher is active below.
        block.release.set_value();
        require(waiting.waitFor(10s), "unblocked");
        require(algo::meshPrimitive(shape, o, primitive).published, "direct primitive");
        sizing::ControlFieldSpecification f;
        f.nodes = {sizing::field::distanceBall({.5, .5, .5}, .1),
                   sizing::field::threshold(0, 0, .35, .2, .6)};
        f.output = 1;
        auto built = sizing::compileControlField(f);
        require(built.compiled, "field compiled");
        auto field = std::make_shared<sizing::ControlField>(std::move(built.field));
        algo::AdaptiveBoxRequest ar;
        ar.targetSize = .6;
        ar.field = field.get();
        ar.threads = 4;
        a = adaptiveBoxTask(ar, field);
        require(algo::meshAdaptiveBox(ar, adaptive).published, "direct adaptive");
        bool refused = false;
        try {
            (void)adaptiveBoxTask(ar);
        } catch (const std::invalid_argument &) {
            refused = true;
        }
        require(refused, "unowned raw field rejected");
        algo::MappedHexBlockRequest hr;
        hr.geometry.cellsU = 2;
        hr.geometry.cellsV = 3;
        hr.executionThreads = 4;
        std::vector<double> fractions{0, .45, 1};
        hr.growth.u.boundaryFractions = fractions;
        h = mappedHexTask(hr);
        require(algo::generateMappedHexBlock(hr, mapped).published(), "direct mapped");
        auto kernel = geo::makeBuiltinKernel();
        std::shared_ptr<const geo::Model> model = kernel->makeBox({0, 0, 0}, {1, 1, 1});
        GeometryJobInput gr;
        gr.model = model;
        gr.options.targetSize = .5;
        gr.options.threads = 4;
        g = geometryTask(gr);
        require(algo::meshGeometryModel(*model, gr.options, geometry).published, "direct geometry");
        auto source = std::make_shared<mesh::Mesh>();
        source->addNode({0, 0, 0}, {2, 1});
        source->addNode({1, 0, 0}, {2, 1});
        source->addNode({0, 1, 0}, {2, 1});
        source->addCell(mesh::CellType::Tri3, {0, 1, 2}, {2, 1});
        algo::LayeredPatchRequest lr;
        lr.threads = 4;
        lr.sizeProfile = &profile.profile;
        l = layeredPatchTask(source, lr);
        require(algo::meshLayeredPatch(*source, lr, layers).published, "direct layers");
    }
    same(primitive, *complete(service.submit(std::move(p)).job));
    same(adaptive, *complete(service.submit(std::move(a)).job));
    same(mapped, *complete(service.submit(std::move(h)).job));
    same(geometry, *complete(service.submit(std::move(g)).job));
    same(layers, *complete(service.submit(std::move(l)).job));
    // Explicit axial arrays (separate from profile mode) also survive their owner.
    MeshTask axial;
    {
        std::vector<double> positions{0, .25, .5, .75, 1};
        algo::PrimitiveShape shape;
        shape.height = 1;
        algo::PrimitiveMeshingOptions o;
        o.axialPositions = positions;
        o.targetSize = .5;
        o.maximumChordError = .04;
        axial = primitiveTask(shape, o);
    }
    (void)complete(service.submit(std::move(axial)).job);
    std::cout << "{\"adapters\":6,\"archive_identity_pairs\":5}\n";
}
void cancellationAndFailure() {
    JobService service({6, 8 * 1024 * 1024, 16 * 1024 * 1024});
    Blocker block;
    auto running = service.submit(block.task()).job;
    block.wait();
    std::atomic<bool> executed = false;
    auto queued = service
                      .submit([&](JobContext &, mesh::Mesh &) {
                          executed = true;
                          return TaskOutcome{};
                      })
                      .job;
    require(queued.cancel(), "queued cancellation accepted");
    require(queued.waitFor(0ms), "queued cancellation immediate");
    require(queued.snapshot().state == JobState::Cancelled && !executed, "queued task never ran");
    require(running.cancel(), "running cancellation accepted");
    block.release.set_value();
    require(running.waitFor(10s), "cancelled running job exits");
    require(running.snapshot().state == JobState::Cancelled && !running.result(),
            "running cancellation suppresses publication");
    auto throwing = service
                        .submit([](JobContext &, mesh::Mesh &) -> TaskOutcome {
                            throw std::runtime_error("failure");
                        })
                        .job;
    require(throwing.waitFor(10s) && throwing.snapshot().error == JobError::Exception,
            "exception contained");
    auto allocation =
        service.submit([](JobContext &, mesh::Mesh &) -> TaskOutcome { throw std::bad_alloc(); })
            .job;
    require(allocation.waitFor(10s) && allocation.snapshot().error == JobError::AllocationFailed,
            "allocation failure contained");
    algo::PrimitiveMeshingOptions bad;
    bad.targetSize = -1;
    auto rejected = service.submit(primitiveTask({}, bad)).job;
    require(rejected.waitFor(10s) && rejected.snapshot().state == JobState::Failed,
            "engine refusal returned");
    require(std::get<algo::PrimitiveMeshingResult>(rejected.report()).reason ==
                algo::PrimitiveReason::InvalidRequest,
            "full refusal report retained");
    service.shutdown();
    require(service.submit(mappedHexTask({})).reason == SubmitReason::Stopped,
            "submission after shutdown rejected");
    service.shutdown();
}
void budgetsAndLifetime() {
    std::shared_ptr<const JobResult> surviving;
    {
        JobService service({2, 680, 680});
        auto j = service.submit(mappedHexTask({})).job;
        surviving = complete(j);
        require(!surviving->summary.hashAvailable && surviving->summary.hash == 0,
                "hash opt-in defaults off");
        require(surviving->summary.logicalBytes == 680 && service.retainedResultBytes() == 680,
                "exact payload accounting");
        j = {};
        auto refused = service.submit(mappedHexTask({})).job;
        require(refused.waitFor(10s) && refused.snapshot().error == JobError::ResultBudget,
                "external result still charged");
        surviving.reset();
        require(service.retainedResultBytes() == 0, "last result reference releases charge");
        surviving = complete(service.submit(mappedHexTask({})).job);
    }
    require(surviving->mesh.nodeCount() == 8, "result survives service destruction");
    {
        JobService service({2, 1024, 2048});
        Blocker block;
        auto a = service.submit(block.task()).job;
        block.wait();
        auto b = service.submit(mappedHexTask({})).job;
        require(service.submit(mappedHexTask({})).reason == SubmitReason::JobLimit,
                "queue/handle admission bound");
        require(b.cancel(), "cancel blocked job");
        block.release.set_value();
        require(a.waitFor(10s), "block released");
    }
    {
        JobService service({2, 128, 128});
        auto j = service.submit(mappedHexTask({})).job;
        require(j.waitFor(10s) && !j.result(), "mapped preflight output bound");
    }
}
void concurrentServices() {
    JobService a({40, 4096, 65536}), b({40, 4096, 65536});
    std::atomic<int> active = 0, collisions = 0;
    auto task = [&](JobContext &c, mesh::Mesh &m) {
        if (active.fetch_add(1) != 0)
            ++collisions;
        std::vector<std::uint64_t> values(5000);
        core::parallelForWithThreads(0, values.size(), 4,
                                     [&](std::size_t i) { values[i] = i * i; });
        auto out = mappedHexTask({})(c, m);
        --active;
        return out;
    };
    std::vector<JobHandle> first, second;
    std::thread t1([&] {
        for (int i = 0; i < 24; ++i)
            first.push_back(a.submit(task).job);
    });
    std::thread t2([&] {
        for (int i = 0; i < 24; ++i)
            second.push_back(b.submit(task).job);
    });
    t1.join();
    t2.join();
    std::atomic<bool> stop = false;
    std::thread poll([&] {
        while (!stop) {
            for (auto &j : first)
                (void)j.snapshot();
            std::this_thread::yield();
        }
    });
    for (auto &j : first)
        (void)complete(j);
    for (auto &j : second)
        (void)complete(j);
    stop = true;
    poll.join();
    require(collisions == 0, "all service dispatches serialized");
    // Waiting on a different service's gate can be cancelled immediately.
    Blocker block;
    auto running = a.submit(block.task()).job;
    block.wait();
    auto waiting = b.submit(mappedHexTask({})).job;
    require(waiting.cancel() && waiting.waitFor(0ms), "cross-service wait cancellation");
    block.release.set_value();
    require(running.waitFor(10s), "cross-service blocker exits");
    std::cout << "{\"concurrent_jobs\":48,\"dispatch_collisions\":" << collisions << "}\n";
}
void publicationRace() {
    JobService service({4, 4096, 16384});
    for (int i = 0; i < 100; ++i) {
        auto job = service.submit(mappedHexTask({})).job;
        require(bool(job), "race job admitted");
        bool accepted = false;
        std::thread cancel([&] { accepted = job.cancel(); });
        require(job.waitFor(10s), "race terminates");
        cancel.join();
        if (accepted)
            require(job.snapshot().state == JobState::Cancelled && !job.result(),
                    "accepted cancellation cannot publish");
        else
            require(job.snapshot().state == JobState::Succeeded && bool(job.result()),
                    "publication wins late cancel");
    }
    std::cout << "{\"publication_races\":100}\n";
}
void shutdownWhileActive() {
    JobService service;
    std::promise<void> entered;
    auto runningOwner = std::make_shared<int>(1), queuedOwner = std::make_shared<int>(2);
    std::weak_ptr<int> runningWeak = runningOwner, queuedWeak = queuedOwner;
    auto running = service
                       .submit([owner = runningOwner, &entered](JobContext &c, mesh::Mesh &) {
                           (void)owner;
                           entered.set_value();
                           while (!c.cancelled())
                               std::this_thread::yield();
                           return TaskOutcome{};
                       })
                       .job;
    require(entered.get_future().wait_for(10s) == std::future_status::ready,
            "shutdown task entered");
    auto queued = service
                      .submit([owner = queuedOwner](JobContext &, mesh::Mesh &) {
                          (void)owner;
                          return TaskOutcome{};
                      })
                      .job;
    runningOwner.reset();
    queuedOwner.reset();
    service.shutdown();
    require(running.snapshot().state == JobState::Cancelled, "shutdown cancels running work");
    require(queued.snapshot().state == JobState::Cancelled, "shutdown cancels queued work");
    require(runningWeak.expired() && queuedWeak.expired(),
            "shutdown releases owned input despite retained handles");
    require(!running.result() && !queued.result(), "shutdown publishes no partial result");
}
} // namespace
int main() {
    try {
        blockerLifetime();
        adapters();
        cancellationAndFailure();
        budgetsAndLifetime();
        concurrentServices();
        publicationRace();
        shutdownWhileActive();
        std::cout << "{\"checks\":" << checks << ",\"passed\":true}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAILED: " << e.what() << '\n';
        return 1;
    }
}
