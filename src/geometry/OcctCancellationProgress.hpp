// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private bridge from the existing call-scoped cancellation flag into OCCT's
// cooperative progress channel. No deadline, geometry tolerance or policy is
// invented here. The caller must keep the flag alive through the kernel call.
#include <Message_ProgressIndicator.hxx>
#include <atomic>
#include <cstdint>

namespace paacs::geo::detail {
class OcctCancellationProgress : public Message_ProgressIndicator {
public:
    explicit OcctCancellationProgress(const std::atomic_bool* cancellation)
        : cancellation_(cancellation) {}
    // Caller-owned indicator, handle and initial root scope payload only.
    // OCCT's subsequently created scopes/workspaces and allocator bookkeeping
    // remain kernel allocations; this is not a process-memory bound.
    static constexpr std::uint64_t retainedPayloadBytes() noexcept {
        return sizeof(OcctCancellationProgress) + sizeof(Message_ProgressScope) +
               sizeof(occ::handle<Message_ProgressIndicator>);
    }

protected:
    bool UserBreak() override {
        // OCCT may poll this concurrently. Do not touch mutable report/budget
        // state or invoke Show/GetPosition from this callback.
        return cancellation_ && cancellation_->load(std::memory_order_relaxed);
    }
    void Show(const Message_ProgressScope&, bool) override {}

private:
    const std::atomic_bool* cancellation_;
};

inline occ::handle<Message_ProgressIndicator> makeOcctCancellationProgress(
    const std::atomic_bool* cancellation) {
    if (!cancellation) return {};
    return new OcctCancellationProgress(cancellation);
}
} // namespace paacs::geo::detail
