// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
// Private optional adapter helper; never installed as public API.
#include "paacs/geometry/OcctModel.hpp"
#include <array>
#include <atomic>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace paacs::geo::detail {
enum class OcctNativeValidationReason : std::uint8_t {
    None, InvalidRequest, BudgetExceeded, Cancelled, AllocationFailed,
    NativeQueryFailed, OutsideMaterial, CopyGeometryMismatch,
    ReferenceStateMismatch, SourceChanged, ThreadLaunchFailed, InternalError
};
struct OcctNativeValidationOptions {
    unsigned maximumThreads = 1; // Explicit admission, hard cap8.
    std::uint64_t maximumPoints = 1'000'000, maximumWork = 20'000'000;
    std::uint64_t maximumWorkingBytes = 128ULL*1024*1024;
    std::uint64_t maximumSnapshotBytes = 16ULL*1024*1024;
    const std::atomic_bool* cancellation = nullptr;
    // Explicit caller policy, in model units. Zero retains strict inclusion.
    // Only OUT points can use this measured band; UNKNOWN always refuses.
    double maximumOutsideDistance = 0;
    bool useTrimmedRayWitnesses = false;
};
struct OcctNativeValidationOutcome {
    bool completed = false, sourceModelUnmodified = false;
    OcctNativeValidationReason reason = OcctNativeValidationReason::InternalError;
    unsigned threadsUsed = 0;
    std::uint64_t points = 0, nativeQueries = 0, originalReferenceQueries = 0;
    std::uint64_t strictAdmissionReuses = 0;
    std::uint64_t trimmedRayAdmissions = 0;
    std::uint64_t copiedGeometryValuePairs = 0, workUsed = 0;
    std::uint64_t maximumWorkingPayloadEstimate = 0, sourceSnapshotBytes = 0;
    std::uint64_t firstFailedPoint = (std::numeric_limits<std::uint64_t>::max)();
    std::uint64_t pointsWithinDistanceBand = 0, nativeDistanceQueries = 0;
    double maximumOutsideDistance = 0;
};
// Evidence from a successful serial strict classification. A later pass
// can reuse a positive result only on the same model, solid, BREP and binary
// point coordinates.
class OcctNativeStrictAdmissions final {
public:
    explicit OcctNativeStrictAdmissions(std::uint64_t maximumRetainedBytes) noexcept
        : maximumRetainedBytes_(maximumRetainedBytes) {}
    std::uint64_t retainedPayloadBytes() const noexcept;

private:
    friend struct OcctNativeInteriorCheck;
    using Key = std::array<std::uint64_t,3>;
    std::uint64_t maximumRetainedBytes_;
    const OcctModel* model_ = nullptr;
    int solidTag_ = -1;
    std::uint64_t tolerance_ = 0;
    bool includesRayAdmissions_ = false;
    std::string source_;
    std::vector<Key> admitted_;
};
// Every input point must classify IN/ON with the original body's native
// tolerance, unless an explicit positive distance allowance is supplied.
// OUT distances are then measured serially against the original solid after
// copied classifier state checks. Each worker owns a deep copy and one local
// classifier; no mutable native adaptor is shared. Copied face/curve values
// and tolerances are checked, up to256 spread queries also use the original
// serial classifier, and original BREP bytes are compared after all queries.
// Work bounds PAACS checks/query calls and logical temporary payload. Native
// copy/classifier allocations/internals, process RSS and time remain outside
// these bounds. Cancellation is checked between cooperative native calls.
// No mesh or externally visible geometry is modified or published.
OcctNativeValidationOutcome validateOcctNativeInteriorPoints(const OcctModel&,
    int solidTag, std::span<const Vec3>, const OcctNativeValidationOptions&,
    const OcctNativeStrictAdmissions* prior = nullptr) noexcept;
// Return every decision after the same copy, state and source checks. An
// outside point marks its row false. On any operational failure the mask is
// unchanged. Completed refers to classification, not admission of all points.
// A prior strict receipt may cover positive rows. Unseen and previously
// refused points still query the original solid.
OcctNativeValidationOutcome classifyOcctNativeInteriorPoints(const OcctModel&,
    int solidTag, std::span<const Vec3>, const OcctNativeValidationOptions&,
    std::vector<unsigned char>& admitted, OcctNativeStrictAdmissions* retained = nullptr,
    const OcctNativeStrictAdmissions* prior = nullptr) noexcept;
const char* occtNativeValidationReasonName(OcctNativeValidationReason) noexcept;
} // namespace paacs::geo::detail
