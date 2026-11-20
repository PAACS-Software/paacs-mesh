// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/sizing/BoundedSizeField.hpp"
#include <array>
#include <atomic>
#include <limits>
#include <span>

namespace paacs::sizing {
using FieldId = std::uint32_t;
inline constexpr FieldId invalidFieldId = UINT32_MAX;
inline constexpr std::size_t maximumControlFieldNodes = 256;
enum class FieldOperation : std::uint8_t {
    Constant,
    DistancePoint,
    DistanceSegment,
    DistanceTriangle,
    DistanceBox,
    DistanceBall,
    ThresholdLinear,
    ThresholdSmooth,
    Minimum,
    Maximum,
    Scale,
    Clamp
};
// Use the named factories below. Dependencies must precede their consumer.
struct FieldNode {
    FieldOperation operation = FieldOperation::Constant;
    FieldId first = invalidFieldId, second = invalidFieldId;
    Vec3 a{}, b{}, c{};
    std::array<double, 4> parameters{1, 0, 0, 0};
};
namespace field {
FieldNode constant(double value) noexcept;
FieldNode distancePoint(Vec3 point) noexcept;
FieldNode distanceSegment(Vec3 start, Vec3 end) noexcept;
FieldNode distanceTriangle(Vec3 a, Vec3 b, Vec3 c) noexcept;
FieldNode distanceBox(Vec3 lower, Vec3 upper) noexcept;
FieldNode distanceBall(Vec3 center, double radius) noexcept;
FieldNode threshold(FieldId distance, double nearDistance, double farDistance, double nearSize,
                    double farSize, bool smooth = false) noexcept;
FieldNode minimum(FieldId first, FieldId second) noexcept;
FieldNode maximum(FieldId first, FieldId second) noexcept;
FieldNode scale(FieldId input, double factor) noexcept;
FieldNode clamp(FieldId input, double minimum, double maximum) noexcept;
} // namespace field
struct ControlFieldSpecification {
    std::vector<FieldNode> nodes;
    FieldId output = invalidFieldId;
    // h(x) = sizeFactor * clamp(graph(x), minimumSize, maximumSize).
    // Limits are applied BEFORE the factor, matching this documented order.
    double minimumSize = .001, maximumSize = 1, sizeFactor = 1;
    std::uint64_t maximumNodes = maximumControlFieldNodes;
    std::uint64_t maximumBytes = 1024 * 1024;
    const std::atomic<bool> *cancellation = nullptr;
};
enum class ControlFieldReason : std::uint8_t {
    None,
    InvalidEnvelope,
    InvalidNode,
    InvalidDependency,
    InvalidOutput,
    DegenerateGeometry,
    BudgetExceeded,
    Cancelled,
    AllocationFailed,
    InternalError
};
struct FieldRange {
    double lower, upper;
};
struct ControlFieldOutcome;
class ControlField final : public BoundedSizeField {
  public:
    bool valid() const noexcept override { return valid_ && !nodes_.empty(); }
    double eval(const Vec3 &point) const noexcept override;
    double minimumOnBox(const Vec3 &lower, const Vec3 &upper) const noexcept override;
    // Conservative enclosure. It may be wider than the exact attainable range
    // for correlated Min/Max branches. Degenerate boxes (points) are accepted.
    FieldRange rangeOnBox(const Vec3 &lower, const Vec3 &upper) const noexcept;
    std::uint64_t boxQueryWorkBound() const noexcept override { return queryWork_; }
    std::size_t nodeCount() const noexcept { return nodes_.size(); }
    std::uint64_t logicalBytes() const noexcept { return nodes_.size() * sizeof(FieldNode); }
    double minimumSize() const noexcept { return minimum_ * factor_; }
    double maximumSize() const noexcept { return maximum_ * factor_; }
    double sizeFactor() const noexcept { return factor_; }
    std::span<const FieldNode> nodes() const noexcept { return nodes_; }

  private:
    friend ControlFieldOutcome compileControlField(const ControlFieldSpecification &) noexcept;
    bool valid_ = false;
    double minimum_ = .001, maximum_ = 1, factor_ = 1;
    std::uint64_t queryWork_ = 0;
    std::vector<FieldNode> nodes_;
};
struct ControlFieldOutcome {
    bool compiled = false;
    ControlFieldReason reason = ControlFieldReason::InternalError;
    FieldId rejectedNode = invalidFieldId;
    std::uint64_t inputNodes = 0, retainedNodes = 0, logicalBytes = 0;
    ControlField field;
};
// Compilation validates every supplied node, removes unreachable nodes and
// retains a flat topological program. Queries allocate nothing and are thread
// safe. No sampled grid can erase a tiny control region. This is a scalar field;
// smooth thresholds do not promise a mesh's adjacent-element growth ratio.
ControlFieldOutcome compileControlField(const ControlFieldSpecification &) noexcept;
const char *controlFieldReasonName(ControlFieldReason) noexcept;
} // namespace paacs::sizing
