// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/CanonicalCurveMesher.hpp"

#include "paacs/mesh/MeshAudit.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <map>
#include <new>
#include <stdexcept>

namespace paacs::algo {
namespace {

std::uint64_t mix(std::uint64_t h, std::uint64_t value) noexcept {
    value ^= value >> 30U;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27U;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6U) + (h >> 2U);
    return h;
}

bool validClassification(const mesh::Classification& value) noexcept {
    return value.dim >= 0 && value.dim <= 3 && value.tag >= 0;
}

}  // namespace

const CanonicalCurveRecord* CanonicalCurveSet::curve(
        int edgeTag) const noexcept {
    const auto found = std::lower_bound(
        curves_.begin(), curves_.end(), edgeTag,
        [](const CanonicalCurveRecord& value, int tag) {
            return value.edgeTag < tag;
        });
    return found != curves_.end() && found->edgeTag == edgeTag
               ? &*found
               : nullptr;
}

std::vector<Index> CanonicalCurveSet::orientedNodes(
        int edgeTag, bool reversed) const {
    const CanonicalCurveRecord* value = curve(edgeTag);
    if (!value) return {};
    std::vector<Index> out = value->naturalNodeOrder;
    if (reversed) std::reverse(out.begin(), out.end());
    return out;
}

CanonicalCurveOutcome buildCanonicalCurves(
        const geo::Model& model, const sizing::MeshSpec& spec) noexcept {
    CanonicalCurveOutcome out;
    try {
        auto edges = model.edges();
        if (edges.empty()) {
            out.status = CanonicalCurveStatus::Rejected;
            out.reason = CanonicalCurveReason::EmptyModel;
            return out;
        }
        std::sort(edges.begin(), edges.end(),
                  [](const geo::GEdge* a, const geo::GEdge* b) {
                      return a->tag() < b->tag();
                  });
        for (std::size_t i = 0; i < edges.size(); ++i) {
            if (!edges[i] || edges[i]->tag() < 0 ||
                (i > 0 && edges[i - 1]->tag() == edges[i]->tag())) {
                out.status = CanonicalCurveStatus::Rejected;
                out.reason = CanonicalCurveReason::InvalidEdge;
                return out;
            }
        }

        // The established edge implementation already evaluates sizing,
        // curvature and endpoint classifications.  V0.10 records its output
        // as an immutable reusable artifact rather than introducing a second
        // discretization algorithm.
        meshEdges(model, spec, out.result.mesh_);

        std::map<int, std::vector<Index>> lineCellsByTag;
        for (Index ci = 0; ci < out.result.mesh_.cellCount(); ++ci) {
            const auto& cell = out.result.mesh_.cell(ci);
            if (cell.type != mesh::CellType::Line2) continue;
            if (cell.onEntity.dim != 1 || cell.onEntity.tag < 0) {
                out.status = CanonicalCurveStatus::Rejected;
                out.reason = CanonicalCurveReason::InvalidClassification;
                return out;
            }
            lineCellsByTag[cell.onEntity.tag].push_back(ci);
        }

        out.result.curves_.reserve(edges.size());
        for (const geo::GEdge* edge : edges) {
            const auto found = lineCellsByTag.find(edge->tag());
            if (found == lineCellsByTag.end() || found->second.empty()) {
                out.status = CanonicalCurveStatus::Rejected;
                out.reason = CanonicalCurveReason::MissingLineCells;
                return out;
            }

            std::map<Index, std::vector<std::pair<Index, Index>>> adjacency;
            for (const Index ci : found->second) {
                const Index* nodes = out.result.mesh_.cellNodes(ci);
                adjacency[nodes[0]].push_back({nodes[1], ci});
                adjacency[nodes[1]].push_back({nodes[0], ci});
            }
            for (const auto& [node, links] : adjacency) {
                (void)node;
                if (links.size() > 2U) {
                    out.status = CanonicalCurveStatus::Rejected;
                    out.reason = CanonicalCurveReason::BranchingCurve;
                    return out;
                }
            }

            Index start = kInvalidIndex;
            Index end = kInvalidIndex;
            for (const auto& [node, links] : adjacency) {
                if (links.size() != 1U) continue;
                const auto cls = out.result.mesh_.nodeClass(node);
                if (cls.dim != 0) continue;
                if (cls.tag == edge->startVertexTag()) start = node;
                if (cls.tag == edge->endVertexTag()) end = node;
            }
            if (start == kInvalidIndex || end == kInvalidIndex) {
                out.status = CanonicalCurveStatus::Rejected;
                out.reason = CanonicalCurveReason::EndpointMismatch;
                return out;
            }

            CanonicalCurveRecord record;
            record.edgeTag = edge->tag();
            record.startVertexTag = edge->startVertexTag();
            record.endVertexTag = edge->endVertexTag();
            record.naturalNodeOrder.reserve(found->second.size() + 1U);
            record.naturalLineOrder.reserve(found->second.size());
            record.naturalNodeOrder.push_back(start);

            Index previous = kInvalidIndex;
            Index current = start;
            while (current != end) {
                const auto at = adjacency.find(current);
                if (at == adjacency.end()) {
                    out.status = CanonicalCurveStatus::Rejected;
                    out.reason = CanonicalCurveReason::DisconnectedCurve;
                    return out;
                }
                std::pair<Index, Index> selected{kInvalidIndex, kInvalidIndex};
                for (const auto& link : at->second) {
                    if (link.first == previous) continue;
                    if (selected.first != kInvalidIndex) {
                        out.status = CanonicalCurveStatus::Rejected;
                        out.reason = CanonicalCurveReason::BranchingCurve;
                        return out;
                    }
                    selected = link;
                }
                if (selected.first == kInvalidIndex) {
                    out.status = CanonicalCurveStatus::Rejected;
                    out.reason = CanonicalCurveReason::DisconnectedCurve;
                    return out;
                }
                previous = current;
                current = selected.first;
                record.naturalLineOrder.push_back(selected.second);
                record.naturalNodeOrder.push_back(current);
                if (record.naturalLineOrder.size() > found->second.size()) {
                    out.status = CanonicalCurveStatus::Rejected;
                    out.reason = CanonicalCurveReason::DisconnectedCurve;
                    return out;
                }
            }
            if (record.naturalLineOrder.size() != found->second.size()) {
                out.status = CanonicalCurveStatus::Rejected;
                out.reason = CanonicalCurveReason::DisconnectedCurve;
                return out;
            }
            out.result.curves_.push_back(std::move(record));
        }

        std::sort(out.result.curves_.begin(), out.result.curves_.end(),
                  [](const auto& a, const auto& b) {
                      return a.edgeTag < b.edgeTag;
                  });
        std::uint64_t hash = 0x5041414353435552ULL;
        for (const auto& record : out.result.curves_) {
            hash = mix(hash, static_cast<std::uint32_t>(record.edgeTag));
            hash = mix(hash,
                       static_cast<std::uint32_t>(record.startVertexTag));
            hash = mix(hash,
                       static_cast<std::uint32_t>(record.endVertexTag));
            hash = mix(hash, record.naturalNodeOrder.size());
            for (const Index node : record.naturalNodeOrder) {
                hash = mix(hash, static_cast<std::uint64_t>(node));
                const Vec3& point = out.result.mesh_.node(node);
                hash = mix(hash, std::bit_cast<std::uint64_t>(point.x));
                hash = mix(hash, std::bit_cast<std::uint64_t>(point.y));
                hash = mix(hash, std::bit_cast<std::uint64_t>(point.z));
            }
        }
        out.result.hash_ = hash;
        out.status = CanonicalCurveStatus::Built;
        out.reason = CanonicalCurveReason::None;
        return out;
    } catch (const std::bad_alloc&) {
        out.status = CanonicalCurveStatus::Failed;
        out.reason = CanonicalCurveReason::AllocationFailed;
    } catch (...) {
        out.status = CanonicalCurveStatus::Failed;
        out.reason = CanonicalCurveReason::InternalError;
    }
    return out;
}

bool appendCanonicalCurves(const CanonicalCurveSet& curves,
                           mesh::Mesh& destination) noexcept {
    try {
        mesh::Mesh staged = destination;
        staged.reserveAppend(
            static_cast<std::size_t>(curves.mesh().nodeCount()),
            static_cast<std::size_t>(curves.mesh().cellCount()),
            static_cast<std::size_t>(curves.mesh().connectivityCount()));
        std::vector<Index> mapping(
            static_cast<std::size_t>(curves.mesh().nodeCount()),
            kInvalidIndex);
        for (Index ni = 0; ni < curves.mesh().nodeCount(); ++ni) {
            const auto cls = curves.mesh().nodeClass(ni);
            if (!validClassification(cls)) return false;
            mapping[static_cast<std::size_t>(ni)] = staged.addNodeUnique(
                curves.mesh().node(ni), cls);
        }
        for (Index ci = 0; ci < curves.mesh().cellCount(); ++ci) {
            const auto& cell = curves.mesh().cell(ci);
            if (cell.type != mesh::CellType::Line2 ||
                !validClassification(cell.onEntity))
                return false;
            const Index* nodes = curves.mesh().cellNodes(ci);
            staged.addCell(mesh::CellType::Line2,
                           {mapping[static_cast<std::size_t>(nodes[0])],
                            mapping[static_cast<std::size_t>(nodes[1])]},
                           cell.onEntity);
        }
        destination.swap(staged);
        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace paacs::algo
