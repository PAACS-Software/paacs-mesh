// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPartition.hpp"
#include "OcctModelInternal.hpp"
#include "OcctTrimTolerance.hpp"
#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <gp_Circ.hxx>
#include <set>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include "OcctSolidClassifier.hpp"
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <Standard_OutOfMemory.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Pln.hxx>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numbers>
#include <ostream>
#include <streambuf>
#include <string>

namespace paacs::geo::detail {
namespace {
using Clock = std::chrono::steady_clock;
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
using Class = mesh::Classification;
struct Stop { OcctPartitionReason reason; };
void require(bool condition, OcctPartitionReason reason) {
    if (!condition) throw Stop{reason};
}
double seconds(Clock::time_point begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}
struct Timer {
    double& destination;
    Clock::time_point begin = Clock::now();
    ~Timer() { destination += seconds(begin); }
};
struct Checkpoints {
    const OcctPartitionOptions& options;
    OcctPartitionOutcome& outcome;
    void check() const {
        require(!options.cancellation || !options.cancellation->load(std::memory_order_relaxed),
                OcctPartitionReason::Cancelled);
    }
    void work() {
        check();
        require(outcome.historyWork < options.maximumHistoryWork,
                OcctPartitionReason::BudgetExceeded);
        ++outcome.historyWork;
    }
};
// Explicit retry only: native span-aware quadrature at two accuracies, with
// the original strict conservation limit retained. Source model metadata and
// the native topology/ordering are never rewritten.
void retrySpanAwareVolume(const ShapeMap& sourceSolids,const ShapeMap& pieces,
                         const std::vector<int>& selected,Checkpoints& checkpoints) {
    auto& out=checkpoints.outcome;
    out.legacySourceVolume=out.sourceVolume;
    out.legacyPartitionVolume=out.partitionVolume;
    std::array<double,2> sourceMass{},pieceMass{},sourceError{},pieceError{};
    for(int pass=0;pass<2;++pass) {
        const double epsilon=pass?1e-10:1e-8;
        auto integrate=[&](const TopoDS_Shape& shape,double& sum,double& maximumError) {
            checkpoints.work();GProp_GProps properties;
            const auto error=BRepGProp::VolumePropertiesGK(shape,properties,epsilon,
                                                          false,true,false,false,false);
            checkpoints.check();
            require(std::isfinite(error)&&error>=0&&std::isfinite(properties.Mass())&&
                    properties.Mass()>0,OcctPartitionReason::InvalidResult);
            sum+=properties.Mass();maximumError=std::max(maximumError,error);
        };
        for(int tag:selected)integrate(sourceSolids(tag),sourceMass[pass],sourceError[pass]);
        for(int tag=1;tag<=pieces.Extent();++tag)integrate(pieces(tag),pieceMass[pass],pieceError[pass]);
    }
    out.sourceReferenceStability=std::abs(sourceMass[1]-sourceMass[0])/sourceMass[1];
    out.partitionReferenceStability=std::abs(pieceMass[1]-pieceMass[0])/pieceMass[1];
    require(std::isfinite(out.sourceReferenceStability)&&std::isfinite(out.partitionReferenceStability)&&
            out.sourceReferenceStability<=checkpoints.options.maximumRelativeVolumeChange&&
            out.partitionReferenceStability<=checkpoints.options.maximumRelativeVolumeChange,
            OcctPartitionReason::VolumeMismatch);
    out.sourceVolume=sourceMass[1];out.partitionVolume=pieceMass[1];
    out.relativeVolumeChange=std::abs(out.partitionVolume-out.sourceVolume)/out.sourceVolume;
    require(std::isfinite(out.relativeVolumeChange),OcctPartitionReason::VolumeMismatch);
    out.sourceIntegrationError=sourceError[1];out.maximumPieceIntegrationError=pieceError[1];
    out.spanAwareVolumeReferenceUsed=true;
}
// Refuse while serializing, before an unbounded string allocation. OCCT's own
// transient allocations remain outside this receipt's storage bound.
class AuditBuffer final : public std::streambuf {
  public:
    explicit AuditBuffer(std::uint64_t maximum) : maximum_(maximum) {}
    std::string value;
  protected:
    std::streamsize xsputn(const char* data, std::streamsize count) override {
        if (count < 0 || static_cast<std::uint64_t>(count) > maximum_ - value.size())
            throw Stop{OcctPartitionReason::BudgetExceeded};
        value.append(data, static_cast<std::size_t>(count));
        return count;
    }
    int_type overflow(int_type character) override {
        if (traits_type::eq_int_type(character, traits_type::eof()))
            return traits_type::not_eof(character);
        const char byte = traits_type::to_char_type(character);
        xsputn(&byte, 1);
        return character;
    }
  private:
    std::uint64_t maximum_;
};
std::string serialize(const TopoDS_Shape& shape, std::uint64_t maximum) {
    AuditBuffer buffer(maximum);
    std::ostream stream(&buffer);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    BRepTools::Write(shape, stream);
    return std::move(buffer.value);
}
void assign(Class& existing, Class candidate) {
    if (existing.dim < 0 || candidate.dim < existing.dim) {
        existing = candidate;
    } else if (candidate.dim == existing.dim && candidate.tag != existing.tag) {
        // Preserve ambiguity until a lower-dimensional native history wins.
        existing.tag = -2;
    }
}
struct Provenance {
    const ShapeMap& vertices;
    const ShapeMap& edges;
    const ShapeMap& faces;
    OcctPartitionOutcome& outcome;
    Checkpoints& checkpoints;
    void apply(const TopoDS_Shape& shape, Class parent, bool allowFace = true) {
        if (shape.IsNull()) return;
        const std::array<TopAbs_ShapeEnum,3> kinds{TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE};
        const std::array<const ShapeMap*,3> maps{&vertices, &edges, &faces};
        const std::array<std::vector<Class>*,3> classes{
            &outcome.parentByVertex, &outcome.parentByEdge, &outcome.parentByFace};
        for (int dim = 0; dim <= std::min(parent.dim, allowFace ? 2 : 1); ++dim) {
            ShapeMap children;
            TopExp::MapShapes(shape, kinds[dim], children);
            for (int i = 1; i <= children.Extent(); ++i) {
                checkpoints.work();
                const int tag = maps[dim]->FindIndex(children(i));
                if (tag) assign((*classes[dim])[tag], parent);
            }
        }
    }
};
void fillInterior(const ShapeMap& entities, TopAbs_ShapeEnum kind,
                  const ShapeMap& pieces, const std::vector<int>& parents,
                  std::vector<Class>& classifications, Checkpoints& checkpoints) {
    std::vector<int> owner(entities.Extent() + 1, -1);
    for (int i = 1; i <= pieces.Extent(); ++i) {
        ShapeMap children;
        TopExp::MapShapes(pieces(i), kind, children);
        for (int j = 1; j <= children.Extent(); ++j) {
            checkpoints.work();
            const int tag = entities.FindIndex(children(j));
            require(tag > 0, OcctPartitionReason::InvalidResult);
            require(owner[tag] < 0 || owner[tag] == parents[i],
                    OcctPartitionReason::InterfaceOwnershipInvalid);
            owner[tag] = parents[i];
        }
    }
    for (int i = 1; i <= entities.Extent(); ++i) {
        require(owner[i] > 0, OcctPartitionReason::InvalidResult);
        require(classifications[i].tag != -2, OcctPartitionReason::ProvenanceLost);
        if (classifications[i].dim < 0) classifications[i] = {3, owner[i]};
    }
}
} // namespace

const char* occtPartitionReasonName(OcctPartitionReason reason) noexcept {
    switch (reason) {
    case OcctPartitionReason::None: return "none";
    case OcctPartitionReason::InvalidRequest: return "invalid_request";
    case OcctPartitionReason::UnsupportedPolicy: return "unsupported_policy";
    case OcctPartitionReason::BudgetExceeded: return "budget_exceeded";
    case OcctPartitionReason::Cancelled: return "cancelled";
    case OcctPartitionReason::AllocationFailed: return "allocation_failed";
    case OcctPartitionReason::KernelFailure: return "kernel_failure";
    case OcctPartitionReason::InvalidResult: return "invalid_result";
    case OcctPartitionReason::ProvenanceLost: return "provenance_lost";
    case OcctPartitionReason::InterfaceOwnershipInvalid: return "interface_ownership_invalid";
    case OcctPartitionReason::VolumeMismatch: return "volume_mismatch";
    case OcctPartitionReason::SourceMutation: return "source_mutation";
    }
    return "unknown";
}

OcctPartitionOutcome partitionOcctModel(const OcctModel& original,
                                       const OcctPartitionOptions& options) noexcept {
    OcctPartitionOutcome out;
    const auto started = Clock::now();
    const auto& source = OcctModelAccess::shape(original);
    std::string before;
    bool sourceSnapshotReady = false;
    Checkpoints checkpoints{options, out};
    try {
        require(std::isfinite(options.integrationRelativeError) &&
                    options.integrationRelativeError > 0 && options.integrationRelativeError < 1 &&
                    std::isfinite(options.maximumRelativeVolumeChange) &&
                    options.maximumRelativeVolumeChange >= 0 &&
                    options.maximumRelativeVolumeChange <= (options.meshConstructionCandidate?1e-4:1e-6) &&
                    options.maximumSnapshotAuditBytes > 0 && options.maximumPieces > 0 &&
                    options.maximumPlanes > 0 && options.maximumHistoryWork > 0,
                OcctPartitionReason::InvalidRequest);
        if(options.meshConstructionCandidate)
            require((options.policy==OcctPartitionPolicy::NativePrincipalFeatureAxes||
                     options.policy==OcctPartitionPolicy::NativeConcaveFeatureAxes)&&
                    !options.allowSampledTrimTolerance&&!options.verifyNativeCommon,
                    OcctPartitionReason::InvalidRequest);
        require(options.policy == OcctPartitionPolicy::CartesianMidplanes ||
                    options.policy == OcctPartitionPolicy::CoaxialAngular ||
                    options.policy == OcctPartitionPolicy::CartesianGrid ||
                    options.policy == OcctPartitionPolicy::NativeFeatureAxes ||
                    options.policy == OcctPartitionPolicy::NativeCircularFeatureAxes ||
                    options.policy == OcctPartitionPolicy::NativeConcaveFeatureAxes ||
                    options.policy == OcctPartitionPolicy::NativePrincipalFeatureAxes,
                OcctPartitionReason::InvalidRequest);
        if(options.allowSampledTrimTolerance)
            require(options.policy==OcctPartitionPolicy::NativeConcaveFeatureAxes&&
                    std::isfinite(options.maximumPhysicalTrimDrift)&&options.maximumPhysicalTrimDrift>0,
                    OcctPartitionReason::InvalidRequest);
        if (options.policy == OcctPartitionPolicy::CartesianGrid)
            require(options.cartesianDivisions >= 2 && options.cartesianDivisions <= 5,
                    OcctPartitionReason::InvalidRequest);
        if (options.policy == OcctPartitionPolicy::CoaxialAngular)
            require(options.angularPlanes >= 2 && options.angularPlanes <= 16,
                    OcctPartitionReason::InvalidRequest);
        if (options.policy == OcctPartitionPolicy::NativeFeatureAxes ||
            options.policy == OcctPartitionPolicy::NativeCircularFeatureAxes)
            require(options.featureAxisPlanes >= 2 && options.featureAxisPlanes <= 8,
                    OcctPartitionReason::InvalidRequest);
        if (options.policy == OcctPartitionPolicy::NativeConcaveFeatureAxes ||
            options.policy == OcctPartitionPolicy::NativePrincipalFeatureAxes)
            require(options.featureAxisPlanes >= 2 && options.featureAxisPlanes <= 16,
                    OcctPartitionReason::InvalidRequest);
        require(options.extraPlanes.size() <= options.maximumPlanes,
                OcctPartitionReason::BudgetExceeded);
        for (const auto& plane : options.extraPlanes) {
            checkpoints.work();
            const auto p = plane.origin, n = plane.normal;
            const double length = std::hypot(n.x,n.y,n.z);
            require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z) &&
                        std::isfinite(n.x) && std::isfinite(n.y) && std::isfinite(n.z) &&
                        std::isfinite(length) && length > 0,
                    OcctPartitionReason::InvalidRequest);
        }
        checkpoints.check();
        BRepBuilderAPI_Copy copier;
        NCollection_List<TopoDS_Shape> arguments, tools;
        ShapeMap sourceSolids, sourceFaces, sourceEdges, sourceVertices;
        TopoDS_Shape copied;
        std::vector<int> selected;
        TopoDS_Compound selectedSource;
        {
            Timer timer{out.preparationSeconds};
            before = serialize(source, options.maximumSnapshotAuditBytes);
            sourceSnapshotReady = true;
            checkpoints.check();
            TopExp::MapShapes(source, TopAbs_SOLID, sourceSolids);
            TopExp::MapShapes(source, TopAbs_FACE, sourceFaces);
            TopExp::MapShapes(source, TopAbs_EDGE, sourceEdges);
            TopExp::MapShapes(source, TopAbs_VERTEX, sourceVertices);
            selected = options.sourceBodyTags;
            if (selected.empty())
                for (int i = 1; i <= sourceSolids.Extent(); ++i) selected.push_back(i);
            std::sort(selected.begin(), selected.end());
            require(!selected.empty() && std::adjacent_find(selected.begin(), selected.end()) == selected.end(),
                    OcctPartitionReason::InvalidRequest);
            for (int tag : selected)
                require(tag > 0 && tag <= sourceSolids.Extent(), OcctPartitionReason::InvalidRequest);
            copier.Perform(source, true, false);
            checkpoints.check();
            copied = copier.Shape();
            require(!copied.IsNull(), OcctPartitionReason::KernelFailure);
            BRep_Builder builder;
            builder.MakeCompound(selectedSource);
            Vec3 low{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::infinity()};
            Vec3 high{-low.x, -low.y, -low.z};
            std::array<std::vector<double>,3> coordinates;
            for (int tag : selected) {
                checkpoints.work();
                const auto body = original.solid(tag);
                Vec3 a, b;
                require(body != nullptr, OcctPartitionReason::InvalidRequest);
                body->bbox(a, b);
                const double lo[]{a.x,a.y,a.z}, hi[]{b.x,b.y,b.z};
                require(std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z) &&
                            std::isfinite(b.x) && std::isfinite(b.y) && std::isfinite(b.z) &&
                            b.x > a.x && b.y > a.y && b.z > a.z,
                        OcctPartitionReason::InvalidResult);
                low.x = std::min(low.x,a.x); low.y = std::min(low.y,a.y); low.z = std::min(low.z,a.z);
                high.x = std::max(high.x,b.x); high.y = std::max(high.y,b.y); high.z = std::max(high.z,b.z);
                const unsigned divisions = options.policy == OcctPartitionPolicy::CartesianGrid ?
                    options.cartesianDivisions : 2;
                const bool useGrid = (options.policy != OcctPartitionPolicy::NativeFeatureAxes &&
                    options.policy != OcctPartitionPolicy::NativeCircularFeatureAxes &&
                    options.policy != OcctPartitionPolicy::NativeConcaveFeatureAxes &&
                    options.policy != OcctPartitionPolicy::NativePrincipalFeatureAxes) ||
                    options.includeCartesianMidplanes;
                for (int axis = 0; useGrid && axis != 3; ++axis)
                  for (unsigned subdivision = 1; subdivision != divisions; ++subdivision) {
                    const double coordinate = lo[axis] + (hi[axis] - lo[axis]) * subdivision/divisions;
                    // Reusing a nearly identical cutting plane only affects the
                    // subdivision choice; entity identity is never welded.
                    if (std::none_of(coordinates[axis].begin(), coordinates[axis].end(),
                        [&](double value) { return std::abs(value - coordinate) <= 1e-8; }))
                        coordinates[axis].push_back(coordinate);
                }
                const auto copy = copier.ModifiedShape(sourceSolids(tag));
                require(!copy.IsNull(), OcctPartitionReason::KernelFailure);
                arguments.Append(copy);
                // Adding a shape to a compound may change its native Free
                // flag. Even this bookkeeping compound must own only copies.
                builder.Add(selectedSource, copy);
                GProp_GProps properties;
                const double error = BRepGProp::VolumeProperties(sourceSolids(tag), properties,
                    options.integrationRelativeError, false, false);
                require(std::isfinite(properties.Mass()) && properties.Mass() > 0 &&
                            std::isfinite(error) && error >= 0,
                        OcctPartitionReason::InvalidResult);
                out.sourceVolume += properties.Mass();
                out.sourceIntegrationError = std::max(out.sourceIntegrationError, error);
            }
            const double extent = (high - low).norm() * 2;
            require(std::isfinite(extent) && extent > 0, OcctPartitionReason::InvalidResult);
            const gp_Pnt center(low.x+(high.x-low.x)*.5, low.y+(high.y-low.y)*.5,
                                low.z+(high.z-low.z)*.5);
            std::vector<std::pair<Vec3,double>> planeEquations;
            auto addPlane = [&](const gp_Pnt& origin, const gp_Dir& normal) {
                checkpoints.work();
                const Vec3 n{normal.X(),normal.Y(),normal.Z()};
                const double offset = n.dot({origin.X(),origin.Y(),origin.Z()});
                for (const auto& [prior,constant] : planeEquations) {
                    checkpoints.work();
                    const double dot = n.dot(prior);
                    if (std::abs(dot) > 1-1e-10 &&
                        std::abs(offset-(dot > 0 ? constant : -constant)) <= 1e-8)
                        return;
                }
                require(out.planes < options.maximumPlanes, OcctPartitionReason::BudgetExceeded);
                BRepBuilderAPI_MakeFace face(gp_Pln(origin,normal),-extent,extent,-extent,extent);
                require(face.IsDone(), OcctPartitionReason::KernelFailure);
                tools.Append(face.Face());
                planeEquations.emplace_back(n,offset);
                ++out.planes;
            };
            if (options.policy != OcctPartitionPolicy::CoaxialAngular) {
                for (int axis = 0; axis != 3; ++axis) {
                    std::sort(coordinates[axis].begin(), coordinates[axis].end());
                    for (double coordinate : coordinates[axis]) {
                        auto origin = center;
                        if (axis == 0) origin.SetX(coordinate);
                        if (axis == 1) origin.SetY(coordinate);
                        if (axis == 2) origin.SetZ(coordinate);
                        addPlane(origin, gp_Dir(axis==0?1:0,axis==1?1:0,axis==2?1:0));
                    }
                }
                if (options.policy == OcctPartitionPolicy::NativeFeatureAxes ||
            options.policy == OcctPartitionPolicy::NativeCircularFeatureAxes) {
                    const Vec3 midpoint{center.X(),center.Y(),center.Z()};
                    auto nativePlane = [&](Vec3 normal,Vec3 origin) {
                        // Center the finite cutting face on the nearest point
                        // of its infinite support plane to the selected bbox.
                        const double offset = normal.dot(origin);
                        const auto p = midpoint-normal*(normal.dot(midpoint)-offset);
                        addPlane(gp_Pnt(p.x,p.y,p.z),gp_Dir(normal.x,normal.y,normal.z));
                    };
                    for (int tag : selected)
                        for (const auto& use : original.faceUses(tag)) {
                            checkpoints.work();
                            const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                            require(descriptor != nullptr,OcctPartitionReason::InvalidResult);
                            if (descriptor->kind == SurfaceKind::Cylinder || descriptor->kind == SurfaceKind::Cone) {
                                const auto x = descriptor->xDirection;
                                const auto y = descriptor->axis.cross(x);
                                for (unsigned i = 0; i != options.featureAxisPlanes; ++i) {
                                    const double angle = std::numbers::pi*i/options.featureAxisPlanes;
                                    nativePlane(x*std::cos(angle)+y*std::sin(angle),descriptor->origin);
                                }
                            } else if (options.includePlanarSupportPlanes && descriptor->kind == SurfaceKind::Plane) {
                                const double halfspan = std::abs(descriptor->axis.x)*(high.x-low.x)*.5 +
                                    std::abs(descriptor->axis.y)*(high.y-low.y)*.5 +
                                    std::abs(descriptor->axis.z)*(high.z-low.z)*.5;
                                if (std::abs(descriptor->axis.dot(descriptor->origin-midpoint)) < halfspan-1e-7)
                                    nativePlane(descriptor->axis,descriptor->origin);
                            }
                        }
                    if (options.policy == OcctPartitionPolicy::NativeCircularFeatureAxes) {
                        // Preserve the proven policy order: source native
                        // face descriptors first, then circular edges by tag.
                        std::set<int> circularCandidates;
                        for (int tag : selected)
                            for (const auto& use : original.faceUses(tag))
                                for (const auto& wire : original.boundaryUses(use.faceTag))
                                    for (const auto& edge : wire.edges)
                                        circularCandidates.insert(edge.edgeTag);
                        for (int tag : circularCandidates) {
                            checkpoints.work();
                            const auto edge = OcctModelAccess::edge(original,tag);
                            if (BRep_Tool::Degenerated(edge)) continue;
                            BRepAdaptor_Curve curve(edge);
                            if (curve.GetType() != GeomAbs_Circle) continue;
                            const auto circle = curve.Circle();
                            const auto position = circle.Position();
                            const Vec3 origin{circle.Location().X(),circle.Location().Y(),circle.Location().Z()};
                            const auto nativeAxis = position.Direction(), nativeX = position.XDirection();
                            const Vec3 axis{nativeAxis.X(),nativeAxis.Y(),nativeAxis.Z()};
                            const Vec3 x{nativeX.X(),nativeX.Y(),nativeX.Z()};
                            const auto y = axis.cross(x);
                            for (unsigned i = 0; i != options.featureAxisPlanes; ++i) {
                                const double angle = std::numbers::pi*i/options.featureAxisPlanes;
                                nativePlane(x*std::cos(angle)+y*std::sin(angle),origin);
                            }
                            const double halfspan = std::abs(axis.x)*(high.x-low.x)*.5 +
                                std::abs(axis.y)*(high.y-low.y)*.5 + std::abs(axis.z)*(high.z-low.z)*.5;
                            if (std::abs(axis.dot(origin-midpoint)) < halfspan-1e-7)
                                nativePlane(axis,origin);
                        }
                    }
                }
                if (options.policy == OcctPartitionPolicy::NativeConcaveFeatureAxes) {
                    // Work from each selected body's actual native faces. A
                    // strict area maximum retains the first native face use
                    // on ties; no fitted axes or source-specific rules enter.
                    for (int tag : selected) {
                        const OcctSurfaceDescriptor* reference = nullptr;
                        double largestArea = -1;
                        for (const auto& use : original.faceUses(tag)) {
                            checkpoints.work();
                            const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                            require(descriptor != nullptr, OcctPartitionReason::InvalidResult);
                            if (descriptor->kind != SurfaceKind::Cylinder) continue;
                            GProp_GProps properties;
                            BRepGProp::SurfaceProperties(sourceFaces(use.faceTag), properties);
                            const double area = properties.Mass();
                            require(std::isfinite(area) && area > 0, OcctPartitionReason::InvalidResult);
                            if (area > largestArea) { largestArea = area; reference = descriptor; }
                        }
                        require(reference != nullptr, OcctPartitionReason::UnsupportedPolicy);
                        Vec3 bodyLow, bodyHigh;
                        original.solid(tag)->bbox(bodyLow, bodyHigh);
                        const auto midpoint = (bodyLow + bodyHigh)*.5;
                        auto nativePlane = [&](Vec3 normal, Vec3 origin) {
                            const double halfspan = std::abs(normal.x)*(bodyHigh.x-bodyLow.x)*.5 +
                                std::abs(normal.y)*(bodyHigh.y-bodyLow.y)*.5 +
                                std::abs(normal.z)*(bodyHigh.z-bodyLow.z)*.5;
                            if (std::abs(normal.dot(origin-midpoint)) >= halfspan-1e-7) return;
                            const double offset = normal.dot(origin);
                            const auto p = midpoint-normal*(normal.dot(midpoint)-offset);
                            addPlane(gp_Pnt(p.x,p.y,p.z),gp_Dir(normal.x,normal.y,normal.z));
                        };
                        const auto x = reference->xDirection;
                        const auto y = reference->axis.cross(x);
                        for (unsigned i = 0; i != options.featureAxisPlanes; ++i) {
                            const double angle = std::numbers::pi*i/options.featureAxisPlanes;
                            nativePlane(x*std::cos(angle)+y*std::sin(angle),reference->origin);
                        }
                        for (const auto& use : original.faceUses(tag)) {
                            checkpoints.work();
                            const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                            require(descriptor != nullptr, OcctPartitionReason::InvalidResult);
                            if ((descriptor->kind == SurfaceKind::Cylinder || descriptor->kind == SurfaceKind::Cone) &&
                                descriptor->orientation == CadOrientation::Reversed) {
                                nativePlane(descriptor->xDirection,descriptor->origin);
                                nativePlane(descriptor->axis.cross(descriptor->xDirection),descriptor->origin);
                            }
                            if (options.includePlanarSupportPlanes && descriptor->kind == SurfaceKind::Plane &&
                                std::abs(descriptor->axis.dot(reference->axis)) > 1-1e-10)
                                nativePlane(descriptor->axis,descriptor->origin);
                        }
                    }
                }
                if (options.policy == OcctPartitionPolicy::NativePrincipalFeatureAxes) {
                    // Retain the independently proven proposal order: native
                    // principal radial planes, then original planar face uses.
                    // Proposal normalization/deduplication is separate from
                    // native entity identity; no coordinates are welded.
                    std::vector<OcctPartitionPlane> proposals;
                    for (int tag : selected) {
                        const OcctSurfaceDescriptor* reference = nullptr;
                        double largestArea = -1;
                        for (const auto& use : original.faceUses(tag)) {
                            checkpoints.work();
                            const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                            require(descriptor != nullptr, OcctPartitionReason::InvalidResult);
                            if (descriptor->kind != SurfaceKind::Cylinder) continue;
                            GProp_GProps properties;
                            BRepGProp::SurfaceProperties(sourceFaces(use.faceTag), properties);
                            const auto area = properties.Mass();
                            require(std::isfinite(area) && area > 0, OcctPartitionReason::InvalidResult);
                            if (area > largestArea) { largestArea = area; reference = descriptor; }
                        }
                        require(reference != nullptr, OcctPartitionReason::UnsupportedPolicy);
                        Vec3 bodyLow, bodyHigh;
                        original.solid(tag)->bbox(bodyLow,bodyHigh);
                        const auto midpoint = (bodyLow+bodyHigh)*.5;
                        auto propose = [&](Vec3 origin,Vec3 normal) {
                            checkpoints.work();
                            normal = normal/normal.norm();
                            const auto offset = normal.dot(origin);
                            const auto halfspan = std::abs(normal.x)*(bodyHigh.x-bodyLow.x)*.5 +
                                std::abs(normal.y)*(bodyHigh.y-bodyLow.y)*.5 +
                                std::abs(normal.z)*(bodyHigh.z-bodyLow.z)*.5;
                            if (std::abs(normal.dot(origin-midpoint)) >= halfspan-1e-7) return;
                            for (const auto& previous : proposals) {
                                checkpoints.work();
                                const auto dot = previous.normal.dot(normal);
                                if (std::abs(dot)>1-1e-10 && std::abs(offset-
                                    (dot>0 ? previous.normal.dot(previous.origin) :
                                             -previous.normal.dot(previous.origin)))<=1e-8) return;
                            }
                            require(proposals.size()<options.maximumPlanes,OcctPartitionReason::BudgetExceeded);
                            proposals.push_back({origin,normal});
                        };
                        const auto x = reference->xDirection;
                        const auto y = reference->axis.cross(x);
                        for (unsigned i=0;i!=options.featureAxisPlanes;++i) {
                            const auto angle = std::numbers::pi*i/options.featureAxisPlanes;
                            propose(reference->origin,x*std::cos(angle)+y*std::sin(angle));
                        }
                        if (options.includePlanarSupportPlanes)
                            for (const auto& use : original.faceUses(tag)) {
                                checkpoints.work();
                                const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                                require(descriptor != nullptr,OcctPartitionReason::InvalidResult);
                                if (descriptor->kind==SurfaceKind::Plane)
                                    propose(descriptor->origin,descriptor->axis);
                            }
                    }
                    // Use precisely the same renormalization and selected-
                    // bbox centering as caller-provided native support cuts.
                    for (const auto& plane : proposals) {
                        const auto length = std::hypot(plane.normal.x,plane.normal.y,plane.normal.z);
                        const auto n = plane.normal/length;
                        const Vec3 midpoint{center.X(),center.Y(),center.Z()};
                        const auto p = midpoint-n*(n.dot(midpoint)-n.dot(plane.origin));
                        require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z),
                                OcctPartitionReason::InvalidRequest);
                        addPlane(gp_Pnt(p.x,p.y,p.z),gp_Dir(n.x,n.y,n.z));
                    }
                }
            } else {
                // This policy is deliberately limited to one body whose
                // cylindrical/conical faces share one axis and whose planar
                // faces are normal to it. The splitter still preserves the
                // actual native body, including hollow profiles and trims.
                require(selected.size() == 1, OcctPartitionReason::UnsupportedPolicy);
                const OcctSurfaceDescriptor* reference = nullptr;
                for (const auto& use : original.faceUses(selected.front())) {
                    const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                    require(descriptor != nullptr, OcctPartitionReason::InvalidResult);
                    if (descriptor->kind == SurfaceKind::Cylinder || descriptor->kind == SurfaceKind::Cone) {
                        if (!reference) reference = descriptor;
                        else require(std::abs(descriptor->axis.dot(reference->axis)) >= 1-1e-10 &&
                                (descriptor->origin-reference->origin).cross(reference->axis).norm() <= 1e-7,
                            OcctPartitionReason::UnsupportedPolicy);
                    } else require(descriptor->kind == SurfaceKind::Plane,
                                   OcctPartitionReason::UnsupportedPolicy);
                }
                require(reference != nullptr, OcctPartitionReason::UnsupportedPolicy);
                for (const auto& use : original.faceUses(selected.front())) {
                    const auto* descriptor = original.surfaceDescriptor(use.faceTag);
                    if (descriptor->kind == SurfaceKind::Plane)
                        require(std::abs(descriptor->axis.dot(reference->axis)) >= 1-1e-10,
                                OcctPartitionReason::UnsupportedPolicy);
                }
                const auto x = reference->xDirection;
                const auto y = reference->axis.cross(x);
                const auto origin = reference->origin;
                for (unsigned i = 0; i != options.angularPlanes; ++i) {
                    const double angle = std::numbers::pi * i / options.angularPlanes;
                    const auto normal = x*std::cos(angle) + y*std::sin(angle);
                    addPlane(gp_Pnt(origin.x,origin.y,origin.z),gp_Dir(normal.x,normal.y,normal.z));
                }
            }
            for (const auto& plane : options.extraPlanes) {
                const double length = std::hypot(plane.normal.x,plane.normal.y,plane.normal.z);
                const auto n = plane.normal / length;
                const Vec3 midpoint{center.X(),center.Y(),center.Z()};
                const auto p = midpoint-n*(n.dot(midpoint)-n.dot(plane.origin));
                require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z),
                        OcctPartitionReason::InvalidRequest);
                addPlane(gp_Pnt(p.x,p.y,p.z),gp_Dir(n.x,n.y,n.z));
            }
            require(!tools.IsEmpty(),OcctPartitionReason::UnsupportedPolicy);
        }
        BRepAlgoAPI_Splitter splitter;
        {
            Timer timer{out.splittingSeconds};
            checkpoints.check();
            splitter.SetArguments(arguments);
            splitter.SetTools(tools);
            splitter.SetNonDestructive(true);
            splitter.SetRunParallel(false);
            splitter.SetToFillHistory(true);
            splitter.Build();
            checkpoints.check();
            require(splitter.IsDone() && !splitter.HasErrors(), OcctPartitionReason::KernelFailure);
        }
        {
            Timer timer{out.validationSeconds};
            const auto result = splitter.Shape();
            require(!result.IsNull(), OcctPartitionReason::InvalidResult);
            ShapeMap pieces, faces, edges, vertices;
            TopExp::MapShapes(result, TopAbs_SOLID, pieces);
            TopExp::MapShapes(result, TopAbs_FACE, faces);
            TopExp::MapShapes(result, TopAbs_EDGE, edges);
            TopExp::MapShapes(result, TopAbs_VERTEX, vertices);
            out.pieces = pieces.Extent();
            require(out.pieces > 0, OcctPartitionReason::InvalidResult);
            require(out.pieces <= options.maximumPieces && out.pieces <= options.snapshotLimits.maximumSolids &&
                        static_cast<std::uint64_t>(faces.Extent()) <= options.snapshotLimits.maximumFaces &&
                        static_cast<std::uint64_t>(edges.Extent()) <= options.snapshotLimits.maximumEdges &&
                        static_cast<std::uint64_t>(vertices.Extent()) <= options.snapshotLimits.maximumVertices,
                    OcctPartitionReason::BudgetExceeded);
            out.parentBodyByPiece.resize(pieces.Extent()+1,-1);
            out.parentByFace.resize(faces.Extent()+1);
            out.parentByEdge.resize(edges.Extent()+1);
            out.parentByVertex.resize(vertices.Extent()+1);
            for (int originalTag : selected) {
                const auto copy = copier.ModifiedShape(sourceSolids(originalTag));
                auto assignPieces = [&](const TopoDS_Shape& shape) {
                    ShapeMap children;
                    TopExp::MapShapes(shape, TopAbs_SOLID, children);
                    for (int i = 1; i <= children.Extent(); ++i) {
                        checkpoints.work();
                        const int tag = pieces.FindIndex(children(i));
                        if (!tag) continue;
                        require(out.parentBodyByPiece[tag] < 0 || out.parentBodyByPiece[tag] == originalTag,
                                OcctPartitionReason::ProvenanceLost);
                        out.parentBodyByPiece[tag] = originalTag;
                    }
                };
                assignPieces(copy);
                for (const auto& modified : splitter.Modified(copy)) assignPieces(modified);
            }
            Provenance provenance{vertices,edges,faces,out,checkpoints};
            const std::array<const ShapeMap*,3> originalMaps{&sourceVertices,&sourceEdges,&sourceFaces};
            for (int dimension = 0; dimension != 3; ++dimension)
                for (int originalTag = 1; originalTag <= originalMaps[dimension]->Extent(); ++originalTag) {
                    checkpoints.work();
                    const auto copy = copier.ModifiedShape((*originalMaps[dimension])(originalTag));
                    if (copy.IsNull()) continue;
                    const Class parent{dimension,originalTag};
                    provenance.apply(copy,parent);
                    for (const auto& modified : splitter.Modified(copy)) provenance.apply(modified,parent);
                    // An original face generates intersection edges/vertices;
                    // it cannot authorize a generated material-interior face.
                    for (const auto& generated : splitter.Generated(copy))
                        provenance.apply(generated,parent,false);
                }
            // Support planes through existing corners can replace vertices
            // without public Modified/Generated entries. The Boolean data
            // structure retains their exact same-domain native representative.
            // This restores source identity; it never compares coordinates or
            // joins distinct result vertices. Multiple original origins remain
            // ambiguous and are refused by the normal provenance gate.
            require(splitter.DSFiller() != nullptr,OcctPartitionReason::KernelFailure);
            const auto& booleanData = splitter.DSFiller()->DS();
            for (int originalTag = 1; originalTag <= sourceVertices.Extent(); ++originalTag) {
                checkpoints.work();
                const auto copy = copier.ModifiedShape(sourceVertices(originalTag));
                if (copy.IsNull()) continue;
                const int sourceIndex = booleanData.Index(copy);
                if (sourceIndex < 0) continue; // The original body was not selected.
                const int representative = booleanData.GetSameDomainIndex(sourceIndex);
                require(representative >= 0 && representative < booleanData.NbShapes(),
                        OcctPartitionReason::KernelFailure);
                const auto& vertex = booleanData.Shape(representative);
                require(vertex.ShapeType() == TopAbs_VERTEX,OcctPartitionReason::KernelFailure);
                const int tag = vertices.FindIndex(vertex);
                if (tag) assign(out.parentByVertex[tag],{0,originalTag});
            }
            for (int i = 1; i <= pieces.Extent(); ++i) {
                checkpoints.work();
                require(out.parentBodyByPiece[i] > 0, OcctPartitionReason::ProvenanceLost);
                const auto piece = TopoDS::Solid(pieces(i));
                require(BRepCheck_Analyzer(piece).IsValid(), OcctPartitionReason::InvalidResult);
                paacs::geo::detail::OcctSolidClassifier classifier(piece);
                classifier.PerformInfinitePoint(1e-7);
                require(classifier.State() == TopAbs_OUT, OcctPartitionReason::InvalidResult);
                GProp_GProps properties;
                const double error = BRepGProp::VolumeProperties(piece, properties,
                    options.integrationRelativeError, false, false);
                require(std::isfinite(properties.Mass()) && properties.Mass() > 0 &&
                            std::isfinite(error) && error >= 0,
                        OcctPartitionReason::InvalidResult);
                out.partitionVolume += properties.Mass();
                out.maximumPieceIntegrationError = std::max(out.maximumPieceIntegrationError,error);
            }
            out.relativeVolumeChange = std::abs(out.partitionVolume-out.sourceVolume)/out.sourceVolume;
            require(std::isfinite(out.relativeVolumeChange),OcctPartitionReason::VolumeMismatch);
            if(out.relativeVolumeChange>options.maximumRelativeVolumeChange&&
               (options.allowSampledTrimTolerance||
                (options.retrySpanAwareVolumeReference&&out.relativeVolumeChange<=1e-4&&pieces.Extent()<=128)))
                retrySpanAwareVolume(sourceSolids,pieces,selected,checkpoints);
            // Preserve the original early refusal and work receipts for every
            // default policy. Only explicit new-policy attempts reach the
            // independent native trimming witness below.
            require(out.relativeVolumeChange<=options.maximumRelativeVolumeChange||
                    options.allowSampledTrimTolerance,OcctPartitionReason::VolumeMismatch);
            std::vector<std::vector<std::pair<int,int>>> faceOwners(faces.Extent()+1);
            for (int i = 1; i <= pieces.Extent(); ++i)
                for (TopExp_Explorer use(pieces(i),TopAbs_FACE); use.More(); use.Next()) {
                    checkpoints.work();
                    const int tag = faces.FindIndex(use.Current());
                    require(tag > 0, OcctPartitionReason::InvalidResult);
                    const auto orientation = use.Current().Orientation();
                    require(orientation == TopAbs_FORWARD || orientation == TopAbs_REVERSED,
                            OcctPartitionReason::InterfaceOwnershipInvalid);
                    const int sense = orientation == faces(tag).Orientation() ? 1 : -1;
                    faceOwners[tag].emplace_back(i,sense);
                }
            for (int tag = 1; tag <= faces.Extent(); ++tag) {
                const auto& owners = faceOwners[tag];
                const auto parent = out.parentByFace[tag];
                if (owners.size() == 1) {
                    require(parent.dim == 2 && parent.tag > 0, OcctPartitionReason::ProvenanceLost);
                    ++out.mappedExteriorFaces;
                } else {
                    require(owners.size() == 2 && owners[0].first != owners[1].first &&
                                owners[0].second == -owners[1].second && parent.dim < 0 &&
                                out.parentBodyByPiece[owners[0].first] == out.parentBodyByPiece[owners[1].first],
                            OcctPartitionReason::InterfaceOwnershipInvalid);
                    out.interfaces.push_back({tag,{owners[0].first,owners[1].first},
                        {owners[0].second,owners[1].second},out.parentBodyByPiece[owners[0].first]});
                    ++out.sharedCutFaces;
                }
            }
            if(out.relativeVolumeChange>options.maximumRelativeVolumeChange) {
                require(options.allowSampledTrimTolerance&&out.spanAwareVolumeReferenceUsed,
                        OcctPartitionReason::VolumeMismatch);
                OcctTrimToleranceInput witness{source,copier,sourceSolids,sourceFaces,sourceEdges,sourceVertices,
                    faces,edges,vertices,selected,out.parentBodyByPiece,out.parentByFace,out.parentByEdge,
                    out.parentByVertex,faceOwners,out.sourceVolume,out.partitionVolume};
                OcctTrimToleranceOptions limits;
                limits.maximumPhysicalDrift=options.maximumPhysicalTrimDrift;
                limits.maximumWork=options.maximumHistoryWork-out.historyWork;
                limits.maximumSnapshotBytes=options.maximumSnapshotAuditBytes;
                limits.cancellation=options.cancellation;
                out.trimTolerance=auditOcctTrimTolerance(witness,limits);
                require(out.trimTolerance.work<=options.maximumHistoryWork-out.historyWork,
                        OcctPartitionReason::BudgetExceeded);
                out.historyWork+=out.trimTolerance.work;checkpoints.check();
                if(!out.trimTolerance.admitted) {
                    switch(out.trimTolerance.reason) {
                    case OcctTrimToleranceReason::BudgetExceeded:require(false,OcctPartitionReason::BudgetExceeded);break;
                    case OcctTrimToleranceReason::Cancelled:require(false,OcctPartitionReason::Cancelled);break;
                    case OcctTrimToleranceReason::AllocationFailed:require(false,OcctPartitionReason::AllocationFailed);break;
                    case OcctTrimToleranceReason::SourceMutation:require(false,OcctPartitionReason::SourceMutation);break;
                    case OcctTrimToleranceReason::SupportMismatch:
                    case OcctTrimToleranceReason::SourceCoverageLost:require(false,OcctPartitionReason::ProvenanceLost);break;
                    case OcctTrimToleranceReason::NativeQueryFailed:require(false,OcctPartitionReason::KernelFailure);break;
                    default:require(false,OcctPartitionReason::VolumeMismatch);
                    }
                }
            }
            fillInterior(vertices,TopAbs_VERTEX,pieces,out.parentBodyByPiece,out.parentByVertex,checkpoints);
            fillInterior(edges,TopAbs_EDGE,pieces,out.parentBodyByPiece,out.parentByEdge,checkpoints);
            fillInterior(faces,TopAbs_FACE,pieces,out.parentBodyByPiece,out.parentByFace,checkpoints);
            if (options.verifyNativeCommon) {
                checkpoints.check();
                BRepAlgoAPI_Common common;
                NCollection_List<TopoDS_Shape> commonArguments, commonTools;
                // The arguments are deep copies, so even a future kernel-side
                // defect cannot alter the source snapshot.
                BRepBuilderAPI_Copy sourceCopy(selectedSource,true,false);
                commonArguments.Append(sourceCopy.Shape());
                commonTools.Append(result);
                common.SetArguments(commonArguments);
                common.SetTools(commonTools);
                common.SetNonDestructive(true);
                common.SetRunParallel(false);
                common.Build();
                checkpoints.check();
                require(common.IsDone() && !common.HasErrors(), OcctPartitionReason::KernelFailure);
                GProp_GProps properties;
                out.nativeCommonIntegrationError = BRepGProp::VolumeProperties(common.Shape(),properties,
                    options.integrationRelativeError,false,false);
                out.nativeCommonVolume = properties.Mass();
                out.nativeCommonChecked = true;
                require(std::isfinite(out.nativeCommonVolume) &&
                            std::isfinite(out.nativeCommonIntegrationError) && out.nativeCommonIntegrationError >= 0 &&
                            std::abs(out.nativeCommonVolume-out.sourceVolume)/out.sourceVolume <=
                                options.maximumRelativeVolumeChange,
                        OcctPartitionReason::VolumeMismatch);
            }
            checkpoints.check();
            auto snapshot = OcctModelAccess::fromShape(result, options.snapshotLimits);
            require(snapshot.imported(), snapshot.reason == CadImportReason::TopologyBudget ?
                OcctPartitionReason::BudgetExceeded : snapshot.reason == CadImportReason::AllocationFailed ?
                OcctPartitionReason::AllocationFailed : OcctPartitionReason::InvalidResult);
            checkpoints.check();
            require(snapshot.solids == out.pieces && snapshot.faces == static_cast<std::uint64_t>(faces.Extent()) &&
                        snapshot.edges == static_cast<std::uint64_t>(edges.Extent()) &&
                        snapshot.vertices == static_cast<std::uint64_t>(vertices.Extent()),
                    OcctPartitionReason::InvalidResult);
            out.sourceImmutabilityVerified = before == serialize(source,options.maximumSnapshotAuditBytes);
            require(out.sourceImmutabilityVerified, OcctPartitionReason::SourceMutation);
            checkpoints.check();
            out.model = std::move(snapshot.model);
            out.reason = OcctPartitionReason::None;
        }
    } catch (const Stop& stop) {
        out.reason = stop.reason;
    } catch (const Standard_OutOfMemory&) {
        out.reason = OcctPartitionReason::AllocationFailed;
    } catch (const std::bad_alloc&) {
        out.reason = OcctPartitionReason::AllocationFailed;
    } catch (const Standard_Failure&) {
        out.reason = OcctPartitionReason::KernelFailure;
    } catch (...) {
        out.reason = OcctPartitionReason::KernelFailure;
    }
    if (!out.model && sourceSnapshotReady) {
        try {
            out.sourceImmutabilityVerified = before == serialize(source,options.maximumSnapshotAuditBytes);
            if (!out.sourceImmutabilityVerified) out.reason = OcctPartitionReason::SourceMutation;
        } catch (...) { /* Keep the original refusal when the audit cannot allocate. */ }
    }
    out.totalSeconds = seconds(started);
    return out;
}
const char* occtPartitionPolicyName(OcctPartitionPolicy policy) noexcept {
    switch (policy) {
        case OcctPartitionPolicy::CartesianMidplanes: return "cartesian_midplanes";
        case OcctPartitionPolicy::CoaxialAngular: return "coaxial_angular";
        case OcctPartitionPolicy::CartesianGrid: return "cartesian_grid";
        case OcctPartitionPolicy::NativeFeatureAxes: return "native_feature_axes";
        case OcctPartitionPolicy::NativeCircularFeatureAxes: return "native_circular_feature_axes";
        case OcctPartitionPolicy::NativeConcaveFeatureAxes: return "native_concave_feature_axes";
        case OcctPartitionPolicy::NativePrincipalFeatureAxes: return "native_principal_feature_axes";
    }
    return "unknown";
}
} // namespace paacs::geo::detail
