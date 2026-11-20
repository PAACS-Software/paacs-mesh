// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctBoxContact.hpp"
#include "OcctModelInternal.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace paacs;
namespace {
void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
std::string snapshot(const geo::OcctModel& model) {
    std::ostringstream stream;
    BRepTools::Write(geo::detail::OcctModelAccess::shape(model), stream);
    return stream.str();
}
std::vector<std::uint8_t> bytes(const mesh::Mesh& value) {
    std::vector<std::uint8_t> out;
    require(sdk::serializeMeshArchive(value, out).completed, "mesh witness serializes");
    return out;
}
bool same(Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
bool sameTopology(const mesh::Mesh& a, const mesh::Mesh& b) {
    if (a.nodeCount() != b.nodeCount() || a.cellCount() != b.cellCount()) return false;
    for (Index i = 0; i < a.nodeCount(); ++i) if (a.nodeClass(i) != b.nodeClass(i)) return false;
    for (Index i = 0; i < a.cellCount(); ++i) {
        if (a.cell(i).type != b.cell(i).type || a.cell(i).onEntity != b.cell(i).onEntity) return false;
        for (int j = 0; j < mesh::nodesPerCell(a.cell(i).type); ++j) if (a.cellNodes(i)[j] != b.cellNodes(i)[j]) return false;
    }
    return true;
}
geo::NativeCadImportOutcome pair(Vec3 at, const TopoDS_Shape& first = BRepPrimAPI_MakeBox(1,1,1).Shape()) {
    TopoDS_Compound compound;
    BRep_Builder builder;
    builder.MakeCompound(compound);
    builder.Add(compound, first);
    builder.Add(compound, BRepPrimAPI_MakeBox(gp_Pnt(at.x,at.y,at.z),1,1,1).Shape());
    return geo::detail::OcctModelAccess::fromShape(compound);
}
struct View final : geo::Model {
    const geo::OcctModel& model; const geo::GSolid* body;
    std::vector<const geo::GFace*> fs; std::vector<const geo::GEdge*> es; std::vector<const geo::GVertex*> vs;
    View(const geo::OcctModel& m, int tag) : model(m), body(m.solid(tag)) {
        std::set<int> edges, vertices;
        for (int face : body->boundaryFaceTags()) {
            fs.push_back(m.face(face));
            for (const auto& loop:m.boundaryUses(face)) for (const auto& use:loop.edges) {
                edges.insert(use.edgeTag); vertices.insert(use.startVertexTag); vertices.insert(use.endVertexTag);
            }
        }
        for (int edge:edges) es.push_back(m.edge(edge));
        for (int vertex:vertices) vs.push_back(m.vertex(vertex));
    }
    std::vector<const geo::GVertex*> vertices() const override { return vs; }
    std::vector<const geo::GEdge*> edges() const override { return es; }
    std::vector<const geo::GFace*> faces() const override { return fs; }
    std::vector<const geo::GSolid*> solids() const override { return {body}; }
    const geo::GVertex* vertex(int tag) const override { return model.vertex(tag); }
    const geo::GEdge* edge(int tag) const override { return model.edge(tag); }
    const geo::GFace* face(int tag) const override { return model.face(tag); }
    const geo::GSolid* solid(int tag) const override { return model.solid(tag); }
};
mesh::Mesh coarse(const geo::OcctModel& model, int tag, double h = 5) {
    View view(model, tag);
    algo::GeometryMeshingOptions options;
    options.targetSize = h;
    options.bodyControls.push_back({});
    options.bodyControls.back().solidTag = tag;
    options.bodyControls.back().elements = algo::GeometryElements::Tet4;
    mesh::Mesh result;
    require(algo::meshGeometryModel(view, options, result).published, "actual planar-core box meshing");
    return result;
}
}
int main() {
    try {
        using namespace geo::detail;
        auto full = pair({1,0,0});
        require(full.imported(), "full face native import");
        const auto sourceBefore = snapshot(*full.model);
        OcctBoxContactOptions options;
        const auto admitted = admitOcctBoxFaceContact(*full.model,1,2,options);
        require(admitted.admitted && admitted.reason == OcctBoxContactReason::None && admitted.sourceImmutabilityVerified &&
                std::abs(admitted.nativeSkinArea-1)<1e-12 && admitted.expectedFullFaceArea==1,
                "complete native full-face contact admitted");
        auto body = coarse(*full.model,1);
        auto original = body;
        body.setNode(0,body.node(0)+Vec3{1e-15,0,0});
        const auto evaluated = evaluateOcctBoxContactCorners(*full.model,1,options,body);
        require(evaluated.published && !evaluated.destinationUnchanged && evaluated.sourceImmutabilityVerified &&
                evaluated.evaluatedVertices==8 && evaluated.nativeInclusionSamples==14 && evaluated.audit.valid && evaluated.quality.valid &&
                evaluated.maximumEndpointCorrection>0 && sameTopology(original,body),
                "owned native endpoint evaluation preserves all topology and classes");
        for (Index i=0;i<body.nodeCount();++i)
            require(same(body.node(i),full.model->vertex(body.nodeClass(i).tag)->position()), "all coordinates exact original native vertex evaluations");
        require(snapshot(*full.model)==sourceBefore, "native source remains byte identical");

        auto partial = pair({1,.5,0});
        const auto partialBefore=snapshot(*partial.model);
        const auto partialResult=admitOcctBoxFaceContact(*partial.model,1,2,options);
        require(!partialResult.admitted && partialResult.reason==OcctBoxContactReason::PartialFaceContact &&
                partialResult.sourceImmutabilityVerified && snapshot(*partial.model)==partialBefore, "partial faces refused without changing source");
        auto overlapping = pair({.5,0,0});
        const auto overlapResult=admitOcctBoxFaceContact(*overlapping.model,1,2,options);
        require(!overlapResult.admitted && overlapResult.reason==OcctBoxContactReason::MaterialOverlap && overlapResult.sourceImmutabilityVerified,
                "positive material overlap refused");
        auto edgeOnly=pair({1,1,0});
        require(admitOcctBoxFaceContact(*edgeOnly.model,1,2,options).reason==OcctBoxContactReason::PartialFaceContact,
                "edge-only contact cannot be claimed as full face");
        auto curved=pair({2,0,0},BRepPrimAPI_MakeCylinder(.5,1).Shape());
        require(admitOcctBoxFaceContact(*curved.model,1,2,options).reason==OcctBoxContactReason::UnsupportedBody,
                "curved source not misidentified as box");
        gp_Trsf rotation;rotation.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),.125);
        auto rotated=pair({2,0,0},BRepBuilderAPI_Transform(BRepPrimAPI_MakeBox(1,1,1).Shape(),rotation,true).Shape());
        require(admitOcctBoxFaceContact(*rotated.model,1,2,options).reason==OcctBoxContactReason::UnsupportedBody,
                "non-axis-aligned boxes remain unsupported");

        auto fine=coarse(*full.model,1,.2);
        const auto fineBefore=bytes(fine);
        const auto fineResult=evaluateOcctBoxContactCorners(*full.model,1,options,fine);
        require(!fineResult.published && fineResult.destinationUnchanged && fineResult.reason==OcctBoxContactReason::UnsupportedResolution && bytes(fine)==fineBefore,
                "fine request is refused transactionally");
        auto wrongBody=coarse(*full.model,2);
        const auto wrongBefore=bytes(wrongBody);
        const auto wrong=evaluateOcctBoxContactCorners(*full.model,1,options,wrongBody);
        require(!wrong.published && wrong.reason==OcctBoxContactReason::InvalidMesh && wrong.destinationUnchanged && bytes(wrongBody)==wrongBefore,
                "foreign original vertex identities refused");
        auto moved=original;moved.setNode(0,moved.node(0)+Vec3{1e-3,0,0});
        const auto movedBefore=bytes(moved);
        const auto movedResult=evaluateOcctBoxContactCorners(*full.model,1,options,moved);
        require(!movedResult.published && movedResult.reason==OcctBoxContactReason::GeometryError && bytes(moved)==movedBefore,
                "out-of-tolerance endpoint cannot be repaired by lookup");
        auto inverted=original;
        mesh::Mesh bad;
        for(Index i=0;i<inverted.nodeCount();++i)bad.addNode(inverted.node(i),inverted.nodeClass(i));
        bool reversed=false;
        for(Index i=0;i<inverted.cellCount();++i){const auto& c=inverted.cell(i);const auto* in=inverted.cellNodes(i);std::array<Index,4> ids{};int count=mesh::nodesPerCell(c.type);for(int j=0;j<count;++j)ids[j]=in[j];if(!reversed&&c.type==mesh::CellType::Tet4){std::swap(ids[0],ids[1]);reversed=true;}bad.addCell(c.type,{ids.data(),static_cast<std::size_t>(count)},c.onEntity);}
        const auto badBefore=bytes(bad);
        const auto badResult=evaluateOcctBoxContactCorners(*full.model,1,options,bad);
        require(!badResult.published && badResult.reason==OcctBoxContactReason::InvalidMesh && bytes(bad)==badBefore, "exact inverted Tet refused");
        for(std::uint64_t bound:{0ULL,1ULL,10ULL}){
            auto limited=options;limited.maximumWork=bound;
            require(admitOcctBoxFaceContact(*full.model,1,2,limited).reason==OcctBoxContactReason::BudgetExceeded,"native pair work budget charged");
            auto candidate=original;const auto before=bytes(candidate);
            require(evaluateOcctBoxContactCorners(*full.model,1,limited,candidate).reason==OcctBoxContactReason::BudgetExceeded&&bytes(candidate)==before,
                    "corner work exhaustion preserves prior mesh");
        }
        auto snapshotLimit=options;snapshotLimit.maximumSnapshotBytes=1;
        require(admitOcctBoxFaceContact(*full.model,1,2,snapshotLimit).reason==OcctBoxContactReason::BudgetExceeded,"native snapshot limit classified");
        auto workingLimit=options;workingLimit.maximumWorkingBytes=1;
        require(admitOcctBoxFaceContact(*full.model,1,2,workingLimit).reason==OcctBoxContactReason::BudgetExceeded,"logical working limit classified");
        std::atomic<bool> cancelled{true};auto stopped=options;stopped.cancellation=&cancelled;
        require(admitOcctBoxFaceContact(*full.model,1,2,stopped).reason==OcctBoxContactReason::Cancelled,"pair cancellation distinct");
        auto candidate=original;const auto candidateBefore=bytes(candidate);
        require(evaluateOcctBoxContactCorners(*full.model,1,stopped,candidate).reason==OcctBoxContactReason::Cancelled&&bytes(candidate)==candidateBefore,"corner cancellation preserves prior mesh");
        require(snapshot(*full.model)==sourceBefore,"all failures preserve source snapshot");
        std::cout<<"Native box full-face admission and transactional endpoint tests passed\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
