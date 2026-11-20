// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctPartition.hpp"
#include "OcctBoundaryDiagonalClosure.hpp"
#include "OcctModelInternal.hpp"
#include "OcctTrimTolerance.hpp"
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <limits>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace paacs;
namespace {
void require(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
std::string serialize(const TopoDS_Shape& shape) {
    std::ostringstream stream;
    BRepTools::Write(shape,stream);
    return stream.str();
}
void inspectDiagonalSelection(const geo::detail::OcctPartitionOutcome& partition) {
    using geo::detail::closeOcctBoundaryDiagonals;
    const auto isExterior = [&](int tag) {
        return partition.parentByFace[static_cast<std::size_t>(tag)].dim == 2 &&
               partition.model->faceOwners(tag).size() == 1;
    };
    const auto check = [&](int tag,int piece,bool exterior) {
        const auto* face = partition.model->face(tag);
        const mesh::Classification faceClass{2,tag};
        Vec2 lo,hi;
        face->paramBounds(lo,hi);
        const auto point = [&](double u,double v) {
            return face->pointAt({lo.u + (hi.u-lo.u)*u,lo.v + (hi.v-lo.v)*v});
        };
        const auto a = point(.25,.25), b = point(.75,.75);
        const auto c = point(.25,.75), d = point(.75,.25);
        mesh::Mesh shell,partial;
        for (const auto node : {a,b,c,d}) {
            shell.addNode(node,faceClass);
            partial.addNode(node,faceClass);
        }
        shell.addCell(mesh::CellType::Tri3,{0,1,2},faceClass);
        shell.addCell(mesh::CellType::Tri3,{1,0,3},faceClass);
        const auto center = (a+b)*.5;
        const auto offset = (b-a).cross(c-a)*.125;
        auto apex = center + offset;
        if (predicates::orient3d(c,d,a,apex) == predicates::Sign::Negative)
            apex = center - offset;
        require(predicates::orient3d(c,d,a,apex) == predicates::Sign::Positive &&
                predicates::orient3d(c,b,d,apex) == predicates::Sign::Positive,
                "planar diagonal fixture tetrahedra are positive");
        partial.addNode(apex,{3,piece});
        partial.addCell(mesh::CellType::Tet4,{2,3,0,4},{3,piece});
        partial.addCell(mesh::CellType::Tet4,{2,1,3,4},{3,piece});
        const auto charge = [](std::uint64_t) {};
        mesh::Mesh unrestricted,selected;
        require(closeOcctBoundaryDiagonals(*partition.model,shell,partial,1,unrestricted,charge).flips == 1,
                "opposite diagonal has both outward tetrahedral owners");
        selected.addNode({91,92,93});
        const auto receipt = closeOcctBoundaryDiagonals(*partition.model,shell,partial,1,selected,charge,
            [&](Index,Index,mesh::Classification cls) { return isExterior(cls.tag); });
        require(exterior ? receipt.flips == 1 && selected.cellCount() == 2
                         : receipt.flips == 0 && selected.cellCount() == 0 &&
                           selected.nodeCount() == 1 && selected.node(0).x == 91,
                "exterior selection permits the CAD face and preserves the shared cut");
    };
    const auto& cut = partition.interfaces.front();
    check(cut.faceTag,cut.pieceTags[0],false);
    for (int tag = 1; tag < static_cast<int>(partition.parentByFace.size()); ++tag) {
        if (!isExterior(tag)) continue;
        check(tag,partition.model->faceOwners(tag).front().solidTag,true);
        return;
    }
    require(false,"partition cube has a derived exterior face");
}
void inspect(const geo::OcctModel& source,const geo::detail::OcctPartitionOutcome& result) {
    using namespace geo::detail;
    require(result.partitioned(),occtPartitionReasonName(result.reason));
    require(result.reason == OcctPartitionReason::None && result.sourceImmutabilityVerified,
            "partition successful and source unchanged");
    require(result.relativeVolumeChange <= 1e-6,"partition conserves native volume");
    require(result.parentBodyByPiece.size() == result.model->solids().size()+1,
            "body provenance dense");
    require(result.parentByFace.size() == result.model->faces().size()+1,
            "face provenance dense");
    auto classification = [&](mesh::Classification parent) {
        require(parent.dim >= 0 && parent.dim <= 3 && parent.tag > 0,
                "every derived entity maps to original topology");
        require(parent.dim == 0 ? source.vertex(parent.tag) != nullptr :
                parent.dim == 1 ? source.edge(parent.tag) != nullptr :
                parent.dim == 2 ? source.face(parent.tag) != nullptr :
                                 source.solid(parent.tag) != nullptr,
                "parent identity resolves on source model");
    };
    for (std::size_t i=1;i<result.parentByVertex.size();++i) classification(result.parentByVertex[i]);
    for (std::size_t i=1;i<result.parentByEdge.size();++i) classification(result.parentByEdge[i]);
    for (std::size_t i=1;i<result.parentByFace.size();++i) classification(result.parentByFace[i]);
    for (const auto& interface : result.interfaces) {
        require(interface.pieceTags[0] != interface.pieceTags[1] &&
                interface.senses[0] == -interface.senses[1] &&
                result.parentBodyByPiece[interface.pieceTags[0]] == interface.originalBodyTag &&
                result.parentBodyByPiece[interface.pieceTags[1]] == interface.originalBodyTag,
                "shared cut has opposite same-source-body owners");
        require(result.model->faceOwners(interface.faceTag).size() == 2 &&
                result.parentByFace[interface.faceTag] == mesh::Classification{3,interface.originalBodyTag},
                "internal shared face is not claimed as original boundary");
    }
}
void inspectTrimWitness(const geo::OcctModel& model) {
    using namespace geo::detail;
    const auto& original=OcctModelAccess::shape(model);
    const auto before=serialize(original);
    BRepBuilderAPI_Copy copy(original,true,false);
    OcctWitnessShapeMap ss,sf,se,sv,df,de,dv;
    TopExp::MapShapes(original,TopAbs_SOLID,ss);TopExp::MapShapes(original,TopAbs_FACE,sf);
    TopExp::MapShapes(original,TopAbs_EDGE,se);TopExp::MapShapes(original,TopAbs_VERTEX,sv);
    TopExp::MapShapes(copy.Shape(),TopAbs_FACE,df);TopExp::MapShapes(copy.Shape(),TopAbs_EDGE,de);
    TopExp::MapShapes(copy.Shape(),TopAbs_VERTEX,dv);
    std::vector<mesh::Classification> pf(df.Extent()+1),pe(de.Extent()+1),pv(dv.Extent()+1);
    auto parents=[&](const OcctWitnessShapeMap&s,const OcctWitnessShapeMap&d,
                     std::vector<mesh::Classification>&p,int dimension){
        for(int tag=1;tag<=s.Extent();++tag){const auto id=d.FindIndex(copy.ModifiedShape(s(tag)));
            require(id>0,"copied native entity identity");p[id]={dimension,tag};}
    };
    parents(sf,df,pf,2);parents(se,de,pe,1);parents(sv,dv,pv,0);
    std::vector<int> selected{1},bodies{-1,1};
    std::vector<std::vector<std::pair<int,int>>> owners(df.Extent()+1);
    for(int tag=1;tag<=df.Extent();++tag)owners[tag].push_back({1,1});
    OcctTrimToleranceInput input{original,copy,ss,sf,se,sv,df,de,dv,selected,bodies,pf,pe,pv,owners,24,24};
    OcctTrimToleranceOptions options;options.maximumPhysicalDrift=.01;
    const auto exact=auditOcctTrimTolerance(input,options);
    require(exact.admitted&&exact.completed&&exact.sourceUnmodified&&exact.supportsUnchanged&&
            exact.originalFacesRetained&&exact.originalEdgesRetained&&exact.originalVerticesRetained&&
            exact.samples>0&&exact.work>exact.samples&&exact.nativeVolumeDrift==0,
            "exact native copy has complete sampled trim evidence");
    auto refused=[&](OcctTrimToleranceOptions limits,OcctTrimToleranceReason reason){
        const auto result=auditOcctTrimTolerance(input,limits);
        require(!result.admitted&&result.reason==reason,"typed sampled-witness refusal");
        require(serialize(original)==before,"trim-witness refusal preserves original geometry");
    };
    auto limited=options;limited.maximumWork=0;refused(limited,OcctTrimToleranceReason::BudgetExceeded);
    limited=options;limited.maximumWork=1;refused(limited,OcctTrimToleranceReason::BudgetExceeded);
    limited=options;limited.maximumLogicalBytes=1;refused(limited,OcctTrimToleranceReason::BudgetExceeded);
    limited=options;limited.maximumSnapshotBytes=1;refused(limited,OcctTrimToleranceReason::BudgetExceeded);
    limited=options;limited.maximumPhysicalDrift=std::numeric_limits<double>::quiet_NaN();
    refused(limited,OcctTrimToleranceReason::InvalidRequest);
    std::atomic_bool cancelled{true};limited=options;limited.cancellation=&cancelled;
    refused(limited,OcctTrimToleranceReason::Cancelled);
    const auto bodyParent=bodies[1];bodies[1]=-1;
    refused(options,OcctTrimToleranceReason::SourceCoverageLost);bodies[1]=bodyParent;
    const auto originalOwners=owners[1];owners[1]={{2,1}};
    refused(options,OcctTrimToleranceReason::SourceCoverageLost);owners[1]=originalOwners;
    // The lower-level witness must refuse missing original identity, even when
    // the native coordinates and mass still match exactly.
    const auto parent=pv[1];pv[1]={2,1};refused(options,OcctTrimToleranceReason::SourceCoverageLost);pv[1]=parent;
    // Reject an independently copied support rather than treating numeric
    // coincidence as original surface identity.
    BRepBuilderAPI_Copy other(original,true,false);
    OcctTrimToleranceInput support{original,other,ss,sf,se,sv,df,de,dv,selected,bodies,pf,pe,pv,owners,24,24};
    require(auditOcctTrimTolerance(support,options).reason==OcctTrimToleranceReason::SupportMismatch,
            "unrelated copied supports cannot claim original history");
    auto largeVolume=input;largeVolume.derivedVolume=25;
    const auto drift=auditOcctTrimTolerance(largeVolume,options);
    require(drift.completed&&!drift.admitted&&drift.reason==OcctTrimToleranceReason::DriftEnvelopeExceeded,
            "finite native drift must fit computed source tolerance envelope");
    // Inject a native endpoint defect in the owned derived copy. The original
    // native vertex/tolerance is authoritative and stays byte-for-byte fixed.
    BRep_Builder builder;const auto vertex=TopoDS::Vertex(dv(1));const auto oldPoint=BRep_Tool::Pnt(vertex);
    builder.UpdateVertex(vertex,oldPoint.Translated(gp_Vec(.001,0,0)),BRep_Tool::Tolerance(vertex));
    const auto badPoint=auditOcctTrimTolerance(input,options);
    require(badPoint.completed&&!badPoint.admitted&&badPoint.outsideOriginalTolerance>0&&
            badPoint.reason==OcctTrimToleranceReason::OriginalToleranceExceeded,
            "derived endpoint outside its original native tolerance is refused");
    require(serialize(original)==before,"native trim witness never changes original tolerances");
}
}
int main(int argc,char** argv) {
    try {
        using namespace geo::detail;
        auto cube = OcctModelAccess::fromShape(BRepPrimAPI_MakeBox(2,3,4).Shape());
        require(cube.imported(),"native cube snapshot");
        inspectTrimWitness(*cube.model);
        const auto before = serialize(OcctModelAccess::shape(*cube.model));
        OcctPartitionOptions options;
        options.verifyNativeCommon = true;
        auto partition = partitionOcctModel(*cube.model,options);
        inspect(*cube.model,partition);
        require(partition.pieces == 8 && partition.planes == 3 && partition.sharedCutFaces == 12 &&
                partition.mappedExteriorFaces == 24 && partition.nativeCommonChecked,
                "cube eight conforming native regions");
        inspectDiagonalSelection(partition);
        require(std::abs(partition.sourceVolume-24) < 1e-10 &&
                std::abs(partition.nativeCommonVolume-24) < 1e-10,"cube native volume proof");
        require(serialize(OcctModelAccess::shape(*cube.model)) == before,"cube serialized immutable");
        auto gridOptions = options;
        gridOptions.policy = OcctPartitionPolicy::CartesianGrid;
        gridOptions.cartesianDivisions = 3;
        gridOptions.verifyNativeCommon = false;
        auto thirds = partitionOcctModel(*cube.model,gridOptions);
        inspect(*cube.model,thirds);
        require(thirds.pieces == 27 && thirds.planes == 6 && thirds.sharedCutFaces == 54,
                "thirds grid native shared ownership");
        gridOptions.cartesianDivisions = 5;
        auto fifths = partitionOcctModel(*cube.model,gridOptions);
        inspect(*cube.model,fifths);
        require(fifths.pieces == 125 && fifths.planes == 12 && fifths.sharedCutFaces == 300,
                "fifths grid bounded native shared ownership");
        auto refusal = [&](const OcctPartitionOptions& request,OcctPartitionReason reason) {
            const auto result = partitionOcctModel(*cube.model,request);
            require(!result.partitioned() && result.reason == reason,"explicit bounded refusal");
            require(serialize(OcctModelAccess::shape(*cube.model)) == before,"refusal leaves source unchanged");
        };
        auto limited = options; limited.maximumPieces = 1;
        refusal(limited,OcctPartitionReason::BudgetExceeded);
        limited = options; limited.maximumPlanes = 1;
        refusal(limited,OcctPartitionReason::BudgetExceeded);
        limited = options; limited.maximumHistoryWork = 1;
        refusal(limited,OcctPartitionReason::BudgetExceeded);
        limited = options; limited.maximumSnapshotAuditBytes = 1;
        refusal(limited,OcctPartitionReason::BudgetExceeded);
        limited = options; limited.snapshotLimits.maximumUses = 1;
        refusal(limited,OcctPartitionReason::BudgetExceeded);
        limited = options; limited.sourceBodyTags = {1,1};
        refusal(limited,OcctPartitionReason::InvalidRequest);
        limited = options; limited.sourceBodyTags = {2};
        refusal(limited,OcctPartitionReason::InvalidRequest);
        limited = options; limited.maximumRelativeVolumeChange = 1e-2;
        refusal(limited,OcctPartitionReason::InvalidRequest);
        limited=options;limited.allowSampledTrimTolerance=true;limited.maximumPhysicalTrimDrift=.05;
        refusal(limited,OcctPartitionReason::InvalidRequest);
        limited.policy=OcctPartitionPolicy::NativeConcaveFeatureAxes;limited.maximumPhysicalTrimDrift=0;
        refusal(limited,OcctPartitionReason::InvalidRequest);
        limited = gridOptions; limited.cartesianDivisions = 6;
        refusal(limited,OcctPartitionReason::InvalidRequest);
        std::atomic_bool cancel{true};
        limited = options; limited.cancellation = &cancel;
        refusal(limited,OcctPartitionReason::Cancelled);
        limited = options; limited.policy = OcctPartitionPolicy::CoaxialAngular;
        refusal(limited,OcctPartitionReason::UnsupportedPolicy);
        const auto tubeShape = BRepAlgoAPI_Cut(BRepPrimAPI_MakeCylinder(4,5).Shape(),
                                              BRepPrimAPI_MakeCylinder(2,5).Shape()).Shape();
        auto tube = OcctModelAccess::fromShape(tubeShape);
        require(tube.imported(),"native hollow cylinder");
        options.policy = OcctPartitionPolicy::CoaxialAngular;
        options.verifyNativeCommon = false;
        const auto sectors = partitionOcctModel(*tube.model,options);
        inspect(*tube.model,sectors);
        require(sectors.pieces == 16 && sectors.sharedCutFaces == 16,
                "angular partition preserves hollow material");
        options.policy = OcctPartitionPolicy::NativeFeatureAxes;
        options.includeCartesianMidplanes = false;
        options.includePlanarSupportPlanes = false;
        const auto quarters = partitionOcctModel(*tube.model,options);
        inspect(*tube.model,quarters);
        require(quarters.pieces == 4 && quarters.planes == 2 && quarters.sharedCutFaces == 4,
                "native feature axes open hollow trim with shared cuts");
        for (unsigned planes : {4U,8U}) {
            options.featureAxisPlanes = planes;
            const auto fineSectors = partitionOcctModel(*tube.model,options);
            inspect(*tube.model,fineSectors);
            require(fineSectors.pieces == 2*planes && fineSectors.sharedCutFaces == 2*planes,
                    "additional native radial cuts preserve hollow provenance");
        }
        // Circular cap support opens a blind bore at its interior floor.
        // Ordinary feature-axis requests remain unchanged and reproducible.
        const auto blindShape = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(8,8,8).Shape(),
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(4,4,3),gp_Dir(0,0,1)),2,8).Shape()).Shape();
        auto blind = OcctModelAccess::fromShape(blindShape);
        require(blind.imported(),"native blind bore");
        options = {}; options.policy = OcctPartitionPolicy::NativeFeatureAxes;
        options.includeCartesianMidplanes = false; options.includePlanarSupportPlanes = false;
        const auto ordinaryBlind = partitionOcctModel(*blind.model,options);
        inspect(*blind.model,ordinaryBlind);
        options.policy = OcctPartitionPolicy::NativeCircularFeatureAxes;
        const auto circularBlind = partitionOcctModel(*blind.model,options);
        inspect(*blind.model,circularBlind);
        require(circularBlind.planes == ordinaryBlind.planes+1 &&
                circularBlind.pieces > ordinaryBlind.pieces,"native circular floor support contributes bounded cut");
        options.policy = OcctPartitionPolicy::NativeFeatureAxes;
        const auto unchangedBlind = partitionOcctModel(*blind.model,options);
        inspect(*blind.model,unchangedBlind);
        require(ordinaryBlind.planes == unchangedBlind.planes && ordinaryBlind.pieces == unchangedBlind.pieces &&
                ordinaryBlind.historyWork == unchangedBlind.historyWork &&
                serialize(OcctModelAccess::shape(*ordinaryBlind.model)) ==
                    serialize(OcctModelAccess::shape(*unchangedBlind.model)),"ordinary policy remains byte reproducible");
        options = {}; options.policy = OcctPartitionPolicy::NativeConcaveFeatureAxes;
        options.featureAxisPlanes = 8; options.includeCartesianMidplanes = false;
        const auto concaveBlind = partitionOcctModel(*blind.model,options);
        inspect(*blind.model,concaveBlind);
        require(concaveBlind.planes == 9 && concaveBlind.pieces > 16,
                "principal radial cuts and interior axial bore floor preserve native body");
        options.maximumPlanes = 1;
        const auto concaveBudget = partitionOcctModel(*blind.model,options);
        require(!concaveBudget.partitioned() && concaveBudget.reason == OcctPartitionReason::BudgetExceeded &&
                concaveBudget.sourceImmutabilityVerified,"concave native planes obey work admission");
        options = {}; options.policy = OcctPartitionPolicy::NativeConcaveFeatureAxes;
        const auto noCylinder = partitionOcctModel(*cube.model,options);
        require(!noCylinder.partitioned() && noCylinder.reason == OcctPartitionReason::UnsupportedPolicy &&
                noCylinder.sourceImmutabilityVerified,"no native cylinder gives explicit policy refusal");
        options.featureAxisPlanes = 17;
        refusal(options,OcctPartitionReason::InvalidRequest);
        options.featureAxisPlanes = 8; options.cancellation = &cancel;
        refusal(options,OcctPartitionReason::Cancelled);
        // A bored L-section has transverse interior planar supports that the
        // axial-only concave policy deliberately leaves uncut. The principal
        // policy adds those exact supports without adding small fillet axes.
        const auto lSection = BRepAlgoAPI_Cut(
            BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(8,3,8).Shape(),
                            BRepPrimAPI_MakeBox(3,8,8).Shape()).Shape(),
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1.5,1.5,0),gp_Dir(0,0,1)),.5,8).Shape()).Shape();
        auto principalSource = OcctModelAccess::fromShape(lSection);
        require(principalSource.imported(),"bored native L-section");
        const auto principalBefore = serialize(OcctModelAccess::shape(*principalSource.model));
        options = {}; options.policy = OcctPartitionPolicy::NativeConcaveFeatureAxes;
        options.featureAxisPlanes = 4; options.includeCartesianMidplanes = false;
        const auto axialOnly = partitionOcctModel(*principalSource.model,options);
        inspect(*principalSource.model,axialOnly);
        options.policy = OcctPartitionPolicy::NativePrincipalFeatureAxes;
        const auto planarSupports = partitionOcctModel(*principalSource.model,options);
        inspect(*principalSource.model,planarSupports);
        {
            auto candidate=options;candidate.meshConstructionCandidate=true;candidate.maximumRelativeVolumeChange=1e-4;
            inspect(*principalSource.model,partitionOcctModel(*principalSource.model,candidate));
            candidate.maximumRelativeVolumeChange=1e-3;auto rejected=partitionOcctModel(*principalSource.model,candidate);
            require(!rejected.partitioned()&&rejected.reason==OcctPartitionReason::InvalidRequest,"construction candidate drift cap");
            candidate.maximumRelativeVolumeChange=1e-4;candidate.meshConstructionCandidate=false;rejected=partitionOcctModel(*principalSource.model,candidate);
            require(!rejected.partitioned()&&rejected.reason==OcctPartitionReason::InvalidRequest,"default native conservation admission stays strict");
            candidate.meshConstructionCandidate=true;candidate.verifyNativeCommon=true;rejected=partitionOcctModel(*principalSource.model,candidate);
            require(!rejected.partitioned()&&rejected.reason==OcctPartitionReason::InvalidRequest,"mesh construction is not a certified native Common");
            candidate.verifyNativeCommon=false;candidate.policy=OcctPartitionPolicy::CartesianMidplanes;rejected=partitionOcctModel(*principalSource.model,candidate);
            require(!rejected.partitioned()&&rejected.reason==OcctPartitionReason::InvalidRequest,"construction candidate requires native feature policy");
        }
        require(planarSupports.planes > axialOnly.planes && planarSupports.pieces > axialOnly.pieces,
                "principal policy includes transverse original planar supports");
        require(serialize(OcctModelAccess::shape(*principalSource.model)) == principalBefore,
                "principal planar-support preparation preserves original CAD");
        options.includePlanarSupportPlanes = false;
        const auto principalAxesOnly = partitionOcctModel(*principalSource.model,options);
        inspect(*principalSource.model,principalAxesOnly);
        require(principalAxesOnly.planes < planarSupports.planes,
                "principal planar-support option is respected");
        gp_Trsf principalRotation;
        principalRotation.SetRotation(gp_Ax1(gp_Pnt(4,-3,2),gp_Dir(1,2,3)),.41);
        auto rotatedPrincipal = OcctModelAccess::fromShape(
            BRepBuilderAPI_Transform(lSection,principalRotation,true).Shape());
        require(rotatedPrincipal.imported(),"rotated bored L-section");
        options.includePlanarSupportPlanes = true;
        inspect(*rotatedPrincipal.model,partitionOcctModel(*rotatedPrincipal.model,options));
        options.maximumPlanes = 1;
        const auto principalBudget = partitionOcctModel(*principalSource.model,options);
        require(!principalBudget.partitioned() && principalBudget.reason == OcctPartitionReason::BudgetExceeded &&
                principalBudget.sourceImmutabilityVerified,"principal support cuts respect plane cap");
        options = {}; options.policy = OcctPartitionPolicy::NativePrincipalFeatureAxes;
        options.maximumHistoryWork = 1;
        const auto principalWork = partitionOcctModel(*principalSource.model,options);
        require(!principalWork.partitioned() && principalWork.reason == OcctPartitionReason::BudgetExceeded,
                "principal proposals spend bounded history work");
        options = {}; options.policy = OcctPartitionPolicy::NativePrincipalFeatureAxes;
        const auto principalNoAxis = partitionOcctModel(*cube.model,options);
        require(!principalNoAxis.partitioned() && principalNoAxis.reason == OcctPartitionReason::UnsupportedPolicy &&
                principalNoAxis.sourceImmutabilityVerified,"principal policy requires an actual native cylinder");
        options.featureAxisPlanes = 17; refusal(options,OcctPartitionReason::InvalidRequest);
        options.featureAxisPlanes = 8; options.cancellation = &cancel; refusal(options,OcctPartitionReason::Cancelled);
        require(std::string_view(occtPartitionPolicyName(OcctPartitionPolicy::NativePrincipalFeatureAxes)) ==
                    "native_principal_feature_axes","principal policy has a named receipt");
        // A rotated native cone and torus retain the same source topology gates.
        auto cone = OcctModelAccess::fromShape(BRepPrimAPI_MakeCone(
            gp_Ax2(gp_Pnt(7,-2,5),gp_Dir(1,2,3)),4,2,5).Shape());
        require(cone.imported(),"rotated native cone");
        options = {}; options.policy = OcctPartitionPolicy::NativeCircularFeatureAxes;
        options.includeCartesianMidplanes = false; options.includePlanarSupportPlanes = false;
        options.featureAxisPlanes = 4;
        inspect(*cone.model,partitionOcctModel(*cone.model,options));
        auto torus = OcctModelAccess::fromShape(BRepPrimAPI_MakeTorus(
            gp_Ax2(gp_Pnt(-3,4,5),gp_Dir(2,1,3)),5,1).Shape());
        require(torus.imported(),"rotated native ring torus");
        inspect(*torus.model,partitionOcctModel(*torus.model,options));
        // Explicit cuts add no new route and preserve identity by original history.
        options = {}; options.extraPlanes = {{{.75,0,0},{1,0,0}}};
        const auto extra = partitionOcctModel(*cube.model,options);
        inspect(*cube.model,extra);
        require(extra.planes == 4 && extra.pieces == 12,"explicit off-midplane cut");
        options.extraPlanes.push_back({{.75,3,4},{-2,0,0}});
        const auto duplicate = partitionOcctModel(*cube.model,options);
        inspect(*cube.model,duplicate);
        require(duplicate.planes == 4 && duplicate.pieces == 12,"equivalent additional planes deduplicate without identity welding");
        options.maximumPlanes = 3; refusal(options,OcctPartitionReason::BudgetExceeded);
        options = {}; options.extraPlanes = {{{0,0,0},{0,0,0}}};
        refusal(options,OcctPartitionReason::InvalidRequest);
        options.extraPlanes = {{{std::numeric_limits<double>::quiet_NaN(),0,0},{1,0,0}}};
        refusal(options,OcctPartitionReason::InvalidRequest);
        options.extraPlanes = {{{0,0,0},{std::numeric_limits<double>::infinity(),0,0}}};
        refusal(options,OcctPartitionReason::InvalidRequest);
        // Native support planes coincide with existing hole-wall corners.
        // OCCT may reconstruct these vertices while omitting public vertex
        // Modified/Generated history; exact same-domain origins must survive.
        const auto squareHole = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(4,4,4).Shape(),
            BRepPrimAPI_MakeBox(gp_Pnt(1,1,0),1,1,4).Shape()).Shape();
        auto holed = OcctModelAccess::fromShape(squareHole);
        require(holed.imported(),"planar through-hole source");
        options = {}; options.policy = OcctPartitionPolicy::NativeFeatureAxes;
        const auto supported = partitionOcctModel(*holed.model,options);
        inspect(*holed.model,supported);
        require(supported.pieces == 16 && std::abs(supported.sourceVolume-60) < 1e-9,
                "coincident support planes preserve square-hole material");
        for (const auto* vertex : holed.model->vertices())
            require(std::find(supported.parentByVertex.begin(),supported.parentByVertex.end(),
                              mesh::Classification{0,vertex->tag()}) != supported.parentByVertex.end(),
                    "every surviving original hole corner retains vertex identity");
        // A selected body keeps original identity even when preceding original
        // bodies are absent from the derived snapshot.
        BRep_Builder builder; TopoDS_Compound assembly; builder.MakeCompound(assembly);
        builder.Add(assembly,BRepPrimAPI_MakeBox(2,3,4).Shape());
        builder.Add(assembly,BRepPrimAPI_MakeBox(gp_Pnt(10,0,0),2,3,4).Shape());
        auto sourceAssembly = OcctModelAccess::fromShape(assembly);
        require(sourceAssembly.imported(),"separated assembly");
        options = {}; options.sourceBodyTags = {2};
        const auto selected = partitionOcctModel(*sourceAssembly.model,options);
        inspect(*sourceAssembly.model,selected);
        for (std::size_t i=1;i<selected.parentBodyByPiece.size();++i)
            require(selected.parentBodyByPiece[i] == 2,"selected original body tag preserved");
        // Optional real corpus witnesses use the same helper and receipts.
        options = {};
        std::filesystem::path outputDirectory;
        for (int i=1;i<argc;++i) {
            const std::string argument = argv[i];
            if (argument == "--angular") { options.policy = OcctPartitionPolicy::CoaxialAngular; continue; }
            if (argument == "--cartesian") { options.policy = OcctPartitionPolicy::CartesianMidplanes; continue; }
            if (argument == "--circular") { options.policy = OcctPartitionPolicy::NativeCircularFeatureAxes; continue; }
            if (argument == "--concave") { options.policy = OcctPartitionPolicy::NativeConcaveFeatureAxes; continue; }
            if (argument == "--principal") { options.policy = OcctPartitionPolicy::NativePrincipalFeatureAxes; continue; }
            if (argument == "--features") { options.policy = OcctPartitionPolicy::NativeFeatureAxes;
                options.includeCartesianMidplanes = true; options.includePlanarSupportPlanes = true; continue; }
            if (argument == "--axes-only") { options.policy = OcctPartitionPolicy::NativeFeatureAxes;
                options.includeCartesianMidplanes = false; options.includePlanarSupportPlanes = false; continue; }
            if (argument.starts_with("--grid=")) { options.policy = OcctPartitionPolicy::CartesianGrid;
                options.cartesianDivisions = static_cast<unsigned>(std::stoul(argument.substr(7))); continue; }
            if (argument.starts_with("--axis-planes=")) { options.featureAxisPlanes =
                static_cast<unsigned>(std::stoul(argument.substr(14))); continue; }
            if (argument.starts_with("--bodies=")) {
                options.sourceBodyTags.clear();
                std::stringstream tags(argument.substr(9)); std::string tag;
                while (std::getline(tags,tag,',')) options.sourceBodyTags.push_back(std::stoi(tag));
                continue;
            }
            if (argument.starts_with("--output-dir=")) {
                outputDirectory = argument.substr(13);
                std::filesystem::create_directories(outputDirectory);
                continue;
            }
            auto imported = geo::importOcctCad(argv[i]);
            require(imported.imported(),"real corpus native import");
            const auto real = partitionOcctModel(*imported.model,options);
            inspect(*imported.model,real);
            if (!outputDirectory.empty())
                require(BRepTools::Write(OcctModelAccess::shape(*real.model),
                    (outputDirectory/(std::filesystem::path(argv[i]).stem().string()+"-partition.brep")).string().c_str()),
                    "private real corpus partition export");
            std::cout << argv[i] << " planes=" << real.planes << " pieces=" << real.pieces << " interfaces=" << real.sharedCutFaces
                      << " volume_error=" << real.relativeVolumeChange << " history=" << real.historyWork
                      << " seconds=" << real.totalSeconds << '\n';
        }
        std::cout << "OCCT partition provenance, conservation, identity, and refusal checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
