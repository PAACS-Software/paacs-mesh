// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

#include "paacs/algo/TetCandidate.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/topo/Workspace.hpp"
#include "../src/algo/TetCandidateDetail.hpp"
#include "test_harness.hpp"

namespace {
using namespace paacs;

std::vector<Vec3> difficultPoints() {
    return {{0, 0, 0}, {3, 0, 0}, {0, 4, 0}, {0, 0, 5},
            {1, 1, 1}, {2, 1, 1}, {1, 2, 1}, {1, 1, 3},
            {2, 2, 2}, {1, 3, 2}, {2, 1, 3}};
}

bool sameGeometry(const algo::TetCandidate& a,
                  const algo::TetCandidate& b) {
    if (a.points.size() != b.points.size() ||
        a.cells.size() != b.cells.size() ||
        a.boundaryFaces.size() != b.boundaryFaces.size())
        return false;
    for (std::size_t i = 0; i < a.points.size(); ++i)
        if (a.points[i].x != b.points[i].x ||
            a.points[i].y != b.points[i].y ||
            a.points[i].z != b.points[i].z)
            return false;
    for (std::size_t i = 0; i < a.cells.size(); ++i)
        if (a.cells[i].vertices != b.cells[i].vertices ||
            a.cells[i].neighbours != b.cells[i].neighbours)
            return false;
    for (std::size_t i = 0; i < a.boundaryFaces.size(); ++i)
        if (a.boundaryFaces[i].vertices != b.boundaryFaces[i].vertices ||
            a.boundaryFaces[i].owner != b.boundaryFaces[i].owner ||
            a.boundaryFaces[i].ownerLocalFace !=
                b.boundaryFaces[i].ownerLocalFace)
            return false;
    return true;
}

bool sameMesh(const mesh::Mesh& a, const mesh::Mesh& b) {
    if (a.nodeCount() != b.nodeCount() || a.cellCount() != b.cellCount() ||
        a.dedupTolerance() != b.dedupTolerance())
        return false;
    for (Index i = 0; i < a.nodeCount(); ++i) {
        const auto& ap = a.node(i);
        const auto& bp = b.node(i);
        const auto ac = a.nodeClass(i);
        const auto bc = b.nodeClass(i);
        if (ap.x != bp.x || ap.y != bp.y || ap.z != bp.z ||
            ac.dim != bc.dim || ac.tag != bc.tag)
            return false;
    }
    for (Index i = 0; i < a.cellCount(); ++i) {
        const auto& ac = a.cell(i);
        const auto& bc = b.cell(i);
        if (ac.type != bc.type || ac.onEntity.dim != bc.onEntity.dim ||
            ac.onEntity.tag != bc.onEntity.tag)
            return false;
        const int count = mesh::nodesPerCell(ac.type);
        for (int j = 0; j < count; ++j)
            if (a.cellNodes(i)[j] != b.cellNodes(i)[j]) return false;
    }
    return true;
}

mesh::Mesh globalForCandidate(
    const algo::TetCandidate& candidate,
    const std::vector<mesh::Classification>& classes) {
    mesh::Mesh global;
    for (std::size_t i = 0; i < candidate.points.size(); ++i)
        global.addNode(candidate.points[i], classes[i]);
    return global;
}

struct TestCountingIndexLess {
    std::uint64_t* comparisons=nullptr;
    bool operator()(Index a,Index b) const noexcept {
        ++*comparisons;
        return a<b;
    }
};

void requireMeasuredHandoffWork(
    const algo::TetHandoffMetrics& metrics,
    const algo::TetCandidate& candidate,
    const std::vector<Index>& selected,
    const std::vector<Index>& globals,
    Index globalNodeCount) {
    const auto points=static_cast<std::uint64_t>(candidate.points.size());
    const auto inputs=static_cast<std::uint64_t>(
        candidate.inputToPoint.size());
    const auto cells=static_cast<std::uint64_t>(candidate.cells.size());
    const auto boundary=static_cast<std::uint64_t>(
        candidate.boundaryFaces.size());
    CHECK(metrics.importCandidateValidationCalls==1);
    CHECK(metrics.importCandidateValidationPointRows==3*points+inputs);
    CHECK(metrics.importCandidateValidationCellRows==cells);
    CHECK(metrics.importCandidateValidationFaceUses==4*cells);
    CHECK(metrics.importCandidateValidationFaceRows==(4*cells+boundary)/2);
    CHECK(metrics.importCandidateValidationBoundaryRows==3*boundary);
    CHECK(metrics.importCandidateValidationDualNeighbourVisits==4*cells);
    CHECK(metrics.importCandidateValidationBoundaryEdgeUses==3*boundary);
    CHECK(metrics.importCandidateValidationBoundaryNeighbourVisits==
          3*boundary);
    CHECK(metrics.importCandidateValidationHullTests==
          boundary*(points-3));
    CHECK(metrics.importCandidateValidationSlotsInitialized==
          2*(points+cells+boundary));

    std::uint64_t expectedSelectionSortComparisons=0;
    auto canonicalSelection=selected;
    std::sort(canonicalSelection.begin(),canonicalSelection.end(),
              [&](Index a,Index b) {
        ++expectedSelectionSortComparisons;
        return a<b;
    });
    CHECK(metrics.importSelectionRowsCopied==selected.size());
    CHECK(metrics.importSelectionSortComparisons==
          expectedSelectionSortComparisons);
    CHECK(metrics.importSelectionDuplicateComparisons==
          (selected.empty()?0:selected.size()-1));
    CHECK(metrics.importSelectionCellVisits==selected.size());
    CHECK(metrics.importSelectionVertexVisits==4*selected.size());
    CHECK(metrics.importSelectionPointSlotsInitialized==points);

    std::vector<bool> used(candidate.points.size(),false);
    for (const Index cell:canonicalSelection)
        for (const Index point:
             candidate.cells[static_cast<std::size_t>(cell)].vertices)
            used[static_cast<std::size_t>(point)]=true;
    std::uint64_t expectedGlobalComparisons=0;
    std::set<Index,TestCountingIndexLess> uniqueGlobals(
        TestCountingIndexLess{&expectedGlobalComparisons});
    std::uint64_t usedPoints=0;
    for (std::size_t point=0;point<used.size();++point) {
        if (!used[point]) continue;
        ++usedPoints;
        (void)uniqueGlobals.insert(globals[point]);
    }
    CHECK(metrics.importGlobalMapRows==usedPoints);
    CHECK(metrics.importGlobalMapComparisons==expectedGlobalComparisons);
    CHECK(metrics.importGlobalMapInsertions==usedPoints);
    CHECK(metrics.importLocalPointSlotsInitialized==points);
    CHECK(metrics.importLocalGlobalSlotsInitialized==
          static_cast<std::uint64_t>(globalNodeCount));
    CHECK(metrics.importLocalCellSlotsInitialized==cells);
    CHECK(metrics.importPersistentNodeSlotsInitialized==usedPoints);
}

void candidateConstructionAndMetrics() {
    TEST("ck3.tet_candidate_complete");
    const auto points = difficultPoints();
    const auto built = algo::buildTetCandidate(points);
    CHECK(built.built());
    if (!built.built()) return;
    const auto& candidate = *built.candidate;
    const auto validation = algo::validateTetCandidate(candidate);
    CHECK(validation.valid);
    CHECK(candidate.metrics.runComplete);
    CHECK(candidate.metrics.identityComplete);
    CHECK(candidate.metrics.adjacencyComplete);
    CHECK(candidate.metrics.boundaryComplete);
    CHECK(candidate.metrics.inputPoints == points.size());
    CHECK(candidate.metrics.uniquePoints == points.size());
    CHECK(candidate.metrics.candidateCells == candidate.cells.size());
    CHECK(candidate.metrics.boundaryFaces == candidate.boundaryFaces.size());
    CHECK(candidate.metrics.validationFaceUses == 4 * candidate.cells.size());
    CHECK(candidate.metrics.construction.runComplete);
    CHECK(candidate.metrics.construction.auditPassed);
    CHECK(candidate.metrics.construction.outputCells == candidate.cells.size());
    CHECK(candidate.metrics.construction.walkCellVisits > 0);
    CHECK(candidate.metrics.construction.conflictQueries > 0);
    CHECK(candidate.metrics.construction.cavityCellVisits > 0);
    CHECK(candidate.metrics.construction.shellFacesEmitted > 0);
    CHECK(candidate.metrics.candidateAccountedBytes >= sizeof(candidate));
    CHECK(candidate.metrics.constructionCandidateConservativeBytes >=
          candidate.metrics.candidateAccountedBytes);

    std::vector<Index> selected(candidate.cells.size());
    for (std::size_t i = 0; i < selected.size(); ++i)
        selected[i] = static_cast<Index>(i);
    std::vector<Index> globals(candidate.points.size());
    std::vector<mesh::Classification> classes(candidate.points.size());
    for (std::size_t i = 0; i < candidate.points.size(); ++i) {
        globals[i] = static_cast<Index>(i);
        classes[i] = {i < 4 ? 2 : 3, i < 4 ? 70 : 9};
    }
    const auto global = globalForCandidate(candidate, classes);
    auto imported = algo::importTetCandidate(
        candidate, selected, globals, classes, global, {3, 9});
    CHECK(imported.imported());
    if (!imported.imported()) return;
    const auto& result = *imported.result;
    CHECK(result.metrics.runComplete);
    CHECK(result.metrics.workspaceValid);
    CHECK(result.metrics.workAccountingComplete);
    CHECK(result.metrics.importedNodes == candidate.points.size());
    CHECK(result.metrics.importedCells == candidate.cells.size());
    requireMeasuredHandoffWork(
        result.metrics,candidate,selected,globals,global.nodeCount());
    CHECK(result.metrics.modeledOverlapPayloadBytes >=
          result.metrics.candidateAccountedBytes);
    CHECK(result.metrics.completeTetRepresentationsAtOverlap == 2);
    CHECK(result.metrics.completeTetRepresentationsBeforePublication == 2);
    CHECK(result.workspace->adjacencyReady());
    CHECK(result.workspace->nodes().activeBySequence().size() ==
          candidate.points.size());
    CHECK(result.workspace->cells().activeBySequence().size() ==
          candidate.cells.size());
    CHECK(result.nodeBindings.size() == candidate.points.size());
    for (const auto& binding : result.nodeBindings) {
        const auto point = static_cast<std::size_t>(binding.candidatePoint);
        CHECK(point < candidate.points.size());
        if (point >= candidate.points.size()) continue;
        CHECK(binding.existingGlobalIndex == globals[point]);
        const auto node = result.workspace->nodes().get(binding.local);
        CHECK(node.coord.x == candidate.points[point].x);
        CHECK(node.coord.y == candidate.points[point].y);
        CHECK(node.coord.z == candidate.points[point].z);
        CHECK(node.cls.dim == classes[point].dim);
        CHECK(node.cls.tag == classes[point].tag);
        CHECK(result.existingGlobalByLocalNodeSlot[
                  static_cast<std::size_t>(binding.local.slot)] ==
              globals[point]);
    }
    for (const auto cellHandle :
         result.workspace->cells().activeBySequence()) {
        const auto cell = result.workspace->cells().get(cellHandle);
        CHECK(cell.type == topo::ElementType::Tet4);
        CHECK(cell.cls.dim == 3);
        CHECK(cell.cls.tag == 9);
        const auto a = result.workspace->nodes().get(cell.nodes[0]).coord;
        const auto b = result.workspace->nodes().get(cell.nodes[1]).coord;
        const auto c = result.workspace->nodes().get(cell.nodes[2]).coord;
        const auto d = result.workspace->nodes().get(cell.nodes[3]).coord;
        CHECK(predicates::orient3d(a, b, c, d) ==
              predicates::Sign::Positive);
    }

    std::array<Index,3> protectedFace{};
    bool foundInteriorFace=false;
    for (const auto& cell:candidate.cells) {
        for (std::size_t localFace=0;localFace<4;++localFace) {
            if (cell.neighbours[localFace]==kInvalidIndex) continue;
            static constexpr std::array<std::array<std::size_t,3>,4> faces{{
                {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
            for (std::size_t i=0;i<3;++i)
                protectedFace[i]=globals[static_cast<std::size_t>(
                    cell.vertices[faces[localFace][i]])];
            foundInteriorFace=true;
            break;
        }
        if (foundInteriorFace) break;
    }
    CHECK(foundInteriorFace);
    if (!foundInteriorFace) return;
    const algo::TetImportProtectedEdge protectedEdge{
        {protectedFace[0],protectedFace[1]},{1,701}};
    const algo::TetImportProtectedFace protectedFaceRow{
        protectedFace,{2,702}};
    auto protectedImport=algo::importTetCandidateWithProtectedBoundary(
        candidate,selected,globals,classes,global,{3,9},
        std::span<const algo::TetImportProtectedEdge>(&protectedEdge,1),
        std::span<const algo::TetImportProtectedFace>(&protectedFaceRow,1));
    CHECK(protectedImport.imported());
    if (protectedImport.imported()) {
        std::vector<topo::NodeHandle> byGlobal(candidate.points.size());
        for (const auto& binding:protectedImport.result->nodeBindings)
            byGlobal[static_cast<std::size_t>(binding.existingGlobalIndex)]=
                binding.local;
        const auto edge=protectedImport.result->workspace->edges().find(
            protectedImport.result->workspace->nodes(),
            byGlobal[static_cast<std::size_t>(protectedFace[0])],
            byGlobal[static_cast<std::size_t>(protectedFace[1])]);
        CHECK(edge.has_value());
        CHECK(edge&&protectedImport.result->workspace->edges().get(*edge).
              requiredConstraint);
        CHECK(edge&&protectedImport.result->workspace->edges().get(*edge).
              cls.dim==1);
        CHECK(edge&&protectedImport.result->workspace->edges().get(*edge).
              cls.tag==701);
        const std::array<topo::NodeHandle,3> localFace{
            byGlobal[static_cast<std::size_t>(protectedFace[0])],
            byGlobal[static_cast<std::size_t>(protectedFace[1])],
            byGlobal[static_cast<std::size_t>(protectedFace[2])]};
        const auto face=protectedImport.result->workspace->faces().find(
            protectedImport.result->workspace->nodes(),localFace);
        CHECK(face.has_value());
        CHECK(face&&protectedImport.result->workspace->faces().get(*face).
              cls.dim==2);
        CHECK(face&&protectedImport.result->workspace->faces().get(*face).
              cls.tag==702);
    }
}

void partialProtectedImportMarksOnlyExactPresent() {
    TEST("ck3.tet_partial_protected_import");
    const auto built=algo::buildTetCandidate(difficultPoints());
    CHECK(built.built());
    if (!built.built()) return;
    const auto& candidate=*built.candidate;

    std::vector<Index> selected(candidate.cells.size());
    std::vector<Index> globals(candidate.points.size());
    std::vector<mesh::Classification> classes(
        candidate.points.size(),mesh::Classification{3,81});
    for (std::size_t i=0;i<selected.size();++i)
        selected[i]=static_cast<Index>(i);
    for (std::size_t i=0;i<globals.size();++i)
        globals[i]=static_cast<Index>(i);
    const auto global=globalForCandidate(candidate,classes);
    const auto globalBefore=global;

    static constexpr std::array<std::array<std::size_t,3>,4> faces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    std::array<Index,3> interiorGlobal{};
    bool foundInterior=false;
    for (const auto& cell:candidate.cells) {
        for (std::size_t localFace=0;localFace<4;++localFace) {
            if (cell.neighbours[localFace]==kInvalidIndex) continue;
            for (std::size_t i=0;i<3;++i)
                interiorGlobal[i]=globals[static_cast<std::size_t>(
                    cell.vertices[faces[localFace][i]])];
            foundInterior=true;
            break;
        }
        if (foundInterior) break;
    }
    CHECK(foundInterior);
    CHECK(!candidate.boundaryFaces.empty());
    if (!foundInterior||candidate.boundaryFaces.empty()) return;
    std::array<Index,3> boundaryGlobal{};
    for (std::size_t i=0;i<3;++i)
        boundaryGlobal[i]=globals[static_cast<std::size_t>(
            candidate.boundaryFaces[0].vertices[i])];

    std::array<Index,2> missingGlobal{kInvalidIndex,kInvalidIndex};
    for (std::size_t a=0;a<candidate.points.size()&&
                         missingGlobal[0]==kInvalidIndex;++a)
        for (std::size_t b=a+1;b<candidate.points.size();++b) {
            bool present=false;
            for (const auto& cell:candidate.cells) {
                bool hasA=false,hasB=false;
                for (const Index point:cell.vertices) {
                    hasA=hasA||point==static_cast<Index>(a);
                    hasB=hasB||point==static_cast<Index>(b);
                }
                if (hasA&&hasB) { present=true; break; }
            }
            if (!present) {
                missingGlobal={globals[a],globals[b]};
                break;
            }
        }
    CHECK(missingGlobal[0]!=kInvalidIndex);
    if (missingGlobal[0]==kInvalidIndex) return;

    const std::array<algo::TetImportProtectedEdge,2> edgeRows{{
        {{interiorGlobal[1],interiorGlobal[0]},{2,701}},
        {{missingGlobal[1],missingGlobal[0]},{1,704}}}};
    const std::array<algo::TetImportProtectedFace,2> faceRows{{
        {{interiorGlobal[2],interiorGlobal[0],interiorGlobal[1]},{2,702}},
        {{boundaryGlobal[1],boundaryGlobal[2],boundaryGlobal[0]},{2,703}}}};

    auto imported=algo::tet_candidate_detail::
        importTetCandidateWithPartialProtectedBoundary(
        candidate,selected,globals,classes,global,{3,81},edgeRows,faceRows);
    CHECK(imported.imported());
    if (!imported.imported()) return;
    std::uint64_t expectedSortComparisons=0;
    auto countedEdges=edgeRows;
    auto countedFaces=faceRows;
    for (auto& row:countedEdges)
        std::sort(row.existingGlobalNodes.begin(),
                  row.existingGlobalNodes.end(),
                  [&](Index a,Index b) {
            ++expectedSortComparisons;
            return a<b;
        });
    std::sort(countedEdges.begin(),countedEdges.end(),
              [&](const auto& a,const auto& b) {
        ++expectedSortComparisons;
        return a.existingGlobalNodes<b.existingGlobalNodes;
    });
    for (auto& row:countedFaces)
        std::sort(row.existingGlobalNodes.begin(),
                  row.existingGlobalNodes.end(),
                  [&](Index a,Index b) {
            ++expectedSortComparisons;
            return a<b;
        });
    std::sort(countedFaces.begin(),countedFaces.end(),
              [&](const auto& a,const auto& b) {
        ++expectedSortComparisons;
        return a.existingGlobalNodes<b.existingGlobalNodes;
    });
    std::vector<Index> protectedGlobals;
    for (const auto& row:edgeRows)
        protectedGlobals.insert(protectedGlobals.end(),
                                row.existingGlobalNodes.begin(),
                                row.existingGlobalNodes.end());
    for (const auto& row:faceRows)
        protectedGlobals.insert(protectedGlobals.end(),
                                row.existingGlobalNodes.begin(),
                                row.existingGlobalNodes.end());
    std::sort(protectedGlobals.begin(),protectedGlobals.end());
    protectedGlobals.erase(
        std::unique(protectedGlobals.begin(),protectedGlobals.end()),
        protectedGlobals.end());
    const auto& partialMetrics=imported.result->metrics;
    requireMeasuredHandoffWork(
        partialMetrics,candidate,selected,globals,global.nodeCount());
    CHECK(partialMetrics.partialImportRowsCopied==4);
    CHECK(partialMetrics.partialImportGlobalSlotsInitialized==
          static_cast<std::uint64_t>(global.nodeCount()));
    CHECK(partialMetrics.partialImportCandidateMappingVisits==
          static_cast<std::uint64_t>(globals.size()));
    CHECK(partialMetrics.partialImportEdgeValidationRows==2);
    CHECK(partialMetrics.partialImportFaceValidationRows==2);
    CHECK(partialMetrics.partialImportProtectedEndpointVisits==20);
    CHECK(partialMetrics.partialImportAuthorityPointValidations==
          static_cast<std::uint64_t>(protectedGlobals.size()));
    CHECK(partialMetrics.partialImportDuplicateKeyComparisons==2);
    CHECK(partialMetrics.partialImportResolutionNodeVisits==10);
    CHECK(partialMetrics.partialImportEdgeResolutionLookups==2);
    CHECK(partialMetrics.partialImportFaceResolutionLookups==2);
    CHECK(partialMetrics.partialImportEdgeDataReads==1);
    CHECK(partialMetrics.partialImportFaceDataCopies==2);
    CHECK(partialMetrics.partialImportFaceIncidenceRowsCopied==3);
    CHECK(partialMetrics.partialImportSortComparisons==
          expectedSortComparisons);
    const auto& workspace=*imported.result->workspace;
    std::vector<topo::NodeHandle> byGlobal(
        static_cast<std::size_t>(global.nodeCount()));
    for (const auto& binding:imported.result->nodeBindings)
        byGlobal[static_cast<std::size_t>(binding.existingGlobalIndex)]=
            binding.local;
    const auto exactEdge=workspace.edges().find(
        workspace.nodes(),
        byGlobal[static_cast<std::size_t>(interiorGlobal[0])],
        byGlobal[static_cast<std::size_t>(interiorGlobal[1])]);
    CHECK(exactEdge.has_value());
    if (exactEdge) {
        const auto data=workspace.edges().get(*exactEdge);
        CHECK(data.requiredConstraint);
        CHECK(data.cls.dim==2);
        CHECK(data.cls.tag==701);
    }
    const auto missingEdge=workspace.edges().find(
        workspace.nodes(),
        byGlobal[static_cast<std::size_t>(missingGlobal[0])],
        byGlobal[static_cast<std::size_t>(missingGlobal[1])]);
    CHECK(!missingEdge.has_value());

    const std::array<topo::NodeHandle,3> exactLocalFace{
        byGlobal[static_cast<std::size_t>(interiorGlobal[0])],
        byGlobal[static_cast<std::size_t>(interiorGlobal[1])],
        byGlobal[static_cast<std::size_t>(interiorGlobal[2])]};
    const auto exactFace=workspace.faces().find(
        workspace.nodes(),exactLocalFace);
    CHECK(exactFace.has_value());
    if (exactFace) {
        const auto data=workspace.faces().get(*exactFace);
        CHECK(data.incident.size()==2);
        CHECK(data.incident.size()==2&&
              data.incident[0].sense*data.incident[1].sense==-1);
        CHECK(data.cls.dim==2);
        CHECK(data.cls.tag==702);
    }
    const std::array<topo::NodeHandle,3> boundaryLocalFace{
        byGlobal[static_cast<std::size_t>(boundaryGlobal[0])],
        byGlobal[static_cast<std::size_t>(boundaryGlobal[1])],
        byGlobal[static_cast<std::size_t>(boundaryGlobal[2])]};
    const auto boundaryFace=workspace.faces().find(
        workspace.nodes(),boundaryLocalFace);
    CHECK(boundaryFace.has_value());
    if (boundaryFace) {
        const auto data=workspace.faces().get(*boundaryFace);
        CHECK(data.incident.size()==1);
        CHECK(data.cls.dim==-1);
        CHECK(data.cls.tag==-1);
    }
    CHECK(sameMesh(global,globalBefore));
    CHECK(algo::validateTetCandidate(candidate).valid);

    // The accepted APIs retain their exact established behavior: ordinary
    // import applies no protected metadata, while the all-or-nothing entry
    // rejects this deliberately mixed present/missing row set.
    auto ordinary=algo::importTetCandidate(
        candidate,selected,globals,classes,global,{3,81});
    CHECK(ordinary.imported());
    if (ordinary.imported()) {
        CHECK(ordinary.result->metrics.partialImportRowsCopied==0);
        CHECK(ordinary.result->metrics.partialImportGlobalSlotsInitialized==0);
        CHECK(ordinary.result->metrics.partialImportCandidateMappingVisits==0);
        CHECK(ordinary.result->metrics.partialImportEdgeValidationRows==0);
        CHECK(ordinary.result->metrics.partialImportFaceValidationRows==0);
        CHECK(ordinary.result->metrics.partialImportProtectedEndpointVisits==0);
        CHECK(ordinary.result->metrics.partialImportAuthorityPointValidations==0);
        CHECK(ordinary.result->metrics.partialImportDuplicateKeyComparisons==0);
        CHECK(ordinary.result->metrics.partialImportResolutionNodeVisits==0);
        CHECK(ordinary.result->metrics.partialImportEdgeResolutionLookups==0);
        CHECK(ordinary.result->metrics.partialImportFaceResolutionLookups==0);
        CHECK(ordinary.result->metrics.partialImportEdgeDataReads==0);
        CHECK(ordinary.result->metrics.partialImportFaceDataCopies==0);
        CHECK(ordinary.result->metrics.
                  partialImportFaceIncidenceRowsCopied==0);
        CHECK(ordinary.result->metrics.partialImportSortComparisons==0);
        std::vector<topo::NodeHandle> ordinaryByGlobal(
            static_cast<std::size_t>(global.nodeCount()));
        for (const auto& binding:ordinary.result->nodeBindings)
            ordinaryByGlobal[static_cast<std::size_t>(
                binding.existingGlobalIndex)]=binding.local;
        const auto edge=ordinary.result->workspace->edges().find(
            ordinary.result->workspace->nodes(),
            ordinaryByGlobal[static_cast<std::size_t>(interiorGlobal[0])],
            ordinaryByGlobal[static_cast<std::size_t>(interiorGlobal[1])]);
        CHECK(edge.has_value());
        CHECK(edge&&!ordinary.result->workspace->edges().get(*edge).
              requiredConstraint);
        CHECK(edge&&ordinary.result->workspace->edges().get(*edge).cls.dim==-1);
        CHECK(edge&&ordinary.result->workspace->edges().get(*edge).cls.tag==-1);
    }
    auto complete=algo::importTetCandidateWithProtectedBoundary(
        candidate,selected,globals,classes,global,{3,81},edgeRows,faceRows);
    CHECK(!complete.imported());
    CHECK(complete.result==nullptr);
    CHECK(complete.reason==algo::TetCandidateReason::InvalidImportSelection);

    auto duplicateEdges=edgeRows;
    duplicateEdges[1]=edgeRows[0];
    duplicateEdges[1].existingGlobalNodes={
        edgeRows[0].existingGlobalNodes[1],
        edgeRows[0].existingGlobalNodes[0]};
    auto duplicate=algo::tet_candidate_detail::
        importTetCandidateWithPartialProtectedBoundary(
        candidate,selected,globals,classes,global,{3,81},
        duplicateEdges,faceRows);
    CHECK(!duplicate.imported());
    CHECK(duplicate.reason==algo::TetCandidateReason::InvalidGlobalMapping);

    auto invalidFaces=faceRows;
    invalidFaces[1].classification={3,703};
    auto invalidClass=algo::tet_candidate_detail::
        importTetCandidateWithPartialProtectedBoundary(
        candidate,selected,globals,classes,global,{3,81},
        edgeRows,invalidFaces);
    CHECK(!invalidClass.imported());
    CHECK(invalidClass.reason==
          algo::TetCandidateReason::InvalidClassification);
    auto invalidEdges=edgeRows;
    invalidEdges[0].classification={0,701};
    auto invalidEdgeClass=algo::tet_candidate_detail::
        importTetCandidateWithPartialProtectedBoundary(
            candidate,selected,globals,classes,global,{3,81},
            invalidEdges,faceRows);
    CHECK(!invalidEdgeClass.imported());
    CHECK(invalidEdgeClass.reason==
          algo::TetCandidateReason::InvalidClassification);
    invalidFaces=faceRows;
    invalidFaces[1].existingGlobalNodes[2]=global.nodeCount();
    auto invalidIdentity=
        algo::tet_candidate_detail::
            importTetCandidateWithPartialProtectedBoundary(
            candidate,selected,globals,classes,global,{3,81},
            edgeRows,invalidFaces);
    CHECK(!invalidIdentity.imported());
    CHECK(invalidIdentity.reason==
          algo::TetCandidateReason::InvalidGlobalMapping);
    CHECK(sameMesh(global,globalBefore));
}

void candidatePermutationAndDifficultInputs() {
    TEST("ck3.tet_candidate_permutation_difficult");
    const auto referencePoints = difficultPoints();
    auto permutedPoints = referencePoints;
    std::rotate(permutedPoints.begin(), permutedPoints.begin() + 4,
                permutedPoints.end());
    std::reverse(permutedPoints.begin(), permutedPoints.end());
    const auto reference = algo::buildTetCandidate(referencePoints);
    const auto permuted = algo::buildTetCandidate(permutedPoints);
    CHECK(reference.built());
    CHECK(permuted.built());
    if (reference.built() && permuted.built())
        CHECK(sameGeometry(*reference.candidate, *permuted.candidate));

    auto duplicatedPoints = referencePoints;
    duplicatedPoints.push_back(referencePoints[4]);
    const auto duplicated = algo::buildTetCandidate(duplicatedPoints);
    CHECK(duplicated.built());
    if (duplicated.built()) {
        CHECK(duplicated.candidate->inputToPoint.size() ==
              duplicatedPoints.size());
        CHECK(duplicated.candidate->points.size() == referencePoints.size());
        CHECK(duplicated.candidate->inputToPoint.back() ==
              duplicated.candidate->inputToPoint[4]);
        CHECK(duplicated.candidate->metrics.construction.pointsDuplicate == 1);
        CHECK(algo::validateTetCandidate(*duplicated.candidate).valid);
    }

    std::vector<Vec3> flattened;
    std::mt19937 rng(778);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int i = 0; i < 80; ++i)
        flattened.push_back(
            {uniform(rng), uniform(rng), uniform(rng) * 1.0e-9});
    const auto flat = algo::buildTetCandidate(flattened);
    if (!flat.built())
        std::printf("  flattened candidate reason=%d\n",
                    static_cast<int>(flat.reason));
    CHECK(flat.built());
    if (flat.built()) {
        CHECK(flat.candidate->metrics.identityComplete);
        CHECK(algo::validateTetCandidate(*flat.candidate).valid);
    }

    std::vector<Vec3> lattice;
    for (int z = 0; z < 3; ++z)
        for (int y = 0; y < 3; ++y)
            for (int x = 0; x < 3; ++x)
                lattice.push_back({double(x), double(y), double(z)});
    const auto grid = algo::buildTetCandidate(lattice);
    CHECK(grid.built());
    if (grid.built()) {
        CHECK(grid.candidate->metrics.identityComplete);
        CHECK(algo::validateTetCandidate(*grid.candidate).valid);
    }
}

void malformedCandidatesAreRejected() {
    TEST("ck3.tet_candidate_malformed");
    const auto built = algo::buildTetCandidate(difficultPoints());
    CHECK(built.built());
    if (!built.built()) return;

    auto repeated = *built.candidate;
    repeated.cells[0].vertices[1] = repeated.cells[0].vertices[0];
    CHECK(algo::validateTetCandidate(repeated).reason ==
          algo::TetCandidateReason::RepeatedCellVertex);

    auto inverted = *built.candidate;
    std::swap(inverted.cells[0].vertices[0], inverted.cells[0].vertices[1]);
    CHECK(algo::validateTetCandidate(inverted).reason ==
          algo::TetCandidateReason::NonPositiveOrientation);

    auto adjacency = *built.candidate;
    bool changed = false;
    for (auto& cell : adjacency.cells)
        for (auto& neighbour : cell.neighbours)
            if (neighbour != kInvalidIndex) {
                neighbour = kInvalidIndex;
                changed = true;
                break;
            }
    CHECK(changed);
    CHECK(algo::validateTetCandidate(adjacency).reason ==
          algo::TetCandidateReason::AdjacencyMismatch);

    auto boundary = *built.candidate;
    boundary.boundaryFaces.pop_back();
    CHECK(algo::validateTetCandidate(boundary).reason ==
          algo::TetCandidateReason::BoundaryExtractionMismatch);

    std::vector<Index> selected{0, 0};
    std::vector<Index> globals(built.candidate->points.size());
    std::vector<mesh::Classification> classes(
        built.candidate->points.size(), mesh::Classification{3, 4});
    for (std::size_t i = 0; i < globals.size(); ++i)
        globals[i] = static_cast<Index>(i);
    const auto global = globalForCandidate(*built.candidate, classes);
    const auto globalBefore = global;
    const auto duplicateSelection = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {3, 4});
    CHECK(!duplicateSelection.imported());
    CHECK(duplicateSelection.result == nullptr);
    CHECK(duplicateSelection.reason ==
          algo::TetCandidateReason::InvalidImportSelection);
    CHECK(sameMesh(global, globalBefore));

    selected.clear();
    const auto emptySelection = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {3, 4});
    CHECK(!emptySelection.imported());
    CHECK(emptySelection.result == nullptr);
    CHECK(emptySelection.reason ==
          algo::TetCandidateReason::InvalidImportSelection);
    CHECK(sameMesh(global, globalBefore));

    selected = {0};
    const Index firstUsed = built.candidate->cells[0].vertices[0];
    const Index secondUsed = built.candidate->cells[0].vertices[1];
    globals[static_cast<std::size_t>(secondUsed)] =
        globals[static_cast<std::size_t>(firstUsed)];
    const auto duplicateGlobal = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {3, 4});
    CHECK(!duplicateGlobal.imported());
    CHECK(duplicateGlobal.result == nullptr);
    CHECK(duplicateGlobal.reason ==
          algo::TetCandidateReason::InvalidGlobalMapping);
    CHECK(sameMesh(global, globalBefore));

    for (std::size_t i = 0; i < globals.size(); ++i)
        globals[i] = static_cast<Index>(i);
    std::swap(globals[static_cast<std::size_t>(firstUsed)],
              globals[static_cast<std::size_t>(secondUsed)]);
    const auto wrongPermutation = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {3, 4});
    CHECK(!wrongPermutation.imported());
    CHECK(wrongPermutation.result == nullptr);
    CHECK(wrongPermutation.reason ==
          algo::TetCandidateReason::InvalidGlobalMapping);
    CHECK(sameMesh(global, globalBefore));

    for (std::size_t i = 0; i < globals.size(); ++i)
        globals[i] = static_cast<Index>(i);
    auto coordinateGlobal = globalForCandidate(*built.candidate, classes);
    const auto coordinateBefore = coordinateGlobal;
    const Vec3& original = coordinateGlobal.node(firstUsed);
    coordinateGlobal.setNode(firstUsed,
                             {original.x + 1.0, original.y, original.z});
    const auto changedCoordinateBefore = coordinateGlobal;
    const auto coordinateMismatch = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, coordinateGlobal,
        {3, 4});
    CHECK(!coordinateMismatch.imported());
    CHECK(coordinateMismatch.result == nullptr);
    CHECK(coordinateMismatch.reason ==
          algo::TetCandidateReason::InvalidGlobalMapping);
    CHECK(sameMesh(coordinateGlobal, changedCoordinateBefore));
    CHECK(!sameMesh(coordinateGlobal, coordinateBefore));

    auto mismatchedClasses = classes;
    mismatchedClasses[static_cast<std::size_t>(firstUsed)] = {2, 91};
    const auto classificationMismatch = algo::importTetCandidate(
        *built.candidate, selected, globals, mismatchedClasses, global,
        {3, 4});
    CHECK(!classificationMismatch.imported());
    CHECK(classificationMismatch.result == nullptr);
    CHECK(classificationMismatch.reason ==
          algo::TetCandidateReason::InvalidGlobalMapping);
    CHECK(sameMesh(global, globalBefore));

    auto invalidNodeClasses = classes;
    invalidNodeClasses[static_cast<std::size_t>(firstUsed)] = {4, 4};
    const auto invalidDimGlobal =
        globalForCandidate(*built.candidate, invalidNodeClasses);
    const auto invalidDimBefore = invalidDimGlobal;
    const auto invalidNodeDim = algo::importTetCandidate(
        *built.candidate, selected, globals, invalidNodeClasses,
        invalidDimGlobal, {3, 4});
    CHECK(!invalidNodeDim.imported());
    CHECK(invalidNodeDim.result == nullptr);
    CHECK(invalidNodeDim.reason ==
          algo::TetCandidateReason::InvalidClassification);
    CHECK(sameMesh(invalidDimGlobal, invalidDimBefore));

    invalidNodeClasses[static_cast<std::size_t>(firstUsed)] = {-1, 4};
    const auto negativeDimGlobal =
        globalForCandidate(*built.candidate, invalidNodeClasses);
    const auto negativeDimBefore = negativeDimGlobal;
    const auto negativeNodeDim = algo::importTetCandidate(
        *built.candidate, selected, globals, invalidNodeClasses,
        negativeDimGlobal, {3, 4});
    CHECK(!negativeNodeDim.imported());
    CHECK(negativeNodeDim.result == nullptr);
    CHECK(negativeNodeDim.reason ==
          algo::TetCandidateReason::InvalidClassification);
    CHECK(sameMesh(negativeDimGlobal, negativeDimBefore));

    invalidNodeClasses[static_cast<std::size_t>(firstUsed)] = {3, -1};
    const auto invalidTagGlobal =
        globalForCandidate(*built.candidate, invalidNodeClasses);
    const auto invalidTagBefore = invalidTagGlobal;
    const auto invalidNodeTag = algo::importTetCandidate(
        *built.candidate, selected, globals, invalidNodeClasses,
        invalidTagGlobal, {3, 4});
    CHECK(!invalidNodeTag.imported());
    CHECK(invalidNodeTag.result == nullptr);
    CHECK(invalidNodeTag.reason ==
          algo::TetCandidateReason::InvalidClassification);
    CHECK(sameMesh(invalidTagGlobal, invalidTagBefore));

    const auto invalidCellDim = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {2, 4});
    CHECK(!invalidCellDim.imported());
    CHECK(invalidCellDim.result == nullptr);
    CHECK(invalidCellDim.reason ==
          algo::TetCandidateReason::InvalidClassification);
    CHECK(sameMesh(global, globalBefore));

    const auto invalidCellTag = algo::importTetCandidate(
        *built.candidate, selected, globals, classes, global, {3, -1});
    CHECK(!invalidCellTag.imported());
    CHECK(invalidCellTag.result == nullptr);
    CHECK(invalidCellTag.reason ==
          algo::TetCandidateReason::InvalidClassification);
    CHECK(sameMesh(global, globalBefore));
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };

const Entry* tetHandoffTests(std::size_t& count) {
    static const Entry tests[] = {
        {"ck3.tet_candidate", candidateConstructionAndMetrics},
        {"ck3.tet_partial_protected_import",
         partialProtectedImportMarksOnlyExactPresent},
        {"ck3.tet_permutation", candidatePermutationAndDifficultInputs},
        {"ck3.tet_malformed", malformedCandidatesAreRejected},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
