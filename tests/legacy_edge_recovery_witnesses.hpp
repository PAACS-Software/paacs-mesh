// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Captured local cavities retained from the old project tests/test_t2b.cpp.
// These are local regression fixtures, never whole-file CAD coverage credit.
#define TEST(name) std::puts(name)
#define CHECK(test) require(bool(test),"legacy cavity witness: " #test)
void abcBoundaryEdgeCavityWitness() {
    TEST("v010 recovery ABC boundary-edge cavity witness");
    std::vector<bnd::TetFacetVertex> vertices{
        {{-4.445, 29.3630452561248, 11.43}, 448},
        {{-4.445, 29.3630452561248, 12.7}, 449},
        {{-3.6030094931986294e-14, 28.347045299999998, 11.43}, 503},
        {{-1.9035991561903372e-15, 32.4110452561248, 12.7}, 508},
        {{-1.7763568394002505e-15, 29.3630452561248, 12.7}, 509},
        {{-1.7581039740940114e-15, 29.3630452561248, 11.43}, 510}};
    std::vector<bnd::TetFacetCell> cells{
        {{{0, 2, 4, 1}}, 3061, {3, 0}},
        {{{0, 4, 3, 1}}, 3062, {3, 0}},
        {{{0, 5, 4, 2}}, 3067, {3, 0}},
        {{{5, 4, 3, 0}}, 3068, {3, 0}}};
    bnd::TetEdgeGeneralCavityInput input;
    input.vertices = vertices;
    input.removedCells = cells;
    input.target = bnd::canonicalTetFacetEdge(1, 5);
    input.region = {3, 0};
    const bnd::TetEdgeGeneralCavityLimits limits{
        64, 32, 20000, 4000000, 40000000};
    const auto solved = bnd::solveTetEdgeGeneralCavity(input, limits);
    std::printf("  ABC edge witness status=%d reason=%d cells=%zu work=%llu\n",
                static_cast<int>(solved.status),
                static_cast<int>(solved.geometryReason), solved.cells.size(),
                static_cast<unsigned long long>(solved.metrics.workUnits));
    CHECK(solved.solved());
    CHECK(solved.cells.size() == 5);
}

void abcNonStarOneSteinerSearchWitness() {
    TEST("v010 recovery ABC non-star one-Steiner search witness");
    std::vector<bnd::TetFacetVertex> vertices{
        {{21.02724708530608, 0.7791722519441544, 11.43}, 791},
        {{23.27511121546057, -3.114242629995042, 11.43}, 820},
        {{24.252975345615067, -0.4908277480558372, 10.16}, 835},
        {{26.95757476177063, -2.540776484029948, 12.7}, 874},
        {{27.00810664298024, 1.9309473684331344, 12.7}, 877},
        {{27.051785368460013, -1.1693120015739824, 11.43}, 878}};
    std::vector<bnd::TetFacetCell> cells{
        {{{2, 1, 4, 0}}, 5236, {3, 0}},
        {{{2, 3, 4, 1}}, 5420, {3, 0}},
        {{{2, 5, 4, 3}}, 5498, {3, 0}}};
    bnd::TetEdgeGeneralCavityInput input;
    input.vertices=vertices;
    input.removedCells=cells;
    input.target=bnd::canonicalTetFacetEdge(0,5);
    input.region={3,0};
    const bnd::TetEdgeGeneralCavityLimits limits{
        64,32,20000,4000000,40000000};
    const auto solved=bnd::solveTetEdgeGeneralCavityOneSteiner(
        input,limits,64);
    std::printf("  ABC non-star status=%d candidate=%llu cells=%zu work=%llu\n",
                static_cast<int>(solved.status),
                static_cast<unsigned long long>(
                    solved.selectedCandidateOrdinal),
                solved.cavity.cells.size(),
                static_cast<unsigned long long>(solved.cumulativeWorkUnits));
    CHECK(solved.solved());
    CHECK(solved.cavity.cells.size()==9);
}

void abcProtectedCandidateOneSteinerWitness() {
    TEST("v010 recovery ABC protected invalid-candidate one-Steiner witness");
    std::vector<bnd::TetFacetVertex> vertices{
        {{-8.759718738016526,26.959642150278892,11.43},392},
        {{-8.501959108476589,25.70764620683979,12.7},394},
        {{-7.439869086551087,22.525972906126572,6.35},404},
        {{-4.494454856011558,26.701428009107048,6.35},446},
        {{-4.434454866976732,27.998046133132192,11.43},452},
        {{-4.154663973216758,26.756403896037806,12.7},457}};
    std::vector<bnd::TetFacetCell> cells{
        {{{2,1,5,0}},2754,{3,0}},
        {{{2,0,4,3}},2756,{3,0}},
        {{{2,0,5,4}},2757,{3,0}}};
    const std::array<bnd::TetFacetEdgeKey,2> protectedEdges{
        bnd::canonicalTetFacetEdge(1,5),
        bnd::canonicalTetFacetEdge(3,5)};
    bnd::TetEdgeGeneralCavityInput input;
    input.vertices=vertices;
    input.removedCells=cells;
    input.target=bnd::canonicalTetFacetEdge(1,3);
    input.protectedEdges=protectedEdges;
    input.region={3,0};
    const bnd::TetEdgeGeneralCavityLimits limits{
        64,32,20000,4000000,40000000};
    const auto existing=bnd::solveTetEdgeGeneralCavity(input,limits);
    CHECK(existing.status==
          bnd::TetEdgeGeneralCavityStatus::InvalidCandidate);
    CHECK(existing.geometryReason==
          bnd::TetFacetGeometryReason::CrossedProtectedEdge);
    const auto solved=bnd::solveTetEdgeGeneralCavityOneSteiner(
        input,limits,96);
    std::printf("  ABC protected candidate status=%d candidate=%llu "
                "attempts=%llu cells=%zu work=%llu\n",
                static_cast<int>(solved.status),
                static_cast<unsigned long long>(
                    solved.selectedCandidateOrdinal),
                static_cast<unsigned long long>(solved.candidateAttempts),
                solved.cavity.cells.size(),
                static_cast<unsigned long long>(solved.cumulativeWorkUnits));
    CHECK(solved.status==
          bnd::TetEdgeGeneralCavityStatus::InvalidCandidate);
    CHECK(solved.geometryReason==
          bnd::TetFacetGeometryReason::CrossedProtectedEdge);
    CHECK(solved.candidateAttempts>1);
}

void abcFailedEdgeOneSteinerSearchWitness() {
    TEST("v010 recovery ABC failed edge one-Steiner search witness");
    std::vector<bnd::TetFacetVertex> vertices{
        {{-16.661975173406542,-22.933241388016565,11.43},282},
        {{-15.814986606702284,-26.307902157751823,12.7},294},
        {{-14.400549777230914,-24.416779954648945,10.16},304},
        {{-13.199474489895191,-23.641917336985557,11.43},329},
        {{-13.155509774610676,-23.666409629438412,12.7},331},
        {{-12.86928926188644,-25.257402303763275,11.43},335}};
    std::vector<bnd::TetFacetCell> cells{
        {{{0,2,4,1}},2162,{3,0}},
        {{{2,3,4,1}},2242,{3,0}},
        {{{5,1,3,2}},2243,{3,0}}};
    bnd::TetEdgeGeneralCavityInput input;
    input.vertices=vertices;
    input.removedCells=cells;
    input.target=bnd::canonicalTetFacetEdge(0,5);
    input.region={3,0};
    const bnd::TetEdgeGeneralCavityLimits limits{
        64,32,20000,4000000,40000000};
    const auto solved=bnd::solveTetEdgeGeneralCavityOneSteiner(
        input,limits,64);
    std::printf("  ABC failed edge status=%d candidate=%llu attempts=%llu "
                "cells=%zu work=%llu\n",static_cast<int>(solved.status),
                static_cast<unsigned long long>(
                    solved.selectedCandidateOrdinal),
                static_cast<unsigned long long>(solved.candidateAttempts),
                solved.cavity.cells.size(),
                static_cast<unsigned long long>(solved.cumulativeWorkUnits));
    CHECK(solved.solved());
    CHECK(solved.selectedCandidateOrdinal<64);
    CHECK(solved.cavity.cells.size()==9);
}

void abcCompletedPathOneSteinerWitness() {
    TEST("v010 recovery ABC completed dual path one-Steiner witness");
    std::vector<bnd::TetFacetVertex> vertices{
        {{-13.199474489895191,-23.641917336985557,11.43},329},
        {{-13.155509774610676,-23.666409629438412,12.7},331},
        {{-12.86928926188644,-25.257402303763275,11.43},335},
        {{-11.712850736856772,-23.85131703969102,10.16},345},
        {{-10.44285073685678,-22.581317039691026,11.43},356},
        {{-9.08072277183882,-25.508956390241842,12.7},379},
        {{-9.016219593816075,-25.531826538936443,11.43},381}};
    std::vector<bnd::TetFacetCell> cells{
        {{{2,0,1,3}},2407,{3,0}},
        {{{2,3,1,4}},2417,{3,0}},
        {{{2,4,5,3}},2454,{3,0}},
        {{{2,3,5,6}},2456,{3,0}}};
    bnd::TetEdgeGeneralCavityInput input;
    input.vertices=vertices;
    input.removedCells=cells;
    input.target=bnd::canonicalTetFacetEdge(0,6);
    input.region={3,0};
    const bnd::TetEdgeGeneralCavityLimits limits{
        64,32,20000,4000000,40000000};
    const auto solved=bnd::solveTetEdgeGeneralCavityOneSteiner(
        input,limits,256);
    std::printf("  ABC completed path status=%d candidate=%llu attempts=%llu "
                "cells=%zu work=%llu\n",static_cast<int>(solved.status),
                static_cast<unsigned long long>(
                    solved.selectedCandidateOrdinal),
                static_cast<unsigned long long>(solved.candidateAttempts),
                solved.cavity.cells.size(),
                static_cast<unsigned long long>(solved.cumulativeWorkUnits));
    CHECK(solved.status==bnd::TetEdgeGeneralCavityStatus::Unavailable);
    CHECK(solved.candidateAttempts==73);
}

#undef TEST
#undef CHECK
