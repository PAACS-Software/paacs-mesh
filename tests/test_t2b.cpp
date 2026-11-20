// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "test_harness.hpp"

#include "paacs/algo/BoundaryConformingMesher.hpp"
#include "paacs/bnd/TetEdgeGeneralCavity.hpp"
#include "paacs/bnd/TetFacetGeometry.hpp"
#include "paacs/predicates/Predicates.hpp"
#include "paacs/geometry/Kernel.hpp"
#include "paacs/sizing/SizeField.hpp"
#include "paacs/mesh/MeshAudit.hpp"
#include "workspace_digest.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

using namespace paacs;

namespace {

struct SevenVertexEdgeWitness {
    std::vector<bnd::TetFacetVertex> vertices;
    std::vector<bnd::TetFacetCell> cells;
    bnd::TetEdgeGeneralCavityInput input{};

    SevenVertexEdgeWitness() {
        vertices={
            {{1,-2.6000000000000001,-2},24},
            {{2,-4,-4},44},
            {{2,-4,-1.3333333333333335},45},
            {{2,-1.9000000000000004,-1},48},
            {{3,-4,-4},59},
            {{3,-4,-1.3333333333333335},60},
            {{3,-1.9000000000000004,-1},63}};
        cells={
            {{{0,4,5,1}},260,{3,0}},
            {{{0,2,6,3}},262,{3,0}},
            {{{0,5,6,2}},263,{3,0}},
            {{{4,5,6,0}},266,{3,0}}};
        input.vertices=vertices;
        input.removedCells=cells;
        input.target={1,3};
        input.region={3,0};
    }
};

struct SevenVertexOracleResult {
    std::vector<std::array<bnd::TetFacetVertexId,4>> candidates;
    std::vector<std::vector<std::array<bnd::TetFacetVertexId,4>>> covers;
    std::size_t cyclicCovers=0;
};

SevenVertexOracleResult exhaustiveSevenVertexOracle(
    const SevenVertexEdgeWitness& witness) {
    using Face=std::array<bnd::TetFacetVertexId,3>;
    struct Candidate {
        std::array<bnd::TetFacetVertexId,4> key{};
        std::array<std::pair<Face,int>,4> faces{};
    };
    const auto faceKey=[](Face face) {
        std::sort(face.begin(),face.end());
        return face;
    };
    const auto faceSense=[&](const Face& directed) {
        std::size_t minimum=0;
        const auto less=[&](auto a,auto b) {
            return std::tie(witness.vertices[a].sequence,a)<
                   std::tie(witness.vertices[b].sequence,b);
        };
        if(less(directed[1],directed[minimum]))minimum=1;
        if(less(directed[2],directed[minimum]))minimum=2;
        return less(directed[(minimum+1)%3],
                    directed[(minimum+2)%3])?+1:-1;
    };
    constexpr std::array<std::array<int,3>,4> localFaces{{
        {{1,2,3}},{{0,3,2}},{{0,1,3}},{{0,2,1}}}};
    const auto rows=[&](const bnd::TetFacetCell& cell) {
        std::array<std::pair<Face,int>,4> result{};
        for(std::size_t f=0;f<localFaces.size();++f) {
            const Face directed{{
                cell.vertices[static_cast<std::size_t>(localFaces[f][0])],
                cell.vertices[static_cast<std::size_t>(localFaces[f][1])],
                cell.vertices[static_cast<std::size_t>(localFaces[f][2])]}};
            result[f]={faceKey(directed),faceSense(directed)};
        }
        return result;
    };
    std::map<Face,std::vector<int>> sourceUses;
    for(const auto& cell:witness.cells)
        for(const auto& [face,sense]:rows(cell))sourceUses[face].push_back(sense);
    std::map<Face,int> boundary;
    for(const auto& [face,senses]:sourceUses)
        if(senses.size()==1)boundary.emplace(face,senses.front());

    std::vector<Candidate> candidates;
    for(std::uint32_t a=0;a<7;++a)
        for(std::uint32_t b=a+1;b<7;++b)
            for(std::uint32_t c=b+1;c<7;++c)
                for(std::uint32_t d=c+1;d<7;++d) {
                    bnd::TetFacetCell cell{{a,b,c,d},0,{3,0}};
                    const auto sign=predicates::orient3d(
                        witness.vertices[a].position,witness.vertices[b].position,
                        witness.vertices[c].position,witness.vertices[d].position);
                    if(sign==predicates::Sign::Zero)continue;
                    if(sign==predicates::Sign::Negative)
                        std::swap(cell.vertices[0],cell.vertices[1]);
                    candidates.push_back({{{a,b,c,d}},rows(cell)});
                }
    std::map<Face,std::vector<std::size_t>> incident;
    for(std::size_t i=0;i<candidates.size();++i)
        for(const auto& [face,sense]:candidates[i].faces) {
            (void)sense;
            incident[face].push_back(i);
        }
    struct Ledger { std::uint8_t count=0; int sense=0; };
    std::map<Face,Ledger> ledger;
    std::vector<std::size_t> selected;
    std::vector<unsigned char> used(candidates.size(),0);
    SevenVertexOracleResult result;
    for(const auto& row:candidates)result.candidates.push_back(row.key);
    const auto unresolved=[&]() -> std::optional<Face> {
        for(const auto& [face,sense]:boundary) {
            (void)sense;
            const auto found=ledger.find(face);
            if(found==ledger.end()||found->second.count==0)return face;
        }
        for(const auto& [face,row]:ledger)
            if(boundary.count(face)==0&&row.count==1)return face;
        return std::nullopt;
    };
    const auto cyclic=[&]() {
        std::map<std::uint32_t,std::set<std::uint32_t>> link;
        for(const auto index:selected) {
            const auto& key=candidates[index].key;
            if(std::find(key.begin(),key.end(),witness.input.target.lo)==key.end()||
               std::find(key.begin(),key.end(),witness.input.target.hi)==key.end())
                continue;
            std::array<std::uint32_t,2> opposite{};
            std::size_t n=0;
            for(const auto vertex:key)
                if(vertex!=witness.input.target.lo&&
                   vertex!=witness.input.target.hi)opposite[n++]=vertex;
            link[opposite[0]].insert(opposite[1]);
            link[opposite[1]].insert(opposite[0]);
        }
        if(link.size()<3)return false;
        for(const auto& [vertex,neighbours]:link) {
            (void)vertex;
            if(neighbours.size()!=2)return false;
        }
        std::set<std::uint32_t> reached{link.begin()->first};
        std::vector<std::uint32_t> pending{link.begin()->first};
        while(!pending.empty()) {
            const auto vertex=pending.back();
            pending.pop_back();
            for(const auto next:link[vertex])
                if(reached.insert(next).second)pending.push_back(next);
        }
        return reached.size()==link.size();
    };
    std::function<void()> visit=[&]() {
        const auto face=unresolved();
        if(!face) {
            std::vector<std::array<std::uint32_t,4>> cover;
            for(const auto index:selected)cover.push_back(candidates[index].key);
            std::sort(cover.begin(),cover.end());
            result.covers.push_back(std::move(cover));
            result.cyclicCovers+=cyclic()?1U:0U;
            return;
        }
        for(const auto index:incident[*face]) {
            if(used[index])continue;
            auto snapshot=ledger;
            bool valid=true;
            for(const auto& [key,sense]:candidates[index].faces) {
                const auto old=ledger.find(key);
                const auto boundaryRow=boundary.find(key);
                if(boundaryRow!=boundary.end()) {
                    if(old!=ledger.end()||sense!=boundaryRow->second) {
                        valid=false;
                        break;
                    }
                } else if(old!=ledger.end()&&
                          (old->second.count>=2||old->second.sense==sense)) {
                    valid=false;
                    break;
                }
                auto& row=ledger[key];
                if(row.count==0)row.sense=sense;
                ++row.count;
            }
            if(valid) {
                used[index]=1;
                selected.push_back(index);
                visit();
                selected.pop_back();
                used[index]=0;
            }
            ledger=std::move(snapshot);
        }
    };
    visit();
    std::sort(result.covers.begin(),result.covers.end());
    return result;
}

void generalEdgeCavityPlainFeasibility() {
    TEST("t2b1 plain seven-vertex solver finds canonical exact cavity");
    SevenVertexEdgeWitness witness;
    const bnd::TetEdgeGeneralCavityLimits limits{
        16,16,128,1000000,10000000};
    const auto solved=bnd::solveTetEdgeGeneralCavity(witness.input,limits);
    std::printf("  t2b1 status=%d reason=%d cells=%zu enum=%llu positive=%llu "
                "states=%llu work=%llu\n",
        static_cast<int>(solved.status),static_cast<int>(solved.geometryReason),
        solved.cells.size(),
        static_cast<unsigned long long>(solved.metrics.enumeratedTets),
        static_cast<unsigned long long>(solved.metrics.positiveTets),
        static_cast<unsigned long long>(solved.metrics.searchStates),
        static_cast<unsigned long long>(solved.metrics.workUnits));
    for(const auto& key:solved.orderedTetKeys)
        std::printf("    tet=(%u,%u,%u,%u)\n",key[0],key[1],key[2],key[3]);
    const auto oracle=exhaustiveSevenVertexOracle(witness);
    CHECK(solved.status==bnd::TetEdgeGeneralCavityStatus::Unavailable);
    CHECK(solved.metrics.runComplete);
    CHECK(solved.metrics.combinationCount==35);
    CHECK(solved.metrics.enumeratedTets==35);
    CHECK(solved.metrics.positiveTets==32);
    CHECK(solved.enumeratedTetKeys==oracle.candidates);
    CHECK(oracle.covers.size()==2);
    CHECK(oracle.cyclicCovers==0);
    const std::vector<std::vector<std::array<std::uint32_t,4>>> expected{{
        {{0,1,4,5},{0,2,3,6},{0,2,5,6},{0,4,5,6}}},
        {{{0,1,4,6},{0,1,5,6},{0,2,3,6},{0,2,5,6},{1,4,5,6}}}};
    CHECK(oracle.covers==expected);
}

void generalEdgeCavityPlainBudgetsAndDeterminism() {
    TEST("t2b1 plain solver pins every hard budget and no-solution order");
    SevenVertexEdgeWitness witness;
    const bnd::TetEdgeGeneralCavityLimits complete{
        16,16,128,1000000,10000000};
    const auto reference=bnd::solveTetEdgeGeneralCavity(witness.input,complete);
    CHECK(reference.status==bnd::TetEdgeGeneralCavityStatus::Unavailable);

    auto limits=complete;
    limits.maxCavityCells=3;
    CHECK(bnd::solveTetEdgeGeneralCavity(witness.input,limits).status==
          bnd::TetEdgeGeneralCavityStatus::BudgetReached);
    limits=complete;
    limits.maxCavityVertices=6;
    CHECK(bnd::solveTetEdgeGeneralCavity(witness.input,limits).status==
          bnd::TetEdgeGeneralCavityStatus::BudgetReached);
    limits=complete;
    limits.maxEnumeratedTets=34;
    const auto enumerationBudget=
        bnd::solveTetEdgeGeneralCavity(witness.input,limits);
    CHECK(enumerationBudget.status==
          bnd::TetEdgeGeneralCavityStatus::BudgetReached);
    CHECK(enumerationBudget.metrics.enumeratedTets==0);
    limits=complete;
    limits.maxSearchStates=1;
    CHECK(bnd::solveTetEdgeGeneralCavity(witness.input,limits).status==
          bnd::TetEdgeGeneralCavityStatus::BudgetReached);
    limits=complete;
    limits.maxWorkUnits=34;
    const auto workBudget=bnd::solveTetEdgeGeneralCavity(witness.input,limits);
    CHECK(workBudget.status==bnd::TetEdgeGeneralCavityStatus::BudgetReached);
    CHECK(workBudget.metrics.workUnits<=limits.maxWorkUnits);
    limits=complete;
    limits.maxCavityCells=0;
    CHECK(bnd::solveTetEdgeGeneralCavity(witness.input,limits).status==
          bnd::TetEdgeGeneralCavityStatus::InvalidInput);

    std::reverse(witness.cells.begin(),witness.cells.end());
    witness.input.removedCells=witness.cells;
    const auto reversed=bnd::solveTetEdgeGeneralCavity(witness.input,complete);
    CHECK(reversed.status==reference.status);
    CHECK(reversed.enumeratedTetKeys==reference.enumeratedTetKeys);
    CHECK(reversed.metrics.combinationCount==reference.metrics.combinationCount);
    CHECK(reversed.metrics.enumeratedTets==reference.metrics.enumeratedTets);
    CHECK(reversed.metrics.positiveTets==reference.metrics.positiveTets);
    CHECK(reversed.metrics.searchStates==reference.metrics.searchStates);
    CHECK(reversed.metrics.workUnits==reference.metrics.workUnits);
}

class Provider final : public algo::VolumeBoundaryOwnershipProvider {
public:
    algo::VolumeBoundaryProviderStamp current{77, 3};
    mutable int ownershipCalls = 0;
    mutable int envelopeCalls = 0;
    mutable int stampCalls = 0;
    int changeAtStampCall = -1;
    bool throwOwnership = false;
    bool allocateOwnershipFailure = false;
    bool wrongOwnershipStamp = false;
    bool wrongOwnershipFace = false;
    bool throwEnvelope = false;
    bool allocateEnvelopeFailure = false;
    std::map<algo::FacePersistentId, std::vector<bnd::VolumeRegionUse>> rows;
    algo::ProviderRecoveryEnvelope supplied{};

    algo::VolumeBoundaryProviderStamp stamp() const noexcept override {
        ++stampCalls;
        if (stampCalls == changeAtStampCall)
            return {current.persistentId,current.revision+1};
        return current;
    }
    algo::VolumeFaceOwnership ownership(
        algo::FacePersistentId face) const override {
        ++ownershipCalls;
        if (allocateOwnershipFailure) throw std::bad_alloc{};
        if (throwOwnership) throw std::runtime_error("ownership failure");
        algo::VolumeFaceOwnership row;
        row.providerStamp = wrongOwnershipStamp
            ? algo::VolumeBoundaryProviderStamp{
                  current.persistentId,current.revision+1}
            : current;
        row.facePersistentId = wrongOwnershipFace ? face+1 : face;
        const auto it = rows.find(face);
        if (it != rows.end()) row.owners = it->second;
        return row;
    }
    algo::ProviderRecoveryEnvelope recoveryEnvelope(
        std::span<const Vec3>) const override {
        ++envelopeCalls;
        if (allocateEnvelopeFailure) throw std::bad_alloc{};
        if (throwEnvelope) throw std::runtime_error("envelope failure");
        auto result = supplied;
        result.providerStamp = current;
        return result;
    }
};

class ReversedFaceModel final : public geo::Model {
public:
    explicit ReversedFaceModel(const geo::Model& source) : source_(source) {}
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        return source_.edges();
    }
    std::vector<const geo::GFace*> faces() const override {
        auto result = source_.faces();
        std::reverse(result.begin(),result.end());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override {
        return source_.edge(tag);
    }
    const geo::GFace* face(int tag) const override {
        return source_.face(tag);
    }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }
private:
    const geo::Model& source_;
};

class DuplicateFaceModel final : public geo::Model {
public:
    explicit DuplicateFaceModel(const geo::Model& source) : source_(source) {}
    std::vector<const geo::GVertex*> vertices() const override {
        return source_.vertices();
    }
    std::vector<const geo::GEdge*> edges() const override {
        return source_.edges();
    }
    std::vector<const geo::GFace*> faces() const override {
        auto result = source_.faces();
        if (!result.empty()) result.push_back(result.front());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override {
        return source_.solids();
    }
    const geo::GVertex* vertex(int tag) const override {
        return source_.vertex(tag);
    }
    const geo::GEdge* edge(int tag) const override { return source_.edge(tag); }
    const geo::GFace* face(int tag) const override { return source_.face(tag); }
    const geo::GSolid* solid(int tag) const override {
        return source_.solid(tag);
    }
private:
    const geo::Model& source_;
};

class EmptyModel final : public geo::Model {
public:
    std::vector<const geo::GVertex*> vertices() const override { return {}; }
    std::vector<const geo::GEdge*> edges() const override { return {}; }
    std::vector<const geo::GFace*> faces() const override { return {}; }
    std::vector<const geo::GSolid*> solids() const override { return {}; }
    const geo::GVertex* vertex(int) const override { return nullptr; }
    const geo::GEdge* edge(int) const override { return nullptr; }
    const geo::GFace* face(int) const override { return nullptr; }
    const geo::GSolid* solid(int) const override { return nullptr; }
};

class InvalidEvaluatedSize final : public sizing::SizeField {
public:
    double eval(const Vec3&) const override { return 0.0; }
};

class ShellEdge final : public geo::GEdge {
public:
    ShellEdge(int tag,const geo::GVertex* a,const geo::GVertex* b)
        : GEdge(tag),a_(a),b_(b) {}
    Vec3 pointAt(double t) const override {
        return a_->position()+(b_->position()-a_->position())*t;
    }
    Vec3 tangentAt(double) const override {
        return (b_->position()-a_->position()).normalized();
    }
    double length() const override {
        return (b_->position()-a_->position()).norm();
    }
    int startVertexTag() const override { return a_->tag(); }
    int endVertexTag() const override { return b_->tag(); }
private:
    const geo::GVertex* a_ = nullptr;
    const geo::GVertex* b_ = nullptr;
};

class ShellFace final : public geo::GFace {
public:
    ShellFace(int tag,Vec3 origin,Vec3 uAxis,Vec3 vAxis,
              std::vector<int> loop)
        : GFace(tag),origin_(origin),uAxis_(uAxis),vAxis_(vAxis),
          loop_(std::move(loop)) {}
    Vec3 pointAt(const Vec2& uv) const override {
        return origin_+uAxis_*uv.u+vAxis_*uv.v;
    }
    Vec2 paramAt(const Vec3& point) const override {
        const Vec3 d = point-origin_;
        return {d.dot(uAxis_)/uAxis_.norm2(),
                d.dot(vAxis_)/vAxis_.norm2()};
    }
    Vec3 normalAt(const Vec2&) const override {
        return uAxis_.cross(vAxis_).normalized();
    }
    void paramBounds(Vec2& lo,Vec2& hi) const override {
        lo = {0,0}; hi = {1,1};
    }
    std::vector<int> boundaryEdgeTags() const override { return loop_; }
    double curvatureAt(const Vec2&) const override { return 0.0; }
    Vec2 metricScale() const override {
        return {uAxis_.norm(),vAxis_.norm()};
    }
private:
    Vec3 origin_{},uAxis_{},vAxis_{};
    std::vector<int> loop_;
};

class SealedCavityModel final : public geo::Model {
public:
    explicit SealedCavityModel(bool disconnected = false) {
        if (disconnected) {
            addShell({-3,-1,-1},{-1,1,1},false);
            addShell({1,-1,-1},{3,1,1},false);
        } else {
            addShell({-2,-2,-2},{2,2,2},false);
            addShell({-0.5,-0.5,-0.5},{0.5,0.5,0.5},true);
        }
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> result;
        for (const auto& value : vertices_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> result;
        for (const auto& value : edges_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GFace*> faces() const override {
        std::vector<const geo::GFace*> result;
        for (const auto& value : faces_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override { return {}; }
    const geo::GVertex* vertex(int tag) const override {
        for (const auto& value : vertices_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for (const auto& value : edges_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        for (const auto& value : faces_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GSolid* solid(int) const override { return nullptr; }
private:
    void addShell(const Vec3& lo,const Vec3& hi,bool inward) {
        const std::size_t first = vertices_.size();
        for (std::size_t i = 0; i < 8; ++i)
            vertices_.push_back(std::make_unique<geo::GVertex>(
                static_cast<int>(vertices_.size()+1),
                Vec3{(i&1U)?hi.x:lo.x,(i&2U)?hi.y:lo.y,
                     (i&4U)?hi.z:lo.z}));
        const std::array<std::array<std::size_t,4>,6> outward{{
            {{0,4,6,2}},{{1,3,7,5}},{{0,1,5,4}},
            {{2,6,7,3}},{{0,2,3,1}},{{4,5,7,6}}}};
        std::map<std::array<std::size_t,2>,int> edgeTag;
        for (auto loop : outward) {
            if (inward) std::reverse(loop.begin()+1,loop.end());
            std::vector<int> faceEdges;
            for (std::size_t i = 0; i < 4; ++i) {
                const std::size_t a = first+loop[i];
                const std::size_t b = first+loop[(i+1)%4];
                auto key = std::array<std::size_t,2>{a,b};
                if (key[1] < key[0]) std::swap(key[0],key[1]);
                auto found = edgeTag.find(key);
                if (found == edgeTag.end()) {
                    const int tag = static_cast<int>(edges_.size()+1);
                    edges_.push_back(std::make_unique<ShellEdge>(
                        tag,vertices_[key[0]].get(),vertices_[key[1]].get()));
                    found = edgeTag.emplace(key,tag).first;
                }
                faceEdges.push_back(found->second);
            }
            const Vec3 p0 = vertices_[first+loop[0]]->position();
            const Vec3 p1 = vertices_[first+loop[1]]->position();
            const Vec3 p3 = vertices_[first+loop[3]]->position();
            faces_.push_back(std::make_unique<ShellFace>(
                static_cast<int>(faces_.size()+1),p0,p1-p0,p3-p0,
                std::move(faceEdges)));
        }
    }
    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<ShellEdge>> edges_;
    std::vector<std::unique_ptr<ShellFace>> faces_;
};

class TwoRegionInterfaceModel final : public geo::Model {
public:
    explicit TwoRegionInterfaceModel(int interfacePresentation = 0) {
        for (int xi = 0; xi < 3; ++xi)
            for (int z = 0; z < 2; ++z)
                for (int y = 0; y < 2; ++y)
                    vertices_.push_back(std::make_unique<geo::GVertex>(
                        static_cast<int>(vertices_.size()+1),
                        Vec3{static_cast<double>(xi-1),
                             static_cast<double>(2*y-1),
                             static_cast<double>(2*z-1)}));
        addFace({v(0,0,0),v(0,0,1),v(0,1,1),v(0,1,0)},{{0,1}});
        addFace({v(2,0,0),v(2,1,0),v(2,1,1),v(2,0,1)},{{1,1}});
        for (int interval = 0; interval < 2; ++interval) {
            const int region = interval;
            const int lo = interval,hi = interval+1;
            addFace({v(lo,0,0),v(hi,0,0),v(hi,0,1),v(lo,0,1)},
                    {{region,1}});
            addFace({v(lo,1,0),v(lo,1,1),v(hi,1,1),v(hi,1,0)},
                    {{region,1}});
            addFace({v(lo,0,0),v(lo,1,0),v(hi,1,0),v(hi,0,0)},
                    {{region,1}});
            addFace({v(lo,0,1),v(hi,0,1),v(hi,1,1),v(lo,1,1)},
                    {{region,1}});
        }
        // Directed +X: region 0 is on the negative side and region 1 on the
        // positive side under the accepted incidence-parity convention.
        std::array<std::size_t,4> interfaceLoop{
            v(1,0,0),v(1,1,0),v(1,1,1),v(1,0,1)};
        std::vector<bnd::VolumeRegionUse> interfaceOwners{{0,1},{1,-1}};
        if (interfacePresentation == 1)
            std::rotate(interfaceLoop.begin(),interfaceLoop.begin()+2,
                        interfaceLoop.end());
        else if (interfacePresentation == 2) {
            std::reverse(interfaceLoop.begin()+1,interfaceLoop.end());
            interfaceOwners = {{0,-1},{1,1}};
        }
        addFace(interfaceLoop,std::move(interfaceOwners));
    }
    const std::vector<bnd::VolumeRegionUse>& owners(int faceTag) const {
        return owners_.at(static_cast<std::size_t>(faceTag-1));
    }
    std::vector<const geo::GVertex*> vertices() const override {
        std::vector<const geo::GVertex*> result;
        for (const auto& value : vertices_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GEdge*> edges() const override {
        std::vector<const geo::GEdge*> result;
        for (const auto& value : edges_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GFace*> faces() const override {
        std::vector<const geo::GFace*> result;
        for (const auto& value : faces_) result.push_back(value.get());
        return result;
    }
    std::vector<const geo::GSolid*> solids() const override { return {}; }
    const geo::GVertex* vertex(int tag) const override {
        for (const auto& value : vertices_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GEdge* edge(int tag) const override {
        for (const auto& value : edges_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GFace* face(int tag) const override {
        for (const auto& value : faces_)
            if (value->tag() == tag) return value.get();
        return nullptr;
    }
    const geo::GSolid* solid(int) const override { return nullptr; }
private:
    std::size_t v(int x,int y,int z) const noexcept {
        return static_cast<std::size_t>(x*4+z*2+y);
    }
    void addFace(std::array<std::size_t,4> loop,
                 std::vector<bnd::VolumeRegionUse> owners) {
        std::vector<int> faceEdges;
        for (std::size_t i = 0; i < 4; ++i) {
            auto key = std::array<std::size_t,2>{loop[i],loop[(i+1)%4]};
            if (key[1] < key[0]) std::swap(key[0],key[1]);
            auto found = edgeTags_.find(key);
            if (found == edgeTags_.end()) {
                const int tag = static_cast<int>(edges_.size()+1);
                edges_.push_back(std::make_unique<ShellEdge>(
                    tag,vertices_[key[0]].get(),vertices_[key[1]].get()));
                found = edgeTags_.emplace(key,tag).first;
            }
            faceEdges.push_back(found->second);
        }
        const Vec3 p0 = vertices_[loop[0]]->position();
        const Vec3 p1 = vertices_[loop[1]]->position();
        const Vec3 p3 = vertices_[loop[3]]->position();
        faces_.push_back(std::make_unique<ShellFace>(
            static_cast<int>(faces_.size()+1),p0,p1-p0,p3-p0,
            std::move(faceEdges)));
        owners_.push_back(std::move(owners));
    }
    std::vector<std::unique_ptr<geo::GVertex>> vertices_;
    std::vector<std::unique_ptr<ShellEdge>> edges_;
    std::vector<std::unique_ptr<ShellFace>> faces_;
    std::vector<std::vector<bnd::VolumeRegionUse>> owners_;
    std::map<std::array<std::size_t,2>,int> edgeTags_;
};

algo::ProviderRecoveryEnvelope tetraEnvelope() {
    algo::ProviderRecoveryEnvelope out;
    out.availability = algo::VolumeProviderAvailability::Available;
    out.guardCoordinates = {{1,1,1},{-1,-1,1},{-1,1,-1},{1,-1,-1}};
    const auto sign = predicates::orient3d(
        out.guardCoordinates[0], out.guardCoordinates[1],
        out.guardCoordinates[2], out.guardCoordinates[3]);
    for (const auto& f : bnd::tet_facet_detail::kTetFaces) {
        std::array<std::size_t, 3> face{
            static_cast<std::size_t>(f[0]),
            static_cast<std::size_t>(f[1]),
            static_cast<std::size_t>(f[2])};
        if (sign == predicates::Sign::Negative) std::swap(face[1], face[2]);
        out.outwardHullFaces.push_back(face);
    }
    return out;
}

algo::ProviderRecoveryEnvelope boxEnvelope(double lo, double hi) {
    algo::ProviderRecoveryEnvelope out;
    out.availability = algo::VolumeProviderAvailability::Available;
    out.guardCoordinates.reserve(8);
    for (std::size_t i = 0; i < 8; ++i)
        out.guardCoordinates.push_back(
            {(i & 1U) ? hi : lo,
             (i & 2U) ? hi : lo,
             (i & 4U) ? hi : lo});
    out.outwardHullFaces = {
        {0,4,2},{4,6,2},{1,3,5},{3,7,5},
        {0,1,4},{1,5,4},{2,6,3},{6,7,3},
        {0,2,1},{2,3,1},{4,5,6},{5,7,6}};
    return out;
}

void configureOneRegionProvider(const geo::Model& model, Provider& provider) {
    const auto solids = model.solids();
    CHECK(solids.size() == 1);
    if (solids.size() != 1) return;
    Vec3 lo{}, hi{};
    solids.front()->bbox(lo,hi);
    const Vec3 center = (lo+hi)*0.5;
    const double span = std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});
    provider.supplied = boxEnvelope(
        std::min({lo.x,lo.y,lo.z})-span,
        std::max({hi.x,hi.y,hi.z})+span);
    for (const geo::GFace* face : model.faces()) {
        Vec2 uvLo{},uvHi{};
        face->paramBounds(uvLo,uvHi);
        const Vec2 uv = (uvLo+uvHi)*0.5;
        const int sense = face->normalAt(uv).dot(
            face->pointAt(uv)-center) > 0.0 ? 1 : -1;
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1] =
            {{0,sense}};
    }
}

bool negativeEnvelopeSides(const algo::TetRecoveryEnvelope& envelope,
                           const Vec3& point) {
    for (const auto& face : envelope.outwardHullFaces) {
        const Vec3& a = envelope.guardCoordinates[face[0]];
        const Vec3& b = envelope.guardCoordinates[face[1]];
        const Vec3& c = envelope.guardCoordinates[face[2]];
        std::size_t witness = envelope.guardCoordinates.size();
        for (std::size_t q = 0; q < envelope.guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            if (predicates::orient3d(a, b, c,
                                     envelope.guardCoordinates[q]) !=
                predicates::Sign::Zero) {
                witness = q;
                break;
            }
        }
        if (witness == envelope.guardCoordinates.size() ||
            predicates::orient3d(a, b, c,
                                 envelope.guardCoordinates[witness]) !=
                predicates::Sign::Negative ||
            predicates::orient3d(a, b, c, point) !=
                predicates::Sign::Negative)
            return false;
    }
    return true;
}

void derivedEnvelopeAllFaces() {
    TEST("t2b derived recovery envelope validates all twelve outward faces");
    Provider provider;
    const std::vector<Vec3> points{{0,0,0},{1,2,3},{-2,4,-1}};
    algo::TetRecoveryEnvelope oneUlp;
    const double inf = std::numeric_limits<double>::infinity();
    const Vec3 oldLo{std::nextafter(-2.0,-inf),
                     std::nextafter(0.0,-inf),
                     std::nextafter(-1.0,-inf)};
    const Vec3 oldHi{std::nextafter(1.0,inf),
                     std::nextafter(4.0,inf),
                     std::nextafter(3.0,inf)};
    for (std::size_t i = 0; i < 8; ++i)
        oneUlp.guardCoordinates.push_back(
            {(i&1U) ? oldHi.x : oldLo.x,
             (i&2U) ? oldHi.y : oldLo.y,
             (i&4U) ? oldHi.z : oldLo.z});
    oneUlp.outwardHullFaces = {
        {0,4,2},{4,6,2},{1,3,5},{3,7,5},
        {0,1,4},{1,5,4},{2,6,3},{6,7,3},
        {0,2,1},{2,3,1},{4,5,6},{5,7,6}};
    CHECK(!algo::boundary_conforming_detail::
              validateTetEnvelopeConstructorFrameForTesting(oneUlp,points));
    auto built = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::DerivedFiniteBox, points,
        provider.current, provider);
    CHECK(built.status == algo::BoundaryConformingV3Status::Published);
    CHECK(built.envelope != nullptr);
    CHECK(provider.envelopeCalls == 0);
    if (!built.envelope) return;
    const int selectedRung = algo::boundary_conforming_detail::
        derivedTetEnvelopeRungForTesting(points,*built.envelope);
    CHECK(selectedRung >= 0 && selectedRung <= 52);
    CHECK(algo::boundary_conforming_detail::
              validateTetEnvelopeConstructorFrameForTesting(
                  *built.envelope,points));
    std::printf("  derived rung=%d lo=(%.17g,%.17g,%.17g) "
                "hi=(%.17g,%.17g,%.17g)\n",
                selectedRung,
                built.envelope->guardCoordinates[0].x,
                built.envelope->guardCoordinates[0].y,
                built.envelope->guardCoordinates[0].z,
                built.envelope->guardCoordinates[7].x,
                built.envelope->guardCoordinates[7].y,
                built.envelope->guardCoordinates[7].z);
    const std::array<std::array<std::size_t,3>,12> expected{{
        {{0,4,2}},{{4,6,2}},{{1,3,5}},{{3,7,5}},
        {{0,1,4}},{{1,5,4}},{{2,6,3}},{{6,7,3}},
        {{0,2,1}},{{2,3,1}},{{4,5,6}},{{5,7,6}}}};
    CHECK(built.envelope->guardCoordinates.size() == 8);
    CHECK(built.envelope->outwardHullFaces.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(built.envelope->outwardHullFaces[i] == expected[i]);
        const auto& face = expected[i];
        const Vec3& a = built.envelope->guardCoordinates[face[0]];
        const Vec3& b = built.envelope->guardCoordinates[face[1]];
        const Vec3& c = built.envelope->guardCoordinates[face[2]];
        std::size_t witness = built.envelope->guardCoordinates.size();
        for (std::size_t q = 0;
             q < built.envelope->guardCoordinates.size(); ++q) {
            if (q == face[0] || q == face[1] || q == face[2]) continue;
            if (predicates::orient3d(
                    a, b, c, built.envelope->guardCoordinates[q]) !=
                predicates::Sign::Zero) {
                witness = q;
                break;
            }
        }
        CHECK(witness < built.envelope->guardCoordinates.size());
        if (witness < built.envelope->guardCoordinates.size())
            CHECK(predicates::orient3d(
                      a, b, c, built.envelope->guardCoordinates[witness]) ==
                  predicates::Sign::Negative);
        for (const Vec3& point : points)
            CHECK(predicates::orient3d(a, b, c, point) ==
                  predicates::Sign::Negative);
    }
}

void derivedEnvelopeRangeAndDeterminism() {
    TEST("t2b derived envelope ladder handles scale range and identity");
    Provider provider;
    const double denorm = std::numeric_limits<double>::denorm_min();
    const std::vector<std::vector<Vec3>> valid{
        {{9e8,-4e8,2e8},{9e8+2.0,-4e8+3.0,2e8+5.0}},
        {{0,0,0},{1e-120,0,0},{0,2e-120,0},{0,0,3e-120}},
        {{4,1,1},{4,2,3},{4,5,2}},
        {{0,0,0}},
        {{denorm,0,0},{2*denorm,denorm,0},{3*denorm,0,denorm}}};
    for (std::size_t caseIndex = 0; caseIndex < valid.size(); ++caseIndex) {
        const auto& points = valid[caseIndex];
        const int callsBefore = provider.envelopeCalls;
        auto built = algo::buildTetRecoveryEnvelope(
            algo::RecoveryEnvelopeMode::DerivedFiniteBox,points,
            provider.current,provider);
        if (built.status != algo::BoundaryConformingV3Status::Published)
            std::printf("  range case=%zu status=%d reason=%d\n",caseIndex,
                        static_cast<int>(built.status),
                        static_cast<int>(built.failure.reason));
        CHECK(built.status == algo::BoundaryConformingV3Status::Published);
        CHECK(built.envelope != nullptr);
        CHECK(provider.envelopeCalls == callsBefore);
        if (built.envelope)
            CHECK(algo::boundary_conforming_detail::
                      validateTetEnvelopeConstructorFrameForTesting(
                          *built.envelope,points));
    }

    std::vector<Vec3> ordered{{0,0,0},{1,2,3},{-2,4,-1}};
    auto a = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::DerivedFiniteBox,ordered,
        provider.current,provider);
    std::reverse(ordered.begin(),ordered.end());
    std::vector<std::uint64_t> allocatorHistory(4096,17);
    CHECK(!allocatorHistory.empty());
    auto b = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::DerivedFiniteBox,ordered,
        provider.current,provider);
    CHECK(a.envelope && b.envelope);
    if (a.envelope && b.envelope) {
        CHECK(a.envelope->guardCoordinates.size() ==
              b.envelope->guardCoordinates.size());
        bool sameCoordinates = a.envelope->guardCoordinates.size() ==
                               b.envelope->guardCoordinates.size();
        for (std::size_t i = 0;
             sameCoordinates && i < a.envelope->guardCoordinates.size();++i) {
            const Vec3& lhs = a.envelope->guardCoordinates[i];
            const Vec3& rhs = b.envelope->guardCoordinates[i];
            sameCoordinates = lhs.x == rhs.x && lhs.y == rhs.y &&
                              lhs.z == rhs.z;
        }
        CHECK(sameCoordinates);
        CHECK(a.envelope->outwardHullFaces == b.envelope->outwardHullFaces);
    }

    const std::vector<Vec3> duplicates{{0,0,0},{0,0,0}};
    auto collision = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::DerivedFiniteBox,duplicates,
        provider.current,provider);
    CHECK(collision.status == algo::BoundaryConformingV3Status::Unsupported);
    CHECK(collision.failure.reason ==
          algo::BoundaryConformingV3Reason::EnvelopeRangeUnrepresentable);
    CHECK(collision.failure.phase ==
          algo::BoundaryConformingV3Phase::RecoveryEnvelope);

    const double maximum = std::numeric_limits<double>::max();
    const std::vector<Vec3> nearLimit{{maximum,maximum,maximum}};
    auto exhausted = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::DerivedFiniteBox,nearLimit,
        provider.current,provider);
    CHECK(exhausted.status == algo::BoundaryConformingV3Status::Unsupported);
    CHECK(exhausted.failure.reason ==
          algo::BoundaryConformingV3Reason::EnvelopeRangeUnrepresentable);
}

void providerEnvelopeNormalizedCollision() {
    TEST("t2b provider envelope normalized collision fails post-envelope");
    Provider provider;
    const std::vector<Vec3> points{{0,0,0},{1,2,3},{-2,4,-1}};
    const double inf = std::numeric_limits<double>::infinity();
    provider.supplied = boxEnvelope(-1.0,1.0);
    const Vec3 lo{std::nextafter(-2.0,-inf),
                  std::nextafter(0.0,-inf),
                  std::nextafter(-1.0,-inf)};
    const Vec3 hi{std::nextafter(1.0,inf),
                  std::nextafter(4.0,inf),
                  std::nextafter(3.0,inf)};
    for (std::size_t i = 0; i < 8; ++i)
        provider.supplied.guardCoordinates[i] =
            {(i&1U) ? hi.x : lo.x,
             (i&2U) ? hi.y : lo.y,
             (i&4U) ? hi.z : lo.z};
    auto failed = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull,points,
        provider.current,provider);
    CHECK(failed.status == algo::BoundaryConformingV3Status::Failed);
    CHECK(failed.failure.reason ==
          algo::BoundaryConformingV3Reason::
              OwnershipProviderContractViolation);
    CHECK(failed.failure.phase ==
          algo::BoundaryConformingV3Phase::RecoveryEnvelope);
    CHECK(failed.failure.providerCheckpoint ==
          algo::VolumeProviderCheckpoint::PostEnvelope);
    CHECK(provider.envelopeCalls == 1);
}

void providerFailureCheckpointsAndInvariance() {
    TEST("t2b provider failures retain exact checkpoint and destination");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0,0,0},{1,1,1});
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.55);
    auto run = [&](Provider& provider,
                   algo::BoundaryConformingV3Request request) {
        mesh::Mesh destination;
        destination.addNode({-9,-8,-7},{0,991});
        const auto before = mesh::meshHash(destination);
        auto outcome = algo::meshModelBoundaryConforming(
            *model,spec,request,destination,provider);
        CHECK(!outcome.published());
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination) == before);
        CHECK(outcome.report != nullptr);
        CHECK(outcome.reportComplete);
        CHECK(outcome.report && outcome.report->reportComplete);
        return outcome;
    };
    auto requestFor = [&](const Provider& provider) {
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {910,1};
        request.expectedOwnershipProviderStamp = provider.current;
        request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
        return request;
    };

    {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        auto request = requestFor(provider);
        request.expectedSurfaceSetStamp = {};
        const auto outcome = run(provider,request);
        CHECK(outcome.status == algo::BoundaryConformingV3Status::Rejected);
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::InvalidRequest);
        CHECK(provider.stampCalls == 0);
        CHECK(provider.ownershipCalls == 0);
        CHECK(provider.envelopeCalls == 0);
    }
    {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        provider.throwOwnership = true;
        const auto outcome = run(provider,requestFor(provider));
        CHECK(outcome.status == algo::BoundaryConformingV3Status::Failed);
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::OwnershipProviderException);
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::Ownership);
        CHECK(provider.envelopeCalls == 0);
    }
    {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        provider.allocateOwnershipFailure = true;
        const auto outcome = run(provider,requestFor(provider));
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::AllocationFailed);
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::Ownership);
    }
    for (int wrong = 0; wrong < 2; ++wrong) {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        provider.wrongOwnershipStamp = wrong == 0;
        provider.wrongOwnershipFace = wrong == 1;
        const auto outcome = run(provider,requestFor(provider));
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::
                  OwnershipProviderContractViolation);
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::Ownership);
    }
    for (int allocation = 0; allocation < 2; ++allocation) {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        provider.allocateEnvelopeFailure = allocation == 0;
        provider.throwEnvelope = allocation == 1;
        auto request = requestFor(provider);
        request.envelopeMode =
            algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull;
        const auto outcome = run(provider,request);
        CHECK(outcome.failure.reason ==
              (allocation == 0
                   ? algo::BoundaryConformingV3Reason::AllocationFailed
                   : algo::BoundaryConformingV3Reason::
                         OwnershipProviderException));
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::RecoveryEnvelope);
        CHECK(provider.envelopeCalls == 1);
    }

    struct StampCase {
        int call = 0;
        algo::BoundaryConformingV3Status status{};
        algo::BoundaryConformingV3Phase phase{};
        algo::VolumeProviderCheckpoint checkpoint{};
    };
    const std::array<StampCase,7> stampCases{{
        {2,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::Ownership,
           algo::VolumeProviderCheckpoint::PreOwnershipTable},
        {3,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::Ownership,
           algo::VolumeProviderCheckpoint::PostOwnershipTable},
        {4,algo::BoundaryConformingV3Status::Rejected,
           algo::BoundaryConformingV3Phase::RecoveryEnvelope,
           algo::VolumeProviderCheckpoint::PreEnvelope},
        {5,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::RecoveryEnvelope,
           algo::VolumeProviderCheckpoint::PostEnvelope},
        {6,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::EdgeRouting,
           algo::VolumeProviderCheckpoint::PreEdgeRecovery},
        {7,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::DomainExtraction,
           algo::VolumeProviderCheckpoint::PreDomainExtraction},
        {8,algo::BoundaryConformingV3Status::Failed,
           algo::BoundaryConformingV3Phase::BoundaryAudit,
           algo::VolumeProviderCheckpoint::PreBoundaryAudit}}};
    for (const auto& expected : stampCases) {
        Provider provider;
        configureOneRegionProvider(*model,provider);
        provider.changeAtStampCall = expected.call;
        const auto outcome = run(provider,requestFor(provider));
        CHECK(outcome.status == expected.status);
        CHECK(outcome.failure.reason ==
              (expected.status == algo::BoundaryConformingV3Status::Rejected
                   ? algo::BoundaryConformingV3Reason::StaleOwnershipProvider
                   : algo::BoundaryConformingV3Reason::
                         OwnershipProviderStampChanged));
        CHECK(outcome.failure.phase == expected.phase);
        CHECK(outcome.failure.providerCheckpoint == expected.checkpoint);
    }
}

void suppliedEnvelopeContract() {
    TEST("t2b supplied envelope is finite bounded stamped and exact");
    Provider provider;
    provider.supplied = tetraEnvelope();
    const std::vector<Vec3> points{{0,0,0}};
    auto built = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull, points,
        provider.current, provider);
    CHECK(built.status == algo::BoundaryConformingV3Status::Published);
    CHECK(built.envelope != nullptr);
    CHECK(provider.envelopeCalls == 1);
    if (built.envelope)
        CHECK(negativeEnvelopeSides(*built.envelope, points[0]));
    provider.envelopeCalls = 0;
    auto overLimit = algo::buildTetRecoveryEnvelope(
        algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull, points,
        provider.current, provider, 3, 508);
    CHECK(overLimit.status == algo::BoundaryConformingV3Status::Failed);
    CHECK(overLimit.failure.reason ==
          algo::BoundaryConformingV3Reason::OwnershipProviderContractViolation);
    CHECK(provider.envelopeCalls == 1);

    auto requireContractFailure = [&](Provider& malformed) {
        const int before = malformed.envelopeCalls;
        const auto failed = algo::buildTetRecoveryEnvelope(
            algo::RecoveryEnvelopeMode::ProviderSuppliedConvexHull,points,
            malformed.current,malformed);
        CHECK(failed.status == algo::BoundaryConformingV3Status::Failed);
        CHECK(failed.failure.reason == algo::BoundaryConformingV3Reason::
                  OwnershipProviderContractViolation);
        CHECK(failed.failure.phase ==
              algo::BoundaryConformingV3Phase::RecoveryEnvelope);
        CHECK(failed.failure.providerCheckpoint ==
              algo::VolumeProviderCheckpoint::PostEnvelope);
        CHECK(malformed.envelopeCalls == before+1);
    };
    {
        Provider malformed;
        malformed.supplied = tetraEnvelope();
        std::swap(malformed.supplied.outwardHullFaces[0][1],
                  malformed.supplied.outwardHullFaces[0][2]);
        requireContractFailure(malformed);
    }
    {
        Provider malformed;
        malformed.supplied = tetraEnvelope();
        malformed.supplied.outwardHullFaces.pop_back();
        requireContractFailure(malformed);
    }
    {
        Provider malformed;
        malformed.supplied = tetraEnvelope();
        malformed.supplied.outwardHullFaces.back() =
            malformed.supplied.outwardHullFaces.front();
        requireContractFailure(malformed);
    }
    {
        Provider malformed;
        malformed.supplied = tetraEnvelope();
        malformed.supplied.guardCoordinates[0] = points.front();
        requireContractFailure(malformed);
    }
}

struct TetraSurface { mesh::Mesh mesh; algo::SurfaceSetResult surface; };

TetraSurface tetraSurface() {
    TetraSurface fixture;
    const std::array<Vec3,4> points{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    for (std::size_t i = 0; i < points.size(); ++i) {
        fixture.mesh.addNode(points[i], {0, static_cast<int>(i)});
        const algo::SharedCurveNodeKey key{500 + i, 1, 0};
        fixture.surface.sharedLedger.emplace(
            key, algo::SharedCurveNodeBinding{
                     key, static_cast<Index>(i), {0, static_cast<int>(i)}});
    }
    fixture.surface.setStamp = {900, 4};
    const auto orientation = predicates::orient3d(
        points[0], points[1], points[2], points[3]);
    std::map<std::array<Index,2>, mesh::Classification> edgeClasses;
    int edgeTag = 100;
    for (Index a = 0; a < 4; ++a)
        for (Index b = a + 1; b < 4; ++b)
            edgeClasses[{a,b}] = {1, edgeTag++};
    for (std::size_t fi = 0;
         fi < bnd::tet_facet_detail::kTetFaces.size(); ++fi) {
        auto face = std::make_unique<algo::SurfaceFaceResult>();
        face->faceStamp = {1000 + fi, 1};
        face->faceClass = {2, static_cast<int>(20 + fi)};
        std::array<Index,3> directed{
            bnd::tet_facet_detail::kTetFaces[fi][0],
            bnd::tet_facet_detail::kTetFaces[fi][1],
            bnd::tet_facet_detail::kTetFaces[fi][2]};
        if (orientation == predicates::Sign::Negative)
            std::swap(directed[1], directed[2]);
        for (Index global : directed) {
            face->localNodes.push_back(
                {fixture.mesh.node(global), fixture.mesh.nodeClass(global)});
            face->publishedGlobalByLocalNode.push_back(global);
        }
        face->localCells.push_back({{0,1,2}, face->faceClass});
        face->publishedGlobalCellByLocalCell.push_back(
            fixture.mesh.addCell(mesh::CellType::Tri3,
                                 {directed[0], directed[1], directed[2]},
                                 face->faceClass));
        for (std::size_t ei = 0; ei < 3; ++ei) {
            Index a = directed[ei], b = directed[(ei + 1) % 3];
            if (b < a) std::swap(a,b);
            algo::PersistentConstraintEdge edge;
            edge.classification = edgeClasses.at({a,b});
            edge.publishedGlobalEndpoints = std::array<Index,2>{a,b};
            face->persistentEdges.push_back(edge);
        }
        fixture.surface.faces.push_back(std::move(face));
    }
    return fixture;
}

void artifactMappingsAndStampPrecedence() {
    TEST("t2b artifact preserves exact mappings and checks stamp first");
    auto fixture = tetraSurface();
    Provider provider;
    for (std::size_t fi = 0; fi < 4; ++fi)
        provider.rows[1000 + fi] = {{7, 1}};
    auto built = algo::buildClosedTri3Artifact(
        fixture.surface, fixture.mesh, fixture.surface.setStamp,
        provider.current, provider);
    CHECK(built.status == algo::BoundaryConformingV3Status::Published);
    CHECK(built.artifact != nullptr);
    CHECK(provider.ownershipCalls == 4);
    if (built.artifact) {
        CHECK(built.artifact->facets.size() == 4);
        CHECK(built.artifact->edges.size() == 6);
        CHECK(built.artifact->regionTags == std::vector<int>{7});
        CHECK(built.metrics.allPairsFallbacks == 0);
        CHECK(built.metrics.treeNodes >= built.artifact->facets.size());
        CHECK(built.metrics.aabbTests > 0);
        CHECK(built.metrics.exactFacetFacetTests > 0);
        CHECK(built.metrics.rejectedContacts == 0);
        CHECK(built.metrics.runComplete);
    }
    provider.ownershipCalls = 0;
    auto stale = algo::buildClosedTri3Artifact(
        fixture.surface, fixture.mesh,
        {fixture.surface.setStamp.persistentId,
         fixture.surface.setStamp.revision + 1},
        provider.current, provider);
    CHECK(stale.status == algo::BoundaryConformingV3Status::Rejected);
    CHECK(stale.failure.reason ==
          algo::BoundaryConformingV3Reason::StaleSurfaceArtifact);
    CHECK(provider.ownershipCalls == 0);

    auto revised = tetraSurface();
    ++revised.surface.setStamp.revision;
    Provider revisedProvider;
    for (std::size_t fi = 0; fi < 4; ++fi)
        revisedProvider.rows[1000+fi] = {{7,1}};
    auto revisedBuilt = algo::buildClosedTri3Artifact(
        revised.surface,revised.mesh,revised.surface.setStamp,
        revisedProvider.current,revisedProvider);
    CHECK(revisedBuilt.artifact != nullptr);
    if (built.artifact && revisedBuilt.artifact) {
        CHECK(revisedBuilt.artifact->surfaceSetStamp ==
              revised.surface.setStamp);
        std::vector<bnd::VolumeFacetKey> originalKeys,revisedKeys;
        for (const auto& facet : built.artifact->facets)
            originalKeys.push_back(facet.key);
        for (const auto& facet : revisedBuilt.artifact->facets)
            revisedKeys.push_back(facet.key);
        CHECK(originalKeys == revisedKeys);
    }

    auto requireMappingFailure = [](TetraSurface broken,
                                    algo::BoundaryConformingV3Reason reason) {
        Provider localProvider;
        for (std::size_t fi = 0; fi < 4; ++fi)
            localProvider.rows[1000 + fi] = {{7,1}};
        const auto before = mesh::meshHash(broken.mesh);
        const auto failed = algo::buildClosedTri3Artifact(
            broken.surface,broken.mesh,broken.surface.setStamp,
            localProvider.current,localProvider);
        CHECK(failed.status == algo::BoundaryConformingV3Status::Rejected);
        CHECK(failed.failure.reason == reason);
        CHECK(failed.failure.phase ==
              algo::BoundaryConformingV3Phase::SurfaceMapping);
        CHECK(localProvider.ownershipCalls == 0);
        CHECK(localProvider.envelopeCalls == 0);
        CHECK(mesh::meshHash(broken.mesh) == before);
    };
    {
        auto broken = tetraSurface();
        broken.surface.faces.front()->publishedGlobalByLocalNode.pop_back();
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodeMappingIncomplete);
    }
    {
        auto broken = tetraSurface();
        auto& mapping =
            broken.surface.faces.front()->publishedGlobalByLocalNode;
        mapping[1] = mapping[0];
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodeMappingAlias);
    }
    {
        auto broken = tetraSurface();
        broken.surface.faces.front()->localNodes.front().coord.x += 0.25;
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodePayloadMismatch);
    }
    {
        auto broken = tetraSurface();
        broken.surface.faces.front()->publishedGlobalCellByLocalCell.clear();
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceCellMappingIncomplete);
    }
    {
        auto broken = tetraSurface();
        broken.surface.faces[1]->publishedGlobalCellByLocalCell[0] =
            broken.surface.faces[0]->publishedGlobalCellByLocalCell[0];
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceCellMappingAlias);
    }
    {
        auto broken = tetraSurface();
        std::swap(broken.surface.faces.front()->localCells.front().nodes[1],
                  broken.surface.faces.front()->localCells.front().nodes[2]);
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceCellPayloadMismatch);
    }
    {
        auto broken = tetraSurface();
        broken.surface.sharedLedger.clear();
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodeMappingAlias);
    }
    {
        auto broken = tetraSurface();
        auto& binding = broken.surface.sharedLedger.begin()->second;
        ++binding.key.nodeOrdinal;
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodePayloadMismatch);
    }
    {
        auto broken = tetraSurface();
        auto& binding = broken.surface.sharedLedger.begin()->second;
        ++binding.classification.tag;
        requireMappingFailure(
            std::move(broken),algo::BoundaryConformingV3Reason::
                                  SurfaceNodePayloadMismatch);
    }

    for (const algo::ArtifactStamp wrong :
         std::array<algo::ArtifactStamp,2>{{{0,0},{901,4}}}) {
        auto broken = tetraSurface();
        Provider localProvider;
        const auto failed = algo::buildClosedTri3Artifact(
            broken.surface,broken.mesh,wrong,localProvider.current,
            localProvider);
        CHECK(failed.status == algo::BoundaryConformingV3Status::Rejected);
        CHECK(failed.failure.reason ==
              algo::BoundaryConformingV3Reason::StaleSurfaceArtifact);
        CHECK(localProvider.stampCalls == 0);
        CHECK(localProvider.ownershipCalls == 0);
    }
}

void artifactBroadPhaseRejectsIntersectingComponents() {
    TEST("t2b artifact BVH rejects intersecting identity-distinct shells");
    auto fixture = tetraSurface();
    auto second = tetraSurface();
    const Index nodeOffset = fixture.mesh.nodeCount();
    const Vec3 shift{0.2,0.2,0.2};
    for (Index i = 0;i < second.mesh.nodeCount();++i)
        fixture.mesh.addNode(second.mesh.node(i)+shift,
                             {0,static_cast<int>(i+4)});
    for (const auto& [key,binding] : second.surface.sharedLedger) {
        auto newKey = key;
        newKey.curvePersistentId += 100;
        fixture.surface.sharedLedger.emplace(
            newKey,algo::SharedCurveNodeBinding{
                newKey,binding.existingGlobalIndex+nodeOffset,
                {0,binding.classification.tag+4}});
    }
    for (auto& face : second.surface.faces) {
        face->faceStamp.persistentId += 100;
        face->faceClass.tag += 100;
        for (auto& cell : face->localCells) cell.classification = face->faceClass;
        for (auto& node : face->localNodes) {
            node.coord = node.coord+shift;
            node.cls.tag += 4;
        }
        for (Index& global : face->publishedGlobalByLocalNode)
            global += nodeOffset;
        face->publishedGlobalCellByLocalCell.clear();
        const auto& v = face->publishedGlobalByLocalNode;
        face->publishedGlobalCellByLocalCell.push_back(
            fixture.mesh.addCell(mesh::CellType::Tri3,
                {v[0],v[1],v[2]},face->faceClass));
        for (auto& edge : face->persistentEdges)
            if (edge.publishedGlobalEndpoints) {
                (*edge.publishedGlobalEndpoints)[0] += nodeOffset;
                (*edge.publishedGlobalEndpoints)[1] += nodeOffset;
                edge.classification.tag += 100;
            }
        fixture.surface.faces.push_back(std::move(face));
    }
    Provider provider;
    for (std::size_t fi = 0;fi < 4;++fi) {
        provider.rows[1000+fi] = {{7,1}};
        provider.rows[1100+fi] = {{7,1}};
    }
    auto built = algo::buildClosedTri3Artifact(
        fixture.surface,fixture.mesh,fixture.surface.setStamp,
        provider.current,provider);
    if (built.failure.reason !=
        algo::BoundaryConformingV3Reason::ArtifactContactConflict)
        std::printf("  broad-phase diagnostic status=%d reason=%d phase=%d\n",
                    static_cast<int>(built.status),
                    static_cast<int>(built.failure.reason),
                    static_cast<int>(built.failure.phase));
    CHECK(built.status == algo::BoundaryConformingV3Status::Rejected);
    CHECK(built.failure.reason ==
          algo::BoundaryConformingV3Reason::ArtifactContactConflict);
    CHECK(built.metrics.allPairsFallbacks == 0);
    CHECK(built.metrics.rejectedContacts == 1);
    CHECK(built.metrics.aabbTests > 0);
}

void artifactShellTaxonomy() {
    TEST("t2b artifact shell audit separates closure orientation manifold");
    auto fixture = tetraSurface();
    Provider provider;
    for (std::size_t fi = 0; fi < 4; ++fi)
        provider.rows[1000+fi] = {{7,1}};
    auto built = algo::buildClosedTri3Artifact(
        fixture.surface,fixture.mesh,fixture.surface.setStamp,
        provider.current,provider);
    CHECK(built.artifact != nullptr);
    if (!built.artifact) return;
    CHECK(algo::boundary_conforming_detail::
              auditClosedTri3RegionShellsForTesting(*built.artifact) ==
          algo::BoundaryConformingV3Reason::None);
    {
        auto open = *built.artifact;
        open.facets.pop_back();
        CHECK(algo::boundary_conforming_detail::
                  auditClosedTri3RegionShellsForTesting(open) ==
              algo::BoundaryConformingV3Reason::NonClosedShellComponent);
    }
    {
        auto reversed = *built.artifact;
        std::swap(reversed.facets.front().directedGlobalNodes[1],
                  reversed.facets.front().directedGlobalNodes[2]);
        CHECK(algo::boundary_conforming_detail::
                  auditClosedTri3RegionShellsForTesting(reversed) ==
              algo::BoundaryConformingV3Reason::NonOrientableShellComponent);
    }
    {
        auto nonmanifold = *built.artifact;
        auto duplicate = nonmanifold.facets.front();
        duplicate.key.facePersistentId += 10000;
        nonmanifold.facets.push_back(std::move(duplicate));
        CHECK(algo::boundary_conforming_detail::
                  auditClosedTri3RegionShellsForTesting(nonmanifold) ==
              algo::BoundaryConformingV3Reason::NonManifoldShellComponent);
    }
}

void artifactBroadPhaseScaling() {
    TEST("t2b artifact AABB broad phase has counted nonquadratic growth");
    const std::array<std::size_t,3> counts{{384,1536,6144}};
    std::array<std::uint64_t,3> work{};
    for (std::size_t level = 0; level < counts.size(); ++level) {
        mesh::Mesh surface;
        bnd::ClosedTri3Artifact artifact;
        artifact.surfaceSetStamp = {930,1};
        artifact.facets.reserve(counts[level]);
        for (std::size_t i = 0; i < counts[level]; ++i) {
            const double x = static_cast<double>(i%128)*4.0;
            const double y = static_cast<double>(i/128)*4.0;
            const Index a = surface.addNode({x,y,0},{2,1});
            const Index b = surface.addNode({x+1,y,0},{2,1});
            const Index c = surface.addNode({x,y+1,0},{2,1});
            artifact.facets.push_back(
                {{static_cast<algo::FacePersistentId>(i+1),0},
                 {a,b,c},{2,1},{{0,1}}});
        }
        algo::ArtifactPreflightMetrics metrics;
        CHECK(algo::boundary_conforming_detail::
                  auditClosedTri3ContactsForTesting(
                      artifact,surface,metrics));
        CHECK(metrics.runComplete);
        CHECK(metrics.facetEntries == counts[level]);
        CHECK(metrics.allPairsFallbacks == 0);
        work[level] = metrics.axisSelectionVisits+metrics.sortComparisons+
                      metrics.nodePairVisits+
                      metrics.aabbTests+metrics.exactEdgeEdgeTests+
                      metrics.exactEdgeFacetTests+
                      metrics.exactFacetFacetTests;
        CHECK(work[level] > 0);
        if (level != 0) CHECK(work[level] <= 8*work[level-1]);
    }
}

std::array<topo::NodeHandle,4> positiveHandles(
    const topo::NodePool& nodes, std::array<topo::NodeHandle,4> value) {
    if (predicates::orient3d(nodes.get(value[0]).coord,
                             nodes.get(value[1]).coord,
                             nodes.get(value[2]).coord,
                             nodes.get(value[3]).coord) ==
        predicates::Sign::Negative)
        std::swap(value[0], value[1]);
    return value;
}

struct TargetEdgeFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::vector<topo::NodeHandle> nodes;
    std::vector<topo::CellHandle> cells;
};

TargetEdgeFixture makeTargetEdgeFixture(bool churn) {
    TargetEdgeFixture fixture;
    topo::NodePool nodes;
    topo::CellPool cells;
    topo::EdgePool edges;
    topo::FacePool faces;

    if (churn) {
        std::vector<topo::NodeHandle> junk;
        for (const Vec3& p : std::array<Vec3,4>{
                 Vec3{20,20,20},Vec3{21,20,20},
                 Vec3{20,21,20},Vec3{20,20,21}})
            junk.push_back(nodes.add(p,{3,99}));
        const auto junkCell = cells.add(
            topo::ElementType::Tet4,positiveHandles(
                nodes,{junk[0],junk[1],junk[2],junk[3]}),
            {3,99},nodes);
        (void)edges.insertCellEdges(nodes,cells,junkCell);
        (void)faces.insertCellFaces(nodes,cells,edges,junkCell);
        faces.removeCellFaces(cells,junkCell);
        edges.removeCellEdges(cells,junkCell);
        cells.remove(junkCell);
        for (const auto node : junk) nodes.remove(node);
    }

    for (const Vec3& p : std::array<Vec3,5>{
             Vec3{0,0,0},Vec3{4,0,0},Vec3{0,4,0},
             Vec3{1,1,3},Vec3{1,1,-3}})
        fixture.nodes.push_back(nodes.add(p,{3,9}));
    fixture.cells.push_back(cells.add(
        topo::ElementType::Tet4,positiveHandles(
            nodes,{fixture.nodes[3],fixture.nodes[0],fixture.nodes[1],
                   fixture.nodes[2]}),
        {3,9},nodes));
    fixture.cells.push_back(cells.add(
        topo::ElementType::Tet4,positiveHandles(
            nodes,{fixture.nodes[4],fixture.nodes[0],fixture.nodes[2],
                   fixture.nodes[1]}),
        {3,9},nodes));
    for (const auto cell : fixture.cells)
        (void)edges.insertCellEdges(nodes,cells,cell);
    for (const auto cell : fixture.cells)
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    fixture.workspace = std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    return fixture;
}

topo::TetEdgeCavitySpec targetEdgeSpec(
    const TargetEdgeFixture& fixture,mesh::Classification target,
    bool reverse=false) {
    topo::TetEdgeCavitySpec spec;
    spec.stamp=topo::WorkspaceStamp::of(*fixture.workspace);
    spec.removed=fixture.cells;
    spec.target=reverse
        ? std::array<topo::NodeHandle,2>{fixture.nodes[4],fixture.nodes[3]}
        : std::array<topo::NodeHandle,2>{fixture.nodes[3],fixture.nodes[4]};
    spec.targetClassification=target;
    spec.region={3,9};
    constexpr std::array<std::array<std::size_t,4>,3> added{{
        {{3,4,0,1}},{{3,4,1,2}},{{3,4,2,0}}}};
    for (const auto ids : added) {
        topo::TetCavityCell cell;
        cell.nodes=positiveHandles(fixture.workspace->nodes(),{
            fixture.nodes[ids[0]],fixture.nodes[ids[1]],
            fixture.nodes[ids[2]],fixture.nodes[ids[3]]});
        cell.cls={3,9};
        spec.added.push_back(cell);
    }
    return spec;
}

void targetDimension2EdgeTransaction() {
    TEST("t2b dimension-2 absent edge transaction and indexed delta");
    const auto run=[](bool reverse,bool churn) {
        auto fixture=makeTargetEdgeFixture(churn);
        auto indexed=bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
        CHECK(indexed.built());
        if (!indexed.context) return std::uint64_t{0};
        const auto spec=targetEdgeSpec(fixture,{2,812},reverse);
        const auto revision=fixture.workspace->topologyRevision();
        auto transaction=fixture.workspace->beginTetCavityTransaction();
        const auto& preview=transaction.stage(spec);
        CHECK(preview.cellsToRetire==2);
        CHECK(preview.cellsToCreate==3);
        CHECK(preview.target==spec.target);
        transaction.prepare();
        CHECK(transaction.prepared());
        const auto predicted=transaction.predictedTargetEdge();
        CHECK(predicted.isSpecified());
        const auto committed=transaction.commitEdge();
        CHECK(committed.targetEdge==predicted);
        CHECK(committed.resultCells.size()==3);
        CHECK(committed.revision==revision+1);
        const auto actual=fixture.workspace->edges().find(
            fixture.workspace->nodes(),fixture.nodes[3],fixture.nodes[4]);
        CHECK(actual.has_value());
        if (actual) {
            CHECK(*actual==predicted);
            const auto data=fixture.workspace->edges().get(*actual);
            CHECK(data.cls.dim==2&&data.cls.tag==812);
            CHECK(data.requiredConstraint);
            CHECK(data.useCount==3);
        }

        const auto updated=bnd::updateTetRecoveryContext(
            *indexed.context,*fixture.workspace,fixture.cells,
            committed.resultCells);
        CHECK(updated.updated());
        CHECK(indexed.context->currentFor(*fixture.workspace));
        CHECK(indexed.context->cells.size()==3);

        bnd::TetFacetSetState state;
        state.workspaceStamp=topo::WorkspaceStamp::of(*fixture.workspace);
        state.artifactStamp={8120,1};
        bnd::TetEdgeRecoveryRequest request;
        request.stamp=state.workspaceStamp;
        request.artifactStamp=state.artifactStamp;
        request.key=bnd::VolumeArtifactEdgeKey::of(3,4);
        request.vertices={fixture.nodes[3],fixture.nodes[4]};
        request.classification={2,812};
        request.region={3,9};
        const auto assessed=bnd::assessTetEdgeIndexed(
            request,state,*fixture.workspace,*indexed.context);
        CHECK(assessed.status==bnd::TetEdgeRecoveryStatus::AlreadyPresent);
        CHECK(assessed.edge==predicted);
        return fixture.workspace->semanticTopologyHashV2();
    };

    const auto forward=run(false,false);
    const auto reversedChurned=run(true,true);
    CHECK(forward!=0);
    CHECK(forward==reversedChurned);
}

void targetDimension2EdgeDomainRejections() {
    TEST("t2b dimension-2 edge transaction rejects invalid domains before prepare");
    struct Invalid {
        mesh::Classification target;
        mesh::Classification region;
        bool required;
    };
    constexpr std::array<Invalid,6> invalid{{
        {{0,812},{3,9},true},
        {{3,812},{3,9},true},
        {{2,-1},{3,9},true},
        {{2,812},{2,9},true},
        {{2,812},{3,-1},true},
        {{2,812},{3,9},false}}};
    for (const auto& row : invalid) {
        auto fixture=makeTargetEdgeFixture(false);
        auto spec=targetEdgeSpec(fixture,row.target);
        spec.region=row.region;
        spec.requiredConstraint=row.required;
        const auto revision=fixture.workspace->topologyRevision();
        const auto digest=paacs_tests::observableDigest(*fixture.workspace);
        {
            auto transaction=fixture.workspace->beginTetCavityTransaction();
            bool rejected=false;
            try {
                (void)transaction.stage(spec);
            } catch (const std::invalid_argument&) {
                rejected=true;
            }
            CHECK(rejected);
            CHECK(!transaction.prepared());
            transaction.abandon();
        }
        CHECK(fixture.workspace->topologyRevision()==revision);
        CHECK(paacs_tests::observableDigest(*fixture.workspace)==digest);
        CHECK(!fixture.workspace->edges().find(
            fixture.workspace->nodes(),fixture.nodes[3],fixture.nodes[4]));
    }
}

void targetEdgeTransaction() {
    TEST("t2b target-edge transaction predicts and atomically marks edge");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::vector<Vec3>{{0,0,0},{4,0,0},{0,4,0},
                                           {1,1,3},{1,1,-3}})
        n.push_back(nodes.add(p, {3,9}));
    topo::CellPool cells;
    std::vector<topo::CellHandle> current;
    current.push_back(cells.add(
        topo::ElementType::Tet4, positiveHandles(nodes, {n[3],n[0],n[1],n[2]}),
        {3,9}, nodes));
    current.push_back(cells.add(
        topo::ElementType::Tet4, positiveHandles(nodes, {n[4],n[0],n[2],n[1]}),
        {3,9}, nodes));
    topo::EdgePool edges;
    for (const auto cell : current)
        (void)edges.insertCellEdges(nodes, cells, cell);
    CHECK(!edges.find(nodes, n[3], n[4]).has_value());
    topo::FacePool faces;
    for (const auto cell : current)
        (void)faces.insertCellFaces(nodes, cells, edges, cell);
    topo::TopologyWorkspace workspace(std::move(nodes), std::move(cells),
                                      std::move(edges), std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(workspace, {3,9});
    CHECK(indexed.built());
    if (indexed.context) {
        CHECK(indexed.context->metrics.nodes == 5);
        CHECK(indexed.context->metrics.cells == 2);
        CHECK(indexed.context->metrics.runComplete);
        CHECK(indexed.context->index.covers(indexed.context->cells));
        CHECK(indexed.context->querySegment(
                  workspace.nodes().get(n[3]).coord,
                  workspace.nodes().get(n[4]).coord).size() == 2);
    }
    topo::TetEdgeCavitySpec spec;
    spec.stamp = topo::WorkspaceStamp::of(workspace);
    spec.removed = current;
    spec.target = {n[3], n[4]};
    spec.targetClassification = {1,55};
    spec.region = {3,9};
    const std::array<std::array<std::size_t,4>,3> added{{
        {{3,4,0,1}},{{3,4,1,2}},{{3,4,2,0}}}};
    for (auto ids : added) {
        topo::TetCavityCell cell;
        cell.nodes = positiveHandles(
            workspace.nodes(), {n[ids[0]], n[ids[1]], n[ids[2]], n[ids[3]]});
        cell.cls = {3,9};
        spec.added.push_back(cell);
    }
    const auto revision = workspace.topologyRevision();
    auto transaction = workspace.beginTetCavityTransaction();
    const auto& preview = transaction.stage(spec);
    CHECK(preview.cellsToRetire == 2);
    CHECK(preview.cellsToCreate == 3);
    transaction.prepare();
    const auto predicted = transaction.predictedTargetEdge();
    CHECK(predicted.isSpecified());
    const auto committed = transaction.commitEdge();
    CHECK(committed.targetEdge == predicted);
    CHECK(committed.resultCells.size() == 3);
    CHECK(committed.revision == revision + 1);
    if (indexed.context) CHECK(!indexed.context->currentFor(workspace));
    const auto actual = workspace.edges().find(workspace.nodes(), n[3], n[4]);
    CHECK(actual.has_value());
    if (actual) {
        CHECK(*actual == predicted);
        const auto data = workspace.edges().get(*actual);
        CHECK(data.requiredConstraint);
        CHECK(data.cls.dim == 1 && data.cls.tag == 55);
        CHECK(data.useCount == 3);
    }
}

void indexedQueryEquivalence() {
    TEST("t2b recovery index equals deterministic exact full AABB scan");
    topo::NodePool nodes;
    topo::CellPool cells;
    for (std::size_t i = 0; i < 64; ++i) {
        const Vec3 base{3.0*static_cast<double>(i%8),
                        3.0*static_cast<double>((i/8)%4),
                        3.0*static_cast<double>(i/32)};
        std::array<topo::NodeHandle,4> cellNodes{
            nodes.add(base,{3,9}),
            nodes.add(base+Vec3{1,0,0},{3,9}),
            nodes.add(base+Vec3{0,1,0},{3,9}),
            nodes.add(base+Vec3{0,0,1},{3,9})};
        (void)cells.add(topo::ElementType::Tet4,
                        positiveHandles(nodes,cellNodes),{3,9},nodes);
    }
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    topo::TopologyWorkspace workspace(std::move(nodes),std::move(cells),
                                      std::move(edges),std::move(faces));
    auto built = bnd::buildTetRecoveryContext(workspace,{3,9});
    CHECK(built.built());
    if (!built.context) return;
    auto overlaps = [](const bnd::TetAabb& a,const bnd::TetAabb& b) {
        return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x &&
               a.lo.y <= b.hi.y && b.lo.y <= a.hi.y &&
               a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
    };
    auto exact = [&](const bnd::TetAabb& query) {
        std::vector<topo::CellHandle> result;
        for (const auto& [handle,record] : built.context->cells)
            if (overlaps(record.bounds,query)) result.push_back(handle);
        std::sort(result.begin(),result.end(),[&](const auto a,const auto b) {
            return built.context->cells.at(a).sequence <
                   built.context->cells.at(b).sequence;
        });
        return result;
    };
    std::uint64_t state = 0x9e3779b97f4a7c15ULL;
    auto sample = [&]() {
        state ^= state << 7U;
        state ^= state >> 9U;
        return static_cast<double>(state%3000U)/100.0-2.0;
    };
    for (std::size_t q = 0; q < 64; ++q) {
        const Vec3 a{sample(),sample(),sample()};
        const Vec3 b{sample(),sample(),sample()};
        const bnd::TetAabb bounds{
            {std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)},
            {std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)}};
        CHECK(built.context->queryBounds(bounds) == exact(bounds));
        CHECK(built.context->querySegment(a,b) == exact(bounds));
        const Vec3 c{sample(),sample(),sample()};
        const bnd::TetAabb triBounds{
            {std::min({a.x,b.x,c.x}),std::min({a.y,b.y,c.y}),
             std::min({a.z,b.z,c.z})},
            {std::max({a.x,b.x,c.x}),std::max({a.y,b.y,c.y}),
             std::max({a.z,b.z,c.z})}};
        CHECK(built.context->queryTriangle({a,b,c}) == exact(triBounds));
    }
    CHECK(built.context->index.covers(built.context->cells));
}

void indexedEdgeRecovery() {
    TEST("t2b indexed edge recovery updates topology context and set state");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::vector<Vec3>{{0,0,0},{4,0,0},{0,4,0},
                                           {1,1,3},{1,1,-3}})
        n.push_back(nodes.add(p, {3,9}));
    topo::CellPool cells;
    std::vector<topo::CellHandle> current;
    current.push_back(cells.add(
        topo::ElementType::Tet4, positiveHandles(nodes, {n[3],n[0],n[1],n[2]}),
        {3,9}, nodes));
    current.push_back(cells.add(
        topo::ElementType::Tet4, positiveHandles(nodes, {n[4],n[0],n[2],n[1]}),
        {3,9}, nodes));
    topo::EdgePool edges;
    for (const auto cell : current)
        (void)edges.insertCellEdges(nodes, cells, cell);
    topo::FacePool faces;
    for (const auto cell : current)
        (void)faces.insertCellFaces(nodes, cells, edges, cell);
    topo::TopologyWorkspace workspace(std::move(nodes), std::move(cells),
                                      std::move(edges), std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(workspace, {3,9});
    CHECK(indexed.built());
    if (!indexed.context) return;
    bnd::TetFacetSetState state;
    state.workspaceStamp = topo::WorkspaceStamp::of(workspace);
    state.artifactStamp = {700,1};
    bnd::TetEdgeRecoveryRequest request;
    request.stamp = state.workspaceStamp;
    request.artifactStamp = state.artifactStamp;
    request.key = bnd::VolumeArtifactEdgeKey::of(30,40);
    request.vertices = {n[3],n[4]};
    request.classification = {1,55};
    request.region = {3,9};
    const auto recovered = bnd::recoverTetEdgeIndexed(
        request,state,workspace,*indexed.context);
    CHECK(recovered.status == bnd::TetEdgeRecoveryStatus::Recovered);
    CHECK(recovered.edge.isSpecified());
    CHECK(state.edgeBindings.count(request.key) == 1);
    CHECK(indexed.context->currentFor(workspace));
    CHECK(indexed.context->cells.size() == 3);
    const auto actual = workspace.edges().find(workspace.nodes(),n[3],n[4]);
    CHECK(actual.has_value() && actual == recovered.edge);
    if (actual) {
        const auto data = workspace.edges().get(*actual);
        CHECK(data.requiredConstraint);
        CHECK(data.cls.dim == 1 && data.cls.tag == 55);
    }
    request.stamp = topo::WorkspaceStamp::of(workspace);
    const auto retry = bnd::recoverTetEdgeIndexed(
        request,state,workspace,*indexed.context);
    CHECK(retry.status == bnd::TetEdgeRecoveryStatus::AlreadyPresent);
    CHECK(retry.edge == recovered.edge);
}

void indexedFacetRecovery() {
    TEST("t2b indexed facet recovery preserves protected edges and retries");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::vector<Vec3>{{0,0,0},{4,0,0},{0,4,0},
                                           {1,1,3},{1,1,-3}})
        n.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    std::vector<topo::CellHandle> current;
    const std::array<std::array<std::size_t,4>,3> initial{{
        {{3,4,0,1}},{{3,4,1,2}},{{3,4,2,0}}}};
    for (const auto ids : initial)
        current.push_back(cells.add(
            topo::ElementType::Tet4,
            positiveHandles(nodes,{n[ids[0]],n[ids[1]],n[ids[2]],n[ids[3]]}),
            {3,9},nodes));
    topo::EdgePool edges;
    for (const auto cell : current)
        (void)edges.insertCellEdges(nodes,cells,cell);
    for (const std::array<std::size_t,2> edge :
         std::array<std::array<std::size_t,2>,3>{{{{0,1}},{{1,2}},{{2,0}}}}) {
        const auto found = edges.find(nodes,n[edge[0]],n[edge[1]]);
        CHECK(found.has_value());
        if (found) {
            edges.setClassification(*found,{1,41});
            edges.setRequiredConstraint(*found,true);
        }
    }
    topo::FacePool faces;
    for (const auto cell : current)
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    topo::TopologyWorkspace workspace(std::move(nodes),std::move(cells),
                                      std::move(edges),std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(workspace,{3,9});
    CHECK(indexed.built());
    if (!indexed.context) return;
    bnd::TetFacetSetState state;
    state.workspaceStamp = topo::WorkspaceStamp::of(workspace);
    state.artifactStamp = {701,1};
    bnd::TetFacetIndexedRequest request;
    request.stamp = state.workspaceStamp;
    request.artifactStamp = state.artifactStamp;
    request.key = {800,0};
    request.vertices = {n[0],n[1],n[2]};
    request.classification = {2,77};
    request.region = {3,9};
    bnd::ClosedTri3Artifact artifact;
    artifact.surfaceSetStamp = state.artifactStamp;
    artifact.regionTags = {9};
    for (const std::array<Index,2> edge :
         std::array<std::array<Index,2>,3>{{{{0,1}},{{1,2}},{{2,0}}}})
        artifact.edges.push_back(
            {bnd::VolumeArtifactEdgeKey::of(edge[0],edge[1]),{1,41},{request.key}});
    artifact.facets.push_back(
        {request.key,{0,1,2},request.classification,{{9,1}}});
    const std::map<Index,topo::NodeHandle> localByGlobal{
        {0,n[0]},{1,n[1]},{2,n[2]}};
    const auto scheduled = bnd::recoverTetFacetSetIndexed(
        artifact,localByGlobal,{3,9},state,workspace,*indexed.context);
    CHECK(scheduled.status == bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(scheduled.reportComplete);
    CHECK(scheduled.processed.size() == 4);
    CHECK(scheduled.schedulerMetrics.routeAssessments <= 8);
    CHECK(scheduled.schedulerMetrics.selectedRebuilds <= 4);
    CHECK(scheduled.schedulerMetrics.fullStageRescans == 0);
    const auto recoveredFace = scheduled.processed.empty()
        ? topo::FaceHandle{}
        : scheduled.processed.back().facetOutcome.face;
    CHECK(!scheduled.processed.empty() &&
          scheduled.processed.back().facetOutcome.status ==
              bnd::TetFacetRecoveryStatus::Recovered);
    CHECK(recoveredFace.isSpecified());
    CHECK(indexed.context->currentFor(workspace));
    CHECK(indexed.context->cells.size() == 2);
    CHECK(state.bindings.count(request.key) == 1);
    if (recoveredFace.isSpecified()) {
        const auto data = workspace.faces().get(recoveredFace);
        CHECK(data.cls.dim == 2 && data.cls.tag == 77);
        CHECK(data.incident.size() == 2);
        if (data.incident.size() == 2)
            CHECK(data.incident[0].sense * data.incident[1].sense == -1);
    }
    request.stamp = topo::WorkspaceStamp::of(workspace);
    const auto retry = bnd::recoverTetFacetIndexed(
        request,state,workspace,*indexed.context);
    CHECK(retry.recovery.status ==
          bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(retry.recovery.face == recoveredFace);
}

void interactingFacetSetRecovery() {
    TEST("t2b set recovers one edge then two interacting absent facets");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::array<Vec3,10>{
             Vec3{-1,-1,0},Vec3{1,-1,0},Vec3{1,1,0},Vec3{-1,1,0},
             Vec3{0,0,2},Vec3{0,0,-2},
             Vec3{10,0,0},Vec3{11,0,0},Vec3{10,1,0},Vec3{10,0,1}})
        n.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    for (std::size_t i = 0;i < 4;++i)
        (void)cells.add(topo::ElementType::Tet4,
            positiveHandles(nodes,{n[4],n[5],n[i],n[(i+1)%4]}),
            {3,9},nodes);
    const auto outside = cells.add(
        topo::ElementType::Tet4,
        positiveHandles(nodes,{n[6],n[7],n[8],n[9]}),{3,9},nodes);
    const auto outsideData = cells.get(outside);
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    for (std::size_t i = 0;i < 4;++i) {
        const auto edge = edges.find(nodes,n[i],n[(i+1)%4]);
        CHECK(edge.has_value());
        if (edge) {
            edges.setClassification(*edge,{1,41});
            edges.setRequiredConstraint(*edge,true);
        }
    }
    CHECK(!edges.find(nodes,n[0],n[2]).has_value());
    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    topo::TopologyWorkspace workspace(std::move(nodes),std::move(cells),
                                      std::move(edges),std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(workspace,{3,9});
    CHECK(indexed.built());
    if (!indexed.context) return;
    bnd::TetFacetSetState state;
    state.workspaceStamp = topo::WorkspaceStamp::of(workspace);
    state.artifactStamp = {703,1};
    bnd::ClosedTri3Artifact artifact;
    artifact.surfaceSetStamp = state.artifactStamp;
    artifact.regionTags = {9};
    artifact.facets = {
        {{810,0},{0,1,2},{2,77},{{9,1}}},
        {{811,0},{0,2,3},{2,78},{{9,1}}}};
    std::map<bnd::VolumeArtifactEdgeKey,bnd::VolumeArtifactEdge> edgeRows;
    for (const auto& facet : artifact.facets)
        for (std::size_t i = 0;i < 3;++i) {
            const auto key = bnd::VolumeArtifactEdgeKey::of(
                facet.directedGlobalNodes[i],
                facet.directedGlobalNodes[(i+1)%3]);
            auto [it,inserted] = edgeRows.emplace(
                key,bnd::VolumeArtifactEdge{key,{1,41},{}});
            (void)inserted;
            it->second.incidentFacets.push_back(facet.key);
        }
    for (auto& [key,row] : edgeRows) {
        (void)key;
        std::sort(row.incidentFacets.begin(),row.incidentFacets.end());
        artifact.edges.push_back(std::move(row));
    }
    const std::map<Index,topo::NodeHandle> localByGlobal{
        {0,n[0]},{1,n[1]},{2,n[2]},{3,n[3]}};
    const auto recovered = bnd::recoverTetFacetSetIndexed(
        artifact,localByGlobal,{3,9},state,workspace,*indexed.context);
    if (recovered.status != bnd::TetFacetSetRecoveryStatus::Complete)
        std::printf("  interacting status=%d reason=%d rows=%zu edges=%zu "
                    "facets=%zu\n",static_cast<int>(recovered.status),
                    static_cast<int>(recovered.reason),
                    recovered.processed.size(),
                    recovered.unattemptedEdges.size(),
                    recovered.unattemptedFacets.size());
    CHECK(recovered.status == bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(recovered.reportComplete);
    std::size_t edgeRecovered = 0,facetsMarked = 0;
    for (const auto& row : recovered.processed) {
        edgeRecovered += row.kind == bnd::TetFacetSetWorkKind::Edge &&
            row.edgeOutcome.status == bnd::TetEdgeRecoveryStatus::Recovered;
        facetsMarked += row.kind == bnd::TetFacetSetWorkKind::Facet &&
            row.facetOutcome.status == bnd::TetFacetRecoveryStatus::Marked;
    }
    CHECK(edgeRecovered == 1);
    CHECK(facetsMarked == 2);
    CHECK(state.edgeBindings.size() == artifact.edges.size());
    CHECK(state.bindings.size() == 2);
    CHECK(workspace.cells().isCurrent(outside));
    if (workspace.cells().isCurrent(outside))
        CHECK(workspace.cells().get(outside).nodes == outsideData.nodes);
    const auto first = state.bindings.find({810,0});
    const auto second = state.bindings.find({811,0});
    CHECK(first != state.bindings.end() && second != state.bindings.end());
    if (first != state.bindings.end() && second != state.bindings.end()) {
        CHECK(first->second.recoveredFace != second->second.recoveredFace);
        CHECK(workspace.faces().isCurrent(first->second.recoveredFace));
        CHECK(workspace.faces().isCurrent(second->second.recoveredFace));
        CHECK(workspace.faces().get(first->second.recoveredFace).cls.tag == 77);
        CHECK(workspace.faces().get(second->second.recoveredFace).cls.tag == 78);
    }
    CHECK(indexed.context->currentFor(workspace));
    const auto retry = bnd::recoverTetFacetSetIndexed(
        artifact,localByGlobal,{3,9},state,workspace,*indexed.context);
    CHECK(retry.status == bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(retry.reportComplete);
    for (const auto& row : retry.processed)
        if (row.kind == bnd::TetFacetSetWorkKind::Facet)
            CHECK(row.facetOutcome.status ==
                  bnd::TetFacetRecoveryStatus::AlreadyPresent);
}

struct MixedPresentSchedulerFixture {
    std::unique_ptr<topo::TopologyWorkspace> workspace;
    std::unique_ptr<bnd::TetRecoveryContext> context;
    bnd::TetFacetSetState state;
    bnd::ClosedTri3Artifact artifact;
    std::map<Index,topo::NodeHandle> localByGlobal;
};

MixedPresentSchedulerFixture makeMixedPresentSchedulerFixture() {
    MixedPresentSchedulerFixture fixture;
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::array<Vec3,11>{
             Vec3{-1,-1,0},Vec3{1,-1,0},Vec3{1,1,0},Vec3{-1,1,0},
             Vec3{0,0,2},Vec3{0,0,-2},
             Vec3{10,0,0},Vec3{11,0,0},Vec3{10,1,0},
             Vec3{10,0,1},Vec3{10,0,-1}})
        n.push_back(nodes.add(p,{3,9}));

    topo::CellPool cells;
    for (std::size_t i=0;i<4;++i)
        (void)cells.add(topo::ElementType::Tet4,
            positiveHandles(nodes,{n[4],n[5],n[i],n[(i+1)%4]}),
            {3,9},nodes);
    (void)cells.add(topo::ElementType::Tet4,
        positiveHandles(nodes,{n[6],n[7],n[8],n[9]}),{3,9},nodes);
    (void)cells.add(topo::ElementType::Tet4,
        positiveHandles(nodes,{n[6],n[8],n[7],n[10]}),{3,9},nodes);

    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    for (std::size_t i=0;i<4;++i) {
        const auto edge=edges.find(nodes,n[i],n[(i+1)%4]);
        if (!edge) throw std::logic_error("missing ring edge");
        edges.setClassification(*edge,{1,41});
        edges.setRequiredConstraint(*edge,true);
    }
    for (const std::array<std::size_t,2> ids :
         std::array<std::array<std::size_t,2>,3>{
             {{{6,7}},{{7,8}},{{8,6}}}}) {
        const auto edge=edges.find(nodes,n[ids[0]],n[ids[1]]);
        if (!edge) throw std::logic_error("missing protected face edge");
        edges.setClassification(*edge,{2,79});
        edges.setRequiredConstraint(*edge,true);
    }

    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    const auto protectedFace=faces.find(
        nodes,std::array<topo::NodeHandle,3>{n[6],n[7],n[8]});
    if (!protectedFace) throw std::logic_error("missing protected face");
    faces.setClassification(*protectedFace,{2,79});

    fixture.workspace=std::make_unique<topo::TopologyWorkspace>(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed=bnd::buildTetRecoveryContext(*fixture.workspace,{3,9});
    if (!indexed.built() || !indexed.context)
        throw std::logic_error("mixed scheduler context build failed");
    fixture.context=std::move(indexed.context);
    fixture.state.workspaceStamp=topo::WorkspaceStamp::of(*fixture.workspace);
    fixture.state.artifactStamp={704,1};
    fixture.artifact.surfaceSetStamp=fixture.state.artifactStamp;
    fixture.artifact.regionTags={9};
    fixture.artifact.facets={
        {{809,0},{6,7,8},{2,79},{{9,1}}},
        {{810,0},{0,1,2},{2,77},{{9,1}}},
        {{811,0},{0,2,3},{2,78},{{9,1}}}};
    std::map<bnd::VolumeArtifactEdgeKey,bnd::VolumeArtifactEdge> edgeRows;
    for (const auto& facet : fixture.artifact.facets)
        for (std::size_t i=0;i<3;++i) {
            const auto key=bnd::VolumeArtifactEdgeKey::of(
                facet.directedGlobalNodes[i],
                facet.directedGlobalNodes[(i+1)%3]);
            const mesh::Classification classification=
                key.canonicalGlobalNodes[0]>=6
                    ? mesh::Classification{2,79}
                    : mesh::Classification{1,41};
            auto [it,inserted]=edgeRows.emplace(
                key,bnd::VolumeArtifactEdge{key,classification,{}});
            if (!inserted &&
                (it->second.classification.dim!=classification.dim ||
                 it->second.classification.tag!=classification.tag))
                throw std::logic_error("mixed scheduler edge class conflict");
            it->second.incidentFacets.push_back(facet.key);
        }
    for (auto& [key,row] : edgeRows) {
        (void)key;
        std::sort(row.incidentFacets.begin(),row.incidentFacets.end());
        fixture.artifact.edges.push_back(std::move(row));
    }
    for (std::size_t i=0;i<n.size();++i)
        fixture.localByGlobal.emplace(static_cast<Index>(i),n[i]);
    return fixture;
}

template<class Handle>
bool sameHandleWithoutPool(Handle a,Handle b) {
    return a.slot==b.slot && a.gen==b.gen;
}

void requireMixedSchedulerReportEquivalent(
    const bnd::TetFacetSetRecoveryOutcome& oracle,
    const bnd::TetFacetSetRecoveryOutcome& optimized) {
    CHECK(optimized.status==oracle.status);
    CHECK(optimized.reason==oracle.reason);
    CHECK(optimized.reportComplete==oracle.reportComplete);
    CHECK(optimized.processed.size()==oracle.processed.size());
    // The accepted scheduler lets exact-present rows participate in cavity
    // dependencies, whereas the opt-in policy is required to remove them
    // from dependency testing.  Raw execution order can therefore differ.
    // Compare every complete report row exactly by its authoritative work
    // key; canonical-prefix order is pinned separately below.
    for (const auto& a:oracle.processed) {
        const auto matches=static_cast<std::size_t>(std::count_if(
            optimized.processed.begin(),optimized.processed.end(),
            [&](const auto& row) {
                return row.kind==a.kind&&
                    (a.kind==bnd::TetFacetSetWorkKind::Edge
                        ? row.edgeKey==a.edgeKey
                        : row.facetKey==a.facetKey);
            }));
        CHECK(matches==1);
        const auto found=std::find_if(
            optimized.processed.begin(),optimized.processed.end(),
            [&](const auto& row) {
                return row.kind==a.kind&&
                    (a.kind==bnd::TetFacetSetWorkKind::Edge
                        ? row.edgeKey==a.edgeKey
                        : row.facetKey==a.facetKey);
            });
        if(found==optimized.processed.end())continue;
        const auto& b=*found;
        CHECK(b.kind==a.kind);
        CHECK(b.edgeKey==a.edgeKey);
        CHECK(b.facetKey==a.facetKey);
        CHECK(b.edgeOutcome.status==a.edgeOutcome.status);
        CHECK(b.edgeOutcome.reason==a.edgeOutcome.reason);
        CHECK(b.edgeOutcome.phase==a.edgeOutcome.phase);
        CHECK(b.edgeOutcome.geometryReason==a.edgeOutcome.geometryReason);
        CHECK(b.edgeOutcome.counters==a.edgeOutcome.counters);
        CHECK(b.edgeOutcome.capacityResiduePermitted==
              a.edgeOutcome.capacityResiduePermitted);
        CHECK(sameHandleWithoutPool(b.edgeOutcome.edge,a.edgeOutcome.edge));
        CHECK(b.facetOutcome.status==a.facetOutcome.status);
        CHECK(b.facetOutcome.reason==a.facetOutcome.reason);
        CHECK(b.facetOutcome.phase==a.facetOutcome.phase);
        CHECK(b.facetOutcome.geometryReason==a.facetOutcome.geometryReason);
        CHECK(b.facetOutcome.counters==a.facetOutcome.counters);
        CHECK(b.facetOutcome.capacityResiduePermitted==
              a.facetOutcome.capacityResiduePermitted);
        CHECK(sameHandleWithoutPool(b.facetOutcome.face,a.facetOutcome.face));
    }
    CHECK(optimized.unattemptedEdges==oracle.unattemptedEdges);
    CHECK(optimized.unattemptedFacets==oracle.unattemptedFacets);
    CHECK(optimized.conflictEdgeVertices==oracle.conflictEdgeVertices);
    CHECK(optimized.conflictEdgeDependencies==oracle.conflictEdgeDependencies);
    CHECK(optimized.conflictVertices==oracle.conflictVertices);
    CHECK(optimized.conflictEdges==oracle.conflictEdges);
}

void requireMixedSchedulerStateEquivalent(
    const bnd::TetFacetSetState& oracle,
    const bnd::TetFacetSetState& optimized) {
    CHECK(optimized.artifactStamp==oracle.artifactStamp);
    CHECK(optimized.workspaceStamp.revision==oracle.workspaceStamp.revision);
    CHECK(optimized.canonicalEdgeSchedule==oracle.canonicalEdgeSchedule);
    CHECK(optimized.canonicalSchedule==oracle.canonicalSchedule);
    CHECK(optimized.edgeBindings.size()==oracle.edgeBindings.size());
    auto oa=oracle.edgeBindings.begin();
    auto ob=optimized.edgeBindings.begin();
    for (;oa!=oracle.edgeBindings.end()&&ob!=optimized.edgeBindings.end();
         ++oa,++ob) {
        CHECK(ob->first==oa->first);
        CHECK(ob->second.key==oa->second.key);
        CHECK(ob->second.classification.dim==oa->second.classification.dim);
        CHECK(ob->second.classification.tag==oa->second.classification.tag);
        CHECK(ob->second.required==oa->second.required);
        CHECK(sameHandleWithoutPool(ob->second.edge,oa->second.edge));
    }
    CHECK(optimized.bindings.size()==oracle.bindings.size());
    auto fa=oracle.bindings.begin();
    auto fb=optimized.bindings.begin();
    for (;fa!=oracle.bindings.end()&&fb!=optimized.bindings.end();++fa,++fb) {
        CHECK(fb->first==fa->first);
        CHECK(fb->second.key==fa->second.key);
        CHECK(fb->second.classification.dim==fa->second.classification.dim);
        CHECK(fb->second.classification.tag==fa->second.classification.tag);
        for (std::size_t i=0;i<3;++i)
            CHECK(sameHandleWithoutPool(
                fb->second.canonicalVertices[i],
                fa->second.canonicalVertices[i]));
        CHECK(sameHandleWithoutPool(
            fb->second.recoveredFace,fa->second.recoveredFace));
    }
}

void requireMixedSchedulerOrderedReportEquivalent(
    const bnd::TetFacetSetRecoveryOutcome& oracle,
    const bnd::TetFacetSetRecoveryOutcome& optimized) {
    requireMixedSchedulerReportEquivalent(oracle,optimized);
    CHECK(optimized.processed.size()==oracle.processed.size());
    const auto count=std::min(
        optimized.processed.size(),oracle.processed.size());
    for(std::size_t i=0;i<count;++i) {
        CHECK(optimized.processed[i].kind==oracle.processed[i].kind);
        if(oracle.processed[i].kind==bnd::TetFacetSetWorkKind::Edge)
            CHECK(optimized.processed[i].edgeKey==oracle.processed[i].edgeKey);
        else
            CHECK(optimized.processed[i].facetKey==oracle.processed[i].facetKey);
    }
}

void mixedPresentSchedulerEquivalence() {
    TEST("t2b mixed exact-present batching matches unoptimized scheduler");
    auto oracle=makeMixedPresentSchedulerFixture();
    const auto oracleOutcome=bnd::recoverTetFacetSetIndexed(
        oracle.artifact,oracle.localByGlobal,{3,9},oracle.state,
        *oracle.workspace,*oracle.context);
    auto optimized=makeMixedPresentSchedulerFixture();
    bnd::TetFacetSetMixedPresentBatchSink sink;
    const auto optimizedOutcome=bnd::recoverTetFacetSetIndexed(
        optimized.artifact,optimized.localByGlobal,{3,9},optimized.state,
        *optimized.workspace,*optimized.context,nullptr,nullptr,&sink);
#if defined(PAACS_TOPO_TESTING)
    auto capOne=makeMixedPresentSchedulerFixture();
    bnd::TetFacetSetMixedPresentBatchSink capOneSink;
    capOneSink.maximumBatchEntitiesForTesting=1;
    const auto capOneOutcome=bnd::recoverTetFacetSetIndexed(
        capOne.artifact,capOne.localByGlobal,{3,9},capOne.state,
        *capOne.workspace,*capOne.context,nullptr,nullptr,&capOneSink);
#endif

    CHECK(oracleOutcome.status==bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(optimizedOutcome.status==bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(sink.metrics.complete);
    requireMixedSchedulerReportEquivalent(oracleOutcome,optimizedOutcome);
    requireMixedSchedulerStateEquivalent(oracle.state,optimized.state);
    CHECK(optimized.workspace->semanticTopologyHashV2()==
          oracle.workspace->semanticTopologyHashV2());
#if defined(PAACS_TOPO_TESTING)
    CHECK(capOneOutcome.status==bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(capOneSink.metrics.complete);
    requireMixedSchedulerOrderedReportEquivalent(
        capOneOutcome,optimizedOutcome);
    requireMixedSchedulerStateEquivalent(capOne.state,optimized.state);
    CHECK(capOne.workspace->semanticTopologyHashV2()==
          optimized.workspace->semanticTopologyHashV2());
    CHECK(capOneSink.metrics.exactPresentSelectedRebuilds==0);
    CHECK(capOneSink.metrics.exactPresentDependencyTests==0);
#endif
    CHECK(optimizedOutcome.processed.size()>=2);
    if(optimizedOutcome.processed.size()>=2) {
        CHECK(optimizedOutcome.processed[0].kind==
              bnd::TetFacetSetWorkKind::Edge);
        CHECK(optimizedOutcome.processed[0].edgeKey==
              bnd::VolumeArtifactEdgeKey::of(0,1));
        CHECK(optimizedOutcome.processed[0].edgeOutcome.status==
              bnd::TetEdgeRecoveryStatus::AlreadyPresent);
        CHECK(optimizedOutcome.processed[1].kind==
              bnd::TetFacetSetWorkKind::Edge);
        CHECK(optimizedOutcome.processed[1].edgeKey==
              bnd::VolumeArtifactEdgeKey::of(0,2));
        CHECK(optimizedOutcome.processed[1].edgeOutcome.status==
              bnd::TetEdgeRecoveryStatus::Recovered);
    }

    const auto edgePresent=static_cast<std::size_t>(
        bnd::TetEdgeRecoveryStatus::AlreadyPresent);
    const auto edgeRecovered=static_cast<std::size_t>(
        bnd::TetEdgeRecoveryStatus::Recovered);
    const auto facetPresent=static_cast<std::size_t>(
        bnd::TetFacetRecoveryStatus::AlreadyPresent);
    const auto facetMarked=static_cast<std::size_t>(
        bnd::TetFacetRecoveryStatus::Marked);
    CHECK(sink.metrics.edgeAssessmentEventStatusCounts[edgePresent]>0);
    CHECK(sink.metrics.edgeAssessmentEventStatusCounts[edgeRecovered]>0);
    CHECK(sink.metrics.facetAssessmentEventStatusCounts[facetPresent]>0);
    CHECK(sink.metrics.facetAssessmentEventStatusCounts[facetMarked]>0);
    CHECK(sink.metrics.presentBindingBatches>0);
    const auto presentReportEntities=static_cast<std::uint64_t>(std::count_if(
        optimizedOutcome.processed.begin(),optimizedOutcome.processed.end(),
        [](const auto& row) {
            return row.kind==bnd::TetFacetSetWorkKind::Edge
                ? row.edgeOutcome.status==
                    bnd::TetEdgeRecoveryStatus::AlreadyPresent
                : row.facetOutcome.status==
                    bnd::TetFacetRecoveryStatus::AlreadyPresent;
        }));
    CHECK(sink.metrics.presentBindingEntities==presentReportEntities);
    CHECK(sink.metrics.edgeAssessmentEventStatusCounts[edgePresent]+
          sink.metrics.facetAssessmentEventStatusCounts[facetPresent]>=
          sink.metrics.presentBindingEntities);
    CHECK(sink.metrics.exactPresentSelectedRebuilds==0);
    CHECK(sink.metrics.exactPresentDependencyTests==0);
    const auto& assessmentWork=sink.metrics.assessmentBindingWork;
    const auto& recoveryWork=sink.metrics.recoveryBindingWork;
    CHECK(assessmentWork.complete);
    CHECK(recoveryWork.complete);
    CHECK(assessmentWork.edgeBindingAuditCalls>0);
    CHECK(assessmentWork.facetBindingAuditCalls>0);
    CHECK(assessmentWork.edgeBindingRowsVisited>0);
    CHECK(assessmentWork.facetBindingRowsVisited>0);
    CHECK(assessmentWork.edgeDataCopyRows>=
          assessmentWork.edgeBindingRowsVisited);
    CHECK(assessmentWork.faceDataCopyRows>=
          assessmentWork.facetBindingRowsVisited);
    CHECK(assessmentWork.faceIncidenceRowsCopied>=
          assessmentWork.faceDataCopyRows);
    CHECK(assessmentWork.facetBindingReceiptBuilds>0);
    CHECK(assessmentWork.facetBindingReceiptReuses>0);
    CHECK(assessmentWork.facetBindingReceiptRowsCopied>0);
    CHECK(assessmentWork.facetBindingReceiptLookups>0);
    CHECK(assessmentWork.facetBindingRequestRowsVisited==0);
    CHECK(recoveryWork.stateCopies>0);
    CHECK(recoveryWork.stateCopyRows>=recoveryWork.stateCopies);

    std::uint64_t recoveredOperations=0;
    for (const auto& row : optimizedOutcome.processed) {
        if (row.kind==bnd::TetFacetSetWorkKind::Edge) {
            recoveredOperations+=row.edgeOutcome.status==
                bnd::TetEdgeRecoveryStatus::Recovered;
            if (row.edgeOutcome.status==
                bnd::TetEdgeRecoveryStatus::AlreadyPresent) {
                const auto binding=optimized.state.edgeBindings.find(
                    row.edgeKey);
                CHECK(binding!=optimized.state.edgeBindings.end());
                if (binding!=optimized.state.edgeBindings.end()) {
                    CHECK(optimized.workspace->edges().isCurrent(
                        binding->second.edge));
                    const auto data=optimized.workspace->edges().get(
                        binding->second.edge);
                    CHECK(data.requiredConstraint);
                    CHECK(data.cls.dim==binding->second.classification.dim);
                    CHECK(data.cls.tag==binding->second.classification.tag);
                }
            }
        } else {
            recoveredOperations+=row.facetOutcome.status==
                bnd::TetFacetRecoveryStatus::Recovered;
            if (row.facetOutcome.status==
                bnd::TetFacetRecoveryStatus::AlreadyPresent) {
                const auto binding=optimized.state.bindings.find(row.facetKey);
                CHECK(binding!=optimized.state.bindings.end());
                if (binding!=optimized.state.bindings.end()) {
                    CHECK(optimized.workspace->faces().isCurrent(
                        binding->second.recoveredFace));
                    const auto data=optimized.workspace->faces().get(
                        binding->second.recoveredFace);
                    CHECK(data.cls.dim==binding->second.classification.dim);
                    CHECK(data.cls.tag==binding->second.classification.tag);
                }
            }
        }
    }
    CHECK(recoveredOperations==1);
    CHECK(sink.metrics.nontrivialRecoveries==recoveredOperations);
    CHECK(recoveryWork.stateCopies>=recoveredOperations);
}

void indexedFacetBindingIdentityAudit() {
    TEST("t2b indexed facet bindings retain exact face identity");
    auto fixture=makeMixedPresentSchedulerFixture();
    bnd::TetFacetSetMixedPresentBatchSink sink;
    const auto built=bnd::recoverTetFacetSetIndexed(
        fixture.artifact,fixture.localByGlobal,{3,9},fixture.state,
        *fixture.workspace,*fixture.context,nullptr,nullptr,&sink);
    CHECK(built.status==bnd::TetFacetSetRecoveryStatus::Complete);
    CHECK(built.reportComplete);
    CHECK(fixture.state.bindings.size()>=2);
    if(fixture.state.bindings.size()<2)return;

    const auto& facet=fixture.artifact.facets.front();
    bnd::TetFacetIndexedRequest request;
    request.stamp=topo::WorkspaceStamp::of(*fixture.workspace);
    request.artifactStamp=fixture.state.artifactStamp;
    request.key=facet.key;
    request.classification=facet.classification;
    request.region={3,9};
    for(std::size_t i=0;i<3;++i)
        request.vertices[i]=fixture.localByGlobal.at(
            facet.directedGlobalNodes[i]);

    const auto workspaceBefore=
        fixture.workspace->semanticTopologyHashV2();
    bnd::TetRecoveryBindingWorkCounters exactWork;
    const auto exactFirst=bnd::assessTetFacetIndexed(
        request,fixture.state,*fixture.workspace,*fixture.context,&exactWork);
    const auto exactSecond=bnd::assessTetFacetIndexed(
        request,fixture.state,*fixture.workspace,*fixture.context,&exactWork);
    CHECK(exactFirst.status==bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(exactSecond.status==bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(exactFirst.face==exactSecond.face);
    CHECK(exactWork.complete);
    CHECK(exactWork.facetBindingAuditCalls==2);
    CHECK(exactWork.facetBindingReceiptBuilds==0);
    CHECK(exactWork.facetBindingReceiptReuses==0);
    CHECK(exactWork.facetBindingDuplicateRowsCompared>0);
    CHECK(fixture.workspace->semanticTopologyHashV2()==workspaceBefore);

    const auto exactStateBefore=fixture.state;
    const auto recoveredFirst=bnd::recoverTetFacetIndexed(
        request,fixture.state,*fixture.workspace,*fixture.context);
    const auto recoveredSecond=bnd::recoverTetFacetIndexed(
        request,fixture.state,*fixture.workspace,*fixture.context);
    CHECK(recoveredFirst.recovery.status==
          bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(recoveredSecond.recovery.status==
          bnd::TetFacetRecoveryStatus::AlreadyPresent);
    CHECK(recoveredFirst.recovery.face==recoveredSecond.recovery.face);
    requireMixedSchedulerStateEquivalent(
        exactStateBefore,fixture.state);
    CHECK(fixture.workspace->semanticTopologyHashV2()==workspaceBefore);

    const auto requireRejected=[&](
        const bnd::TetFacetSetState& state,
        const bnd::TetFacetIndexedRequest& candidate) {
        bnd::TetRecoveryBindingWorkCounters work;
        const auto outcome=bnd::assessTetFacetIndexed(
            candidate,state,*fixture.workspace,*fixture.context,&work);
        CHECK(outcome.status==bnd::TetFacetRecoveryStatus::Rejected);
        CHECK(outcome.reason==bnd::TetFacetRecoveryReason::StateKeyMismatch);
        CHECK(outcome.phase==bnd::TetFacetRecoveryPhase::Preflight);
        CHECK(fixture.workspace->semanticTopologyHashV2()==workspaceBefore);
        return work;
    };

    auto wrongTriple=fixture.state;
    auto wrong=wrongTriple.bindings.find(request.key);
    CHECK(wrong!=wrongTriple.bindings.end());
    if(wrong==wrongTriple.bindings.end())return;
    wrong->second.canonicalVertices={
        fixture.localByGlobal.at(0),fixture.localByGlobal.at(1),
        fixture.localByGlobal.at(3)};
    (void)requireRejected(wrongTriple,request);

    auto staleVertex=fixture.state;
    auto stale=staleVertex.bindings.find(request.key);
    CHECK(stale!=staleVertex.bindings.end());
    if(stale==staleVertex.bindings.end())return;
    ++stale->second.canonicalVertices[0].gen;
    (void)requireRejected(staleVertex,request);

    auto sameKeyOtherFace=fixture.state;
    auto source=sameKeyOtherFace.bindings.begin();
    if(source->first==request.key)++source;
    CHECK(source!=sameKeyOtherFace.bindings.end());
    if(source==sameKeyOtherFace.bindings.end())return;
    const auto sourceKey=source->first;
    auto replacement=source->second;
    replacement.key=request.key;
    sameKeyOtherFace.bindings.erase(sourceKey);
    sameKeyOtherFace.bindings[request.key]=replacement;
    auto sameClassRequest=request;
    sameClassRequest.classification=replacement.classification;
    (void)requireRejected(sameKeyOtherFace,sameClassRequest);

    auto duplicateTriple=fixture.state;
    auto first=duplicateTriple.bindings.begin();
    auto second=std::next(first);
    second->second.canonicalVertices=first->second.canonicalVertices;
    const auto duplicateTripleWork=requireRejected(duplicateTriple,request);
    CHECK(duplicateTripleWork.facetBindingDuplicateRowsCompared==0);
    CHECK(duplicateTripleWork.facetBindingRowsVisited==2);

    auto duplicateHandle=fixture.state;
    first=duplicateHandle.bindings.begin();
    second=std::next(first);
    second->second.recoveredFace=first->second.recoveredFace;
    const auto duplicateHandleWork=requireRejected(duplicateHandle,request);
    CHECK(duplicateHandleWork.facetBindingDuplicateRowsCompared==0);
    CHECK(duplicateHandleWork.facetBindingRowsVisited==2);
}

void mixedPresentEntryBindingValidationPrecedesSchedule() {
    TEST("t2b mixed entry binding validation precedes schedule");
    auto fixture=makeMixedPresentSchedulerFixture();
    CHECK(fixture.artifact.facets.size()>=2);
    if(fixture.artifact.facets.size()<2)return;

    // Facet 810 is deliberately not an exact-present candidate in this
    // workspace.  A corrupt retained row for it therefore used to survive
    // the fast preclassification pass until after canonical schedules had
    // already been installed.
    const auto& facet=fixture.artifact.facets[1];
    std::array<topo::NodeHandle,3> target{
        fixture.localByGlobal.at(facet.directedGlobalNodes[0]),
        fixture.localByGlobal.at(facet.directedGlobalNodes[1]),
        fixture.localByGlobal.at(facet.directedGlobalNodes[2])};
    CHECK(!fixture.workspace->faces().find(
        fixture.workspace->nodes(),target).has_value());
    std::sort(target.begin(),target.end(),[&](const auto a,const auto b) {
        return fixture.workspace->nodes().seqOf(a)<
               fixture.workspace->nodes().seqOf(b);
    });
    const std::array<topo::NodeHandle,3> other{
        fixture.localByGlobal.at(6),fixture.localByGlobal.at(7),
        fixture.localByGlobal.at(8)};
    const auto wrongFace=fixture.workspace->faces().find(
        fixture.workspace->nodes(),other);
    CHECK(wrongFace.has_value());
    if(!wrongFace)return;
    fixture.state.bindings.emplace(facet.key,
        bnd::PersistentTetFacetBinding{
            facet.key,target,facet.classification,*wrongFace});

    const bnd::VolumeArtifactEdgeKey sentinelEdge=
        bnd::VolumeArtifactEdgeKey::of(90,91);
    const bnd::VolumeFacetKey sentinelFacet{9999,7};
    fixture.state.canonicalEdgeSchedule={sentinelEdge};
    fixture.state.canonicalSchedule={sentinelFacet};
    const auto stateBefore=fixture.state;
    const auto digestBefore=
        paacs_tests::observableDigest(*fixture.workspace);
    const auto revisionBefore=fixture.workspace->topologyRevision();

    std::array<bnd::TetFacetSetFrontierRow,1> frontierRows{};
    std::array<topo::CellHandle,1> frontierCells{};
    std::array<std::size_t,1> frontierDependencies{};
    bnd::TetFacetSetFrontierSink frontier{
        frontierRows,frontierCells,frontierDependencies};
    std::array<bnd::VolumeFacetKey,1> cycleMembers{};
    std::array<std::array<bnd::VolumeFacetKey,2>,1> cycleDependencies{};
    bnd::TetFacetConflictCycleSink cycle{
        cycleMembers,cycleDependencies};
    bnd::TetFacetSetMixedPresentBatchSink sink;
    const auto outcome=bnd::recoverTetFacetSetIndexed(
        fixture.artifact,fixture.localByGlobal,{3,9},fixture.state,
        *fixture.workspace,*fixture.context,&frontier,&cycle,&sink);

    CHECK(outcome.status==bnd::TetFacetSetRecoveryStatus::Rejected);
    CHECK(outcome.reason==bnd::TetFacetSetRecoveryReason::None);
    CHECK(!outcome.reportComplete);
    CHECK(outcome.processed.empty());
    CHECK(outcome.unattemptedEdges.empty());
    CHECK(outcome.unattemptedFacets.empty());
    CHECK(outcome.schedulerMetrics.routeAssessments==0);
    CHECK(outcome.schedulerMetrics.selectedRebuilds==0);
    CHECK(sink.metrics.complete);
    CHECK(sink.metrics.presentBindingBatches==0);
    CHECK(sink.metrics.presentBindingEntities==0);
    CHECK(sink.metrics.stateCopies==0);
    CHECK(sink.metrics.nontrivialRecoveries==0);
    CHECK(sink.metrics.assessmentCounters==bnd::TetFacetGeometryCounters{});
    CHECK(sink.metrics.recoveryCounters==bnd::TetFacetGeometryCounters{});
    CHECK(sink.metrics.assessmentBindingWork.facetBindingAuditCalls==1);
    CHECK(sink.metrics.assessmentBindingWork.facetBindingReceiptBuilds==1);
    CHECK(sink.metrics.assessmentBindingWork.facetBindingRowsVisited==1);
    CHECK(sink.metrics.recoveryBindingWork==
          bnd::TetRecoveryBindingWorkCounters{});
    CHECK(frontier.rowsUsed==0);
    CHECK(frontier.supportCellsUsed==0);
    CHECK(frontier.dependenciesUsed==0);
    CHECK(!cycle.captured);
    CHECK(cycle.membersUsed==0);
    CHECK(cycle.dependenciesUsed==0);
    CHECK(fixture.state.workspaceStamp==stateBefore.workspaceStamp);
    CHECK(fixture.state.artifactStamp==stateBefore.artifactStamp);
    requireMixedSchedulerStateEquivalent(stateBefore,fixture.state);
    CHECK(fixture.state.canonicalEdgeSchedule==
          std::vector<bnd::VolumeArtifactEdgeKey>{sentinelEdge});
    CHECK(fixture.state.canonicalSchedule==
          std::vector<bnd::VolumeFacetKey>{sentinelFacet});
    CHECK(fixture.workspace->topologyRevision()==revisionBefore);
    CHECK(fixture.context->currentFor(*fixture.workspace));
    CHECK(paacs_tests::observableDigest(*fixture.workspace)==digestBefore);
}

void mixedPresentSchedulerPermanentPrefix() {
    TEST("t2b mixed scheduler permanent stop precedes exact-present batch");
    topo::NodePool nodes;
    std::array<topo::NodeHandle,4> n{
        nodes.add({0,0,0},{3,9}),nodes.add({1,0,0},{3,9}),
        nodes.add({0,1,0},{3,9}),nodes.add({0,0,1},{3,9})};
    topo::CellPool cells;
    (void)cells.add(topo::ElementType::Tet4,positiveHandles(nodes,n),
                    {3,9},nodes);
    topo::EdgePool edges;
    for (const auto cell : cells.activeBySequence())
        (void)edges.insertCellEdges(nodes,cells,cell);
    const auto exact=edges.find(nodes,n[0],n[1]);
    const auto conflict=edges.find(nodes,n[1],n[2]);
    CHECK(exact.has_value());
    CHECK(conflict.has_value());
    if (!exact || !conflict) return;
    edges.setClassification(*exact,{1,41});
    edges.setRequiredConstraint(*exact,true);
    edges.setClassification(*conflict,{1,42});
    edges.setRequiredConstraint(*conflict,true);
    topo::FacePool faces;
    for (const auto cell : cells.activeBySequence())
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    topo::TopologyWorkspace workspace(
        std::move(nodes),std::move(cells),std::move(edges),std::move(faces));
    auto indexed=bnd::buildTetRecoveryContext(workspace,{3,9});
    CHECK(indexed.built());
    if (!indexed.context) return;
    bnd::TetFacetSetState state;
    state.workspaceStamp=topo::WorkspaceStamp::of(workspace);
    state.artifactStamp={705,1};
    bnd::ClosedTri3Artifact artifact;
    artifact.surfaceSetStamp=state.artifactStamp;
    artifact.regionTags={9};
    const auto exactKey=bnd::VolumeArtifactEdgeKey::of(0,1);
    const auto conflictKey=bnd::VolumeArtifactEdgeKey::of(1,2);
    artifact.edges={
        {exactKey,{1,41},{}},{conflictKey,{1,99},{}}};
    const std::map<Index,topo::NodeHandle> localByGlobal{
        {0,n[0]},{1,n[1]},{2,n[2]},{3,n[3]}};
    const auto digest=workspace.semanticTopologyHashV2();
    bnd::TetFacetSetMixedPresentBatchSink sink;
    const auto outcome=bnd::recoverTetFacetSetIndexed(
        artifact,localByGlobal,{3,9},state,workspace,*indexed.context,
        nullptr,nullptr,&sink);
    CHECK(outcome.status==bnd::TetFacetSetRecoveryStatus::Refused);
    CHECK(outcome.reportComplete);
    CHECK(outcome.processed.size()==1);
    if (outcome.processed.size()==1) {
        CHECK(outcome.processed[0].kind==bnd::TetFacetSetWorkKind::Edge);
        CHECK(outcome.processed[0].edgeKey==conflictKey);
        CHECK(outcome.processed[0].edgeOutcome.reason==
              bnd::TetEdgeRecoveryReason::ProtectedConstraintConflict);
    }
    CHECK(outcome.unattemptedEdges==
          std::vector<bnd::VolumeArtifactEdgeKey>{exactKey});
    CHECK(state.edgeBindings.empty());
    CHECK(state.bindings.empty());
    CHECK(sink.metrics.presentBindingBatches==0);
    CHECK(sink.metrics.presentBindingEntities==0);
    CHECK(sink.metrics.stateCopies==0);
    CHECK(sink.metrics.nontrivialRecoveries==0);
    CHECK(sink.metrics.exactPresentSelectedRebuilds==0);
    CHECK(sink.metrics.exactPresentDependencyTests==0);
    CHECK(workspace.semanticTopologyHashV2()==digest);
}

void schedulerPermanentStopPrecedesReadyCommit() {
    TEST("t2b scheduler stops on later permanent facet before ready commit");
    topo::NodePool nodes;
    std::vector<topo::NodeHandle> n;
    for (const Vec3& p : std::vector<Vec3>{{0,0,0},{4,0,0},{0,4,0},
                                           {1,1,3},{1,1,-3}})
        n.push_back(nodes.add(p,{3,9}));
    topo::CellPool cells;
    std::vector<topo::CellHandle> current;
    const std::array<std::array<std::size_t,4>,3> initial{{
        {{3,4,0,1}},{{3,4,1,2}},{{3,4,2,0}}}};
    for (const auto ids : initial)
        current.push_back(cells.add(
            topo::ElementType::Tet4,
            positiveHandles(nodes,{n[ids[0]],n[ids[1]],n[ids[2]],n[ids[3]]}),
            {3,9},nodes));
    topo::EdgePool edges;
    for (const auto cell : current)
        (void)edges.insertCellEdges(nodes,cells,cell);
    topo::FacePool faces;
    for (const auto cell : current)
        (void)faces.insertCellFaces(nodes,cells,edges,cell);
    topo::TopologyWorkspace workspace(std::move(nodes),std::move(cells),
                                      std::move(edges),std::move(faces));
    auto indexed = bnd::buildTetRecoveryContext(workspace,{3,9});
    CHECK(indexed.built());
    if (!indexed.context) return;
    bnd::TetFacetSetState state;
    state.workspaceStamp = topo::WorkspaceStamp::of(workspace);
    state.artifactStamp = {702,1};
    bnd::ClosedTri3Artifact artifact;
    artifact.surfaceSetStamp = state.artifactStamp;
    artifact.regionTags = {9};
    const bnd::VolumeFacetKey readyKey{800,0};
    const bnd::VolumeFacetKey permanentKey{801,0};
    artifact.facets.push_back(
        {readyKey,{0,1,2},{2,77},{{9,1}}});
    // This exact face is one-sided in the current complex and therefore a
    // permanent T2-A refusal. Its key sorts after the recoverable facet.
    artifact.facets.push_back(
        {permanentKey,{0,1,3},{2,77},{{9,1}}});
    std::map<bnd::VolumeArtifactEdgeKey,bnd::VolumeArtifactEdge> edgeRows;
    for (const auto& facet : artifact.facets)
        for (std::size_t i = 0; i < 3; ++i) {
            const auto key = bnd::VolumeArtifactEdgeKey::of(
                facet.directedGlobalNodes[i],
                facet.directedGlobalNodes[(i+1)%3]);
            auto [it,inserted] = edgeRows.emplace(
                key,bnd::VolumeArtifactEdge{key,{2,77},{}});
            (void)inserted;
            it->second.incidentFacets.push_back(facet.key);
        }
    for (auto& [key,row] : edgeRows) {
        (void)key;
        artifact.edges.push_back(std::move(row));
    }
    const std::map<Index,topo::NodeHandle> localByGlobal{
        {0,n[0]},{1,n[1]},{2,n[2]},{3,n[3]},{4,n[4]}};
    const auto cellsBefore = workspace.cells().activeBySequence();
    CHECK(!workspace.faces().find(
        workspace.nodes(),std::array<topo::NodeHandle,3>{n[0],n[1],n[2]}));
    const auto scheduled = bnd::recoverTetFacetSetIndexed(
        artifact,localByGlobal,{3,9},state,workspace,*indexed.context);
    CHECK(scheduled.status == bnd::TetFacetSetRecoveryStatus::Refused);
    CHECK(scheduled.reportComplete);
    CHECK(!scheduled.processed.empty());
    if (!scheduled.processed.empty()) {
        const auto& stop = scheduled.processed.back();
        CHECK(stop.kind == bnd::TetFacetSetWorkKind::Facet);
        CHECK(stop.facetKey == permanentKey);
        CHECK(stop.facetOutcome.reason ==
              bnd::TetFacetRecoveryReason::OneSidedTargetUnsupported);
    }
    CHECK(scheduled.unattemptedFacets.size() == 1);
    CHECK(!scheduled.unattemptedFacets.empty() &&
          scheduled.unattemptedFacets[0] == readyKey);
    CHECK(workspace.cells().activeBySequence() == cellsBefore);
    CHECK(!workspace.faces().find(
        workspace.nodes(),std::array<topo::NodeHandle,3>{n[0],n[1],n[2]}));
    CHECK(scheduled.schedulerMetrics.fullStageRescans == 0);
}

void canonicalConflictScc() {
    TEST("t2b scheduler canonical SCC and acyclic discriminator");
    std::vector<std::vector<std::size_t>> graph{{1},{0},{2},{2}};
    std::vector<bool> pending(graph.size(),true);
    std::uint64_t visits = 0;
    const auto first = bnd::canonicalTetRecoveryConflictCycleForTesting(
        graph,pending,visits);
    CHECK(first == std::vector<std::size_t>({0,1}));
    CHECK(visits == 4);
    std::reverse(graph[0].begin(),graph[0].end());
    std::uint64_t repeatVisits = 0;
    const auto repeat = bnd::canonicalTetRecoveryConflictCycleForTesting(
        graph,pending,repeatVisits);
    CHECK(repeat == first);
    CHECK(repeatVisits == visits);
    pending[0] = false;
    pending[1] = false;
    std::uint64_t reducedVisits = 0;
    const auto reduced = bnd::canonicalTetRecoveryConflictCycleForTesting(
        graph,pending,reducedVisits);
    CHECK(reduced == std::vector<std::size_t>({2}));
    CHECK(reducedVisits == 2);
    const std::vector<std::vector<std::size_t>> acyclic{{1},{2},{}};
    const std::vector<bool> allPending(3,true);
    std::uint64_t acyclicVisits = 0;
    CHECK(bnd::canonicalTetRecoveryConflictCycleForTesting(
              acyclic,allPending,acyclicVisits).empty());
    CHECK(acyclicVisits == 3);
    CHECK(bnd::tetRecoveryBlockerGraphValidForTesting(acyclic,allPending));
    auto invalid = acyclic;
    invalid[0] = {3};
    CHECK(!bnd::tetRecoveryBlockerGraphValidForTesting(invalid,allPending));
    auto completedBlocker = allPending;
    completedBlocker[1] = false;
    invalid = {{1},{},{}};
    CHECK(!bnd::tetRecoveryBlockerGraphValidForTesting(
        invalid,completedBlocker));

    const topo::WorkspaceStamp workspaceStamp{1,2,3,4,5};
    const algo::ArtifactStamp artifactStamp{6,7};
    const std::array<bnd::VolumeFacetKey,2> sourceMembers{{
        {41,4},{40,4}}};
    const std::array<std::array<bnd::VolumeFacetKey,2>,2>
        sourceDependencies{{
            {{{41,4},{40,4}}},
            {{{40,4},{41,4}}}}};
    std::array<bnd::VolumeFacetKey,2> memberStorage{};
    std::array<std::array<bnd::VolumeFacetKey,2>,2> dependencyStorage{};
    bnd::TetFacetConflictCycleSink sink{
        memberStorage,dependencyStorage};
    bnd::fillTetFacetConflictCycleSinkForTesting(
        sourceMembers,sourceDependencies,workspaceStamp,artifactStamp,sink);
    CHECK(sink.captured && sink.complete);
    CHECK(sink.workspaceStamp == workspaceStamp);
    CHECK(sink.artifactStamp == artifactStamp);
    CHECK(sink.membersUsed == 2 && sink.membersRequired == 2);
    CHECK(sink.dependenciesUsed == 2 && sink.dependenciesRequired == 2);
    CHECK((memberStorage[0] == bnd::VolumeFacetKey{40,4}));
    CHECK((memberStorage[1] == bnd::VolumeFacetKey{41,4}));
    CHECK((dependencyStorage[0] ==
           std::array<bnd::VolumeFacetKey,2>{{{40,4},{41,4}}}));
    CHECK((dependencyStorage[1] ==
           std::array<bnd::VolumeFacetKey,2>{{{41,4},{40,4}}}));

    std::array<bnd::VolumeFacetKey,1> shortMemberStorage{};
    std::array<std::array<bnd::VolumeFacetKey,2>,1>
        shortDependencyStorage{};
    bnd::TetFacetConflictCycleSink shortSink{
        shortMemberStorage,shortDependencyStorage};
    bnd::fillTetFacetConflictCycleSinkForTesting(
        sourceMembers,sourceDependencies,workspaceStamp,artifactStamp,
        shortSink);
    CHECK(shortSink.captured && !shortSink.complete);
    CHECK(shortSink.workspaceStamp == workspaceStamp);
    CHECK(shortSink.artifactStamp == artifactStamp);
    CHECK(shortSink.membersUsed == 1 && shortSink.membersRequired == 2);
    CHECK(shortSink.dependenciesUsed == 1 &&
          shortSink.dependenciesRequired == 2);

    bnd::TetFacetConflictCycleSink zeroSink{};
    bnd::fillTetFacetConflictCycleSinkForTesting(
        sourceMembers,sourceDependencies,workspaceStamp,artifactStamp,
        zeroSink);
    CHECK(zeroSink.captured && !zeroSink.complete);
    CHECK(zeroSink.membersUsed == 0 && zeroSink.membersRequired == 2);
    CHECK(zeroSink.dependenciesUsed == 0 &&
          zeroSink.dependenciesRequired == 2);
}

void boundaryConformingV3Box() {
    TEST("t2b boundary-conforming-v3 box publishes zero-mismatch Tet4");
    auto kernel = geo::makeBuiltinKernel();
    auto model = kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    configureOneRegionProvider(*model,provider);
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(0.55);
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp = {901,1};
    request.expectedOwnershipProviderStamp = provider.current;
    request.envelopeMode =
        algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    mesh::Mesh destination;
    const auto outcome = algo::meshModelBoundaryConforming(
        *model,spec,request,destination,provider);
    if (!outcome.published()) {
        std::printf("  v3 box diagnostic status=%d reason=%d phase=%d "
                    "witness=%llu rows=%zu\n",
                    static_cast<int>(outcome.status),
                    static_cast<int>(outcome.failure.reason),
                    static_cast<int>(outcome.failure.phase),
                    static_cast<unsigned long long>(
                        outcome.failure.exactWitnessDigest),
                    outcome.report ? outcome.report->processed.size() : 0U);
        if (outcome.report)
            for (const auto& row : outcome.report->processed)
                if (row.status == algo::BoundaryConformingV3WorkStatus::Rejected ||
                    row.status == algo::BoundaryConformingV3WorkStatus::Unsupported ||
                    row.status == algo::BoundaryConformingV3WorkStatus::Refused ||
                    row.status == algo::BoundaryConformingV3WorkStatus::Failed)
                std::printf("    row=%zu kind=%d status=%d edge=(%llu,%llu) "
                            "facet=%llu edgeReason=%d facetReason=%d "
                            "facetGeom=%d\n",
                    row.scheduleOrdinal,static_cast<int>(row.key.kind),
                    static_cast<int>(row.status),
                    static_cast<unsigned long long>(
                        row.key.edgeKey.canonicalGlobalNodes[0]),
                    static_cast<unsigned long long>(
                        row.key.edgeKey.canonicalGlobalNodes[1]),
                    static_cast<unsigned long long>(
                        row.key.facetKey.facePersistentId),
                    static_cast<int>(row.failure.edgeOutcome.reason),
                    static_cast<int>(row.failure.facetOutcome.reason),
                    static_cast<int>(row.failure.facetOutcome.geometryReason));
    }
    CHECK(outcome.published());
    CHECK(outcome.retainedResult != nullptr);
    CHECK(outcome.surfaceOutcome != nullptr);
    CHECK(outcome.surfaceOutcome &&
          outcome.surfaceOutcome->setStatus == algo::SurfaceSetStatus::Published);
    if (outcome.retainedResult) {
        CHECK(outcome.retainedResult->constraintState.artifactStamp ==
              request.expectedSurfaceSetStamp);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->retainedArtifact &&
              outcome.surfaceOutcome->retainedArtifact->setStamp ==
                  request.expectedSurfaceSetStamp);
        CHECK(outcome.retainedResult->boundaryAudit.passed);
        CHECK(outcome.retainedResult->boundaryAudit.orphanExterior == 0);
        CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior == 0);
        CHECK(outcome.retainedResult->recoveryMetrics.
              productionGlobalScanFallbacks == 0);
        CHECK(outcome.retainedResult->recoveryMetrics.fullStageRescans == 0);
        CHECK(outcome.retainedResult->recoveryMetrics.runComplete);
        CHECK(outcome.retainedResult->extractionMetrics.runComplete);
        CHECK(outcome.retainedResult->extractionMetrics.
              centroidClassifierCalls == 0);
        CHECK(outcome.retainedResult->extractionMetrics.
              boundaryPeelingCalls == 0);
        CHECK(outcome.retainedResult->extractionMetrics.
              hullCapPeelingCalls == 0);
        const auto workCount =
            outcome.retainedResult->constraintState.edgeBindings.size()+
            outcome.retainedResult->constraintState.bindings.size();
        CHECK(outcome.retainedResult->recoveryMetrics.
              batchPreflightAssessments == workCount);
        CHECK(outcome.retainedResult->recoveryMetrics.
              alreadyPresentBatchEntities == workCount);
        CHECK(outcome.retainedResult->recoveryMetrics.selectedRebuilds == 0);
        CHECK(outcome.retainedResult->recoveryMetrics.
              protectedImportFallbacks == 0);
        CHECK(outcome.retainedResult->recoveryMetrics.routeAssessments+
              outcome.retainedResult->recoveryMetrics.selectedRebuilds <=
              4*workCount);
    }
    CHECK(destination.countCells(mesh::CellType::Tri3) > 0);
    CHECK(destination.countCells(mesh::CellType::Tet4) > 0);
    const auto audit = mesh::auditMesh(destination);
    CHECK(audit.boundaryOrphans == 0);
    CHECK(audit.boundaryUnexpected == 0);
    CHECK(audit.nInverted == 0);
}

void boundaryConformingV3SurfaceStops() {
    TEST("t2b v3 preserves complete nested surface stop outcomes");
    auto kernel = geo::makeBuiltinKernel();
    auto box = kernel->makeBox({0,0,0},{1,1,1});
    Provider provider;
    sizing::MeshSpec validSpec;
    validSpec.field = std::make_shared<sizing::ConstantField>(1.0);
    auto run = [&](const geo::Model& model,const sizing::MeshSpec& spec,
                   std::uint64_t stamp) {
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {stamp,1};
        request.expectedOwnershipProviderStamp = provider.current;
        mesh::Mesh destination;
        destination.addNode({-3,-4,-5},{0,988});
        const auto before = mesh::meshHash(destination);
        auto outcome = algo::meshModelBoundaryConforming(
            model,spec,request,destination,provider);
        CHECK(!outcome.published());
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination) == before);
        CHECK(outcome.reportComplete);
        CHECK(outcome.report && outcome.report->reportComplete);
        CHECK(outcome.report && outcome.report->processed.empty());
        CHECK(outcome.report && outcome.report->unattempted.empty());
        CHECK(outcome.surfaceOutcome != nullptr);
        CHECK(provider.ownershipCalls == 0);
        CHECK(provider.envelopeCalls == 0);
        return outcome;
    };
    {
        EmptyModel empty;
        auto outcome = run(empty,validSpec,970);
        CHECK(outcome.status == algo::BoundaryConformingV3Status::Refused);
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::SurfaceConstructionRefused);
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::SurfaceConstruction);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->setStatus ==
                  algo::SurfaceSetStatus::Refused);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->reportComplete);
    }
    {
        DuplicateFaceModel duplicate(*box);
        auto outcome = run(duplicate,validSpec,971);
        CHECK(outcome.status == algo::BoundaryConformingV3Status::Rejected);
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::SurfaceConstructionRejected);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->setStatus ==
                  algo::SurfaceSetStatus::Rejected);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->reportComplete);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->processed.size() ==
                  duplicate.faces().size());
    }
    {
        sizing::MeshSpec invalid = validSpec;
        invalid.field = std::make_shared<InvalidEvaluatedSize>();
        auto outcome = run(*box,invalid,972);
        CHECK(outcome.status == algo::BoundaryConformingV3Status::Failed);
        CHECK(outcome.failure.reason ==
              algo::BoundaryConformingV3Reason::SizingEvaluationFailed);
        CHECK(outcome.failure.phase ==
              algo::BoundaryConformingV3Phase::SurfaceConstruction);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->setStatus ==
                  algo::SurfaceSetStatus::Failed);
        CHECK(outcome.surfaceOutcome &&
              outcome.surfaceOutcome->reportComplete);
    }
}

void boundaryConformingV3CylinderAndThinBox() {
    TEST("t2b boundary-conforming-v3 cylinder and thin box finish 0/0");
    auto kernel = geo::makeBuiltinKernel();
    struct Fixture {
        std::unique_ptr<geo::Model> model;
        double h = 0.0;
        std::uint64_t stamp = 0;
    };
    std::array<Fixture,2> fixtures{
        Fixture{kernel->makeCylinder({0,0,0},0.5,1.0),0.15,902},
        Fixture{kernel->makeBox({0,0,0},{1,1,1e-3}),0.2,903}};
    std::uint64_t totalEdgeCommits = 0,totalFacetCommits = 0;
    for (auto& fixture : fixtures) {
        Provider provider;
        configureOneRegionProvider(*fixture.model,provider);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(fixture.h);
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {fixture.stamp,1};
        request.expectedOwnershipProviderStamp = provider.current;
        request.envelopeMode =
            algo::RecoveryEnvelopeMode::DerivedFiniteBox;
        mesh::Mesh destination;
        const auto outcome = algo::meshModelBoundaryConforming(
            *fixture.model,spec,request,destination,provider);
        if (!outcome.published())
            std::printf("  v3 fixture stamp=%llu status=%d reason=%d "
                        "phase=%d rows=%zu\n",
                static_cast<unsigned long long>(fixture.stamp),
                static_cast<int>(outcome.status),
                static_cast<int>(outcome.failure.reason),
                static_cast<int>(outcome.failure.phase),
                outcome.report ? outcome.report->processed.size() : 0U);
        CHECK(outcome.published());
        CHECK(outcome.retainedResult != nullptr);
        if (outcome.retainedResult) {
            CHECK(outcome.retainedResult->boundaryAudit.passed);
            CHECK(outcome.retainedResult->boundaryAudit.orphanExterior == 0);
            CHECK(outcome.retainedResult->boundaryAudit.unexpectedExterior == 0);
            CHECK(outcome.retainedResult->recoveryMetrics.
                  productionGlobalScanFallbacks == 0);
            CHECK(outcome.retainedResult->recoveryMetrics.fullStageRescans == 0);
            CHECK(outcome.retainedResult->extractionMetrics.
                  centroidClassifierCalls == 0);
            CHECK(outcome.retainedResult->extractionMetrics.
                  boundaryPeelingCalls == 0);
            CHECK(outcome.retainedResult->extractionMetrics.
                  hullCapPeelingCalls == 0);
            totalEdgeCommits += outcome.retainedResult->recoveryMetrics.
                                    edgeCavityCommits;
            totalFacetCommits += outcome.retainedResult->recoveryMetrics.
                                     facetCommits;
        }
        const auto audit = mesh::auditMesh(destination);
        CHECK(audit.boundaryOrphans == 0);
        CHECK(audit.boundaryUnexpected == 0);
        CHECK(audit.nInverted == 0);
    }
    std::printf("  v3 cylinder/thin recovered edges=%llu facets=%llu\n",
        static_cast<unsigned long long>(totalEdgeCommits),
        static_cast<unsigned long long>(totalFacetCommits));
    // These two convex artifacts match the frozen hull triangulation exactly;
    // the separate interacting scheduler fixture pins genuine mutation.
}

void boundaryConformingV3Determinism() {
    TEST("t2b v3 face order allocator history and 1/4/auto determinism");
    auto kernel = geo::makeBuiltinKernel();
    auto source = kernel->makeCylinder({0,0,0},0.5,1.0);
    ReversedFaceModel reversed(*source);
    struct Evidence {
        std::uint64_t strict = 0;
        std::uint64_t semantic = 0;
        Index nodes = 0;
        Index cells = 0;
        std::vector<bnd::VolumeArtifactEdgeKey> edgeSchedule;
        std::vector<bnd::VolumeFacetKey> facetSchedule;
        std::uint64_t routeAssessments = 0;
        std::uint64_t dependencyEdges = 0;
        std::uint64_t selectedRebuilds = 0;
    };
    std::array<Evidence,3> evidence{};
    const std::array<int,3> threads{{1,4,0}};
    for (std::size_t run = 0; run < threads.size(); ++run) {
        topo::NodePool churnNodes;
        topo::CellPool churnCells;
        for (std::size_t i = 0; i < run*17; ++i)
            (void)churnNodes.add(
                {static_cast<double>(i),-1,2},{0,static_cast<int>(i)});
        CHECK(churnCells.stats().live == 0);
        const geo::Model& presented = run == 1
            ? static_cast<const geo::Model&>(reversed)
            : static_cast<const geo::Model&>(*source);
        Provider provider;
        configureOneRegionProvider(presented,provider);
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(0.15);
        spec.threads = threads[run];
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {920,1};
        request.expectedOwnershipProviderStamp = provider.current;
        request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
        mesh::Mesh destination;
        const auto outcome = algo::meshModelBoundaryConforming(
            presented,spec,request,destination,provider);
        CHECK(outcome.published());
        CHECK(outcome.retainedResult != nullptr);
        if (!outcome.retainedResult) continue;
        const auto semantic = mesh::semanticHashEx(destination);
        CHECK(semantic.canonical);
        evidence[run].strict = mesh::meshHash(destination);
        evidence[run].semantic = semantic.value;
        evidence[run].nodes = destination.nodeCount();
        evidence[run].cells = destination.cellCount();
        evidence[run].edgeSchedule =
            outcome.retainedResult->constraintState.canonicalEdgeSchedule;
        evidence[run].facetSchedule =
            outcome.retainedResult->constraintState.canonicalSchedule;
        evidence[run].routeAssessments =
            outcome.retainedResult->recoveryMetrics.routeAssessments;
        evidence[run].dependencyEdges =
            outcome.retainedResult->recoveryMetrics.dependencyEdges;
        evidence[run].selectedRebuilds =
            outcome.retainedResult->recoveryMetrics.selectedRebuilds;
    }
    for (std::size_t run = 1; run < evidence.size(); ++run) {
        CHECK(evidence[run].strict == evidence[0].strict);
        CHECK(evidence[run].semantic == evidence[0].semantic);
        CHECK(evidence[run].nodes == evidence[0].nodes);
        CHECK(evidence[run].cells == evidence[0].cells);
        CHECK(evidence[run].edgeSchedule == evidence[0].edgeSchedule);
        CHECK(evidence[run].facetSchedule == evidence[0].facetSchedule);
        CHECK(evidence[run].routeAssessments ==
              evidence[0].routeAssessments);
        CHECK(evidence[run].dependencyEdges == evidence[0].dependencyEdges);
        CHECK(evidence[run].selectedRebuilds ==
              evidence[0].selectedRebuilds);
    }
}

void boundaryConformingV3SealedCavity() {
    TEST("t2b v3 outer plus sealed-cavity components extract topologically");
    SealedCavityModel model;
    Provider provider;
    provider.supplied = boxEnvelope(-4,4);
    for (const geo::GFace* face : model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1] =
            {{0,1}};
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(10.0);
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp = {940,1};
    request.expectedOwnershipProviderStamp = provider.current;
    request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    mesh::Mesh destination;
    const auto outcome = algo::meshModelBoundaryConforming(
        model,spec,request,destination,provider);
    if (!outcome.published())
        std::printf("  sealed cavity status=%d reason=%d phase=%d rows=%zu\n",
            static_cast<int>(outcome.status),
            static_cast<int>(outcome.failure.reason),
            static_cast<int>(outcome.failure.phase),
            outcome.report ? outcome.report->processed.size() : 0U);
    CHECK(outcome.published());
    CHECK(outcome.retainedResult != nullptr);
    if (!outcome.retainedResult) return;
    const auto& result = *outcome.retainedResult;
    CHECK(result.boundaryAudit.passed);
    CHECK(result.boundaryAudit.orphanExterior == 0);
    CHECK(result.boundaryAudit.unexpectedExterior == 0);
    CHECK(result.boundaryAudit.prescribedExterior ==
          result.boundaryAudit.matchedExterior);
    CHECK(result.extractionMetrics.exteriorComponents >= 2);
    CHECK(result.extractionMetrics.retainedRegionComponents == 1);
    CHECK(result.extractionMetrics.centroidClassifierCalls == 0);
    CHECK(result.extractionMetrics.boundaryPeelingCalls == 0);
    CHECK(result.extractionMetrics.hullCapPeelingCalls == 0);
    mesh::AuditOptions options;
    options.boundaryMode = mesh::BoundaryComparisonMode::Required;
    options.boundaryOrphanBudget = 0;
    options.boundaryUnexpectedBudget = 0;
    const auto audit = mesh::auditMesh(destination,options);
    CHECK(audit.boundaryOrphans == 0);
    CHECK(audit.boundaryUnexpected == 0);
    CHECK(audit.nInverted == 0);
}

void boundaryConformingV3DisconnectedRegionRefusal() {
    TEST("t2b v3 disconnected solids sharing one region are unsupported");
    SealedCavityModel model(true);
    Provider provider;
    provider.supplied = boxEnvelope(-5,5);
    for (const geo::GFace* face : model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1] =
            {{0,1}};
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(10.0);
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp = {941,1};
    request.expectedOwnershipProviderStamp = provider.current;
    request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    mesh::Mesh destination;
    destination.addNode({-9,-9,-9},{0,77});
    const auto before = mesh::meshHash(destination);
    const auto outcome = algo::meshModelBoundaryConforming(
        model,spec,request,destination,provider);
    CHECK(outcome.status == algo::BoundaryConformingV3Status::Unsupported);
    CHECK(outcome.failure.reason ==
          algo::BoundaryConformingV3Reason::DisconnectedRegionTagReuse);
    CHECK(outcome.failure.phase ==
          algo::BoundaryConformingV3Phase::DomainExtraction);
    CHECK(outcome.destinationUnchanged);
    CHECK(mesh::meshHash(destination) == before);
    CHECK(outcome.reportComplete);
    CHECK(outcome.report && outcome.report->reportComplete);
}

void boundaryConformingV3TwoRegionInterface() {
    TEST("t2b v3 two-region interface has exact owner-side parity");
    Index referenceTetCount = -1;
    for (int presentation = 0; presentation < 3; ++presentation) {
        TwoRegionInterfaceModel model(presentation);
        Provider provider;
        provider.supplied = boxEnvelope(-4,4);
        for (const geo::GFace* face : model.faces())
            provider.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1] =
                model.owners(face->tag());
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(10.0);
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {942,1};
        request.expectedOwnershipProviderStamp = provider.current;
        request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
        mesh::Mesh destination;
        const auto outcome = algo::meshModelBoundaryConforming(
            model,spec,request,destination,provider);
        if (!outcome.published())
            std::printf("  interface presentation=%d status=%d reason=%d "
                        "phase=%d rows=%zu witness=0x%016llx\n",
                presentation,static_cast<int>(outcome.status),
                static_cast<int>(outcome.failure.reason),
                static_cast<int>(outcome.failure.phase),
                outcome.report ? outcome.report->processed.size() : 0U,
                static_cast<unsigned long long>(
                    outcome.failure.exactWitnessDigest));
        CHECK(outcome.published());
        CHECK(outcome.retainedResult != nullptr);
        if (!outcome.retainedResult) continue;
        const auto& result = *outcome.retainedResult;
        CHECK(result.boundaryAudit.passed);
        CHECK(result.boundaryAudit.prescribedInterface > 0);
        CHECK(result.boundaryAudit.matchedInterface ==
              result.boundaryAudit.prescribedInterface);
        CHECK(result.boundaryAudit.missingInterface == 0);
        CHECK(result.boundaryAudit.unexpectedInterface == 0);
        CHECK(result.boundaryAudit.orphanExterior == 0);
        CHECK(result.boundaryAudit.unexpectedExterior == 0);
        CHECK(result.extractionMetrics.retainedRegionComponents == 2);
        std::set<int> cellRegions;
        Index tetCount = 0;
        for (Index cell = 0; cell < destination.cellCount(); ++cell)
            if (destination.cell(cell).type == mesh::CellType::Tet4) {
                cellRegions.insert(destination.cell(cell).onEntity.tag);
                ++tetCount;
            }
        CHECK(cellRegions == std::set<int>({0,1}));
        if (referenceTetCount < 0) referenceTetCount = tetCount;
        CHECK(tetCount == referenceTetCount);
    }
}

void boundaryConformingV3OwnershipAdversaries() {
    TEST("t2b v3 rejects malformed and inconsistent interface ownership");
    struct Case {
        std::vector<bnd::VolumeRegionUse> interfaceOwners;
        algo::BoundaryConformingV3Status status{};
        algo::BoundaryConformingV3Reason reason{};
        algo::BoundaryConformingV3Phase phase{};
    };
    const std::array<Case,6> cases{{
        {{},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::InvalidOwnerMultiplicity,
         algo::BoundaryConformingV3Phase::Ownership},
        {{{0,1},{1,-1},{2,1}},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::InvalidOwnerMultiplicity,
         algo::BoundaryConformingV3Phase::Ownership},
        {{{0,0},{1,-1}},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::InvalidOwnerSense,
         algo::BoundaryConformingV3Phase::Ownership},
        {{{0,1},{0,-1}},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::InvalidOwnerSide,
         algo::BoundaryConformingV3Phase::Ownership},
        {{{0,1},{1,1}},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::InvalidOwnerSide,
         algo::BoundaryConformingV3Phase::Ownership},
        {{{0,1}},algo::BoundaryConformingV3Status::Rejected,
         algo::BoundaryConformingV3Reason::NonClosedShellComponent,
         algo::BoundaryConformingV3Phase::IndexedPreflight}
    }};
    for (std::size_t ci = 0; ci < cases.size(); ++ci) {
        TwoRegionInterfaceModel model;
        Provider provider;
        for (const geo::GFace* face : model.faces())
            provider.rows[
                static_cast<algo::FacePersistentId>(face->tag())+1] =
                model.owners(face->tag());
        provider.rows[12] = cases[ci].interfaceOwners;
        sizing::MeshSpec spec;
        spec.field = std::make_shared<sizing::ConstantField>(10.0);
        algo::BoundaryConformingV3Request request;
        request.expectedSurfaceSetStamp = {950+ci,1};
        request.expectedOwnershipProviderStamp = provider.current;
        request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
        mesh::Mesh destination;
        destination.addNode({-8,-7,-6},{0,990});
        const auto before = mesh::meshHash(destination);
        const auto outcome = algo::meshModelBoundaryConforming(
            model,spec,request,destination,provider);
        if (outcome.status != cases[ci].status ||
            outcome.failure.reason != cases[ci].reason)
            std::printf("  owner case=%zu status=%d reason=%d phase=%d\n",
                ci,static_cast<int>(outcome.status),
                static_cast<int>(outcome.failure.reason),
                static_cast<int>(outcome.failure.phase));
        CHECK(outcome.status == cases[ci].status);
        CHECK(outcome.failure.reason == cases[ci].reason);
        CHECK(outcome.failure.phase == cases[ci].phase);
        CHECK(outcome.destinationUnchanged);
        CHECK(mesh::meshHash(destination) == before);
        CHECK(outcome.reportComplete);
        CHECK(outcome.report && outcome.report->reportComplete);
    }

    // Reversing only the declared interface ownership, while leaving the
    // authoritative directed surface presentation fixed, breaks both
    // region-shell orientations and must stop before recovery.
    TwoRegionInterfaceModel model;
    Provider provider;
    for (const geo::GFace* face : model.faces())
        provider.rows[static_cast<algo::FacePersistentId>(face->tag())+1] =
            model.owners(face->tag());
    provider.rows[12] = {{0,-1},{1,1}};
    sizing::MeshSpec spec;
    spec.field = std::make_shared<sizing::ConstantField>(10.0);
    algo::BoundaryConformingV3Request request;
    request.expectedSurfaceSetStamp = {960,1};
    request.expectedOwnershipProviderStamp = provider.current;
    request.envelopeMode = algo::RecoveryEnvelopeMode::DerivedFiniteBox;
    mesh::Mesh destination;
    destination.addNode({-5,-6,-7},{0,989});
    const auto before = mesh::meshHash(destination);
    const auto conflict = algo::meshModelBoundaryConforming(
        model,spec,request,destination,provider);
    if (conflict.failure.reason !=
        algo::BoundaryConformingV3Reason::NonOrientableShellComponent)
        std::printf("  owner conflict status=%d reason=%d phase=%d\n",
            static_cast<int>(conflict.status),
            static_cast<int>(conflict.failure.reason),
            static_cast<int>(conflict.failure.phase));
    CHECK(conflict.status == algo::BoundaryConformingV3Status::Rejected);
    CHECK(conflict.failure.reason ==
          algo::BoundaryConformingV3Reason::NonOrientableShellComponent);
    CHECK(conflict.failure.phase ==
          algo::BoundaryConformingV3Phase::IndexedPreflight);
    CHECK(conflict.destinationUnchanged);
    CHECK(mesh::meshHash(destination) == before);
    CHECK(conflict.reportComplete);
}

}  // namespace

namespace paacs_tests {
struct Entry { const char* name; void (*fn)(); };
const Entry* t2bTests(std::size_t& count) {
    static const Entry tests[] = {
        {"t2b1.general_cavity.plain",generalEdgeCavityPlainFeasibility},
        {"t2b1.general_cavity.budgets",
         generalEdgeCavityPlainBudgetsAndDeterminism},
        {"t2b.envelope.derived", derivedEnvelopeAllFaces},
        {"t2b.envelope.derived_range", derivedEnvelopeRangeAndDeterminism},
        {"t2b.envelope.supplied", suppliedEnvelopeContract},
        {"t2b.envelope.provider_normalized_collision",
         providerEnvelopeNormalizedCollision},
        {"t2b.provider.failure_checkpoints",
         providerFailureCheckpointsAndInvariance},
        {"t2b.artifact.mapping", artifactMappingsAndStampPrecedence},
        {"t2b.artifact.broad_phase",
         artifactBroadPhaseRejectsIntersectingComponents},
        {"t2b.artifact.shell_taxonomy", artifactShellTaxonomy},
        {"t2b.artifact.broad_phase_scaling", artifactBroadPhaseScaling},
        {"t2b.edge.transaction", targetEdgeTransaction},
        {"t2b.edge.transaction_dim2", targetDimension2EdgeTransaction},
        {"t2b.edge.transaction_dim2_rejections",
         targetDimension2EdgeDomainRejections},
        {"t2b.index.equivalence", indexedQueryEquivalence},
        {"t2b.edge.indexed", indexedEdgeRecovery},
        {"t2b.facet.indexed", indexedFacetRecovery},
        {"t2b.scheduler.interacting", interactingFacetSetRecovery},
        {"t2b.scheduler.mixed_present_equivalence",
         mixedPresentSchedulerEquivalence},
        {"t2b.facet.binding_identity",indexedFacetBindingIdentityAudit},
        {"t2b.scheduler.mixed_present_entry_binding",
         mixedPresentEntryBindingValidationPrecedesSchedule},
        {"t2b.scheduler.mixed_present_permanent_prefix",
         mixedPresentSchedulerPermanentPrefix},
        {"t2b.scheduler.permanent_precedence",
         schedulerPermanentStopPrecedesReadyCommit},
        {"t2b.scheduler.canonical_scc", canonicalConflictScc},
        {"t2b.v3.box", boundaryConformingV3Box},
        {"t2b.v3.surface_stops", boundaryConformingV3SurfaceStops},
        {"t2b.v3.cylinder_thin", boundaryConformingV3CylinderAndThinBox},
        {"t2b.v3.determinism", boundaryConformingV3Determinism},
        {"t2b.v3.sealed_cavity", boundaryConformingV3SealedCavity},
        {"t2b.v3.disconnected_region",
         boundaryConformingV3DisconnectedRegionRefusal},
        {"t2b.v3.two_region_interface",
         boundaryConformingV3TwoRegionInterface},
        {"t2b.v3.ownership_adversaries",
         boundaryConformingV3OwnershipAdversaries},
    };
    count = sizeof(tests) / sizeof(tests[0]);
    return tests;
}
}  // namespace paacs_tests
