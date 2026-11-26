// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctInitialTopologyQuality.hpp"
#include "OcctInitialDelaunayScorer.hpp"
#include "../src/bnd/TetConstraintRecovery.hpp"
#include "paacs/mesh/LinearVolumeAudit.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include <bit>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <new>
#include <set>
#include <string_view>

namespace allocation_probe {
long failAfter=-1,cancelAfter=-1;
std::atomic_bool* cancel=nullptr;
void before(){
    if(cancelAfter>=0&&cancelAfter--==0&&cancel)cancel->store(true);
    if(failAfter>=0&&failAfter--==0)throw std::bad_alloc{};
}
}
void* operator new(std::size_t n){allocation_probe::before();if(auto p=std::malloc(n?n:1))return p;throw std::bad_alloc{};}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p)noexcept{std::free(p);}
void operator delete[](void* p)noexcept{std::free(p);}
void operator delete(void* p,std::size_t)noexcept{std::free(p);}
void operator delete[](void* p,std::size_t)noexcept{std::free(p);}

using namespace paacs;
using namespace paacs::geo::detail;
namespace {
using Carrier=bnd::detail::TetNumericalTopology;
using Complex=bnd::detail::BistellarComplex;
using Cell=bnd::TetFacetCell;
using Vertex=bnd::TetFacetVertexId;
constexpr std::uint64_t MiB=1024ULL*1024;
constexpr mesh::Classification region{3,17};
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

Cell positive(Complex& c,std::array<Vertex,4> v){
    Cell cell{v,0,region};require(c.positive(cell),"owned cavity has a zero-volume tetrahedron");return cell;
}
// Three cells surround an edge through the interior of a triangular link.
// One link edge approaches that central edge, producing a true positive
// sliver. The alternate two-cell bipyramid uses the same five points/shell.
Carrier fixture(double scale=1,unsigned rotation=0,bool good=false){
    Carrier result;result.region=region;result.sourceNodes=5;result.available=true;
    auto& c=result.complex;c.walkConstraints=true;c.maximumStoredCells=128;
    const std::array<Vec3,5> points{{{0,0,-1},{0,0,1},{1,0,0},{-1,1e-15,0},{0,-1,0}}};
    const auto transform=[&](Vec3 p){
        if(rotation==1)p={p.y,p.z,p.x};if(rotation==2)p={-p.z,p.y,p.x};return p*scale;
    };
    for(auto p:points)c.vertices.push_back({transform(p),c.vertices.size()+1});
    for(unsigned i=0;i<8;++i)c.vertices.push_back({transform({(i&1)?8.:-8.,(i&2)?8.:-8.,(i&4)?8.:-8.}),c.vertices.size()+1});
    if(good){c.insert(positive(c,{0,2,3,4}));c.insert(positive(c,{1,2,4,3}));}
    else {c.insert(positive(c,{0,1,2,3}));c.insert(positive(c,{0,1,3,4}));c.insert(positive(c,{0,1,4,2}));}
    for(const auto& [face,owners]:c.faces)if(owners.size()==1){
        const auto& cell=c.cells[*owners.begin()];
        for(unsigned side=0;side<4;++side){const auto directed=bnd::tet_facet_detail::outwardFace(cell,side);
            if(Complex::key(directed)==face){result.sourceFaces.emplace(face,directed);break;}}
    }
    result.liveCells=c.cells.size();return result;
}
std::vector<std::uint64_t> digest(const Carrier& carrier){
    std::vector<std::uint64_t> values;const auto add=[&](auto v){values.push_back(std::uint64_t(v));};
    const auto& c=carrier.complex;
    add(carrier.available);add(carrier.sourceNodes);add(carrier.liveCells);add(carrier.region.dim);add(carrier.region.tag);
    add(c.vertices.size());add(c.cells.size());add(bool(c.chargeWork));
    for(const auto& v:c.vertices){add(std::bit_cast<std::uint64_t>(v.position.x));
        add(std::bit_cast<std::uint64_t>(v.position.y));add(std::bit_cast<std::uint64_t>(v.position.z));add(v.sequence);}
    for(const auto& cell:c.cells){for(auto v:cell.vertices)add(v);add(cell.sequence);add(cell.classification.dim);add(cell.classification.tag);}
    for(auto x:c.alive)add(x);for(const auto& row:c.neighbours)for(auto x:row)add(x);
    for(const auto& [face,owners]:c.faces){for(auto x:face)add(x);add(owners.size());for(auto x:owners)add(x);}
    for(const auto& [edge,owners]:c.edges){for(auto x:edge)add(x);add(owners.size());for(auto x:owners)add(x);}
    for(const auto& face:c.protectedFaces)for(auto x:face)add(x);
    for(const auto& edge:c.protectedEdges)for(auto x:edge)add(x);
    for(const auto& [key,face]:carrier.sourceFaces){for(auto x:key)add(x);for(auto x:face)add(x);}
    return values;
}
void coherent(const Carrier& topology,bool requireAllSourceFaces=true){
    const auto& c=topology.complex;
    std::map<Complex::Face,std::set<std::size_t>> faces;
    std::map<Complex::Edge,std::set<std::size_t>> edges;
    std::size_t live=0;
    for(std::size_t id=0;id<c.cells.size();++id)if(c.alive[id]){
        ++live;const auto& cell=c.cells[id];const auto& v=cell.vertices;
        require(cell.classification==topology.region&&predicates::orient3d(c.vertices[v[0]].position,
            c.vertices[v[1]].position,c.vertices[v[2]].position,c.vertices[v[3]].position)==predicates::Sign::Positive,
            "cleanup changed ownership or exact orientation");
        for(auto face:Complex::faceKeys(cell))faces[face].insert(id);
        for(auto edge:Complex::edgeKeys(cell))edges[edge].insert(id);
    }
    require(live==topology.liveCells&&faces==c.faces&&edges==c.edges,"cleanup incidence/live count is inconsistent");
    for(const auto& [face,owners]:faces){
        require(owners.size()<=2,"cleanup created a nonmanifold face");
        for(auto id:owners){const auto keys=Complex::faceKeys(c.cells[id]);
            const auto side=std::ranges::find(keys,face)-keys.begin();
            const auto other=owners.size()==1?Complex::noCell:(*owners.begin()==id?*owners.rbegin():*owners.begin());
            require(c.neighbours[id][side]==other,"cleanup did not retain reciprocal adjacency");}
    }
    for(const auto& [key,directed]:topology.sourceFaces){
        if(!c.faces.contains(key)){require(!requireAllSourceFaces,"cleanup removed an original source facet");continue;}
        bool found=false;
        for(auto id:c.faces.at(key))for(unsigned side=0;side<4;++side){
            const auto f=bnd::tet_facet_detail::outwardFace(c.cells[id],side);
            for(unsigned j=0;j<3;++j)found|=f[0]==directed[j]&&f[1]==directed[(j+1)%3]&&f[2]==directed[(j+2)%3];}
        require(found,"cleanup changed a source facet's orientation");
    }
}
void ownedSliverAndSimilarity(){
    for(double scale:{1e-6,1.,1e6})for(unsigned rotation=0;rotation<3;++rotation){
        auto topology=fixture(scale,rotation);const auto vertices=topology.complex.vertices;
        const auto source=topology.sourceFaces;
        OcctInitialTopologyQualityLimits limits;limits.maximumWorkingBytes=64*MiB;
        const auto result=finalizeOcctInitialTopologyQuality(topology,limits);
        require(result.completed&&result.changed&&result.edits==1&&result.flips32==1&&
            result.severeBefore==1&&result.severeAfter==0&&result.minimumMeanRatioBefore<1e-8&&
            result.minimumMeanRatioAfter>.5&&result.minimumScaledJacobianAfter>.3,
            "generic three-to-two cleanup did not remove the owned near-coplanar outlier");
        require(result.work<=limits.maximumWork&&result.maximumModeledBytes<=limits.maximumWorkingBytes&&
            topology.liveCells==2&&!topology.complex.chargeWork&&topology.sourceFaces==source&&
            topology.complex.protectedFaces.empty()&&topology.complex.protectedEdges.empty(),
            "cleanup exceeded resources or leaked temporary source protection/callbacks");
        for(std::size_t i=0;i<vertices.size();++i)require(
            std::bit_cast<std::uint64_t>(vertices[i].position.x)==std::bit_cast<std::uint64_t>(topology.complex.vertices[i].position.x)&&
            std::bit_cast<std::uint64_t>(vertices[i].position.y)==std::bit_cast<std::uint64_t>(topology.complex.vertices[i].position.y)&&
            std::bit_cast<std::uint64_t>(vertices[i].position.z)==std::bit_cast<std::uint64_t>(topology.complex.vertices[i].position.z)&&
            vertices[i].sequence==topology.complex.vertices[i].sequence,"cleanup moved or renumbered a vertex");
        coherent(topology);
    }
}
void noOpAndConstraints(){
    auto good=fixture(1,0,true);const auto before=digest(good);
    auto result=finalizeOcctInitialTopologyQuality(good);
    require(result.completed&&!result.changed&&std::string_view(result.diagnostic)=="no_severe_candidates"&&
        digest(good)==before,"good topology was copied/reordered or edited");
    for(bool face:{false,true}){
        auto topology=fixture();
        if(face)topology.complex.protectedFaces.insert(Complex::key(Complex::Face{0,1,2}));
        else topology.complex.protectedEdges.insert(Complex::Edge{0,1});
        const auto original=digest(topology);result=finalizeOcctInitialTopologyQuality(topology);
        require(result.completed&&!result.changed&&digest(topology)==original,
            "cleanup erased a protected internal constraint");
    }
}
void candidateTruncation(){
    auto topology=fixture();auto& c=topology.complex;
    std::array<Vec3,8> guards{};for(unsigned i=0;i<8;++i)guards[i]=c.vertices[i+5].position;
    c.vertices.resize(5);
    for(unsigned i=0;i<5;++i)c.vertices.push_back({c.vertices[i].position+Vec3{4,0,0},c.vertices.size()+1});
    for(auto p:guards)c.vertices.push_back({p,c.vertices.size()+1});
    for(const std::array<Vertex,4> v:std::array<std::array<Vertex,4>,3>{{{5,6,7,8},{5,6,8,9},{5,6,9,7}}})
        c.insert(positive(c,v));
    topology.sourceNodes=10;topology.liveCells=6;topology.sourceFaces.clear();
    for(const auto& [face,owners]:c.faces)if(owners.size()==1){
        const auto& cell=c.cells[*owners.begin()];
        for(unsigned side=0;side<4;++side){const auto directed=bnd::tet_facet_detail::outwardFace(cell,side);
            if(Complex::key(directed)==face){topology.sourceFaces.emplace(face,directed);break;}}
    }
    OcctInitialTopologyQualityLimits limits;limits.maximumCandidates=1;
    const auto result=finalizeOcctInitialTopologyQuality(topology,limits);
    require(result.completed&&result.changed&&result.edits==1&&result.severeBefore==2&&result.severeAfter==1&&
        result.candidateListTruncated&&result.budgetLimited&&std::string_view(result.diagnostic)=="candidate_budget",
        "truncated candidate search claimed unbounded completion or discarded its certified edit");
    coherent(topology);
}
void boundedAndTransactional(){
    const auto original=fixture();
    auto full=original;const auto completed=finalizeOcctInitialTopologyQuality(full);
    require(completed.changed&&completed.work>32,"owned work-bound control did not complete");
    for(std::uint64_t work:std::array<std::uint64_t,4>{0,1,32,completed.work-1}){
        auto topology=original;const auto before=digest(topology);OcctInitialTopologyQualityLimits limits;
        limits.maximumWork=work;const auto result=finalizeOcctInitialTopologyQuality(topology,limits);
        require(result.work<=work&&!result.changed&&digest(topology)==before,
            "work refusal mutated the original or exceeded the requested cap");
    }
    for(std::uint64_t bytes:std::array<std::uint64_t,3>{0,8*MiB,2*original.retainedModeledBytes()-1}){
        auto topology=original;const auto before=digest(topology);OcctInitialTopologyQualityLimits limits;
        limits.maximumWorkingBytes=bytes;const auto result=finalizeOcctInitialTopologyQuality(topology,limits);
        require(result.budgetLimited&&!result.changed&&digest(topology)==before,
            "memory refusal mutated the original private topology");
    }
    auto disabled=original;const auto before=digest(disabled);OcctInitialTopologyQualityLimits limits;
    limits.maximumEdits=0;const auto result=finalizeOcctInitialTopologyQuality(disabled,limits);
    require(!result.changed&&digest(disabled)==before,"disabled cleanup changed connectivity");
}
void allocationAndCancellation(){
    unsigned allocationStops=0;
    for(long fail=0;fail<24;++fail){
        auto topology=fixture();const auto before=digest(topology);
        allocation_probe::failAfter=fail;
        const auto result=finalizeOcctInitialTopologyQuality(topology);
        allocation_probe::failAfter=-1;
        require(!result.changed&&std::string_view(result.diagnostic)=="allocation_failed"&&digest(topology)==before,
            "allocation refusal leaked a partially edited owner");++allocationStops;
    }
    require(allocationStops==24,"allocation failure coverage did not reach staged-copy allocations");
    auto topology=fixture();const auto before=digest(topology);std::atomic_bool stop=true;
    OcctInitialTopologyQualityLimits limits;limits.cancellation=&stop;
    auto result=finalizeOcctInitialTopologyQuality(topology,limits);
    require(result.cancelled&&!result.changed&&result.work==0&&digest(topology)==before,
        "pre-cancelled cleanup performed work or mutated topology");
    stop=false;allocation_probe::cancel=&stop;allocation_probe::cancelAfter=4;
    result=finalizeOcctInitialTopologyQuality(topology,limits);
    allocation_probe::cancelAfter=-1;allocation_probe::cancel=nullptr;
    require(result.cancelled&&!result.changed&&result.work>0&&result.work<=limits.maximumWork&&digest(topology)==before,
        "cancellation inside staged allocation did not preserve the old owner and work receipt");
}
void archiveDiagnostic(const char* path,double threshold,const char* outputDirectory){
    sdk::MeshArchiveLimits archiveLimits;archiveLimits.maxBytes=128*MiB;
    archiveLimits.maxNodes=100'000;archiveLimits.maxCells=500'000;archiveLimits.maxConnectivity=4'000'000;
    mesh::Mesh saved,boundary;const auto loaded=sdk::readMeshArchive(path,saved,archiveLimits);
    require(loaded.completed,"could not read diagnostic archive");
    int tag=-1;
    for(Index i=0;i<saved.nodeCount();++i)boundary.addNode(saved.node(i),saved.nodeClass(i));
    for(Index i=0;i<saved.cellCount();++i){const auto& cell=saved.cell(i);const auto* v=saved.cellNodes(i);
        if(cell.type==mesh::CellType::Tri3)boundary.addCell(cell.type,std::span<const Index>(v,3),cell.onEntity);
        if(cell.type==mesh::CellType::Tet4){if(tag<0)tag=cell.onEntity.tag;require(tag==cell.onEntity.tag,"diagnostic archive must contain one body");}}
    require(tag>=0&&boundary.cellCount()>0,"diagnostic archive has no body/source triangles");
    OcctInitialDelaunayScorer scorer;std::uint64_t bootstrap=0;
    require(scorer.initialize(boundary,{3,tag},100'000,384*MiB,[&](std::uint64_t n){
        require(n<=50'000'000-bootstrap,"diagnostic bootstrap exceeded its cap");bootstrap+=n;},
        {},0,50'000'000,nullptr,true),"reconstructed diagnostic bootstrap refused");
    Carrier carrier;require(scorer.releaseNumericalTopology(carrier,{3,tag}),"diagnostic carrier release failed");
    OcctInitialTopologyQualityLimits limits;limits.maximumWork=30'000'000;limits.maximumWorkingBytes=384*MiB;
    limits.searchBelowMeanRatio=threshold;
    const auto result=finalizeOcctInitialTopologyQuality(carrier,limits);
    std::cout<<"{\"scope\":\"reconstructed_bootstrap_not_whole_pipeline_trace\",\"source_nodes\":"<<boundary.nodeCount()
        <<",\"source_triangles\":"<<boundary.cellCount()<<",\"bootstrap_work\":"<<bootstrap
        <<",\"changed\":"<<(result.changed?"true":"false")<<",\"diagnostic\":\""<<result.diagnostic
        <<"\",\"severe_before\":"<<result.severeBefore<<",\"severe_after\":"<<result.severeAfter
        <<",\"minimum_mean_before\":"<<result.minimumMeanRatioBefore<<",\"minimum_mean_after\":"<<result.minimumMeanRatioAfter
        <<",\"flips32\":"<<result.flips32<<",\"edits\":"<<result.edits<<",\"work\":"<<result.work
        <<",\"peak_modeled_bytes\":"<<result.maximumModeledBytes<<"}\n";
    require(result.completed,"diagnostic finalization did not complete");coherent(carrier,false);
    mesh::Mesh seed;const auto recovery=bnd::detail::recoverTetConstraints(boundary,{3,tag},50'000'000,384*MiB,
        1'000'000,nullptr,seed,0,false,true,true,true,&carrier);
    require(recovery.initialTopologyUsed&&std::string_view(recovery.diagnostic)=="complete_seed",
        "diagnostic source constraints did not recover completely");
    std::uint64_t severe=0;double minimum=1;std::uint64_t scan=0;
    Complex measured;for(const auto p:seed.nodes())measured.vertices.push_back({p,measured.vertices.size()+1});
    for(Index i=0;i<seed.cellCount();++i)if(seed.cell(i).type==mesh::CellType::Tet4){
        const auto* v=seed.cellNodes(i);Cell c{{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,{3,tag}};
        const auto q=occt_initial_topology_quality_detail::shape(measured,c,[&](std::uint64_t n){scan+=n;});
        severe+=q.mean<limits.searchBelowMeanRatio;minimum=std::min(minimum,q.mean);
    }
    std::cout<<"{\"scope\":\"reconstructed_material_seed_not_Native_CAD_audit\",\"tetrahedra\":"<<seed.cellCount()
        <<",\"severe\":"<<severe<<",\"minimum_mean_ratio\":"<<minimum<<",\"work\":"<<recovery.work<<"}\n";
    require(seed.nodeCount()<=100'000&&seed.cellCount()<=500'000-boundary.cellCount(),"diagnostic output exceeds its count allowance");
    for(Index i=0;i<boundary.nodeCount();++i){
        require(seed.nodeClass(i)==boundary.nodeClass(i),"diagnostic changed source node classification");
        for(unsigned axis=0;axis<3;++axis){
            const auto coordinate=[](Vec3 p,unsigned a){return a==0?p.x:a==1?p.y:p.z;};
            require(std::bit_cast<std::uint64_t>(coordinate(seed.node(i),axis))==
                std::bit_cast<std::uint64_t>(coordinate(boundary.node(i),axis)),"diagnostic changed a source coordinate");
        }
    }
    const auto tetrahedra=seed.cellCount();
    for(Index i=0;i<boundary.cellCount();++i){const auto* v=boundary.cellNodes(i);
        seed.addCell(mesh::CellType::Tri3,std::span<const Index>(v,3),boundary.cell(i).onEntity);}
    require(mesh::auditLinearVolumeBoundary(seed).valid,"diagnostic seed disagrees with the tagged source boundary");
    if(outputDirectory){
        const std::filesystem::path destination(outputDirectory);
        require(std::filesystem::is_directory(destination.parent_path())&&!std::filesystem::exists(destination)&&
            std::filesystem::create_directory(destination),"diagnostic output must be a fresh directory under an existing parent");
        require(sdk::writeMeshArchiveFile(seed,destination/"seed-with-boundary.paacsmesh").completed,
            "could not save diagnostic seed archive");
        std::ofstream receipt(destination/"receipt.json",std::ios::binary);
        require(bool(receipt),"could not create diagnostic receipt");
        receipt<<std::setprecision(17)
            <<"{\n  \"schema\": \"paacs-private-topology-shape-diagnostic-v1\",\n"
            <<"  \"scope\": \"reconstructed_bootstrap_and_constrained_material_seed_not_whole_Native_pipeline\",\n"
            <<"  \"input_archive_filename\": "<<std::quoted(std::filesystem::path(path).filename().string())<<",\n"
            <<"  \"output_archive_filename\": \"seed-with-boundary.paacsmesh\",\n"
            <<"  \"native_CAD_audit_rerun\": false,\n  \"independent_intersection_audit_performed\": false,\n"
            <<"  \"source_node_coordinates_and_tags_unchanged\": true,\n  \"source_triangles_and_tags_unchanged\": true,\n"
            <<"  \"linear_volume_boundary_audit_valid\": true,\n"
            <<"  \"body_tag\": "<<tag<<",\n  \"source_nodes\": "<<boundary.nodeCount()
            <<",\n  \"source_triangles\": "<<boundary.cellCount()<<",\n  \"bootstrap_work\": "<<bootstrap
            <<",\n  \"requested_search_below_mean_ratio\": "<<threshold
            <<",\n  \"finalization\": {\n    \"diagnostic\": "<<std::quoted(result.diagnostic)
            <<",\n    \"changed\": "<<(result.changed?"true":"false")
            <<",\n    \"candidate_list_truncated\": "<<(result.candidateListTruncated?"true":"false")
            <<",\n    \"budget_limited\": "<<(result.budgetLimited?"true":"false")
            <<",\n    \"severe_before\": "<<result.severeBefore<<",\n    \"severe_after\": "<<result.severeAfter
            <<",\n    \"minimum_mean_ratio_before\": "<<result.minimumMeanRatioBefore
            <<",\n    \"minimum_mean_ratio_after\": "<<result.minimumMeanRatioAfter
            <<",\n    \"minimum_scaled_jacobian_before\": "<<result.minimumScaledJacobianBefore
            <<",\n    \"minimum_scaled_jacobian_after\": "<<result.minimumScaledJacobianAfter
            <<",\n    \"edits\": "<<result.edits<<",\n    \"flips23\": "<<result.flips23
            <<",\n    \"flips32\": "<<result.flips32<<",\n    \"flips44\": "<<result.flips44
            <<",\n    \"work\": "<<result.work<<",\n    \"peak_modeled_bytes\": "<<result.maximumModeledBytes
            <<"\n  },\n  \"constrained_material_seed\": {\n    \"diagnostic\": "<<std::quoted(recovery.diagnostic)
            <<",\n    \"tetrahedra\": "<<tetrahedra<<",\n    \"cells_below_requested_trigger\": "<<severe
            <<",\n    \"minimum_mean_ratio\": "<<minimum<<",\n    \"work\": "<<recovery.work<<"\n  }\n}\n";
        receipt.close();require(bool(receipt),"could not finish diagnostic receipt");
    }
}
}
int main(int argc,char** argv){try{
    if(argc>=3&&argc<=5&&std::string_view(argv[1])=="--archive"){
        double threshold=1e-6;
        if(argc>=4){std::size_t consumed=0;threshold=std::stod(argv[3],&consumed);
            require(consumed==std::string_view(argv[3]).size()&&std::isfinite(threshold)&&threshold>0&&threshold<1,
                "diagnostic threshold must be finite and between zero and one");}
        archiveDiagnostic(argv[2],threshold,argc==5?argv[4]:nullptr);return 0;
    }
    require(argc==1,"expected no arguments or --archive PATH [THRESHOLD] [NEW_OUTPUT_DIRECTORY]");
    ownedSliverAndSimilarity();noOpAndConstraints();candidateTruncation();boundedAndTransactional();allocationAndCancellation();
    std::cout<<"initial topology shape finalization, source protection and transactional budgets passed\n";return 0;
}catch(const std::exception& error){allocation_probe::failAfter=-1;allocation_probe::cancelAfter=-1;
    std::cerr<<error.what()<<'\n';return 1;}
catch(...){allocation_probe::failAfter=-1;allocation_probe::cancelAfter=-1;
    std::cerr<<"unexpected private topology-quality exception\n";return 1;}}
