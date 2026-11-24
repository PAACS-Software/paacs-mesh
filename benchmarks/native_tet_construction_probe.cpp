// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// Private, opt-in diagnostics. Numerical snapshots are NOT material meshes.
#include "paacs/geometry/OcctMesher.hpp"
#include "paacs/algo/GeometryMesher.hpp"
#include "paacs/mesh/Quality.hpp"
#include "paacs/sdk/MeshArchive.hpp"
#include "../src/bnd/TetBistellarRecovery.hpp"
#include "../src/bnd/TetConstructionProbe.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {
using namespace paacs;
using Key=std::array<std::array<std::uint64_t,3>,4>;
Key key(const std::array<Vec3,4>& p) {
    Key result{};
    for(unsigned j=0;j<4;++j) {
        const std::array coordinates{p[j].x,p[j].y,p[j].z};
        for(unsigned a=0;a<3;++a)result[j][a]=coordinates[a]==0?0:std::bit_cast<std::uint64_t>(coordinates[a]);
    }
    std::sort(result.begin(),result.end());return result;
}
std::pair<double,double> shape(const std::array<Vec3,4>& p) {
    const auto det=predicates::orient3dDeterminantEstimate(p[0],p[1],p[2],p[3]);
    double sum=0;for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j)sum+=(p[i]-p[j]).norm2();
    const auto mr=det>0&&sum>0?12*std::pow(det/2.,2./3.)/sum:0;
    return {std::isfinite(mr)?mr:0,mesh::tetScaledJacobian(p[0],p[1],p[2],p[3])};
}
void number(std::ostream& output,double value) { if(std::isfinite(value))output<<value;else output<<"null"; }
struct Probe {
    std::ofstream output;
    std::filesystem::path directory;
    std::vector<Key> targets;
    std::uint64_t visits=0,events=0,records=0;
    bool truncated=false,failed=false;
    explicit Probe(const std::filesystem::path& path):output(path/"stages.jsonl"),directory(path) {
        output<<std::setprecision(17);
    }
    void dump(const bnd::detail::TetConstructionObservation& view) {
        if(!view.completeSnapshot||!view.complex||std::string_view(view.stage)!="initial_scoring_final")return;
        const auto& complex=*view.complex;
        if(complex.vertices.size()>100000||complex.cells.size()>500000) { truncated=true;return; }
        mesh::Mesh snapshot;
        for(const auto& v:complex.vertices)snapshot.addNode(v.position,{3,view.bodyTag});
        for(std::size_t i=0;i<complex.cells.size();++i)if(complex.alive[i]) {
            const auto& cell=complex.cells[i];
            const std::array<Index,4> ids{Index(cell.vertices[0]),Index(cell.vertices[1]),Index(cell.vertices[2]),Index(cell.vertices[3])};
            snapshot.addCell(mesh::CellType::Tet4,ids,{3,view.bodyTag});
        }
        if(!sdk::writeMeshArchiveFile(snapshot,directory/"initial_unconstrained_with_guards.paacsmesh").completed)failed=true;
    }
    void observe(const bnd::detail::TetConstructionObservation& view) {
        ++events;if(visits>=2'000'000||records>=8192){truncated=true;return;}
        const auto total=view.complex?view.complex->cells.size():std::size_t(view.mesh->cellCount());
        const auto first=view.completeSnapshot?std::size_t{0}:view.firstCell;
        const auto last=view.completeSnapshot?total:std::min(total,first+view.cellCount);
        std::uint64_t scanned=0,guards=0,severe=0;double minimum=1,minimumJac=1;
        std::array<std::uint64_t,5> relevance{};
        std::vector<std::size_t> found;
        struct Poor { std::size_t id;std::array<std::size_t,4> ids;std::array<Vec3,4> p;double mean,jac;unsigned state; };
        std::vector<Poor> worst;
        for(std::size_t i=first;i<last;++i) {
            if(++visits>2'000'000){truncated=true;break;}
            std::array<Vec3,4> p;std::array<std::size_t,4> ids{};bool guard=false;unsigned state=0;
            if(view.complex) {
                if(!view.complex->alive[i])continue;
                const auto& c=view.complex->cells[i];
                for(unsigned j=0;j<4;++j) {
                    ids[j]=c.vertices[j];p[j]=view.complex->vertices[ids[j]].position;
                    guard|=ids[j]>=view.guardBegin&&ids[j]-view.guardBegin<view.guardCount;
                }
                if(i<view.relevance.size())state=view.relevance[i];
            }else {
                if(view.mesh->cell(Index(i)).type!=mesh::CellType::Tet4)continue;
                const auto* n=view.mesh->cellNodes(Index(i));
                for(unsigned j=0;j<4;++j){ids[j]=std::size_t(n[j]);p[j]=view.mesh->node(n[j]);}
            }
            if(guard){++guards;continue;}
            ++scanned;if(state<relevance.size())++relevance[state];
            const auto q=shape(p);minimum=std::min(minimum,q.first);minimumJac=std::min(minimumJac,q.second);
            const auto signature=key(p);
            for(std::size_t j=0;j<targets.size();++j)if(signature==targets[j])found.push_back(j);
            if(q.first<1e-8){++severe;if(worst.size()<16)worst.push_back({i,ids,p,q.first,q.second,state});}
        }
        if(view.completeSnapshot||severe||!found.empty()) {
            ++records;
            output<<"{\"stage\":\""<<view.stage<<"\",\"body\":"<<view.bodyTag
                <<",\"complete_snapshot\":"<<(view.completeSnapshot?"true":"false")
                <<",\"original_coordinate_retry\":"<<(view.originalCoordinateRetryUsed?"true":"false")
                <<",\"guard_begin\":"<<view.guardBegin<<",\"guard_count\":"<<view.guardCount
                <<",\"scanned_nonguard_tets\":"<<scanned<<",\"guards_excluded\":"<<guards
                <<",\"minimum_mean_ratio\":";number(output,minimum);
            output<<",\"minimum_scaled_jacobian\":";number(output,minimumJac);
            output<<",\"severe_count\":"<<severe<<",\"relevance_counts\":[";
            for(unsigned j=0;j<relevance.size();++j){if(j)output<<',';output<<relevance[j];}
            output<<"],\"targets_present\":[";
            for(std::size_t j=0;j<found.size();++j){if(j)output<<',';output<<found[j];}
            output<<"],\"severe_samples\":[";
            for(std::size_t j=0;j<worst.size();++j){if(j)output<<',';const auto& row=worst[j];
                output<<"{\"cell\":"<<row.id<<",\"relevance\":"<<row.state<<",\"mean_ratio\":";number(output,row.mean);
                output<<",\"scaled_jacobian\":";number(output,row.jac);output<<",\"nodes\":[";
                for(unsigned n=0;n<4;++n){if(n)output<<',';output<<row.ids[n];}
                output<<"],\"points\":[";
                for(unsigned n=0;n<4;++n){if(n)output<<',';output<<'[';number(output,row.p[n].x);output<<',';number(output,row.p[n].y);output<<',';number(output,row.p[n].z);output<<']';}
                output<<"]}";
            }
            output<<"]}\n";output.flush();
        }
        dump(view);
    }
    static void callback(const bnd::detail::TetConstructionObservation& view,void* context) noexcept {
        auto& self=*static_cast<Probe*>(context);try{self.observe(view);}catch(...){self.failed=true;}
    }
};
}
int main(int argc,char** argv) {
    using namespace paacs;
    if(argc==4&&std::string_view(argv[1])=="--analyze-planar") {
        auto cad=geo::importOcctCad(argv[2]);if(!cad.imported())return 2;
        algo::GeometryMeshingOptions options;options.targetSize=std::stod(argv[3]);
        for(const auto* body:cad.model->solids()) {
            algo::GeometryBodyControl control;control.solidTag=body->tag();control.elements=algo::GeometryElements::Tet4;
            options.bodyControls.push_back(control);
        }
        const auto result=algo::analyzeGeometryStrategies(*cad.model,options);
        std::cout<<"{\"complete\":"<<(result.complete?"true":"false")<<",\"reason\":\""<<algo::geometryReasonName(result.reason)<<"\",\"bodies\":[";
        for(std::size_t i=0;i<result.bodies.size();++i){if(i)std::cout<<',';const auto& b=result.bodies[i];
            std::cout<<"{\"body\":"<<b.solidTag<<",\"strategy\":\""<<algo::geometryStrategyName(b.strategy)
                <<"\",\"estimated_nodes\":"<<b.estimatedNodes<<",\"estimated_cells\":"<<b.estimatedCells
                <<",\"exact_count_estimate\":"<<(b.exactCountEstimate?"true":"false")<<"}";
        }
        std::cout<<"]}\n";return result.complete?0:1;
    }
    if(argc<7||argc>8){std::cerr<<"usage: probe source.step global-size wall-size output-directory reference.paacsmesh target-ids-comma-separated [seconds]\n";return 2;}
    try {
        const std::filesystem::path directory=argv[4];
        if(std::filesystem::exists(directory)||!std::filesystem::create_directories(directory))throw std::runtime_error("output directory must be new");
        Probe probe(directory);if(!probe.output)throw std::runtime_error("cannot create probe output");
        mesh::Mesh reference;sdk::MeshArchiveLimits limits;
        limits.maxNodes=100000;limits.maxCells=500000;limits.maxBytes=128ULL*1024*1024;
        if(!sdk::readMeshArchiveFile(argv[5],reference,limits).completed)throw std::runtime_error("cannot load reference");
        std::istringstream requested(argv[6]);std::string token;
        while(std::getline(requested,token,',')) {
            std::size_t parsed=0;const auto raw=std::stoll(token,&parsed);
            if(parsed!=token.size()||raw<0||raw>=reference.cellCount()||probe.targets.size()>=32)
                throw std::runtime_error("reference target outside probe limit");
            const auto id=Index(raw);
            if(reference.cell(id).type!=mesh::CellType::Tet4)throw std::runtime_error("reference target unavailable");
            const auto* ids=reference.cellNodes(id);std::array<Vec3,4> p;
            for(unsigned j=0;j<4;++j)p[j]=reference.node(ids[j]);probe.targets.push_back(key(p));
        }
        auto cad=geo::importOcctCad(argv[1]);if(!cad.imported())throw std::runtime_error("CAD import failed");
        auto options=geo::normalNativeCadMeshOptions();options.geometry.targetSize=std::stod(argv[2]);
        options.geometry.threads=1;options.maximumNativeValidationThreads=1;
        options.geometry.maximumNodes=500000;options.geometry.maximumCells=3000000;
        options.geometry.maximumOutputBytes=256ULL*1024*1024;options.geometry.maximumAuditTemporaryBytes=512ULL*1024*1024;
        options.maximumPrimitiveWorkingBytes=512ULL*1024*1024;
        int wall=-1;for(const auto* face:cad.model->faces())if(face->kind()==geo::SurfaceKind::Cylinder){if(wall!=-1)throw std::runtime_error("wall selection ambiguous");wall=face->tag();}
        if(wall<0)throw std::runtime_error("no cylindrical wall");
        options.faceSizes={{wall,std::stod(argv[3]),1.3}};options.useSpatialSizing=true;
        if(probe.targets.empty())throw std::runtime_error("at least one reference target is required");
        const auto seconds=argc==8?std::stod(argv[7]):300.;if(!(seconds>0&&seconds<=3600))throw std::runtime_error("seconds outside probe limit");
        std::atomic_bool cancellation=false;options.geometry.cancellation=&cancellation;
        std::jthread deadline([&](std::stop_token stop){
            const auto until=std::chrono::steady_clock::now()+std::chrono::duration<double>(seconds);
            while(!stop.stop_requested()&&std::chrono::steady_clock::now()<until)std::this_thread::sleep_for(std::chrono::milliseconds(10));
            if(!stop.stop_requested())cancellation.store(true);
        });
        bnd::detail::ScopedTetConstructionObserver observing(Probe::callback,&probe);
        mesh::Mesh destination;const auto outcome=geo::meshNativeCad(*cad.model,options,destination);
        deadline.request_stop();
        std::ofstream result(directory/"result.json");
        result<<"{\"published\":"<<(outcome.published?"true":"false")<<",\"reason\":\""<<geo::nativeCadMeshReasonName(outcome.reason)
            <<"\",\"probe_failed\":"<<(probe.failed?"true":"false")<<",\"probe_truncated\":"<<(probe.truncated?"true":"false")
            <<",\"probe_visits\":"<<probe.visits<<",\"probe_events\":"<<probe.events<<",\"probe_records\":"<<probe.records
            <<",\"scope\":\"private numerical states; not final material or solver certificates\"}\n";
        if(outcome.published&&!sdk::writeMeshArchiveFile(destination,directory/"published.paacsmesh").completed)return 2;
        return outcome.published&&!probe.failed&&!probe.truncated?0:1;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
}
