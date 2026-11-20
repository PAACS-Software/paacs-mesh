// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/algo/ConvexToTetMesher.hpp>
#include <paacs/algo/LayeredPatchMesher.hpp>
#include <paacs/algo/MixedMeshQuality.hpp>
#include <paacs/geometry/OcctMesher.hpp>
#include <paacs/mesh/LinearVolumeAudit.hpp>
#include <paacs/sdk/MeshArchive.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace paacs;

namespace {
constexpr double halfWidth=5., length=40., tolerance=1e-8;
constexpr int coreInterfaceTag=0;
constexpr std::array groupNames{"inlet","outlet","x_minus","x_plus","y_minus","y_plus"};

void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}

bool close(double a,double b) {
    return std::abs(a-b)<tolerance;
}

struct Duct {
    const geo::OcctModel& model;
    const geo::GSolid& body;
    std::array<int,6> tags{};
    std::vector<const geo::GEdge*> edges;
    std::vector<const geo::GVertex*> vertices;

    explicit Duct(const geo::OcctModel& source)
        :model(source),body(*source.solids().front()),edges(source.edges()),vertices(source.vertices()) {
        Vec3 lo,hi;
        body.bbox(lo,hi);
        require((lo-Vec3{-halfWidth,-halfWidth,0.}).norm()<1e-5&&
                (hi-Vec3{halfWidth,halfWidth,length}).norm()<1e-5,"Unexpected square duct dimensions");
        for (auto tag:body.boundaryFaceTags()) {
            const auto& face=*model.face(tag);
            require(face.kind()==geo::SurfaceKind::Plane,"Expected six original planar faces");
            const auto on=planes(face.pointAt({0.,0.}));
            const auto normal=face.normalAt({0.,0.});
            const auto axis=std::abs(normal.z)>.99 ? 0u : std::abs(normal.x)>.99 ? 2u : 4u;
            require(std::max({std::abs(normal.x),std::abs(normal.y),std::abs(normal.z)})>.99 &&
                    (on[axis]||on[axis+1]),"Unexpected original box face");
            const auto index=axis+!on[axis];
            require(!tags[index],"Duplicate original box face");
            tags[index]=tag;
        }
        require(std::ranges::none_of(tags,[](int tag){return tag<=0;}),"Missing original CAD face");
    }

    static std::array<bool,6> planes(Vec3 point) {
        return {close(point.z,0.),close(point.z,length),close(point.x,-halfWidth),
                close(point.x,halfWidth),close(point.y,-halfWidth),close(point.y,halfWidth)};
    }

    mesh::Classification classify(Vec3 point) const {
        const auto on=planes(point);
        const auto count=std::ranges::count(on,true);
        if (!count) return {3,body.tag()};
        if (count==1) return {2,tags[std::ranges::find(on,true)-on.begin()]};
        if (count==2) {
            for (auto edge:edges)
                if ((edge->closestPoint(point)-point).norm()<tolerance) return {1,edge->tag()};
            throw std::runtime_error("Could not classify original CAD edge point");
        }
        for (auto vertex:vertices)
            if ((vertex->position()-point).norm()<tolerance) return {0,vertex->tag()};
        throw std::runtime_error("Could not classify original CAD vertex");
    }

    int boundaryTag(const mesh::Mesh& mesh,std::span<const Index> ids) const {
        auto on=planes(mesh.node(ids.front()));
        for (const auto id:ids.subspan(1)) {
            const auto next=planes(mesh.node(id));
            for (unsigned side=0;side<on.size();++side) on[side]=on[side]&&next[side];
        }
        require(std::ranges::count(on,true)<=1,"Degenerate exterior facet");
        const auto i=std::ranges::find(on,true);
        return i==on.end() ? 0 : tags[i-on.begin()];
    }
};

mesh::Mesh coreTetrahedra(const Duct& duct,double half,unsigned intervals,unsigned axialLayers) {
    const auto stride=Index(intervals+1);
    const auto id=[&](Index i,Index j){return stride*j+i;};
    mesh::Mesh section, hexes, tets;
    for (Index j=0;j<stride;++j) for (Index i=0;i<stride;++i)
        section.addNode({-half+2*half*i/intervals,-half+2*half*j/intervals,0.},{2,duct.tags[0]});
    for (Index j=0;j+1<stride;++j) for (Index i=0;i+1<stride;++i)
        section.addCell(mesh::CellType::Quad4,{id(i,j),id(i+1,j),id(i+1,j+1),id(i,j+1)},{2,duct.tags[0]});
    algo::LayeredPatchRequest extrusion;
    extrusion.volumeClassification={3,duct.body.tag()};
    extrusion.targetClassification={2,duct.tags[1]};
    extrusion.sideClassification={2,coreInterfaceTag};
    extrusion.schedule.wallLayers=axialLayers;
    extrusion.schedule.firstLayerThickness=length/axialLayers;
    extrusion.schedule.growthRatio=1.;
    extrusion.schedule.totalDepth=length;
    extrusion.schedule.maximumCoreThickness=length/axialLayers;
    const auto swept=algo::meshLayeredPatch(section,extrusion,hexes);
    require(swept.published,"PAACS could not sweep the inner structured core");
    const auto converted=algo::convertConvexCellsToTet4(hexes,{},tets);
    require(converted.published,"PAACS could not convert the core into conforming Tet4");
    return tets;
}

void exteriorFace(const Duct& duct,mesh::Mesh& result,std::span<const Index> vertices) {
    const int tag=duct.boundaryTag(result,vertices);
    if (tag) result.addCell(vertices.size()==3 ? mesh::CellType::Tri3 : mesh::CellType::Quad4,vertices,{2,tag});
}

void sweepWallPrisms(const Duct& duct,const mesh::Mesh& core,const algo::LayerSchedule& walls,
                    mesh::Mesh& result) {
    const auto half=halfWidth-walls.wallThickness;
    std::vector<std::vector<Index>> rings(walls.offsets.size()-1,
        std::vector<Index>(core.nodeCount(),kInvalidIndex));
    for (Index i=0;i<core.nodeCount();++i) {
        const auto point=core.node(i);
        if (!close(std::max(std::abs(point.x),std::abs(point.y)),half)) continue;
        for (std::size_t layer=0;layer<rings.size();++layer) {
            const auto scale=(halfWidth-walls.offsets[layer])/half;
            const Vec3 projected{scale*point.x,scale*point.y,point.z};
            rings[layer][i]=result.addNode(projected,duct.classify(projected));
        }
    }
    const auto vertex=[&](std::size_t layer,Index i) {
        const auto node=layer==rings.size() ? i : rings[layer][i];
        require(node!=kInvalidIndex,"Expected core wall ring node");
        return node;
    };
    for (Index i=0;i<core.cellCount();++i) {
        const auto face=core.cell(i);
        if (face.type!=mesh::CellType::Tri3||face.onEntity.tag!=coreInterfaceTag) continue;
        const auto* src=core.cellNodes(i);
        // Core face is outward from the core. Reverse its order for an inward prism sweep.
        const std::array<Index,3> inward{src[0],src[2],src[1]};
        for (std::size_t layer=0;layer<rings.size();++layer) {
            std::array<Index,6> prism{};
            for (std::size_t j=0;j<inward.size();++j) {
                prism[j]=vertex(layer,inward[j]);
                prism[j+3]=vertex(layer+1,inward[j]);
            }
            result.addCell(mesh::CellType::Prism6,prism,{3,duct.body.tag()});
            if (layer==0) exteriorFace(duct,result,std::array{prism[0],prism[2],prism[1]});
            for (unsigned edge=0;edge<3;++edge) {
                const auto j=(edge+1)%3;
                exteriorFace(duct,result,std::array{prism[edge],prism[j],prism[j+3],prism[edge+3]});
            }
        }
    }
}

mesh::Mesh hybrid(const Duct& duct,const algo::LayerSchedule& walls,unsigned intervals,unsigned axialLayers) {
    const auto core=coreTetrahedra(duct,halfWidth-walls.wallThickness,intervals,axialLayers);
    mesh::Mesh result;
    for (auto p:core.nodes()) result.addNode(p,duct.classify(p));
    for (Index i=0;i<core.cellCount();++i) {
        const auto& cell=core.cell(i);
        if (cell.type==mesh::CellType::Tri3&&cell.onEntity.tag==coreInterfaceTag) continue;
        result.addCell(cell.type,{core.cellNodes(i),std::size_t(mesh::nodesPerCell(cell.type))},cell.onEntity);
    }
    sweepWallPrisms(duct,core,walls,result);
    return result;
}
} // namespace

int main(int argc,char** argv) {
    if (argc!=8) {
        std::cerr<<"Usage: mesh-hybrid-duct input.step output.paacsmesh first-height-mm growth wall-layers core-intervals axial-layers\n";
        return 2;
    }
    try {
        const auto start=std::chrono::steady_clock::now();
        const auto imported=geo::importOcctCad(argv[1],{});
        require(imported.imported()&&imported.solids==1&&imported.faces==6&&imported.curvedFaces==0,
                "Expected one imported straight duct fixture");
        const Duct duct(*imported.model);
        const auto layerCount=std::stoul(argv[5]), intervals=std::stoul(argv[6]), axialLayers=std::stoul(argv[7]);
        require(layerCount&&layerCount<=12&&intervals>=2&&intervals<=48&&axialLayers>=2&&axialLayers<=48,
                "Duct mesh intervals are out of range");
        algo::LayerScheduleRequest sizing;
        sizing.wallLayers=unsigned(layerCount);
        sizing.firstLayerThickness=std::stod(argv[3]);
        sizing.growthRatio=std::stod(argv[4]);
        sizing.totalDepth=0.;
        const auto walls=algo::planLayerSchedule(sizing);
        require(walls.planned&&walls.wallThickness<halfWidth*.9,"Invalid wall-layer schedule for square section");
        const auto mesh=hybrid(duct,walls,unsigned(intervals),unsigned(axialLayers));
        const auto audit=mesh::auditLinearVolumeBoundary(mesh);
        const auto quality=algo::analyzeMixedMeshQuality(mesh);
        require(audit.valid&&quality.valid&&std::abs(quality.totalVolume-4000.)<1e-7,
                "Hybrid duct volume audit failed");
        const auto pri=mesh.countCells(mesh::CellType::Prism6);
        const auto tet=mesh.countCells(mesh::CellType::Tet4);
        const auto interfaces=8*intervals*axialLayers;
        require(pri==Index(interfaces*layerCount)&&tet==Index(6*intervals*intervals*axialLayers),
                "Unexpected hybrid cell inventory");
        require(sdk::writeMeshArchive(mesh,argv[2]).completed,"Cannot save hybrid duct archive");
        std::cout<<std::setprecision(16)
            <<"{\"route\":\"paacs_mitered_prism_walls_tet_core\",\"nodes\":"<<mesh.nodeCount()
            <<",\"prism6\":"<<pri<<",\"tet4\":"<<tet<<",\"prism_wall_layers\":"<<layerCount
            <<",\"core_intervals_per_axis\":"<<intervals<<",\"axial_layers\":"<<axialLayers
            <<",\"tet_prism_interface_triangles\":"<<interfaces
            <<",\"wall_height_mm\":"<<walls.wallThickness
            <<",\"minimum_mean_ratio\":"<<quality.minimumQuality<<",\"linear_volume_audit_valid\":true"
            <<",\"cad_volume_mm3\":"<<imported.model->volume(duct.body.tag())
            <<",\"mesh_volume_mm3\":"<<quality.totalVolume
            <<",\"volume_tag\":"<<duct.body.tag()<<",\"face_tags\":{";
        for (std::size_t i=0;i<groupNames.size();++i)
            std::cout<<(i?",\"":"\"")<<groupNames[i]<<"\":"<<duct.tags[i];
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"},\"generation_seconds\":"<<elapsed<<"}\n";
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
