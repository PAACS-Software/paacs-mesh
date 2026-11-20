// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

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
constexpr double halfWidth=5., length=40.;
constexpr double tolerance=1e-8;

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool close(double a,double b) {
    return std::abs(a-b)<tolerance;
}

class Duct {
public:
    const geo::OcctModel& model;
    const geo::GSolid& solid;
    // Inlet, outlet and four wall sides in order.
    std::array<int,6> faceTags{};

    explicit Duct(const geo::OcctModel& source) : model(source), solid(*source.solids().front()) {
        Vec3 lo,hi;
        solid.bbox(lo,hi);
        require((lo-Vec3{-halfWidth,-halfWidth,0.}).norm()<1e-5&&
                (hi-Vec3{halfWidth,halfWidth,length}).norm()<1e-5,
                "Unexpected axial box dimensions");
        for (auto tag:solid.boundaryFaceTags()) {
            const auto& face=*model.face(tag);
            require(face.kind()==geo::SurfaceKind::Plane,"Expected all original faces planar");
            const auto point=face.pointAt({0.,0.});
            const auto normal=face.normalAt({0.,0.});
            std::size_t side;
            if (std::abs(normal.z)>.99) {
                require(close(point.z,0.)||close(point.z,length),"Unexpected end plane");
                side=close(point.z,0.) ? 0 : 1;
            } else if (std::abs(normal.x)>.99) {
                require(close(std::abs(point.x),halfWidth),"Unexpected x wall plane");
                side=point.x<0. ? 2 : 3;
            } else {
                require(std::abs(normal.y)>.99&&close(std::abs(point.y),halfWidth),"Unexpected y wall plane");
                side=point.y<0. ? 4 : 5;
            }
            require(!faceTags[side],"Duplicate face class");
            faceTags[side]=tag;
        }
        require(std::ranges::none_of(faceTags,[](int tag){return tag<=0;}),"Missing face class");
    }

    std::array<bool,6> planes(Vec3 p) const {
        return {close(p.z,0.),close(p.z,length),close(p.x,-halfWidth),close(p.x,halfWidth),
                close(p.y,-halfWidth),close(p.y,halfWidth)};
    }

    int faceOf(const mesh::Mesh& mesh,Index cell) const {
        auto common=planes(mesh.node(mesh.cellNodes(cell)[0]));
        for (int j=1;j<mesh::nodesPerCell(mesh.cell(cell).type);++j) {
            const auto on=planes(mesh.node(mesh.cellNodes(cell)[j]));
            for (std::size_t side=0;side<common.size();++side) common[side]=common[side]&&on[side];
        }
        require(std::ranges::count(common,true)==1,"Boundary facet does not have one original CAD owner");
        return faceTags[std::ranges::find(common,true)-common.begin()];
    }

    mesh::Mesh classified(const mesh::Mesh& input) const {
        mesh::Mesh result;
        const auto edges=model.edges();
        const auto vertices=model.vertices();
        for (auto point:input.nodes()) {
            const auto on=planes(point);
            const auto count=std::ranges::count(on,true);
            mesh::Classification classification{3,solid.tag()};
            if (count==1) {
                classification={2,faceTags[std::ranges::find(on,true)-on.begin()]};
            } else if (count==2) {
                for (const auto* edge:edges) if ((edge->closestPoint(point)-point).norm()<tolerance) {
                    classification={1,edge->tag()};
                    break;
                }
                require(classification.dim==1,"Missing original CAD edge classification");
            } else if (count==3) {
                for (const auto* vertex:vertices) if ((vertex->position()-point).norm()<tolerance) {
                    classification={0,vertex->tag()};
                    break;
                }
                require(classification.dim==0,"Missing original CAD vertex classification");
            }
            result.addNode(point,classification);
        }
        for (Index i=0;i<input.cellCount();++i) {
            const auto cell=input.cell(i);
            const auto classification=cell.type==mesh::CellType::Prism6 ? cell.onEntity :
                mesh::Classification{2,faceOf(input,i)};
            result.addCell(cell.type,{input.cellNodes(i),std::size_t(mesh::nodesPerCell(cell.type))},classification);
        }
        return result;
    }
};

std::vector<double> mirror(const algo::LayerSchedule& half) {
    std::vector<double> stations;
    stations.reserve(2*half.offsets.size()-1);
    for (auto offset:half.offsets) stations.push_back(offset-halfWidth);
    for (auto i=half.offsets.size()-1;i>0;--i) stations.push_back(halfWidth-half.offsets[i-1]);
    return stations;
}

mesh::Mesh sweep(const Duct& duct,const std::vector<double>& stations,unsigned layers) {
    mesh::Mesh section,swept;
    const auto count=Index(stations.size());
    const auto id=[&](Index i,Index j) {return count*j+i;};
    for (const auto y:stations) for (const auto x:stations)
        section.addNode({x,y,0.},{2,duct.faceTags[0]});
    for (Index j=0;j+1<count;++j) for (Index i=0;i+1<count;++i) {
        const auto a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);
        // Alternate the diagonals across the square section.
        const auto triangle=[&](Index u,Index v,Index w) {
            section.addCell(mesh::CellType::Tri3,{u,v,w},{2,duct.faceTags[0]});
        };
        if ((i+j)%2==0) {
            triangle(a,b,c);
            triangle(a,c,d);
        } else {
            triangle(a,b,d);
            triangle(b,c,d);
        }
    }
    algo::LayeredPatchRequest request;
    request.schedule.wallLayers=layers;
    request.schedule.firstLayerThickness=length/layers;
    request.schedule.growthRatio=1.;
    request.schedule.totalDepth=length;
    request.schedule.maximumCoreThickness=length/layers;
    request.volumeClassification={3,duct.solid.tag()};
    request.targetClassification={2,duct.faceTags[1]};
    request.sideClassification={2,0};
    const auto outcome=algo::meshLayeredPatch(section,request,swept);
    if (!outcome.published) {
        std::cerr<<"PAACS sweep reason "<<static_cast<int>(outcome.reason)<<'\n';
        throw std::runtime_error("PAACS straight prism sweep failed");
    }
    return duct.classified(swept);
}
} // namespace

int main(int argc,char** argv) {
    if (argc!=8) {
        std::cerr<<"Usage: mesh-straight-duct input.step output.paacsmesh first-height-mm growth wall-rows core-cap-mm axial-layers\n";
        return 2;
    }
    try {
        const auto started=std::chrono::steady_clock::now();
        const auto imported=geo::importOcctCad(argv[1],{});
        require(imported.imported()&&imported.solids==1&&imported.faces==6&&imported.curvedFaces==0,
                "Expected one imported straight square duct");
        const Duct duct(*imported.model);
        const auto wallRows=std::stoul(argv[5]), layers=std::stoul(argv[7]);
        require(wallRows>=1&&wallRows<=16&&layers>=1&&layers<=64,"Invalid duct layer count");
        algo::LayerScheduleRequest across;
        across.wallLayers=unsigned(wallRows);
        across.firstLayerThickness=std::stod(argv[3]);
        across.growthRatio=std::stod(argv[4]);
        across.maximumAllowedGrowth=1.3;
        across.totalDepth=halfWidth;
        across.maximumCoreThickness=std::stod(argv[6]);
        across.maximumLayers=48;
        const auto half=algo::planLayerSchedule(across);
        if (!half.planned) {
            std::cerr<<"PAACS layer planner reason "<<static_cast<int>(half.reason)<<'\n';
            throw std::runtime_error("PAACS wall layer planner refused half-section");
        }
        const auto stations=mirror(half);
        const auto mesh=sweep(duct,stations,unsigned(layers));
        const auto audit=mesh::auditLinearVolumeBoundary(mesh);
        const auto quality=algo::analyzeMixedMeshQuality(mesh);
        require(audit.valid&&quality.valid&&close(quality.totalVolume,4000.),"Invalid duct volume or orientation");
        std::array<Index,6> groups{};
        Index prisms=0;
        for (Index i=0;i<mesh.cellCount();++i) {
            const auto cell=mesh.cell(i);
            if (cell.type==mesh::CellType::Prism6) {
                ++prisms;
            } else {
                const auto which=std::ranges::find(duct.faceTags,cell.onEntity.tag);
                require(which!=duct.faceTags.end(),"Unknown CAD face classification");
                const auto side=which-duct.faceTags.begin();
                require((side<2&&cell.type==mesh::CellType::Tri3)||(side>=2&&cell.type==mesh::CellType::Quad4),
                        "Unexpected duct boundary facet type");
                ++groups[side];
            }
        }
        const auto n=Index(stations.size()-1);
        require(prisms==2*n*n*Index(layers)&&groups[0]==2*n*n&&groups[1]==2*n*n,
                "Incorrect prism or end face count");
        for (int side=2;side<6;++side)
            require(groups[side]==n*Index(layers),"Missing wetted wall facets");
        require(sdk::writeMeshArchive(mesh,argv[2]).completed,"Could not write duct archive");
        constexpr std::array groupNames{"inlet","outlet","x_minus","x_plus","y_minus","y_plus"};
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::cout<<std::setprecision(16)
            <<"{\"route\":\"paacs_streamwise_graded_prism_sweep\",\"nodes\":"<<mesh.nodeCount()
            <<",\"prism6\":"<<prisms<<",\"wall_rows\":"<<wallRows<<",\"axial_layers\":"<<layers
            <<",\"core_rows_per_half\":"<<half.coreLayers<<",\"cross_section_intervals_per_axis\":"<<n
            <<",\"dedicated_wall_depth_mm\":"<<half.wallThickness<<",\"cross_section_core_maximum_mm\":"
            <<half.maximumThickness<<",\"cell_volume_total_mm3\":"<<quality.totalVolume
            <<",\"original_cad_volume_mm3\":"<<imported.model->volume(duct.solid.tag())
            <<",\"minimum_mean_ratio\":"<<quality.minimumQuality<<",\"linear_volume_audit_valid\":true"
            <<",\"volume_tag\":"<<duct.solid.tag()<<",\"face_tags\":{";
        for (std::size_t i=0;i<groupNames.size();++i)
            std::cout<<(i?",\"":"\"")<<groupNames[i]<<"\":"<<duct.faceTags[i];
        std::cout<<"},\"face_counts\":{";
        for (std::size_t i=0;i<groupNames.size();++i)
            std::cout<<(i?",\"":"\"")<<groupNames[i]<<"\":"<<groups[i];
        std::cout<<"},\"generation_seconds\":"<<elapsed<<"}\n";
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
