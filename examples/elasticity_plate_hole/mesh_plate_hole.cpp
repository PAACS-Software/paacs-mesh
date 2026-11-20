// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include <paacs/algo/ConvexToTetMesher.hpp>
#include <paacs/algo/Delaunay2DScalable.hpp>
#include <paacs/algo/LayeredPatchMesher.hpp>
#include <paacs/algo/PlanarConstraintRecovery.hpp>
#include <paacs/geometry/OcctMesher.hpp>
#include <paacs/mesh/LinearVolumeAudit.hpp>
#include <paacs/sdk/MeshArchive.hpp>
#include <paacs/sizing/RegionSizeField.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

using namespace paacs;

namespace {
constexpr double tol = 1e-8;
constexpr std::size_t maximumSizingPoints = 75000;

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

double radius(Vec2 point) {
    return std::hypot(point.u, point.v);
}

bool close(double a, double b) {
    return std::abs(a-b) < tol;
}

struct Plate {
    const geo::OcctModel& model;
    const geo::GSolid& solid;
    Vec3 lo, hi;
    // bottom, top, left, right, lower y, upper y, circular bore.
    std::array<int,7> tags{};
    const geo::GEdge* rim = nullptr;
    double holeRadius = 0.;
    double first = 0., last = 0.;

    Plate(const geo::OcctModel& imported, const geo::GSolid& body) : model(imported), solid(body) {
        solid.bbox(lo,hi);
        // The fixture is nominally exact. Snap only the values verified below.
        require((lo-Vec3{-50.,-25.,-1.}).norm()<1e-5 &&
                (hi-Vec3{50.,25.,1.}).norm()<1e-5, "Unexpected plate fixture bounds");
        lo = {-50.,-25.,-1.};
        hi = {50.,25.,1.};
        for (auto tag : solid.boundaryFaceTags()) {
            const auto& face = *model.face(tag);
            if (face.kind()==geo::SurfaceKind::Cylinder) {
                require(tags[6]==0, "Expected one cylindrical wall");
                tags[6] = tag;
                continue;
            }
            require(face.kind()==geo::SurfaceKind::Plane, "Expected planar rectangle caps/walls");
            const auto sample = face.pointAt({0.,0.});
            const auto normal = face.normalAt({0.,0.});
            std::size_t index;
            if (std::abs(normal.z)>.99) index = close(sample.z,lo.z) ? 0 : 1;
            else if (std::abs(normal.x)>.99) index = close(sample.x,lo.x) ? 2 : 3;
            else if (std::abs(normal.y)>.99) index = close(sample.y,lo.y) ? 4 : 5;
            else throw std::runtime_error("Unexpected plane normal");
            require(tags[index]==0, "Duplicate CAD face classification");
            tags[index] = tag;
        }
        require(std::ranges::none_of(tags, [](int tag) {return tag<=0;}), "Missing CAD face classification");
        for (const auto edge : model.edges()) {
            if (edge->kind()!=geo::CurveKind::Circle) continue;
            double a,b;
            edge->parameterBounds(a,b);
            const auto sample = edge->pointAt(a);
            if (!close(sample.z,lo.z)) continue;
            require(rim==nullptr, "Expected one closed circular lower rim");
            rim = edge;
            first = a;
            last = b;
            holeRadius = radius({sample.x,sample.y});
        }
        require(rim && close(holeRadius,5.) && close(rim->length(),2*std::numbers::pi*holeRadius),
                "Expected one native circular radius-5 mm bore");
    }

    int side(Vec2 a,Vec2 b) const {
        if (close(a.u,lo.x)&&close(b.u,lo.x)) return tags[2];
        if (close(a.u,hi.x)&&close(b.u,hi.x)) return tags[3];
        if (close(a.v,lo.y)&&close(b.v,lo.y)) return tags[4];
        if (close(a.v,hi.y)&&close(b.v,hi.y)) return tags[5];
        if (close(radius(a),holeRadius)&&close(radius(b),holeRadius)) return tags[6];
        return 0;
    }
};

struct Section {
    std::vector<Vec2> points;
    std::vector<algo::PlanarSegment> constraints;
    std::size_t holeBegin = 0, holeEnd = 0;
};

void closeRing(Section& section, std::size_t start) {
    const auto end = section.points.size();
    for (auto i=start;i<end;++i)
        section.constraints.push_back({Index(i),Index(i+1==end?start:i+1)});
}

Section sectionPoints(const Plate& plate,double far,double near) {
    Section section;
    const auto lo=plate.lo, hi=plate.hi;
    const auto nx = static_cast<int>(std::ceil((hi.x-lo.x)/far));
    const auto ny = static_cast<int>(std::ceil((hi.y-lo.y)/far));
    const auto x = [&](int i) { return lo.x+(hi.x-lo.x)*i/nx; };
    const auto y = [&](int i) { return lo.y+(hi.y-lo.y)*i/ny; };
    for (int i=0;i<nx;++i) section.points.push_back({x(i),lo.y});
    for (int i=0;i<ny;++i) section.points.push_back({hi.x,y(i)});
    for (int i=nx;i>0;--i) section.points.push_back({x(i),hi.y});
    for (int i=ny;i>0;--i) section.points.push_back({lo.x,y(i)});
    closeRing(section,0);
    section.holeBegin = section.points.size();
    const auto rimCount = 4*std::max(3,static_cast<int>(std::ceil(plate.rim->length()/(4*near))));
    for (int i=0;i<rimCount;++i) {
        const auto p = plate.rim->pointAt(plate.first+(plate.last-plate.first)*i/rimCount);
        require(close(radius({p.x,p.y}),plate.holeRadius)&&close(p.z,lo.z),"Unexpected point on native rim");
        section.points.push_back({p.x,p.y});
    }
    closeRing(section,section.holeBegin);
    section.holeEnd=section.points.size();
    // Both runs use the same outer grid. The hole neighborhood uses rings
    // with controlled growth to their shared nominal far-field spacing.
    double ring = plate.holeRadius;
    double spacing = near;
    unsigned ringIndex = 0;
    while (ring + spacing <= 12.) {
        ring += spacing;
        const auto n=std::max(8,static_cast<int>(std::ceil(2*std::numbers::pi*ring/spacing)));
        const double offset = ringIndex++%2==0 ? .5 : 0.;
        for (int i=0;i<n;++i) {
            const auto angle=2*std::numbers::pi*(i+offset)/n;
            section.points.push_back({ring*std::cos(angle),ring*std::sin(angle)});
        }
        spacing = std::min(far,1.32*spacing);
    }
    const double gridStart=std::max(12.+.35*far,ring+.6*far);
    for (int j=1;j<ny;++j) for (int i=1;i<nx;++i) {
        const Vec2 p{x(i),y(j)};
        if (radius(p)>=gridStart) section.points.push_back(p);
    }
    return section;
}

// The plate is a straight sweep. Both the selected cylindrical wall and
// the through-thickness box have the same distance field on each section.
struct PlateSizeField {
    sizing::RegionSizeField regions;
    double radius, boreSize, slope;

    double minimum(Vec2 lo, Vec2 hi) const {
        const double regional = regions.minimumOnBox({lo.u, lo.v, 0}, {hi.u, hi.v, 0});
        if (boreSize == 0) return regional;
        const double inner = std::hypot(std::max({0., lo.u, -hi.u}), std::max({0., lo.v, -hi.v}));
        const double outer = std::hypot(std::max(std::abs(lo.u), std::abs(hi.u)),
                                        std::max(std::abs(lo.v), std::abs(hi.v)));
        const double distance = std::max({0., inner - radius, radius - outer});
        return std::min(regional, boreSize + slope * distance);
    }
};

// Subdivide the original curve in parameter order. A midpoint square of
// half the interval's physical arclength encloses the entire interval.
template<class Evaluate>
void sampleCurve(Section& section, const PlateSizeField& field, Evaluate pointAt,
                 double length, unsigned seedIntervals, double scale) {
    const auto split = [&](auto&& visit, double start, double end) -> void {
        const double middle = (start + end) / 2;
        const auto point = pointAt(middle);
        const double half = (end - start) * length / 2;
        const Vec2 reach{half, half};
        if (2 * half > scale * field.minimum(point - reach, point + reach)) {
            require(start < middle && middle < end, "Plate size is below curve parameter resolution");
            visit(visit, start, middle);
            visit(visit, middle, end);
        } else {
            require(section.points.size() < maximumSizingPoints, "Plate sizing point budget exceeded");
            section.points.push_back(pointAt(start));
        }
    };
    for (unsigned i = 0; i < seedIntervals; ++i)
        split(split, double(i) / seedIntervals, double(i + 1) / seedIntervals);
}

Section sectionPoints(const Plate& plate, const PlateSizeField& field) {
    Section section;
    const Vec2 lo{plate.lo.x, plate.lo.y}, hi{plate.hi.x, plate.hi.y};
    const std::array corners{lo, Vec2{hi.u, lo.v}, hi, Vec2{lo.u, hi.v}};
    for (std::size_t i = 0; i < corners.size(); ++i) {
        const auto a = corners[i], delta = corners[(i + 1) % corners.size()] - a;
        sampleCurve(section, field, [=](double t) { return a + delta * t; },
                    std::hypot(delta.u, delta.v), 1, .9);
    }
    closeRing(section, 0);
    section.holeBegin = section.points.size();
    sampleCurve(section, field, [&](double t) {
        const auto p = plate.rim->pointAt(plate.first + t * (plate.last - plate.first));
        return Vec2{p.x, p.y};
    }, plate.rim->length(), 24, 1.);
    closeRing(section, section.holeBegin);
    section.holeEnd = section.points.size();

    const auto sampleBlock = [&](auto&& visit, Vec2 lower, Vec2 upper) -> void {
        const auto center = (lower + upper) * .5;
        const auto width = upper - lower;
        const Vec2 furthest{std::max(std::abs(lower.u), std::abs(upper.u)),
                            std::max(std::abs(lower.v), std::abs(upper.v))};
        if (radius(furthest) <= plate.holeRadius) return;
        if (std::max(width.u, width.v) > field.minimum(lower, upper)) {
            // Bisect the longer side, keeping a deterministic sampling grid.
            const bool splitX = width.u >= width.v;
            require(splitX ? lower.u < center.u && center.u < upper.u : lower.v < center.v && center.v < upper.v,
                    "Plate size is below section coordinate resolution");
            visit(visit, lower, splitX ? Vec2{center.u, upper.v} : Vec2{upper.u, center.v});
            visit(visit, splitX ? Vec2{center.u, lower.v} : Vec2{lower.u, center.v}, upper);
            return;
        }
        if (radius(center) > plate.holeRadius + .28 * std::min(width.u, width.v)) {
            require(section.points.size() < maximumSizingPoints, "Plate sizing point budget exceeded");
            section.points.push_back(center);
        }
    };
    sampleBlock(sampleBlock, lo, hi);
    return section;
}

PlateSizeField compilePlateField(const Plate& plate, char** args) {
    const int bodyTag = std::stoi(args[0]), faceTag = std::stoi(args[1]);
    const double global = std::stod(args[2]), floor = std::stod(args[3]), growth = std::stod(args[4]);
    const double bore = std::stod(args[5]), box = std::stod(args[6]);
    const Vec3 lo{std::stod(args[7]), std::stod(args[8]), std::stod(args[9])};
    const Vec3 hi{std::stod(args[10]), std::stod(args[11]), std::stod(args[12])};
    require(bodyTag == plate.solid.tag() && (faceTag == 0 || faceTag == plate.tags[6]),
            "Saved body or cylindrical bore face selection does not match the STEP import");
    const auto* wall = plate.model.surfaceDescriptor(plate.tags[6]);
    require(wall && close(wall->origin.x, 0.) && close(wall->origin.y, 0.) &&
            close(std::abs(wall->axis.z), 1.) && close(wall->radius, plate.holeRadius),
            "The plate sizing field requires the original centered, Z-aligned cylindrical bore");
    require(std::isfinite(global) && global <= 10 && std::isfinite(floor) && floor > 0 && floor <= global &&
            std::isfinite(growth) && growth >= 1 && std::isfinite(bore) && bore >= 0 && bore <= global &&
            std::isfinite(box) && box >= 0 && box <= global && (faceTag != 0) == (bore > 0),
            "Invalid plate size, minimum, growth or local target");
    sizing::RegionSizeSpecification specification;
    specification.baseSize = global;
    specification.minimumSize = floor;
    if (box > 0) {
        require(lo.z <= plate.lo.z && hi.z >= plate.hi.z && lo.x < plate.hi.x && hi.x > plate.lo.x &&
                lo.y < plate.hi.y && hi.y > plate.lo.y, "The sizing box must cross the full thickness of this swept plate");
        sizing::RegionSizeControl control;
        control.lower = lo;
        control.upper = hi;
        control.targetSize = box;
        control.growthRatio = growth;
        specification.controls.push_back(control);
    }
    auto compiled = sizing::compileRegionSizeField(specification);
    require(compiled.compiled, "Invalid plate sizing box or field");
    return {std::move(compiled.field), plate.holeRadius, faceTag ? std::max(floor, bore) : 0.,
            2 * (1 - 2 / (growth + 1))};
}

bool insideHole(const Section& section,Vec2 point) {
    bool inside=false;
    auto b = section.points[section.holeEnd-1];
    for (auto i=section.holeBegin;i<section.holeEnd;++i) {
        const auto a=section.points[i];
        if ((a.v>point.v)!=(b.v>point.v) &&
            point.u<(b.u-a.u)*(point.v-a.v)/(b.v-a.v)+a.u) inside=!inside;
        b=a;
    }
    return inside;
}

mesh::Mesh tetrahedralize(const Plate& plate, const Section& section,unsigned layers) {
    auto triangles=algo::delaunay2DScalable(section.points);
    const auto recovered=algo::recoverPlanarConstraints(section.points,section.constraints,triangles);
    require(recovered.completed,"PAACS failed to recover a planar rim constraint");
    mesh::Mesh cap,swept,classified,tets;
    for (const auto p:section.points) cap.addNode({p.u,p.v,plate.lo.z},{2,plate.tags[0]});
    for (const auto tri:triangles) {
        const auto a=section.points[tri.v[0]], b=section.points[tri.v[1]], c=section.points[tri.v[2]];
        if (!insideHole(section,(a+b+c)*(1./3.))) cap.addCell(mesh::CellType::Tri3,tri.v,{2,plate.tags[0]});
    }
    algo::LayeredPatchRequest request;
    const auto depth=plate.hi.z-plate.lo.z;
    request.volumeClassification={3,plate.solid.tag()};
    request.targetClassification={2,plate.tags[1]};
    request.sideClassification={2,0};
    request.schedule.wallLayers=layers;
    request.schedule.firstLayerThickness=depth/layers;
    request.schedule.growthRatio=1.;
    request.schedule.totalDepth=depth;
    request.schedule.maximumCoreThickness=depth/layers;
    request.maximumSourcePairChecks=200'000'000;
    const auto sweep=algo::meshLayeredPatch(cap,request,swept);
    require(sweep.published,"PAACS failed to sweep the section through plate thickness");
    const auto sectionCount=section.points.size();
    require(swept.nodeCount()==Index(sectionCount*(layers+1)), "Unexpected PAACS sweep layout");
    const auto vertices=plate.model.vertices();
    const auto edges=plate.model.edges();
    for (Index i=0;i<swept.nodeCount();++i) {
        const auto point=swept.node(i);
        const Vec2 xy{point.x,point.y};
        const int wall=plate.side(xy,xy);
        mesh::Classification classification{3,plate.solid.tag()};
        if (close(point.z,plate.lo.z)) classification={2,plate.tags[0]};
        else if (close(point.z,plate.hi.z)) classification={2,plate.tags[1]};
        else if (wall) classification={2,wall};
        if (wall) {
            for (const auto vertex:vertices) if ((vertex->position()-point).norm()<tol) {
                classification={0,vertex->tag()};
                break;
            }
            if (classification.dim!=0) for (const auto edge:edges) {
                if ((edge->closestPoint(point)-point).norm()<tol) {
                    classification={1,edge->tag()};
                    break;
                }
            }
        }
        classified.addNode(point,classification);
    }
    for (Index i=0;i<swept.cellCount();++i) {
        const auto cell=swept.cell(i);
        const auto ids=swept.cellNodes(i);
        const auto count=mesh::nodesPerCell(cell.type);
        auto classification=cell.onEntity;
        if (classification==request.sideClassification) {
            const auto a=section.points[ids[0]%sectionCount];
            for (int j=1;j<count;++j) {
                const auto b=section.points[ids[j]%sectionCount];
                if (close(a.u,b.u)&&close(a.v,b.v)) continue;
                classification={2,plate.side(a,b)};
                break;
            }
            require(classification.tag>0,"Could not preserve an original swept wall tag");
        }
        classified.addCell(cell.type,{ids,static_cast<std::size_t>(count)},classification);
    }
    require(algo::convertConvexCellsToTet4(classified,{},tets).published,"PAACS Tet4 conversion failed");
    require(mesh::auditLinearVolumeBoundary(tets).valid,"PAACS final linear volume audit failed");
    return tets;
}
} // namespace

int main(int argc,char** argv) {
    const bool controlled = argc == 18 && std::string_view(argv[3]) == "--controls";
    if (argc != 6 && !controlled) {
        std::cerr<<"Usage: mesh-plate-hole input.step output.paacsmesh far-mm rim-mm layers\n"
                 <<"  or: mesh-plate-hole input.step output.paacsmesh --controls body face global min growth "
                    "bore-size box-size xmin ymin zmin xmax ymax zmax layers\n";
        return 2;
    }
    try {
        const auto start=std::chrono::steady_clock::now();
        const auto meshOptions=geo::NativeCadImportOptions{};
        const auto input=geo::importOcctCad(argv[1],meshOptions);
        require(input.imported()&&input.solids==1&&input.faces==7&&input.curvedFaces==1,
                "Expected one imported circular-hole plate fixture");
        const Plate plate(*input.model,*input.model->solids().front());
        const auto layers=std::stoul(argv[argc-1]);
        require(layers>=1&&layers<=16,"Invalid through-thickness layer count");
        Section section;
        if (controlled) {
            section = sectionPoints(plate, compilePlateField(plate, argv + 4));
        } else {
            const auto far=std::stod(argv[3]), near=std::stod(argv[4]);
            require(std::isfinite(far)&&std::isfinite(near)&&far>=near&&near>=.15&&far<=10.,"Invalid section size");
            section = sectionPoints(plate,far,near);
        }
        const auto tets=tetrahedralize(plate,section,static_cast<unsigned>(layers));
        const auto geometrySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::uint64_t count=0, centroidRejects=0;
        double maxChord=0., volume=0.;
        std::array<unsigned,7> faceCounts{};
        for (Index i=0;i<tets.cellCount();++i) {
            const auto cell=tets.cell(i);
            const auto ids=tets.cellNodes(i);
            if (cell.type==mesh::CellType::Tet4) {
                ++count;
                const auto a=tets.node(ids[0]),b=tets.node(ids[1]),c=tets.node(ids[2]),d=tets.node(ids[3]);
                volume+=(b-a).dot((c-a).cross(d-a))/6.;
                if (!plate.solid.isInside((a+b+c+d)*.25)) ++centroidRejects;
            } else {
                const auto tag=std::ranges::find(plate.tags,cell.onEntity.tag);
                require(tag!=plate.tags.end(),"Unexpected boundary face tag");
                ++faceCounts[tag-plate.tags.begin()];
                if (cell.onEntity.tag!=plate.tags[6]) continue;
                for (int j=0;j<3;++j) {
                    const auto a=tets.node(ids[j]), b=tets.node(ids[(j+1)%3]);
                    require(close(radius({a.x,a.y}),plate.holeRadius),"Hole node differs from native radius");
                    const auto p=(a+b)*.5;
                    maxChord=std::max(maxChord,plate.holeRadius-radius({p.x,p.y}));
                }
            }
        }
        require(count&&centroidRejects==0&&std::ranges::all_of(faceCounts,[](unsigned n){return n>0;}),
                "Original-CAD material or face validation failed");
        const auto cadVolume=input.model->volume(plate.solid.tag());
        require(std::abs(volume-cadVolume)/cadVolume<.001,"Plate CAD volume mismatch");
        require(sdk::writeMeshArchive(tets,argv[2]).completed,"Cannot write Tet4 mesh archive");
        std::cout<<std::setprecision(16)
            <<"{\"route\":\"paacs_constrained_cap_sweep\",\"nodes\":"<<tets.nodeCount()
            <<",\"tet4\":"<<count<<",\"layers\":"<<layers<<",\"bore_tag\":"<<plate.tags[6]
            <<",\"body_tag\":"<<plate.solid.tag()
            <<",\"plate_radius_mm\":"<<plate.holeRadius<<",\"radius_nodes_per_rim\":"
            <<section.holeEnd-section.holeBegin<<",\"maximum_bore_chord_error_mm\":"<<maxChord
            <<",\"outside_centroids\":"<<centroidRejects<<",\"original_cad_volume_mm3\":"<<cadVolume
            <<",\"tet_mesh_volume_mm3\":"<<volume<<",\"generation_seconds\":"<<geometrySeconds
            <<",\"original_face_tags\":[";
        for (std::size_t i=0;i<plate.tags.size();++i) std::cout<<(i?",":"")<<plate.tags[i];
        std::cout<<"]}\n";
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
