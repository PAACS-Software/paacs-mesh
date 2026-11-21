// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "OcctDelaunayCavity.hpp"
#include "OcctDelaunayQuality.hpp"
#include "OcctSpatialSizing.hpp"
#include "paacs/mesh/Mesh.hpp"
#include "paacs/mesh/Quality.hpp"
#include <queue>
#include <span>

namespace paacs::geo::detail {

struct LocalDelaunayLimits {
    std::uint64_t work, queries, nodes, tets, workingBytes;
    double size;
    const OcctSpatialSizing* field = nullptr;
    bool improveSpatialQuality = false;
};

// Only used on an embedded, oriented Tet4 fill of a single original solid.
// A cavity keeps its old oriented shell and never crosses a protected face.
// Stage the vertex and all replacement cells before changing the live complex.
template<class AdmitNode, class AdmitCentroid, class Charge>
std::size_t insertOcctDelaunaySites(mesh::Mesh& mesh, std::span<const Vec3> supports,
    const LocalDelaunayLimits& limits, AdmitNode&& admitNode, AdmitCentroid&& admitCentroid, Charge&& charge,
    std::span<const unsigned char> material = {}, unsigned* qualityEdits = nullptr) {
    using Complex = bnd::detail::BistellarComplex;
    using Cell = bnd::TetFacetCell;
    using Vertex = bnd::TetFacetVertexId;
    using Sign = predicates::Sign;
    using bnd::tet_facet_detail::outwardFace;
    constexpr std::size_t cavityLimit = 64, maximumRefinementNodes = 1024;
    const bool sizing = supports.empty();
    const double minimumConeShape = sizing && limits.field && limits.improveSpatialQuality ? 1e-6 : 0.;
    const bool optimize = sizing && qualityEdits;
    if (qualityEdits) *qualityEdits=0;
    const auto maximumInsertions = sizing ? maximumRefinementNodes : supports.size();
    auto tetCount = std::size_t(mesh.countCells(mesh::CellType::Tet4));
    const auto originalNodes = std::size_t(mesh.nodeCount());
    if (!tetCount || (!material.empty() && material.size()!=tetCount) ||
        tetCount > limits.tets || originalNodes >= limits.nodes ||
        !(limits.size > 0.) || !limits.queries) return 0;
    // Append-only cell IDs keep stale heap entries harmless and avoid updating
    // unaffected neighbors. Allow for the sparse face/edge and AABB indices.
    const auto maximumStored = std::min<std::uint64_t>(tetCount + maximumInsertions*(2*cavityLimit+2),
        limits.workingBytes/4096);
    if (maximumStored <= tetCount) return 0;

    struct TrialLimit {};
    Complex live;
    live.walkConstraints = true;
    live.maximumStoredCells = std::size_t(maximumStored);
    std::uint64_t work = 0, queries = 0;
    auto workLimit=optimize ? limits.work-limits.work/4 : limits.work;
    auto queryLimit=optimize ? limits.queries-limits.queries/4 : limits.queries;
    std::vector<unsigned char> admitted(material.begin(),material.end());
    live.chargeWork = [&](std::uint64_t amount) {
        if (amount > workLimit-work) throw TrialLimit{};
        charge(amount);
        work += amount;
    };
    const auto point = [&](Vertex v) { return live.vertices[v].position; };
    const auto locateSupport = [&](Vec3 p) {
        auto found=Complex::noCell;
        live.visit({p,p},[&](std::size_t owner) {
            for (unsigned side=0;side<4;++side) {
                live.charge();
                const auto f=outwardFace(live.cells[owner],side);
                if (predicates::orient3d(point(f[0]),point(f[1]),point(f[2]),p)==Sign::Positive)
                    return true;
            }
            found=owner;
            return false;
        });
        return found;
    };
    const auto shape = [&](const Cell& cell) {
        const auto& v=cell.vertices;
        return mesh::tetScaledJacobian(point(v[0]),point(v[1]),point(v[2]),point(v[3]));
    };
    const auto checkNode = [&](Vec3 p) {
        if (queries>=queryLimit) throw TrialLimit{};
        ++queries;
        return admitNode(p);
    };
    const auto checkCentroid = [&](Vec3 p) {
        if (queries>=queryLimit) throw TrialLimit{};
        ++queries;
        return admitCentroid(p);
    };
    const auto commit = [&](Complex::Proposal& patch, Vec3 p, double localSize) {
        if (tetCount-patch.removed.size()+patch.added.size()>limits.tets ||
            patch.added.size()>live.maximumStoredCells-live.cells.size()) return false;
        if (!admitted.empty() && std::ranges::all_of(patch.removed,[&](std::size_t id) {
            return admitted[id]!=0;
        })) return false;
        const auto staged=Vertex(live.vertices.size());
        live.vertices.push_back({p,std::uint64_t(staged)+1});
        const auto acceptable = [&] {
            if (!live.valid(patch,std::nullopt,staged)) return false;
            if (minimumConeShape>0.) {
                for (const auto& cell:patch.added) {
                    live.charge(4);
                    if (!(shape(cell)>=minimumConeShape)) return false;
                }
            }
            if (!admitted.empty()) {
                double minimum=1.;
                for (const auto id:patch.removed) minimum=std::min(minimum,shape(live.cells[id]));
                const auto floor=std::min(minimum,.08);
                for (const auto& cell:patch.added)
                    if (!(shape(cell)>=floor)) return false;
            }
            if (sizing) {
                for (const auto& cell:patch.added) for (const auto v:cell.vertices) {
                    if (v!=staged && (point(v)-p).norm()<(limits.field?.50:.05)*localSize) return false;
                }
            }
            if (!live.preservesSurvivors(patch)) return false;
            if (!checkNode(p)) return false;
            for (const auto& cell:patch.added) {
                Vec3 center{};
                for (const auto v:cell.vertices) center=center+point(v)*.25;
                if (!checkCentroid(center)) return false;
            }
            live.apply(patch);
            if (!admitted.empty()) admitted.resize(live.cells.size(),1);
            return true;
        };
        try {
            if (acceptable()) {
                tetCount=tetCount-patch.removed.size()+patch.added.size();
                return true;
            }
        } catch (const TrialLimit&) {
            live.vertices.pop_back();
            throw;
        }
        live.vertices.pop_back();
        return false;
    };
    // Carry the measured mesh size and the background size independently.
    // The volume background for the native global-size route is constant.
    struct Size { double mesh=0., background=0.; };
    std::vector<Size> sizes(originalNodes);
    std::vector<unsigned char> fixed(originalNodes);
    struct Candidate {
        double ratio;
        std::size_t tet;
        Vec3 center, centroid;
        bool operator<(const Candidate& other) const {
            return ratio==other.ratio ? tet>other.tet : ratio<other.ratio;
        }
    };
    std::priority_queue<Candidate> queue;
    const auto fieldSize = [&](Vec3 p) {
        if (!limits.field) return limits.size;
        live.charge(limits.field->queryWork());
        return limits.field->eval(p);
    };
    const auto enqueue = [&](std::size_t id) {
        live.charge();
        const auto& v=live.cells[id].vertices;
        const auto a=point(v[1])-point(v[0]), b=point(v[2])-point(v[0]), c=point(v[3])-point(v[0]);
        const auto determinant=2.*predicates::orient3dDeterminantEstimate(
            point(v[0]),point(v[1]),point(v[2]),point(v[3]));
        if (!(determinant>0.) || !std::isfinite(determinant)) return;
        const auto offset=(b.cross(c)*a.norm2()+c.cross(a)*b.norm2()+a.cross(b)*c.norm2())/determinant;
        const auto radius=offset.norm();
        const auto circumcenter=point(v[0])+offset;
        if (limits.field) {
            // A sliver may have a very large circumsphere with short edges.
            // Rank actual oversize cells and use the circumcenter only when
            // it is local enough to influence the cell that was queued.
            const auto centroid=point(v[0])+(a+b+c)*.25;
            std::array lengths{a.norm(),b.norm(),c.norm(),(a-b).norm(),(a-c).norm(),(b-c).norm()};
            std::ranges::sort(lengths);
            const auto size=fieldSize(centroid);
            const auto ratio=(lengths[2]+lengths[3])/(2*size);
            if (ratio>1.10) {
                const bool local=std::isfinite(radius)&&(circumcenter-centroid).norm()<.60*lengths.back();
                queue.push({ratio,id,local?circumcenter:centroid,centroid});
            }
            return;
        }
        double meshSize=0.,background=0.;
        for (auto vertex:v) {
            meshSize+=sizes[vertex].mesh*.25;
            background+=sizes[vertex].background*.25;
        }
        const auto size=std::min(meshSize,background);
        if (size>0. && std::isfinite(radius) && radius>=size)
            queue.push({radius/size,id,circumcenter,{}});
    };
    const auto interpolate = [&](std::size_t owner,Vec3 p) {
        const auto& v=live.cells[owner].vertices;
        const auto a=point(v[0]),b=point(v[1]),c=point(v[2]),d=point(v[3]);
        const std::array weights{
            predicates::orient3dDeterminantEstimate(p,b,c,d),
            predicates::orient3dDeterminantEstimate(a,p,c,d),
            predicates::orient3dDeterminantEstimate(a,b,p,d),
            predicates::orient3dDeterminantEstimate(a,b,c,p)};
        double total=0.,weighted=0.;
        for (unsigned i=0;i<4;++i) {
            const auto weight=std::max(0.,weights[i]);
            total+=weight;
            weighted+=weight*sizes[v[i]].mesh;
        }
        return total>0. && std::isfinite(weighted) ? weighted/total : limits.size;
    };
    std::size_t inserted=0;
    bool initialized=false;
    try {
        live.vertices.reserve(std::min<std::uint64_t>(limits.nodes,originalNodes+maximumInsertions));
        for (const auto p:mesh.nodes()) {
            live.charge();
            live.vertices.push_back({p,live.vertices.size()+1});
        }
        for (Index i=0;i<mesh.cellCount();++i) {
            live.charge();
            const auto& item=mesh.cell(i);
            const auto* v=mesh.cellNodes(i);
            if (item.type==mesh::CellType::Tet4) {
                live.insert({{Vertex(v[0]),Vertex(v[1]),Vertex(v[2]),Vertex(v[3])},0,item.onEntity});
                continue;
            }
            if (item.type!=mesh::CellType::Tri3) return 0;
            const Complex::Face face{Vertex(v[0]),Vertex(v[1]),Vertex(v[2])};
            live.protectedFaces.insert(Complex::key(face));
            for (unsigned j=0;j<3;++j) {
                const auto a=face[j], b=face[(j+1)%3];
                live.protectedEdges.insert(Complex::key(Complex::Edge{a,b}));
                if (!sizing) continue;
                fixed[a]=fixed[b]=1;
                if (!limits.field) {
                    const auto length=(point(a)-point(b)).norm();
                    sizes[a].mesh=std::max(sizes[a].mesh,length);
                    sizes[b].mesh=std::max(sizes[b].mesh,length);
                }
            }
        }
        initialized=true;
        if (sizing) {
            if (!limits.field) {
                for (const auto& cell:live.cells) for (const auto& e:Complex::edgeKeys(cell)) {
                    live.charge();
                    if (fixed[e[0]] && fixed[e[1]]) continue;
                    const auto length=(point(e[0])-point(e[1])).norm();
                    for (const auto v:e) if (!fixed[v]) sizes[v].mesh=std::max(sizes[v].mesh,length);
                }
                for (std::size_t i=0;i<originalNodes;++i) {
                    live.charge();
                    if (!(sizes[i].mesh>0.)) sizes[i].mesh=limits.size;
                    sizes[i].background=sizes[i].mesh;
                }
            }
            const auto insertAt = [&](std::size_t start, Vec3 p, bool localPatch, bool updateQueue) {
                std::size_t owner=Complex::noCell;
                auto patch=buildOcctDelaunayCavity(live,start,p,cavityLimit,owner,!localPatch,
                                                   minimumConeShape);
                if (!patch) return false;
                const auto background=fieldSize(p);
                const auto size=limits.field?background:interpolate(owner,p);
                const auto first=live.cells.size();
                if (!commit(*patch,p,limits.field?background:std::min(size,background))) return false;
                ++inserted;
                if (!limits.field) sizes.push_back({size,background});
                if (updateQueue)
                    for (auto id=first;id<live.cells.size();++id) enqueue(id);
                return true;
            };
            if (limits.field) {
                // An enclosed control can be smaller than the initial cells.
                // Sample inside the box before the circumcenter queue so the
                // region does not depend on intersecting an original vertex.
                const auto& field=limits.field->boxes();
                const auto boxes=field.controls();
                const auto siteShare=std::max<std::size_t>(1,maximumInsertions/2/
                    std::max<std::size_t>(1,boxes.size()));
                for (const auto& box:boxes) {
                    live.charge();
                    const auto h=std::max(field.minimumSize(),box.targetSize);
                    if (h>=field.baseSize()) continue;
                    const auto extent=box.upper-box.lower;
                    const std::array length{extent.x,extent.y,extent.z};
                    std::array<unsigned,3> count{};
                    for (unsigned i=0;i<3;++i) count[i]=static_cast<unsigned>(std::clamp(
                        std::ceil(length[i]/h),1.,double(siteShare)));
                    while (std::size_t(count[0])*count[1]*count[2]>siteShare) {
                        auto largest=std::max_element(count.begin(),count.end());
                        *largest=(*largest+1)/2;
                    }
                    for (unsigned i=0;i<count[0];++i)
                        for (unsigned j=0;j<count[1];++j)
                            for (unsigned k=0;k<count[2];++k) {
                                live.charge();
                                if (inserted>=maximumInsertions/2||live.vertices.size()>=limits.nodes) break;
                                const Vec3 site{
                                    std::lerp(box.lower.x,box.upper.x,(i+.5)/count[0]),
                                    std::lerp(box.lower.y,box.upper.y,(j+.5)/count[1]),
                                    std::lerp(box.lower.z,box.upper.z,(k+.5)/count[2])};
                                const auto owner=locateSupport(site);
                                if (owner==Complex::noCell) continue;
                                const auto localSize=fieldSize(site);
                                const auto& vertices=live.cells[owner].vertices;
                                if (std::ranges::any_of(vertices,[&](auto v) {
                                    return (point(v)-site).norm()<.42*localSize;
                                })) continue;
                                if (!insertAt(owner,site,false,false)) insertAt(owner,site,true,false);
                            }
                }
            }
            for (std::size_t i=0;i<live.cells.size();++i)
                if (live.alive[i]) enqueue(i);
            while (!queue.empty() && inserted<maximumInsertions && live.vertices.size()<limits.nodes) {
                live.charge();
                const auto candidate=queue.top();
                queue.pop();
                if (!live.alive[candidate.tet]) continue;
                if (!insertAt(candidate.tet,candidate.center,false,true) && limits.field &&
                    (candidate.center-candidate.centroid).norm2()>0.)
                    insertAt(candidate.tet,candidate.centroid,false,true);
            }
        } else {
            for (const auto p:supports) {
                live.charge();
                if (live.vertices.size()>=limits.nodes) break;
                const auto owner=locateSupport(p);
                if (owner==Complex::noCell) continue;
                // Thin CAD wedges can block an entire in-sphere cavity. Try
                // the smaller visible patch at the same verified support site.
                std::size_t containing=Complex::noCell;
                auto patch=buildOcctDelaunayCavity(live,owner,p,32,containing);
                if (patch && commit(*patch,p,limits.size)) {
                    ++inserted;
                    continue;
                }
                patch=buildOcctDelaunayCavity(live,owner,p,32,containing,false);
                if (patch && commit(*patch,p,limits.size)) ++inserted;
            }
        }
    } catch (const TrialLimit&) {
        // Work stops between transactions. Every committed edit is usable.
    }
    if (initialized&&optimize) {
        workLimit=limits.work;
        queryLimit=limits.queries;
        try {
            optimizeOcctDelaunayQuality(live,fixed,checkCentroid,checkNode,*qualityEdits,limits.tets);
        } catch (const TrialLimit&) {
            // Keep the edits admitted before the limit.
        }
    }
    if (!inserted && !(qualityEdits && *qualityEdits)) return 0;
    mesh::Mesh result;
    const auto region=live.cells.front().classification;
    for (Index i=0;i<mesh.nodeCount();++i) result.addNode(point(Vertex(i)),mesh.nodeClass(i));
    for (std::size_t i=originalNodes;i<live.vertices.size();++i) result.addNode(point(Vertex(i)),region);
    for (Index i=0;i<mesh.cellCount();++i) if (mesh.cell(i).type==mesh::CellType::Tri3)
        result.addCell(mesh::CellType::Tri3,std::span(mesh.cellNodes(i),3),mesh.cell(i).onEntity);
    for (std::size_t i=0;i<live.cells.size();++i) if (live.alive[i]) {
        const auto& tet=live.cells[i];
        const auto& v=tet.vertices;
        result.addCell(mesh::CellType::Tet4,{Index(v[0]),Index(v[1]),Index(v[2]),Index(v[3])},tet.classification);
    }
    charge(0);
    mesh.swap(result);
    return inserted;
}

} // namespace paacs::geo::detail
