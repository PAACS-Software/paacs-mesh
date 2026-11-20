// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "OcctWireRefinement.hpp"
#include "OcctNativeWireAuthority.hpp"
#include <Standard_OutOfMemory.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <new>
#include <set>

namespace paacs::geo::detail {
namespace {
using Edge = std::array<Index,2>;
struct Stop { const char* diagnostic; OcctWireRefinementReason reason; };
void need(bool ok, const char* why,
          OcctWireRefinementReason reason=OcctWireRefinementReason::UnsupportedGeometry) {
    if(!ok) throw Stop{why,reason};
}
Edge key(Index a,Index b) { if(b<a)std::swap(a,b);return {a,b}; }
bool finite(Vec3 p) { return std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z); }
struct Segment { std::set<int> faces; unsigned uses=0; std::vector<Index> chain;std::array<double,2> parameterKey{}; };
}
OcctWireRefinementOutcome refineOcctWireRefinement(const OcctModel& model,
    std::span<const OcctWireBoundary> inputs,const OcctWireRefinementOptions& options) noexcept {
    OcctWireRefinementOutcome out;
    const auto tick=[&](std::uint64_t n=1){
        if(options.limits.cancellation&&options.limits.cancellation->load())
            throw Stop{"cancelled",OcctWireRefinementReason::Cancelled};
        need(n<=options.limits.maximumWork-out.workUsed,"work-budget",OcctWireRefinementReason::BudgetExceeded);
        out.workUsed+=n;
    };
    const auto certificate=[&](const mesh::Mesh& surface){
        auto limits=options.limits;limits.maximumWork-=out.workUsed;
        const auto c=certifyOcctBoundarySurface(surface,limits);tick(c.workUsed);
        if(!c.certified){
            const auto reason=c.reason==OcctStarReason::WorkBudget||c.reason==OcctStarReason::OutputBudget?
                OcctWireRefinementReason::BudgetExceeded:c.reason==OcctStarReason::Cancelled?
                OcctWireRefinementReason::Cancelled:c.reason==OcctStarReason::AllocationFailed?
                OcctWireRefinementReason::AllocationFailed:OcctWireRefinementReason::InvalidBoundary;
            throw Stop{c.diagnostic,reason};
        }
    };
    try {
        tick(0);
        need(!inputs.empty()&&inputs.size()<=64&&options.edgeTag>0&&
             options.segmentsPerExistingSegment>=2&&options.segmentsPerExistingSegment<=32,
             "invalid-refinement-options",OcctWireRefinementReason::InvalidRequest);
        const auto* wire=model.edge(options.edgeTag);
        need(wire!=nullptr,"missing-native-wire",OcctWireRefinementReason::InvalidRequest);
        // A trimmed arc can have a periodic underlying curve without crossing
        // its seam. Distinct endpoints and the finite native parameter interval
        // give an unambiguous chain; every pcurve use is still checked below.
        const bool closedWire=wire->startVertexTag()>=0&&wire->startVertexTag()==wire->endVertexTag();
        need(!closedWire||options.allowClosedWireCandidates,"closed-periodic-wire-requires-coordinated-remeshing");
        need(wire->kind()==CurveKind::Line||options.allowNonlinearWireCandidates,"nonlinear-wire-candidate-not-enabled");
        double lower=0,upper=0;wire->parameterBounds(lower,upper);
        need(std::isfinite(lower)&&std::isfinite(upper)&&upper>lower,"invalid-native-parameter-range");
        const auto start=wire->pointAt(lower),end=wire->pointAt(upper),delta=end-start;
        const auto squared=delta.norm2();
        need(finite(start)&&finite(end)&&std::isfinite(squared)&&(closedWire||squared>0),"invalid-linear-wire");
        const auto rawTolerance=model.tolerance(1,options.edgeTag);
        need(std::isfinite(rawTolerance)&&rawTolerance>0,"invalid-source-wire-tolerance");
        std::set<int> bodies,affectedFaces;
        std::set<const mesh::Mesh*> destinations;
        std::uint64_t estimated=0;
        for(const auto& item:inputs){tick();
            need(item.boundary&&model.solid(item.bodyTag)&&bodies.insert(item.bodyTag).second&&
                 destinations.insert(item.boundary).second,"invalid-or-duplicate-body-boundary",OcctWireRefinementReason::InvalidRequest);
            const auto& surface=*item.boundary;
            need(std::uint64_t(surface.nodeCount())<=options.limits.maximumBoundaryNodes&&
                 std::uint64_t(surface.cellCount())<=options.limits.maximumBoundaryTriangles,
                 "input-count-budget",OcctWireRefinementReason::BudgetExceeded);
            // Conservative payload admission; native query heap is excluded.
            const auto rows=std::uint64_t(surface.nodeCount())+std::uint64_t(surface.cellCount());
            need(rows<=(UINT64_MAX-estimated)/512,"working-overflow",OcctWireRefinementReason::BudgetExceeded);
            estimated+=rows*512;
            for(const int face:model.solid(item.bodyTag)->boundaryFaceTags())
                for(const auto& ring:model.boundaryUses(face))for(const auto& use:ring.edges){tick();
                    if(use.edgeTag==options.edgeTag&&!use.degenerated)affectedFaces.insert(face);}
        }
        need(options.maximumAddedNodes<=(UINT64_MAX-estimated)/512,"working-overflow",OcctWireRefinementReason::BudgetExceeded);
        estimated+=options.maximumAddedNodes*512;
        need(estimated<=options.limits.maximumWorkingBytes,"working-payload-budget",OcctWireRefinementReason::BudgetExceeded);
        need(!affectedFaces.empty(),"wire-not-owned-by-supplied-bodies",OcctWireRefinementReason::InvalidRequest);
        for(const int face:affectedFaces)for(const auto& owner:model.faceOwners(face)){tick();
            need(bodies.contains(owner.solidTag),"missing-shared-face-coowner");}
        std::vector<mesh::Mesh> rebuilt;rebuilt.reserve(inputs.size());
        // Parameter chains must match across coowners before any publication.
        using ParameterSegment=std::array<double,2>;
        std::map<int,std::map<int,std::set<ParameterSegment>>> sharedChains;
        std::uint64_t newNodes=0,newTriangles=0;
        for(const auto& item:inputs){
            const auto& surface=*item.boundary;certificate(surface);
            std::set<int> bodyFaces;
            for(const int face:model.solid(item.bodyTag)->boundaryFaceTags()){tick();bodyFaces.insert(face);}
            std::map<Edge,Segment> segments;
            for(Index ci=0;ci<surface.cellCount();++ci){tick();const auto& cell=surface.cell(ci);
                need(cell.type==mesh::CellType::Tri3&&cell.onEntity.dim==2&&bodyFaces.contains(cell.onEntity.tag),"invalid-native-face-ownership",OcctWireRefinementReason::InvalidBoundary);
                const auto* ids=surface.cellNodes(ci);
                for(unsigned j=0;j<3;++j){tick();auto& s=segments[key(ids[j],ids[(j+1)%3])];++s.uses;s.faces.insert(cell.onEntity.tag);}
            }
            mesh::Mesh candidate;
            for(Index i=0;i<surface.nodeCount();++i){tick();candidate.addNode(surface.node(i),surface.nodeClass(i));}
            std::set<Edge> selectedEdges;
            for(const auto& [edge,segment]:segments){tick();if(actualOcctWire(model,surface.nodeClass(edge[0]),surface.nodeClass(edge[1]),segment.faces,[&]{tick();})==options.edgeTag)selectedEdges.insert(edge);}
            std::map<Edge,std::array<double,2>> closedParameters;
            if(closedWire){
                // Order the whole loop by native parameter, rather than
                // guessing the seam from angular/parameter distance. This
                // also supports nonuniform closed spline parameterizations.
                std::map<Index,unsigned> degree;for(auto edge:selectedEdges)for(auto id:edge)++degree[id];
                need(degree.size()>=3&&selectedEdges.size()==degree.size(),"closed-wire-cycle-incomplete");
                Index anchor=-1;std::vector<std::pair<double,Index>> ordered;
                for(auto [id,count]:degree){tick();need(count==2,"closed-wire-node-degree");const auto cls=surface.nodeClass(id);
                    if(cls.dim==0&&cls.tag==wire->startVertexTag()){need(anchor<0,"duplicate-closed-wire-anchor");anchor=id;continue;}
                    const auto t=wire->closestParameter(surface.node(id));need(std::isfinite(t)&&t>lower&&t<upper,"closed-wire-parameter-ambiguous");ordered.push_back({t,id});}
                need(anchor>=0,"closed-wire-anchor-missing");std::sort(ordered.begin(),ordered.end());ordered.insert(ordered.begin(),{lower,anchor});
                for(std::size_t i=0;i<ordered.size();++i){tick();const auto a=ordered[i],b=i+1<ordered.size()?ordered[i+1]:std::pair<double,Index>{upper,anchor};
                    need(b.first>a.first,"closed-wire-parameter-order");const auto edge=key(a.second,b.second);need(selectedEdges.contains(edge),"closed-wire-cycle-disagrees-with-native-parameters");
                    closedParameters.emplace(edge,edge[0]==a.second?std::array{a.first,b.first}:std::array{b.first,a.first});}
            }
            std::size_t selected=0;
            for(auto& [edge,segment]:segments){tick();
                if(!selectedEdges.contains(edge))continue;
                ++selected;need(segment.uses==2,"selected-wire-is-not-closed",OcctWireRefinementReason::InvalidBoundary);
                std::array<double,2> parameters;
                for(unsigned i=0;i<2;++i){tick();const auto p=surface.node(edge[i]);
                    const auto fraction=squared>0?(p-start).dot(delta)/squared:0.;
                    const auto cls=surface.nodeClass(edge[i]);
                    parameters[i]=closedWire?closedParameters.at(edge)[i]:cls.dim==0&&cls.tag==wire->startVertexTag()?lower:
                        cls.dim==0&&cls.tag==wire->endVertexTag()?upper:wire->closestParameter(p);
                    const auto c=surface.nodeClass(edge[i]);
                    const auto tolerance=c.dim==0?std::max(rawTolerance,model.tolerance(0,c.tag)):rawTolerance;
                    need(finite(p)&&std::isfinite(fraction)&&std::isfinite(parameters[i])&&parameters[i]>=lower&&parameters[i]<=upper&&
                         (p-wire->pointAt(parameters[i])).norm()<=tolerance,"original-node-outside-native-wire");
                }
                need(parameters[0]!=parameters[1],"collapsed-native-wire-segment");
                auto parameterKey=parameters;std::sort(parameterKey.begin(),parameterKey.end());
                segment.parameterKey=parameterKey;
                for(const int face:segment.faces)if(model.faceOwners(face).size()>1)
                    sharedChains[face][item.bodyTag].insert(parameterKey);
                segment.chain.push_back(edge[0]);
                for(unsigned j=1;j<options.segmentsPerExistingSegment;++j){tick();
                    need(newNodes<options.maximumAddedNodes,"added-node-budget",OcctWireRefinementReason::BudgetExceeded);
                    // Evaluate in canonical parameter order, even when bodies
                    // use opposite local endpoint order, for bitwise identity.
                    const unsigned ordinal=parameters[0]<parameters[1]?j:options.segmentsPerExistingSegment-j;
                    const auto t=parameterKey[0]+(parameterKey[1]-parameterKey[0])*(double(ordinal)/options.segmentsPerExistingSegment);
                    const auto point=wire->pointAt(t);need(finite(point),"nonfinite-new-wire-point");
                    // Check every actual use, including both branches of seams.
                    for(const int face:segment.faces){bool found=false;
                        for(const auto& ring:model.boundaryUses(face))for(const auto& use:ring.edges){tick();
                            if(use.edgeTag!=options.edgeTag||use.degenerated)continue;
                            Vec2 uv;tick();++out.nativeUseChecks;
                            need(model.parameterOnUse(use.useId,t,uv)&&model.face(face)->containsParam(uv)&&
                                 (model.face(face)->pointAt(uv)-point).norm()<=std::max(rawTolerance,model.tolerance(2,face)),"native-pcurve-or-trim-disagreement");
                            found=true;
                        }
                        need(found,"missing-incident-native-use");
                    }
                    segment.chain.push_back(candidate.addNode(point,{1,options.edgeTag}));++newNodes;
                }
                segment.chain.push_back(edge[1]);
            }
            need(selected>0,"selected-wire-has-no-surface-segment");
            for(Index ci=0;ci<surface.cellCount();++ci){tick();const auto* ids=surface.cellNodes(ci);int split=-1;std::vector<std::pair<std::array<double,2>,Edge>> splits;
                for(unsigned j=0;j<3;++j)if(!segments.at(key(ids[j],ids[(j+1)%3])).chain.empty()){
                    const auto edge=key(ids[j],ids[(j+1)%3]);splits.push_back({segments.at(edge).parameterKey,edge});split=int(j);
                }
                if(split<0){candidate.addCell(mesh::CellType::Tri3,std::span<const Index>(ids,3),surface.cell(ci).onEntity);continue;}
                if(splits.size()>1){
                    // Canonical parameter order gives coowners the same
                    // diagonals even when their face windings are reversed.
                    std::sort(splits.begin(),splits.end());std::vector<std::array<Index,3>> local{{ids[0],ids[1],ids[2]}};
                    for(const auto& [parameter,edge]:splits){(void)parameter;const auto& chain=segments.at(edge).chain;std::vector<std::array<Index,3>> next;bool found=false;
                        for(auto t:local){bool changed=false;for(unsigned side=0;side<3;++side)if(key(t[side],t[(side+1)%3])==edge){need(!found,"duplicate-local-wire-edge");found=true;changed=true;const bool forward=chain.front()==t[side];
                            for(std::size_t j=0;j+1<chain.size();++j){tick();next.push_back({forward?chain[j]:chain[chain.size()-1-j],forward?chain[j+1]:chain[chain.size()-2-j],t[(side+2)%3]});}break;}
                            if(!changed)next.push_back(t);}
                        need(found,"missing-local-wire-edge");local.swap(next);
                    }
                    for(auto t:local){tick();candidate.addCell(mesh::CellType::Tri3,t,surface.cell(ci).onEntity);}newTriangles+=local.size()-1;continue;
                }
                const auto a=ids[split],b=ids[(split+1)%3],c=ids[(split+2)%3];
                const auto& chain=segments.at(key(a,b)).chain;const bool forward=chain.front()==a;
                for(std::size_t j=0;j+1<chain.size();++j){tick();
                    const auto p=forward?chain[j]:chain[chain.size()-1-j];
                    const auto q=forward?chain[j+1]:chain[chain.size()-2-j];
                    const std::array triangle{p,q,c};candidate.addCell(mesh::CellType::Tri3,triangle,surface.cell(ci).onEntity);
                }
                newTriangles+=chain.size()-2;
            }
            need(std::uint64_t(candidate.nodeCount())<=options.limits.maximumBoundaryNodes&&
                 std::uint64_t(candidate.cellCount())<=options.limits.maximumBoundaryTriangles,"output-count-budget",OcctWireRefinementReason::BudgetExceeded);
            const auto outputBytes=std::uint64_t(candidate.nodeCount())*(sizeof(Vec3)+sizeof(mesh::Classification))+
                std::uint64_t(candidate.cellCount())*sizeof(mesh::Cell)+std::uint64_t(candidate.connectivityCount())*sizeof(Index);
            need(outputBytes<=options.limits.maximumOutputBytes,"output-payload-budget",OcctWireRefinementReason::BudgetExceeded);
            certificate(candidate);rebuilt.push_back(std::move(candidate));
        }
        for(const auto& [face,owners]:sharedChains){const auto& expected=owners.begin()->second;
            need(owners.size()==model.faceOwners(face).size(),"missing-shared-face-chain");
            for(const auto& [body,chain]:owners){tick(chain.size());need(chain==expected,"incompatible-shared-wire-chains");}
        }
        tick(0); // The only mutation is a nonthrowing batch publication.
        for(std::size_t i=0;i<inputs.size();++i)inputs[i].boundary->swap(rebuilt[i]);
        out.published=true;out.reason=OcctWireRefinementReason::None;out.diagnostic="none";
        out.addedNodes=newNodes;out.addedTriangles=newTriangles;
    }catch(const Stop& stop){out.reason=stop.reason;out.diagnostic=stop.diagnostic;}
    catch(const Standard_OutOfMemory&){out.reason=OcctWireRefinementReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(const std::bad_alloc&){out.reason=OcctWireRefinementReason::AllocationFailed;out.diagnostic="allocation-failed";}
    catch(const Standard_Failure&){out.reason=OcctWireRefinementReason::NativeQueryFailed;out.diagnostic="native-query-failed";}
    catch(...){out.reason=OcctWireRefinementReason::InternalError;out.diagnostic="internal-error";}
    return out;
}
} // namespace paacs::geo::detail
