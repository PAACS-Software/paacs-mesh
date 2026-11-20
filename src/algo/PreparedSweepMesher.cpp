// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/SmallSort.hpp"
#include "paacs/algo/PreparedSweepMesher.hpp"

#include "paacs/core/Parallel.hpp"
#include "paacs/mesh/CertifiedMeshStorage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <vector>

namespace paacs::algo {
namespace {

struct EdgeKey {
    Index a = kInvalidIndex, b = kInvalidIndex;
    friend bool operator<(const EdgeKey& x, const EdgeKey& y) noexcept {
        return x.a < y.a || (x.a == y.a && x.b < y.b);
    }
};
struct EdgeUse { std::uint32_t count = 0; Index a = -1, b = -1; };
struct SourceCell {
    mesh::CellType type = mesh::CellType::Tri3;
    std::array<Index,4> nodes{};
    std::uint8_t count = 0;
    std::uint64_t connectivityOffset = 0;
    mesh::Classification classification{};
};
struct SourceAnalysis {
    std::vector<SourceCell> cells;
    std::vector<std::array<Index,2>> boundaryEdges;
    std::vector<mesh::Classification> boundaryClassifications;
    std::uint64_t triangles = 0, quads = 0;
    std::uint64_t connectivityPerLayer = 0;
};

bool finite(const Vec3& p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool validClass(mesh::Classification c, int dim) noexcept {
    return c.dim == dim && c.tag >= 0;
}
EdgeKey key(Index a, Index b) noexcept {
    return a < b ? EdgeKey{a,b} : EdgeKey{b,a};
}
bool checkedAdd(std::uint64_t a, std::uint64_t b,
                std::uint64_t& out) noexcept {
    if (b > (std::numeric_limits<std::uint64_t>::max)() - a) return false;
    out = a + b; return true;
}
bool checkedMul(std::uint64_t a, std::uint64_t b,
                std::uint64_t& out) noexcept {
    if (a && b > (std::numeric_limits<std::uint64_t>::max)() / a) return false;
    out = a * b; return true;
}
PreparedSweepReason analyzeSource(const mesh::Mesh& input,
                                  SourceAnalysis& out) {
    if (input.nodeCount() == 0 || input.cellCount() == 0)
        return PreparedSweepReason::EmptySource;
    std::map<EdgeKey, EdgeUse> edges;
    out.cells.reserve(static_cast<std::size_t>(input.cellCount()));
    for (Index ci=0; ci<input.cellCount(); ++ci) {
        const auto& c=input.cell(ci);
        if (c.type != mesh::CellType::Tri3 &&
            c.type != mesh::CellType::Quad4)
            return PreparedSweepReason::UnsupportedSourceCell;
        if (!validClass(c.onEntity,2))
            return PreparedSweepReason::InvalidClassification;
        SourceCell row;
        row.type=c.type;
        row.count=static_cast<std::uint8_t>(mesh::nodesPerCell(c.type));
        row.connectivityOffset=out.connectivityPerLayer;
        out.connectivityPerLayer += c.type==mesh::CellType::Tri3 ? 6U : 8U;
        row.classification=c.onEntity;
        std::array<Index,4> sorted{};
        if(c.first<0 || c.first>input.connectivityCount() ||
           row.count>input.connectivityCount()-c.first)
            return PreparedSweepReason::InvalidSourceConnectivity;
        const Index* conn=input.cellNodes(ci);
        for (std::uint8_t n=0;n<row.count;++n) {
            row.nodes[n]=conn[n]; sorted[n]=conn[n];
            if (conn[n]<0 || conn[n]>=input.nodeCount() ||
                !finite(input.node(conn[n])))
                return PreparedSweepReason::InvalidSourceConnectivity;
        }
        if (!core::sortSmallPrefix(sorted, row.count))
            return PreparedSweepReason::InvalidSourceConnectivity;
        if (std::adjacent_find(sorted.begin(), sorted.begin()+row.count) !=
            sorted.begin()+row.count)
            return PreparedSweepReason::InvalidSourceConnectivity;
        for (std::uint8_t e=0;e<row.count;++e) {
            const Index a=row.nodes[e];
            const Index b=row.nodes[(e+1U)%row.count];
            EdgeUse& use=edges[key(a,b)];
            if (use.count==0) { use.a=a; use.b=b; }
            else if (use.count==1 && (use.a!=b || use.b!=a))
                return PreparedSweepReason::InconsistentSourceOrientation;
            if (++use.count>2) return PreparedSweepReason::NonmanifoldSource;
        }
        if (c.type==mesh::CellType::Tri3) ++out.triangles; else ++out.quads;
        out.cells.push_back(row);
    }
    for (const auto& [ignored,use]:edges) {
        (void)ignored;
        if (use.count==1) out.boundaryEdges.push_back({use.a,use.b});
        else if (use.count!=2) return PreparedSweepReason::NonmanifoldSource;
    }
    return PreparedSweepReason::None;
}

double signedTet(const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&d) noexcept {
    return (b-a).cross(c-a).dot(d-a)/6.0;
}
double wedgeVolume(const std::array<Vec3,6>& p) noexcept {
    // Three non-overlapping tetrahedra with positive orientation for the
    // canonical Prism6 ordering.
    const std::array<double,3> volumes{
        signedTet(p[0],p[1],p[2],p[3]), signedTet(p[1],p[4],p[2],p[3]),
        signedTet(p[2],p[4],p[5],p[3])};
    double sum=0;
    for(double v:volumes) { if(!(v>0) || !std::isfinite(v)) return 0; sum+=v; }
    return sum;
}
double hexVolume(const std::array<Vec3,8>& p) noexcept {
    const std::array<double,6> volumes{
        signedTet(p[0],p[1],p[2],p[6]), signedTet(p[0],p[2],p[3],p[6]),
        signedTet(p[0],p[3],p[7],p[6]), signedTet(p[0],p[7],p[4],p[6]),
        signedTet(p[0],p[4],p[5],p[6]), signedTet(p[0],p[5],p[1],p[6])};
    double sum=0;
    for(double v:volumes) { if(!(v>0) || !std::isfinite(v)) return 0; sum+=v; }
    return sum;
}
PreparedSweepOutcome rejected(PreparedSweepReason reason) noexcept {
    PreparedSweepOutcome out;
    out.status=PreparedSweepStatus::Rejected;
    out.reason=reason; out.phase=PreparedSweepPhase::Preflight;
    return out;
}
Index layeredNode(Index local,std::uint32_t layer,
                  std::uint64_t sourceNodes) noexcept {
    return static_cast<Index>(static_cast<std::uint64_t>(layer)*sourceNodes+
                              static_cast<std::uint64_t>(local));
}
} // namespace

PreparedSweepAuditReport auditPreparedMixedSweep(const mesh::Mesh& value) noexcept {
    PreparedSweepAuditReport out;
    out.nodes=static_cast<std::uint64_t>(value.nodeCount());
    try {
        double minimum=(std::numeric_limits<double>::infinity)(), maximum=0.0;
        long double total=0.0L;
        for (Index ci=0;ci<value.cellCount();++ci) {
            const auto& c=value.cell(ci);
            const Index width=mesh::nodesPerCell(c.type);
            if (width==0) { ++out.unsupportedCells; continue; }
            if (c.first<0 || c.first>value.connectivityCount() ||
                width>value.connectivityCount()-c.first) {
                ++out.invalidConnectivity; continue;
            }
            const Index* conn=value.cellNodes(ci);
            if (c.type==mesh::CellType::Tri3) { ++out.boundaryTriangles; continue; }
            if (c.type==mesh::CellType::Quad4) { ++out.boundaryQuads; continue; }
            double volume=0.0;
            if (c.type==mesh::CellType::Prism6) {
                ++out.wedges; std::array<Vec3,6> p{};
                for (std::size_t n=0;n<p.size();++n) {
                    if (conn[n]<0 || conn[n]>=value.nodeCount()) { ++out.invalidConnectivity; continue; }
                    p[n]=value.node(conn[n]);
                }
                volume=wedgeVolume(p);
            } else if (c.type==mesh::CellType::Hex8) {
                ++out.hexes; std::array<Vec3,8> p{};
                for (std::size_t n=0;n<p.size();++n) {
                    if (conn[n]<0 || conn[n]>=value.nodeCount()) { ++out.invalidConnectivity; continue; }
                    p[n]=value.node(conn[n]);
                }
                volume=hexVolume(p);
            } else { ++out.unsupportedCells; continue; }
            if (!validClass(c.onEntity,3)) ++out.invalidClassifications;
            if (!std::isfinite(volume) || !(volume>0.0)) ++out.nonPositiveCells;
            else { minimum=std::min(minimum,volume); maximum=std::max(maximum,volume); total+=volume; }
        }
        if (out.wedges+out.hexes==0) minimum=0.0;
        out.minimumCellVolume=minimum; out.maximumCellVolume=maximum;
        out.totalVolume=static_cast<double>(total);
        out.valid=out.wedges+out.hexes!=0 && out.unsupportedCells==0 &&
                  out.invalidConnectivity==0 && out.invalidClassifications==0 &&
                  out.nonPositiveCells==0;
    } catch (const std::bad_alloc&) { out.complete=false; out.allocationFailed=true; }
    catch (...) { out.complete=false; }
    return out;
}

PreparedSweepOutcome generatePreparedMixedSweep(
        const mesh::Mesh& sourceSurface,
        const PreparedMixedSweepRequest& request,
        mesh::Mesh& destination) noexcept {
    if (request.layers==0) return rejected(PreparedSweepReason::InvalidLayerCount);
    if (!std::isfinite(request.cellToCellGrowth) || !(request.cellToCellGrowth>0.0) ||
        !std::isfinite(request.maximumAllowedGrowth) || request.maximumAllowedGrowth<1.0 ||
        std::max(request.cellToCellGrowth,1.0/request.cellToCellGrowth)>request.maximumAllowedGrowth)
        return rejected(PreparedSweepReason::InvalidGrowth);
    if (!validClass(request.volumeClassification,3) ||
        (request.emitBoundaryFaces &&
         (!validClass(request.targetClassification,2) ||
          !validClass(request.sideClassification,2))))
        return rejected(PreparedSweepReason::InvalidClassification);
    SourceAnalysis source;
    try {
        const auto reason=analyzeSource(sourceSurface,source);
        if (reason!=PreparedSweepReason::None) return rejected(reason);
        source.boundaryClassifications.assign(source.boundaryEdges.size(),request.sideClassification);
        if (!request.sideEdgeClassifications.empty()) {
            if (request.sideEdgeClassifications.size()!=source.boundaryEdges.size())
                return rejected(PreparedSweepReason::InvalidClassification);
            std::map<EdgeKey,mesh::Classification> mapped;
            for (const auto& row:request.sideEdgeClassifications)
                if (!validClass(row.classification,2) ||
                    !mapped.emplace(key(row.first,row.second),row.classification).second)
                    return rejected(PreparedSweepReason::InvalidClassification);
            for (std::size_t i=0;i<source.boundaryEdges.size();++i) {
                const auto edge=source.boundaryEdges[i];
                const auto it=mapped.find(key(edge[0],edge[1]));
                if (it==mapped.end())return rejected(PreparedSweepReason::InvalidClassification);
                source.boundaryClassifications[i]=it->second;
            }
        }
    } catch (const std::bad_alloc&) {
        PreparedSweepOutcome out; out.reason=PreparedSweepReason::AllocationFailed;
        out.phase=PreparedSweepPhase::Preflight; return out;
    } catch (...) {
        PreparedSweepOutcome out; out.reason=PreparedSweepReason::InternalError;
        out.phase=PreparedSweepPhase::Preflight; return out;
    }
    const std::uint64_t sourceNodes=static_cast<std::uint64_t>(sourceSurface.nodeCount());
    const bool explicitTarget=!request.targetNodes.empty();
    // Classifications are also consumed by the implicit translation route.
    // Validate them before either route can index the caller's array.
    if (!request.targetNodeClassifications.empty() &&
        request.targetNodeClassifications.size()!=static_cast<std::size_t>(sourceNodes))
        return rejected(PreparedSweepReason::InvalidTarget);
    for (const auto cls:request.targetNodeClassifications)
        if (cls.dim<0 || cls.dim>3 || cls.tag<0)
            return rejected(PreparedSweepReason::InvalidClassification);
    if (!request.intermediateNodeClassifications.empty() &&
        request.intermediateNodeClassifications.size()!=static_cast<std::size_t>(sourceNodes))
        return rejected(PreparedSweepReason::InvalidTarget);
    for (const auto cls:request.intermediateNodeClassifications)
        if (cls.dim<0 || cls.dim>3 || cls.tag<0)
            return rejected(PreparedSweepReason::InvalidClassification);
    Vec3 unit{};
    if (explicitTarget) {
        if (request.targetNodes.size()!=static_cast<std::size_t>(sourceNodes) ||
            (!request.targetNodeClassifications.empty() &&
             request.targetNodeClassifications.size()!=static_cast<std::size_t>(sourceNodes)))
            return rejected(PreparedSweepReason::InvalidTarget);
        for (const auto& p:request.targetNodes) if (!finite(p)) return rejected(PreparedSweepReason::InvalidTarget);
    } else {
        if (!finite(request.direction) || !(request.direction.norm2()>0.0) ||
            !std::isfinite(request.totalThickness) || !(request.totalThickness>0.0))
            return rejected(PreparedSweepReason::InvalidThickness);
        unit=request.direction.normalized();
    }
    std::uint64_t nodes=0,wedges=0,hexes=0,boundaryFaces=0,cells=0,connectivity=0;
    std::uint64_t sourceFaces=0,sideFaces=0,tmp=0;
    if (!checkedMul(sourceNodes,static_cast<std::uint64_t>(request.layers)+1U,nodes) ||
        !checkedMul(source.triangles,request.layers,wedges) ||
        !checkedMul(source.quads,request.layers,hexes) ||
        !checkedAdd(source.triangles,source.quads,sourceFaces))
        return rejected(PreparedSweepReason::CountOverflow);
    if (request.emitBoundaryFaces) {
        if (!checkedMul(source.boundaryEdges.size(),request.layers,sideFaces) ||
            !checkedMul(sourceFaces,2U,boundaryFaces) ||
            !checkedAdd(boundaryFaces,sideFaces,boundaryFaces))
            return rejected(PreparedSweepReason::CountOverflow);
    }
    if (!checkedAdd(wedges,hexes,cells) || !checkedAdd(cells,boundaryFaces,cells))
        return rejected(PreparedSweepReason::CountOverflow);
    if (!checkedMul(wedges,6U,connectivity) || !checkedMul(hexes,8U,tmp) ||
        !checkedAdd(connectivity,tmp,connectivity))
        return rejected(PreparedSweepReason::CountOverflow);
    if (request.emitBoundaryFaces) {
        std::uint64_t triBoundary=0,quadBoundary=0;
        if (!checkedMul(source.triangles,2U,triBoundary) ||
            !checkedMul(triBoundary,3U,tmp) ||
            !checkedAdd(connectivity,tmp,connectivity) ||
            !checkedMul(source.quads,2U,quadBoundary) ||
            !checkedAdd(quadBoundary,sideFaces,quadBoundary) ||
            !checkedMul(quadBoundary,4U,tmp) ||
            !checkedAdd(connectivity,tmp,connectivity))
            return rejected(PreparedSweepReason::CountOverflow);
    }
    if (nodes>request.maximumNodes || cells>request.maximumCells ||
        boundaryFaces>request.maximumBoundaryFaces ||
        nodes>static_cast<std::uint64_t>((std::numeric_limits<Index>::max)()) ||
        connectivity>static_cast<std::uint64_t>((std::numeric_limits<Index>::max)()))
        return rejected(PreparedSweepReason::BudgetExceeded);

    std::vector<double> fractions;
    double measuredGrowth=1;
    try {
        if(!request.layerFractions.empty()) {
            if(request.layerFractions.size()!=static_cast<std::size_t>(request.layers)+1U ||
               request.layerFractions.front()!=0 || request.layerFractions.back()!=1)
                return rejected(PreparedSweepReason::InvalidGrowth);
            fractions=request.layerFractions;
        } else {
            fractions.resize(static_cast<std::size_t>(request.layers)+1U);
            long double width=1.0L,total=0.0L;
            for (std::uint32_t l=0;l<request.layers;++l) {
                total+=width; fractions[l+1U]=static_cast<double>(total);
                width*=request.cellToCellGrowth;
            }
            if (!(total>0.0L) || !std::isfinite(static_cast<double>(total)))
                return rejected(PreparedSweepReason::InvalidGrowth);
            for (double& f:fractions) f/=static_cast<double>(total);
        }
        double previous=0;
        for(std::size_t i=1;i<fractions.size();++i) {
            const double width=fractions[i]-fractions[i-1];
            if(!(width>0) || !std::isfinite(width)) return rejected(PreparedSweepReason::InvalidGrowth);
            if(previous>0) measuredGrowth=std::max({measuredGrowth,width/previous,previous/width});
            previous=width;
        }
        if(!request.layerFractions.empty() && measuredGrowth>request.maximumAllowedGrowth*(1+128*std::numeric_limits<double>::epsilon()))
            return rejected(PreparedSweepReason::InvalidGrowth);
    } catch (const std::bad_alloc&) {
        PreparedSweepOutcome out; out.reason=PreparedSweepReason::AllocationFailed;
        out.phase=PreparedSweepPhase::Schedule; return out;
    }

    PreparedSweepAuditReport certifiedAudit{};
    bool useCertifiedAudit = !request.fullAudit && !explicitTarget;
    if (useCertifiedAudit) {
        double minimumSourceVolume =
            (std::numeric_limits<double>::infinity)();
        double maximumSourceVolume = 0.0;
        long double totalSourceVolume = 0.0L;
        const Vec3 displacement = unit * request.totalThickness;
        double coordinateScale=0, maximumSpan=std::hypot(displacement.x,displacement.y,displacement.z);
        for (const auto& cell : source.cells) {
            for(std::uint8_t i=0;i<cell.count;++i) {
                const auto p=sourceSurface.node(cell.nodes[i]);
                const auto q=p+displacement;
                coordinateScale=std::max({coordinateScale,std::abs(p.x),std::abs(p.y),std::abs(p.z),
                    std::abs(q.x),std::abs(q.y),std::abs(q.z)});
                const auto delta=p-sourceSurface.node(cell.nodes[0]);
                maximumSpan=std::max(maximumSpan,std::hypot(delta.x,delta.y,delta.z));
            }
            double volume = 0.0;
            if (cell.type == mesh::CellType::Tri3) {
                std::array<Vec3, 6> points{};
                for (std::size_t i = 0; i < 3; ++i) {
                    points[i] = sourceSurface.node(cell.nodes[i]);
                    points[i + 3] = points[i] + displacement;
                }
                volume = wedgeVolume(points);
            } else {
                std::array<Vec3, 8> points{};
                for (std::size_t i = 0; i < 4; ++i) {
                    points[i] = sourceSurface.node(cell.nodes[i]);
                    points[i + 4] = points[i] + displacement;
                }
                volume = hexVolume(points);
            }
            if (!std::isfinite(volume) || !(volume > 0.0))
                return rejected(PreparedSweepReason::NonPositiveCell);
            minimumSourceVolume = std::min(minimumSourceVolume, volume);
            maximumSourceVolume = std::max(maximumSourceVolume, volume);
            totalSourceVolume += volume;
        }
        double minimumLayerFraction =
            (std::numeric_limits<double>::infinity)();
        double maximumLayerFraction = 0.0;
        for (std::uint32_t layer = 0; layer < request.layers; ++layer) {
            const double delta = fractions[static_cast<std::size_t>(layer) + 1U] -
                                 fractions[static_cast<std::size_t>(layer)];
            if (!std::isfinite(delta) || !(delta > 0.0))
                return rejected(PreparedSweepReason::InvalidGrowth);
            minimumLayerFraction = std::min(minimumLayerFraction, delta);
            maximumLayerFraction = std::max(maximumLayerFraction, delta);
        }
        certifiedAudit.nodes = nodes;
        certifiedAudit.wedges = wedges;
        certifiedAudit.hexes = hexes;
        if (request.emitBoundaryFaces) {
            certifiedAudit.boundaryTriangles = source.triangles * 2U;
            certifiedAudit.boundaryQuads = source.quads * 2U + sideFaces;
        }
        certifiedAudit.minimumCellVolume =
            minimumSourceVolume * minimumLayerFraction;
        certifiedAudit.maximumCellVolume =
            maximumSourceVolume * maximumLayerFraction;
        certifiedAudit.totalVolume = static_cast<double>(totalSourceVolume);
        certifiedAudit.valid =
            std::isfinite(certifiedAudit.minimumCellVolume) &&
            certifiedAudit.minimumCellVolume > 0.0 &&
            std::isfinite(certifiedAudit.maximumCellVolume) &&
            certifiedAudit.maximumCellVolume > 0.0 &&
            std::isfinite(certifiedAudit.totalVolume) &&
            certifiedAudit.totalVolume > 0.0;
        if (!certifiedAudit.valid)
            return rejected(PreparedSweepReason::AuditFailed);
        // Near coordinate-resolution limits, check actual stored cells rather
        // than relying on the ideal translated source and positive fractions.
        const long double guard = 512.L*std::numeric_limits<double>::epsilon()*
            coordinateScale*static_cast<long double>(maximumSpan)*maximumSpan;
        if (!(certifiedAudit.minimumCellVolume > guard)) useCertifiedAudit=false;
    }

    PreparedSweepOutcome out; out.phase=PreparedSweepPhase::Construction;
    try {
        mesh::detail::CertifiedMeshStorage storage(
            static_cast<std::size_t>(nodes),
            static_cast<std::size_t>(cells),
            static_cast<std::size_t>(connectivity));
        const std::uint64_t sourceCellCount=source.cells.size();

        core::parallelForWithThreads(
            0,static_cast<std::size_t>(nodes),request.executionThreads,
            [&](std::size_t ordinal){
                const std::uint64_t layer=ordinal/sourceNodes;
                const Index node=static_cast<Index>(ordinal%sourceNodes);
                const double f=fractions[static_cast<std::size_t>(layer)];
                const Vec3 p0=sourceSurface.node(node);
                const Vec3 p1=explicitTarget
                    ?request.targetNodes[static_cast<std::size_t>(node)]
                    :p0+unit*request.totalThickness;
                mesh::Classification cls=request.volumeClassification;
                if (!request.intermediateNodeClassifications.empty())
                    cls=request.intermediateNodeClassifications[static_cast<std::size_t>(node)];
                if(layer==0)cls=sourceSurface.nodeClass(node);
                else if(layer==request.layers)
                    cls=request.targetNodeClassifications.empty()
                        ?request.targetClassification
                        :request.targetNodeClassifications[
                            static_cast<std::size_t>(node)];
                storage.setNode(ordinal,p0+(p1-p0)*f,cls);
            });

        const std::uint64_t volumeCells=wedges+hexes;
        const std::uint64_t volumeConnectivity=
            source.connectivityPerLayer*request.layers;
        core::parallelForWithThreads(
            0,static_cast<std::size_t>(volumeCells),request.executionThreads,
            [&](std::size_t ordinal){
                const std::uint64_t layer=ordinal/sourceCellCount;
                const std::uint64_t local=ordinal%sourceCellCount;
                const auto& c=source.cells[static_cast<std::size_t>(local)];
                const std::uint64_t first=layer*source.connectivityPerLayer+
                                          c.connectivityOffset;
                if(c.type==mesh::CellType::Tri3){
                    const std::array<Index,6> conn{
                        layeredNode(c.nodes[0],layer,sourceNodes),
                        layeredNode(c.nodes[1],layer,sourceNodes),
                        layeredNode(c.nodes[2],layer,sourceNodes),
                        layeredNode(c.nodes[0],layer+1U,sourceNodes),
                        layeredNode(c.nodes[1],layer+1U,sourceNodes),
                        layeredNode(c.nodes[2],layer+1U,sourceNodes)};
                    storage.setCell(ordinal,static_cast<std::size_t>(first),
                                    mesh::CellType::Prism6,conn,
                                    request.volumeClassification);
                }else{
                    const std::array<Index,8> conn{
                        layeredNode(c.nodes[0],layer,sourceNodes),
                        layeredNode(c.nodes[1],layer,sourceNodes),
                        layeredNode(c.nodes[2],layer,sourceNodes),
                        layeredNode(c.nodes[3],layer,sourceNodes),
                        layeredNode(c.nodes[0],layer+1U,sourceNodes),
                        layeredNode(c.nodes[1],layer+1U,sourceNodes),
                        layeredNode(c.nodes[2],layer+1U,sourceNodes),
                        layeredNode(c.nodes[3],layer+1U,sourceNodes)};
                    storage.setCell(ordinal,static_cast<std::size_t>(first),
                                    mesh::CellType::Hex8,conn,
                                    request.volumeClassification);
                }
            });

        if(request.emitBoundaryFaces){
            const std::uint64_t boundaryPairCells=2U*sourceCellCount;
            core::parallelForWithThreads(
                0,static_cast<std::size_t>(sourceCellCount),
                request.executionThreads,[&](std::size_t ordinal){
                    const auto& c=source.cells[ordinal];
                    const std::uint64_t firstCell=volumeCells+2U*ordinal;
                    const std::uint64_t firstConnectivity=
                        volumeConnectivity+c.connectivityOffset;
                    if(c.type==mesh::CellType::Tri3){
                        const std::array<Index,3> bottom{
                            c.nodes[2],c.nodes[1],c.nodes[0]};
                        const std::array<Index,3> top{
                            layeredNode(c.nodes[0],request.layers,sourceNodes),
                            layeredNode(c.nodes[1],request.layers,sourceNodes),
                            layeredNode(c.nodes[2],request.layers,sourceNodes)};
                        storage.setCell(firstCell,firstConnectivity,
                            mesh::CellType::Tri3,bottom,c.classification);
                        storage.setCell(firstCell+1U,firstConnectivity+3U,
                            mesh::CellType::Tri3,top,request.targetClassification);
                    }else{
                        const std::array<Index,4> bottom{
                            c.nodes[3],c.nodes[2],c.nodes[1],c.nodes[0]};
                        const std::array<Index,4> top{
                            layeredNode(c.nodes[0],request.layers,sourceNodes),
                            layeredNode(c.nodes[1],request.layers,sourceNodes),
                            layeredNode(c.nodes[2],request.layers,sourceNodes),
                            layeredNode(c.nodes[3],request.layers,sourceNodes)};
                        storage.setCell(firstCell,firstConnectivity,
                            mesh::CellType::Quad4,bottom,c.classification);
                        storage.setCell(firstCell+1U,firstConnectivity+4U,
                            mesh::CellType::Quad4,top,request.targetClassification);
                    }
                });
            const std::uint64_t sideCellBase=volumeCells+boundaryPairCells;
            const std::uint64_t sideConnectivityBase=
                volumeConnectivity+source.connectivityPerLayer;
            core::parallelForWithThreads(
                0,static_cast<std::size_t>(sideFaces),request.executionThreads,
                [&](std::size_t ordinal){
                    const std::uint64_t layer=ordinal/source.boundaryEdges.size();
                    const std::uint64_t edgeOrdinal=
                        ordinal%source.boundaryEdges.size();
                    const auto edge=source.boundaryEdges[
                        static_cast<std::size_t>(edgeOrdinal)];
                    const std::array<Index,4> side{
                        layeredNode(edge[0],layer,sourceNodes),
                        layeredNode(edge[1],layer,sourceNodes),
                        layeredNode(edge[1],layer+1U,sourceNodes),
                        layeredNode(edge[0],layer+1U,sourceNodes)};
                    storage.setCell(sideCellBase+ordinal,
                        sideConnectivityBase+ordinal*4U,
                        mesh::CellType::Quad4,side,
                        source.boundaryClassifications[static_cast<std::size_t>(edgeOrdinal)]);
                });
        }
        mesh::Mesh local;
        storage.publishTo(local);
        out.phase=PreparedSweepPhase::Audit;
        out.metrics.usedFullGeometricAudit = !useCertifiedAudit;
        out.metrics.audit = useCertifiedAudit
            ? certifiedAudit
            : auditPreparedMixedSweep(local);
        if (!out.metrics.audit.complete || !out.metrics.audit.valid) {
            out.reason=PreparedSweepReason::AuditFailed; return out;
        }
        local.swap(destination);
        out.status=PreparedSweepStatus::Published; out.reason=PreparedSweepReason::None;
        out.phase=PreparedSweepPhase::Publication; out.destinationUnchanged=false;
        out.metrics.sourceNodes=sourceNodes; out.metrics.sourceTriangles=source.triangles;
        out.metrics.sourceQuads=source.quads; out.metrics.boundaryEdges=source.boundaryEdges.size();
        out.metrics.nodes=nodes; out.metrics.wedges=wedges; out.metrics.hexes=hexes;
        out.metrics.boundaryFaces=boundaryFaces; out.metrics.connectivityEntries=connectivity;
        out.metrics.maximumAdjacentGrowth=request.layerFractions.empty()
            ?std::max(request.cellToCellGrowth,1.0/request.cellToCellGrowth):measuredGrowth;
        out.metrics.usedExplicitLayerFractions=!request.layerFractions.empty();
        out.metrics.minimumCellVolume=out.metrics.audit.minimumCellVolume;
        out.metrics.totalVolume=out.metrics.audit.totalVolume;
        return out;
    } catch (const std::bad_alloc&) { out.reason=PreparedSweepReason::AllocationFailed; }
    catch (...) { out.reason=PreparedSweepReason::InternalError; }
    return out;
}

}  // namespace paacs::algo
