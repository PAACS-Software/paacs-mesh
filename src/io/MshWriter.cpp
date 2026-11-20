// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/topo/ElementTraits.hpp"
#include "paacs/algo/Meshers.hpp"

#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace paacs::io {

using mesh::CellType;
using mesh::Mesh;

namespace {
static int mshType(mesh::CellType t) {
    // Single source of truth: ElementTraits (M4-B.0). The previous local
    // switch duplicated the mapping (review finding).
    topo::ElementType et;
    if (!topo::elementTypeOf(t, et)) return 0;
    return topo::traits(et).mshId;
}
} // namespace

// MSH 4.1 ASCII. Nodes and elements are grouped into entity blocks by
// their (dim, tag) classification - this is what makes physical-group /
// boundary-condition workflows possible downstream.
void writeMsh(const Mesh& m, const std::string& path) {
    std::ofstream os(path);
    if (!os) throw std::runtime_error("cannot open " + path);

    // Group nodes and cells by classification.
    std::map<std::pair<int, int>, std::vector<Index>> nodeBlocks;
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto& c = m.nodeClass(i);
        nodeBlocks[{c.dim < 0 ? 3 : c.dim, c.tag < 0 ? 1 : c.tag}].push_back(i);
    }
    std::map<std::tuple<int, int, int>, std::vector<Index>> elemBlocks;
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i);
        const int type = mshType(c.type);
        if (!type) continue;
        elemBlocks[{c.onEntity.dim < 0 ? 3 : c.onEntity.dim,
                    c.onEntity.tag < 0 ? 1 : c.onEntity.tag, type}].push_back(i);
    }

    os << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n";

    // Minimal $Entities section with zero-filled bounding data.
    std::set<std::pair<int, int>> ents;
    for (const auto& [k, v] : nodeBlocks) ents.insert(k);
    for (const auto& [k, v] : elemBlocks)
        ents.insert({std::get<0>(k), std::get<1>(k)});
    int nE[4] = {0, 0, 0, 0};
    for (const auto& e : ents) ++nE[e.first];
    os << "$Entities\n" << nE[0] << ' ' << nE[1] << ' ' << nE[2] << ' '
       << nE[3] << '\n';
    for (int dim = 0; dim <= 3; ++dim)
        for (const auto& e : ents) {
            if (e.first != dim) continue;
            if (dim == 0) os << e.second << " 0 0 0 0\n";
            else os << e.second << " 0 0 0 0 0 0 0 0\n";  // bbox + 0 physicals + 0 bounding
        }
    os << "$EndEntities\n";

    // $Nodes: block-structured, 1-based global tags.
    os << "$Nodes\n" << nodeBlocks.size() << ' ' << m.nodeCount() << " 1 "
       << m.nodeCount() << '\n';
    for (const auto& [key, ids] : nodeBlocks) {
        os << key.first << ' ' << key.second << " 0 " << ids.size() << '\n';
        for (Index i : ids) os << (i + 1) << '\n';
        for (Index i : ids) {
            const auto& p = m.node(i);
            os << p.x << ' ' << p.y << ' ' << p.z << '\n';
        }
    }
    os << "$EndNodes\n";

    // $Elements.
    Index nElems = 0;
    for (const auto& [k, v] : elemBlocks) nElems += static_cast<Index>(v.size());
    os << "$Elements\n" << elemBlocks.size() << ' ' << nElems << " 1 "
       << nElems << '\n';
    Index etag = 1;
    for (const auto& [key, ids] : elemBlocks) {
        os << std::get<0>(key) << ' ' << std::get<1>(key) << ' '
           << std::get<2>(key) << ' ' << ids.size() << '\n';
        for (Index ci : ids) {
            os << etag++;
            const auto& c = m.cell(ci);
            const Index* v = m.cellNodes(ci);
            for (int k = 0; k < mesh::nodesPerCell(c.type); ++k)
                os << ' ' << (v[k] + 1);
            os << '\n';
        }
    }
    os << "$EndElements\n";
}

} // namespace paacs::io
