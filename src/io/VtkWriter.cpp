// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/topo/ElementTraits.hpp"
#include "paacs/algo/Meshers.hpp"

#include <fstream>
#include <stdexcept>

namespace paacs::io {

using mesh::CellType;
using mesh::Mesh;

static int vtkType(CellType t) {
    // Single source of truth: ElementTraits (M4-B.0).
    paacs::topo::ElementType et;
    if (!paacs::topo::elementTypeOf(t, et)) return 0;
    return paacs::topo::traits(et).vtkId;
}

void writeVtk(const Mesh& m, const std::string& path, bool volumeOnly) {
    std::ofstream os(path);
    if (!os) throw std::runtime_error("cannot open " + path);

    auto keep = [&](CellType t) {
        return !volumeOnly || t == CellType::Tet4 || t == CellType::Tet10 ||
               t == CellType::Prism6 || t == CellType::Hex8 ||
               t == CellType::Pyramid5;
    };

    os << "# vtk DataFile Version 3.0\nPAACS Mesh output\nASCII\n"
       << "DATASET UNSTRUCTURED_GRID\n"
       << "POINTS " << m.nodeCount() << " double\n";
    for (Index i = 0; i < m.nodeCount(); ++i) {
        const auto& p = m.node(i);
        os << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }

    Index nCells = 0, listLen = 0;
    for (Index i = 0; i < m.cellCount(); ++i)
        if (keep(m.cell(i).type)) {
            ++nCells;
            listLen += 1 + mesh::nodesPerCell(m.cell(i).type);
        }

    os << "CELLS " << nCells << ' ' << listLen << '\n';
    for (Index i = 0; i < m.cellCount(); ++i) {
        const auto& c = m.cell(i);
        if (!keep(c.type)) continue;
        const int nn = mesh::nodesPerCell(c.type);
        os << nn;
        const Index* v = m.cellNodes(i);
        for (int k = 0; k < nn; ++k) os << ' ' << v[k];
        os << '\n';
    }

    os << "CELL_TYPES " << nCells << '\n';
    for (Index i = 0; i < m.cellCount(); ++i)
        if (keep(m.cell(i).type)) os << vtkType(m.cell(i).type) << '\n';

    // Classification tags as cell data - lets ParaView colour by CAD entity.
    os << "CELL_DATA " << nCells << "\nSCALARS entity_tag int 1\n"
       << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < m.cellCount(); ++i)
        if (keep(m.cell(i).type)) os << m.cell(i).onEntity.tag << '\n';
}

} // namespace paacs::io
