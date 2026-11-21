// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/topo/ElementTraits.hpp"
#include "paacs/algo/Meshers.hpp"

#include <fstream>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

namespace paacs::io {

using mesh::CellType;
using mesh::Mesh;

static int vtkType(CellType t) {
    // Single source of truth: ElementTraits (M4-B.0).
    paacs::topo::ElementType et;
    if (!paacs::topo::elementTypeOf(t, et)) return 0;
    return paacs::topo::traits(et).vtkId;
}

static void serializeVtk(const Mesh& m, std::ostream& os, bool volumeOnly) {
    os.imbue(std::locale::classic());
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
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

    // Both parts of the CAD classification are needed: face and body tags can
    // share the same integer. Preserve node classification as well as cells.
    // Keep CELL_DATA/entity_tag first for existing PAACS readers that consume
    // just this original field and ignore the additional classifications.
    os << "CELL_DATA " << nCells << "\nSCALARS entity_tag int 1\n"
       << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < m.cellCount(); ++i)
        if (keep(m.cell(i).type)) os << m.cell(i).onEntity.tag << '\n';
    os << "SCALARS entity_dimension int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < m.cellCount(); ++i)
        if (keep(m.cell(i).type)) os << m.cell(i).onEntity.dim << '\n';
    os << "POINT_DATA " << m.nodeCount() << "\nSCALARS entity_tag int 1\n"
       << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < m.nodeCount(); ++i)
        os << m.nodeClass(i).tag << '\n';
    os << "SCALARS entity_dimension int 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < m.nodeCount(); ++i)
        os << m.nodeClass(i).dim << '\n';
}

void writeVtkFile(const Mesh& m, const std::filesystem::path& path, bool volumeOnly) {
    if (path.empty()) throw std::runtime_error("empty VTK output path");
    const auto destination = std::filesystem::absolute(path);
    struct Temporary {
        std::filesystem::path path;
        bool owned = false;
        ~Temporary() {
            if (owned) {
                std::error_code error;
                std::filesystem::remove(path, error);
            }
        }
    } temporary;
    std::ofstream os;
    static std::atomic<std::uint64_t> sequence{0};
    const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
    for (unsigned attempt = 0; attempt < 32; ++attempt) {
        temporary.path = destination;
        temporary.path += ".paacs-tmp-" + std::to_string(tick) + "-" +
                          std::to_string(sequence.fetch_add(1, std::memory_order_relaxed));
        try {
            os.open(temporary.path, std::ios::out | std::ios::binary | std::ios::noreplace);
        } catch (...) {
            temporary.owned = os.is_open();
            throw;
        }
        temporary.owned = os.is_open();
        if (temporary.owned) break;
        std::error_code error;
        if (!std::filesystem::exists(temporary.path, error) || error)
            throw std::runtime_error("cannot create VTK output file");
        os.clear();
    }
    if (!temporary.owned) throw std::runtime_error("cannot create unique VTK output file");
    serializeVtk(m, os, volumeOnly);
    os.close();
    if (!os) throw std::runtime_error("cannot finish writing VTK output file");
#if defined(_WIN32)
    if (!MoveFileExW(temporary.path.c_str(), destination.c_str(), MOVEFILE_REPLACE_EXISTING))
        throw std::runtime_error("cannot publish VTK output file");
#else
    std::error_code error;
    std::filesystem::rename(temporary.path, destination, error);
    if (error) throw std::runtime_error("cannot publish VTK output file");
#endif
    temporary.owned = false;
}

void writeVtk(const Mesh& m, const std::string& path, bool volumeOnly) {
    writeVtkFile(m, std::filesystem::path(path), volumeOnly);
}

} // namespace paacs::io
