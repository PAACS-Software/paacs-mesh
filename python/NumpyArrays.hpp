// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "paacs/mesh/Mesh.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace paacs::python {
namespace numpy_detail {
namespace py = pybind11;

inline constexpr std::array cellTypes{
    mesh::CellType::Line2, mesh::CellType::Line3, mesh::CellType::Tri3,
    mesh::CellType::Tri6, mesh::CellType::Quad4, mesh::CellType::Tet4,
    mesh::CellType::Tet10, mesh::CellType::Prism6, mesh::CellType::Hex8,
    mesh::CellType::Pyramid5,
};
inline constexpr std::array cellNames{
    "line2", "line3", "tri3", "tri6", "quad4", "tet4", "tet10",
    "prism6", "hex8", "pyramid5",
};

inline std::size_t blockIndex(mesh::CellType type) {
    switch (type) {
        case mesh::CellType::Line2: return 0;
        case mesh::CellType::Line3: return 1;
        case mesh::CellType::Tri3: return 2;
        case mesh::CellType::Tri6: return 3;
        case mesh::CellType::Quad4: return 4;
        case mesh::CellType::Tet4: return 5;
        case mesh::CellType::Tet10: return 6;
        case mesh::CellType::Prism6: return 7;
        case mesh::CellType::Hex8: return 8;
        case mesh::CellType::Pyramid5: return 9;
    }
    throw std::invalid_argument("to_numpy encountered an unsupported native cell type");
}

// Check each array's dimensions and byte size against NumPy's signed size
// range, then accumulate the payload without wrapping the caller's budget.
inline void addArrayBytes(std::uint64_t& total, std::uint64_t rows,
                          std::uint64_t columns, std::uint64_t itemBytes) {
    constexpr auto signedLimit = static_cast<std::uint64_t>(
        (std::numeric_limits<py::ssize_t>::max)());
    if (rows > signedLimit || columns > signedLimit ||
        (rows != 0 && columns > signedLimit / rows))
        throw std::overflow_error("to_numpy array dimensions exceed the NumPy size range");
    const auto items = rows * columns;
    if (items > signedLimit / itemBytes)
        throw std::overflow_error("to_numpy array bytes exceed the NumPy size range");
    const auto bytes = items * itemBytes;
    if (bytes > (std::numeric_limits<std::uint64_t>::max)() - total)
        throw std::overflow_error("to_numpy total array bytes overflow");
    total += bytes;
}

inline void requireNumpy() {
    try {
        py::module_::import("numpy");
    } catch (const py::error_already_set& error) {
        if (error.matches(PyExc_ModuleNotFoundError)) {
            const py::object name = error.value().attr("name");
            if (name.equal(py::str("numpy")))
                throw py::import_error(
                    "Mesh.to_numpy requires NumPy; install paacs-mesh[numpy] or numpy.");
        }
        throw;
    }
}

struct BlockArrays {
    py::array_t<std::int64_t> connectivity;
    py::array_t<std::int32_t> dimensions, tags;
    py::array_t<std::int64_t> indices;
    std::int64_t* connectivityData;
    std::int32_t* dimensionData;
    std::int32_t* tagData;
    std::int64_t* indexData;
    py::ssize_t width;
    py::ssize_t next = 0;

    BlockArrays(py::ssize_t rows, py::ssize_t nodesPerCell)
        : connectivity(std::array{rows, nodesPerCell}), dimensions(rows), tags(rows), indices(rows),
          connectivityData(connectivity.mutable_data()), dimensionData(dimensions.mutable_data()),
          tagData(tags.mutable_data()), indexData(indices.mutable_data()), width(nodesPerCell) {}

    py::dict dictionary() const {
        // These arrays own their allocations. Re-enabling writes can only
        // change the snapshot; it cannot expose or mutate native mesh data.
        connectivity.attr("setflags")(false);
        dimensions.attr("setflags")(false);
        tags.attr("setflags")(false);
        indices.attr("setflags")(false);
        py::dict result;
        result["connectivity"] = connectivity;
        result["entity_dimensions"] = dimensions;
        result["entity_tags"] = tags;
        result["cell_indices"] = indices;
        return result;
    }
};
} // namespace numpy_detail

inline pybind11::dict meshArrays(const mesh::Mesh& native, std::uint64_t maximumBytes) {
    namespace py = pybind11;
    using namespace numpy_detail;
    static_assert(std::numeric_limits<int>::digits <= std::numeric_limits<std::int32_t>::digits,
                  "Native classifications must fit the public int32 array dtype");
    const auto nodeCount = native.nodeCount(), cellCount = native.cellCount();
    if (nodeCount < 0 || cellCount < 0)
        throw std::overflow_error("to_numpy encountered a negative native mesh count");
    std::array<std::uint64_t, cellTypes.size()> counts{};
    {
        py::gil_scoped_release release;
        for (Index cell = 0; cell < cellCount; ++cell)
            ++counts[blockIndex(native.cell(cell).type)];
    }

    std::uint64_t bytes = 0;
    const auto nodes = static_cast<std::uint64_t>(nodeCount);
    addArrayBytes(bytes, nodes, 3, sizeof(double));
    addArrayBytes(bytes, nodes, 1, sizeof(std::int32_t));
    addArrayBytes(bytes, nodes, 1, sizeof(std::int32_t));
    for (std::size_t block = 0; block < counts.size(); ++block) {
        if (counts[block] == 0) continue;
        const auto width = static_cast<std::uint64_t>(mesh::nodesPerCell(cellTypes[block]));
        addArrayBytes(bytes, counts[block], width, sizeof(std::int64_t));
        addArrayBytes(bytes, counts[block], 1, sizeof(std::int32_t));
        addArrayBytes(bytes, counts[block], 1, sizeof(std::int32_t));
        addArrayBytes(bytes, counts[block], 1, sizeof(std::int64_t));
    }
    if (bytes > maximumBytes) {
        const auto message = "Mesh.to_numpy requires " + std::to_string(bytes) +
            " array payload bytes, exceeding max_bytes=" + std::to_string(maximumBytes);
        PyErr_SetString(PyExc_MemoryError, message.c_str());
        throw py::error_already_set();
    }

    requireNumpy();
    const auto rows = static_cast<py::ssize_t>(nodeCount);
    py::array_t<double> points(std::array{rows, py::ssize_t{3}});
    py::array_t<std::int32_t> nodeDimensions(rows), nodeTags(rows);
    auto* pointData = points.mutable_data();
    auto* nodeDimensionData = nodeDimensions.mutable_data();
    auto* nodeTagData = nodeTags.mutable_data();
    // An absent optional constructs no ndarray, including for empty blocks.
    std::array<std::optional<BlockArrays>, cellTypes.size()> blocks;
    for (std::size_t block = 0; block < counts.size(); ++block)
        if (counts[block] != 0)
            blocks[block].emplace(static_cast<py::ssize_t>(counts[block]),
                                  mesh::nodesPerCell(cellTypes[block]));
    {
        // All arrays and pointers are ready before releasing the GIL. They
        // remain local until copying finishes; this scope uses only C++ data.
        py::gil_scoped_release release;
        for (Index node = 0; node < nodeCount; ++node) {
            const auto& point = native.node(node);
            pointData[3 * node] = point.x;
            pointData[3 * node + 1] = point.y;
            pointData[3 * node + 2] = point.z;
            const auto classification = native.nodeClass(node);
            nodeDimensionData[node] = static_cast<std::int32_t>(classification.dim);
            nodeTagData[node] = static_cast<std::int32_t>(classification.tag);
        }
        for (Index cell = 0; cell < cellCount; ++cell) {
            const auto& source = native.cell(cell);
            auto& block = *blocks[blockIndex(source.type)];
            const auto row = block.next++;
            const auto* connectivity = native.cellNodes(cell);
            for (py::ssize_t local = 0; local < block.width; ++local)
                block.connectivityData[row * block.width + local] = connectivity[local];
            block.dimensionData[row] = static_cast<std::int32_t>(source.onEntity.dim);
            block.tagData[row] = static_cast<std::int32_t>(source.onEntity.tag);
            block.indexData[row] = cell;
        }
    }
    points.attr("setflags")(false);
    nodeDimensions.attr("setflags")(false);
    nodeTags.attr("setflags")(false);
    py::dict cellBlocks, result;
    for (std::size_t block = 0; block < blocks.size(); ++block)
        if (blocks[block]) cellBlocks[cellNames[block]] = blocks[block]->dictionary();
    result["points"] = points;
    result["node_entity_dimensions"] = nodeDimensions;
    result["node_entity_tags"] = nodeTags;
    result["cells"] = cellBlocks;
    result["units"] = "mm";
    return result;
}
} // namespace paacs::python
