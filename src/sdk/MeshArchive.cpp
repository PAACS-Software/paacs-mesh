// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshArchive.hpp"
#include "paacs/sdk/Version.hpp"
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <new>
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

namespace paacs::sdk {
namespace {
constexpr std::array<std::uint8_t, 8> kMagic{'P', 'A', 'A', 'C', 'S', 'M', '8', 0};
constexpr std::uint64_t kHeaderBytes = 40, kNodeBytes = 32, kCellBytes = 12;
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
static_assert(std::endian::native == std::endian::little ||
              std::endian::native == std::endian::big);

struct Writer {
    std::span<std::uint8_t> out;
    std::size_t at = 0;
    void u8(std::uint8_t v) { out[at++] = v; }
    template <class T> void integer(T v) {
        if constexpr (std::endian::native == std::endian::big)
            v = std::byteswap(v);
        std::memcpy(out.data() + at, &v, sizeof v);
        at += sizeof v;
    }
    void f64(double v) { integer(std::bit_cast<std::uint64_t>(v)); }
};
struct Reader {
    std::span<const std::uint8_t> in;
    std::size_t at = 0;
    bool ok = true;
    std::uint8_t u8() {
        if (at == in.size()) {
            ok = false;
            return 0;
        }
        return in[at++];
    }
    template <class T> T integer() {
        if (sizeof(T) > in.size() - at) {
            ok = false;
            return 0;
        }
        T v;
        std::memcpy(&v, in.data() + at, sizeof v);
        at += sizeof v;
        if constexpr (std::endian::native == std::endian::big)
            v = std::byteswap(v);
        return v;
    }
    double f64() { return std::bit_cast<double>(integer<std::uint64_t>()); }
    mesh::Classification classification() {
        return {integer<std::int32_t>(), integer<std::int32_t>()};
    }
};
MeshArchiveOutcome failure(MeshArchiveReason reason) noexcept {
    return {false, reason, 0};
}
struct TemporaryArchive {
    std::filesystem::path path;
    bool owned = false;
    ~TemporaryArchive() {
        if (owned) {
            std::error_code error;
            std::filesystem::remove(path, error);
        }
    }
};
bool validClassification(mesh::Classification c) noexcept {
    return c.dim >= 0 && c.dim <= 3 && c.tag >= 0;
}
bool validPoint(const Vec3 &p) noexcept {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
bool encodedSize(std::uint64_t nodes, std::uint64_t cells, std::uint64_t connectivity,
                 std::uint64_t &size) noexcept {
    size = kHeaderBytes;
    for (auto [count, width] : {std::pair{nodes, kNodeBytes}, std::pair{cells, kCellBytes},
                                std::pair{connectivity, std::uint64_t{8}}}) {
        if (count > (std::numeric_limits<std::uint64_t>::max() - size) / width)
            return false;
        size += count * width;
    }
    return size <= std::numeric_limits<std::size_t>::max();
}
bool valid(const mesh::Mesh &m) noexcept {
    for (Index n = 0; n < m.nodeCount(); ++n)
        if (!validPoint(m.node(n)) || !validClassification(m.nodeClass(n)))
            return false;
    std::uint64_t count = 0;
    for (Index ci = 0; ci < m.cellCount(); ++ci) {
        const auto &c = m.cell(ci);
        const int width = mesh::nodesPerCell(c.type);
        if (width <= 0 || !validClassification(c.onEntity) || c.first < 0 ||
            static_cast<std::uint64_t>(c.first) != count)
            return false;
        count += static_cast<std::uint64_t>(width);
        if (count > static_cast<std::uint64_t>(m.connectivityCount()))
            return false;
        for (int k = 0; k < width; ++k)
            if (m.cellNodes(ci)[k] < 0 || m.cellNodes(ci)[k] >= m.nodeCount())
                return false;
    }
    return count == static_cast<std::uint64_t>(m.connectivityCount());
}
} // namespace

MeshArchiveOutcome serializeMeshArchive(const mesh::Mesh &value,
                                        std::vector<std::uint8_t> &bytes) noexcept {
    if (!valid(value))
        return failure(MeshArchiveReason::InvalidMesh);
    try {
        std::uint64_t size;
        if (!encodedSize(value.nodeCount(), value.cellCount(), value.connectivityCount(), size))
            return failure(MeshArchiveReason::CountOverflow);
        // One exact allocation; fixed-width stores replace byte-at-a-time growth.
        std::vector<std::uint8_t> buffer(static_cast<std::size_t>(size));
        Writer w{buffer};
        for (auto b : kMagic)
            w.u8(b);
        w.integer(kMeshArchiveVersion);
        w.integer(std::uint32_t{0});
        w.integer(static_cast<std::uint64_t>(value.nodeCount()));
        w.integer(static_cast<std::uint64_t>(value.cellCount()));
        w.integer(static_cast<std::uint64_t>(value.connectivityCount()));
        for (Index n = 0; n < value.nodeCount(); ++n) {
            const auto &p = value.node(n);
            const auto c = value.nodeClass(n);
            w.f64(p.x);
            w.f64(p.y);
            w.f64(p.z);
            w.integer(static_cast<std::int32_t>(c.dim));
            w.integer(static_cast<std::int32_t>(c.tag));
        }
        for (Index ci = 0; ci < value.cellCount(); ++ci) {
            const auto &c = value.cell(ci);
            w.u8(static_cast<std::uint8_t>(c.type));
            w.u8(0);
            w.u8(0);
            w.u8(0);
            w.integer(static_cast<std::int32_t>(c.onEntity.dim));
            w.integer(static_cast<std::int32_t>(c.onEntity.tag));
            for (int k = 0; k < mesh::nodesPerCell(c.type); ++k)
                w.integer(static_cast<std::int64_t>(value.cellNodes(ci)[k]));
        }
        bytes.swap(buffer);
        return {true, MeshArchiveReason::None, size};
    } catch (const std::bad_alloc &) {
        return failure(MeshArchiveReason::AllocationFailed);
    } catch (const std::length_error &) {
        return failure(MeshArchiveReason::CountOverflow);
    } catch (...) {
        return failure(MeshArchiveReason::InternalError);
    }
}

MeshArchiveOutcome deserializeMeshArchive(std::span<const std::uint8_t> bytes, mesh::Mesh &value,
                                          const MeshArchiveLimits &limits) noexcept {
    if (bytes.size() > limits.maxBytes)
        return failure(MeshArchiveReason::ResourceLimit);
    if (bytes.size() < kHeaderBytes)
        return failure(MeshArchiveReason::Truncated);
    Reader r{bytes};
    for (auto expected : kMagic)
        if (r.u8() != expected)
            return failure(MeshArchiveReason::InvalidHeader);
    if (r.integer<std::uint32_t>() != kMeshArchiveVersion)
        return failure(MeshArchiveReason::UnsupportedVersion);
    if (r.integer<std::uint32_t>() != 0)
        return failure(MeshArchiveReason::InvalidHeader);
    const auto nodes = r.integer<std::uint64_t>();
    const auto cells = r.integer<std::uint64_t>();
    const auto connectivity = r.integer<std::uint64_t>();
    const auto indexLimit = static_cast<std::uint64_t>(std::numeric_limits<Index>::max());
    std::uint64_t expected;
    if (nodes > indexLimit || cells > indexLimit || connectivity > indexLimit ||
        !encodedSize(nodes, cells, connectivity, expected))
        return failure(MeshArchiveReason::CountOverflow);
    if (expected > bytes.size())
        return failure(MeshArchiveReason::Truncated);
    if (expected != bytes.size())
        return failure(MeshArchiveReason::InvalidHeader);
    if (nodes > limits.maxNodes || cells > limits.maxCells || connectivity > limits.maxConnectivity)
        return failure(MeshArchiveReason::ResourceLimit);
    // Complete structural preflight BEFORE any allocation controlled by input.
    for (std::uint64_t n = 0; n < nodes; ++n) {
        const Vec3 point{r.f64(), r.f64(), r.f64()};
        if (!validPoint(point) || !validClassification(r.classification()))
            return failure(MeshArchiveReason::InvalidMesh);
    }
    std::uint64_t readConnectivity = 0;
    for (std::uint64_t ci = 0; ci < cells; ++ci) {
        const int width = mesh::nodesPerCell(static_cast<mesh::CellType>(r.u8()));
        const auto reserved1 = r.u8(), reserved2 = r.u8(), reserved3 = r.u8();
        const auto classification = r.classification();
        if (!r.ok)
            return failure(MeshArchiveReason::Truncated);
        if (reserved1 || reserved2 || reserved3)
            return failure(MeshArchiveReason::InvalidHeader);
        if (width <= 0 || static_cast<std::uint64_t>(width) > connectivity - readConnectivity)
            return failure(MeshArchiveReason::InvalidCell);
        if (!validClassification(classification))
            return failure(MeshArchiveReason::InvalidMesh);
        for (int k = 0; k < width; ++k) {
            const auto index = r.integer<std::int64_t>();
            if (!r.ok)
                return failure(MeshArchiveReason::Truncated);
            if (index < 0 || static_cast<std::uint64_t>(index) >= nodes)
                return failure(MeshArchiveReason::InvalidCell);
        }
        readConnectivity += static_cast<std::uint64_t>(width);
    }
    if (!r.ok || r.at != bytes.size() || readConnectivity != connectivity)
        return failure(MeshArchiveReason::InvalidHeader);
    try {
        mesh::Mesh local;
        local.reserveAppend(static_cast<std::size_t>(nodes), static_cast<std::size_t>(cells),
                            static_cast<std::size_t>(connectivity));
        r.at = static_cast<std::size_t>(kHeaderBytes);
        for (std::uint64_t n = 0; n < nodes; ++n) {
            const Vec3 point{r.f64(), r.f64(), r.f64()};
            (void)local.addNode(point, r.classification());
        }
        std::array<Index, 10> conn{};
        for (std::uint64_t ci = 0; ci < cells; ++ci) {
            const auto type = static_cast<mesh::CellType>(r.u8());
            r.at += 3;
            const auto classification = r.classification();
            const auto width = static_cast<std::size_t>(mesh::nodesPerCell(type));
            for (std::size_t k = 0; k < width; ++k)
                conn[k] = r.integer<std::int64_t>();
            (void)local.addCell(type, std::span<const Index>(conn.data(), width), classification);
        }
        local.swap(value);
        return {true, MeshArchiveReason::None, bytes.size()};
    } catch (const std::bad_alloc &) {
        return failure(MeshArchiveReason::AllocationFailed);
    } catch (const std::length_error &) {
        return failure(MeshArchiveReason::CountOverflow);
    } catch (...) {
        return failure(MeshArchiveReason::InternalError);
    }
}
MeshArchiveOutcome deserializeMeshArchive(const std::vector<std::uint8_t> &bytes,
                                          mesh::Mesh &value) noexcept {
    return deserializeMeshArchive(std::span<const std::uint8_t>(bytes), value);
}
MeshArchiveOutcome writeMeshArchive(const mesh::Mesh &value, const std::string &path) noexcept {
    std::vector<std::uint8_t> bytes;
    auto out = serializeMeshArchive(value, bytes);
    if (!out.completed)
        return out;
    if (bytes.size() > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max()))
        return failure(MeshArchiveReason::CountOverflow);
    try {
        // Build privately next to the destination, so final publication uses a
        // same-filesystem rename. noreplace prevents a colliding temporary name
        // from truncating an existing file, including across concurrent writers.
        if (path.empty())
            return failure(MeshArchiveReason::IoError);
        const auto destination = std::filesystem::absolute(std::filesystem::path(path));
        TemporaryArchive temporary;
        std::ofstream stream;
        static std::atomic<std::uint64_t> sequence{0};
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (unsigned attempt = 0; attempt < 32; ++attempt) {
            temporary.path = destination;
            temporary.path += ".paacs-tmp-" + std::to_string(tick) + "-" +
                              std::to_string(sequence.fetch_add(1, std::memory_order_relaxed));
            try {
                stream.open(temporary.path, std::ios::binary | std::ios::noreplace);
            } catch (...) {
                temporary.owned = stream.is_open();
                throw;
            }
            temporary.owned = stream.is_open();
            if (temporary.owned)
                break;
            std::error_code error;
            if (!std::filesystem::exists(temporary.path, error) || error)
                return failure(MeshArchiveReason::IoError);
            stream.clear();
        }
        if (!temporary.owned)
            return failure(MeshArchiveReason::IoError);
        stream.write(reinterpret_cast<const char *>(bytes.data()),
                     static_cast<std::streamsize>(bytes.size()));
        // close() reports delayed write errors before the destination is touched.
        stream.close();
        if (!stream)
            return failure(MeshArchiveReason::IoError);
#if defined(_WIN32)
        // std::filesystem::rename cannot replace an existing file on Windows.
        if (!MoveFileExW(temporary.path.c_str(), destination.c_str(), MOVEFILE_REPLACE_EXISTING))
            return failure(MeshArchiveReason::IoError);
#else
        std::error_code error;
        std::filesystem::rename(temporary.path, destination, error);
        if (error)
            return failure(MeshArchiveReason::IoError);
#endif
        temporary.owned = false;
        return out;
    } catch (const std::bad_alloc &) {
        return failure(MeshArchiveReason::AllocationFailed);
    } catch (...) {
        return failure(MeshArchiveReason::IoError);
    }
}
MeshArchiveOutcome readMeshArchive(const std::string &path, mesh::Mesh &value,
                                   const MeshArchiveLimits &limits) noexcept {
    try {
        std::ifstream stream(path, std::ios::binary | std::ios::ate);
        if (!stream)
            return failure(MeshArchiveReason::IoError);
        const auto size = stream.tellg();
        if (size < 0)
            return failure(MeshArchiveReason::IoError);
        const auto count = static_cast<std::uint64_t>(size);
        if (count > limits.maxBytes)
            return failure(MeshArchiveReason::ResourceLimit);
        if (count > std::numeric_limits<std::size_t>::max() ||
            count > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max()))
            return failure(MeshArchiveReason::CountOverflow);
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(count));
        stream.seekg(0);
        stream.read(reinterpret_cast<char *>(bytes.data()), static_cast<std::streamsize>(count));
        if (!stream)
            return failure(MeshArchiveReason::IoError);
        return deserializeMeshArchive(std::span<const std::uint8_t>(bytes), value, limits);
    } catch (const std::bad_alloc &) {
        return failure(MeshArchiveReason::AllocationFailed);
    } catch (...) {
        return failure(MeshArchiveReason::IoError);
    }
}
MeshArchiveOutcome readMeshArchive(const std::string &path, mesh::Mesh &value) noexcept {
    return readMeshArchive(path, value, MeshArchiveLimits{});
}
} // namespace paacs::sdk
