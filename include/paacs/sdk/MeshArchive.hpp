// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
#include "paacs/mesh/Mesh.hpp"
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace paacs::sdk {
enum class MeshArchiveReason : std::uint8_t {
    None, InvalidMesh, InvalidHeader, UnsupportedVersion, Truncated,
    CountOverflow, InvalidCell, IoError, AllocationFailed, InternalError,
    ResourceLimit
};
struct MeshArchiveOutcome {
    bool completed = false;
    MeshArchiveReason reason = MeshArchiveReason::InternalError;
    std::uint64_t bytes = 0;
};
// Explicit budgets for untrusted input. Raise them deliberately for large
// meshes. maxBytes bounds encoded bytes, not peak decoded memory consumption.
struct MeshArchiveLimits {
    std::uint64_t maxBytes = 1ULL << 30;
    std::uint64_t maxNodes = 50'000'000;
    std::uint64_t maxCells = 50'000'000;
    std::uint64_t maxConnectivity = 500'000'000;
};
MeshArchiveOutcome serializeMeshArchive(
    const mesh::Mesh& value, std::vector<std::uint8_t>& bytes) noexcept;
MeshArchiveOutcome deserializeMeshArchive(
    const std::vector<std::uint8_t>& bytes, mesh::Mesh& value) noexcept;
// On every failure value remains unchanged. No input ownership is transferred.
MeshArchiveOutcome deserializeMeshArchive(
    std::span<const std::uint8_t> bytes, mesh::Mesh& value,
    const MeshArchiveLimits& limits = {}) noexcept;
// Serializes privately, closes an exclusively created same-directory temporary,
// then replaces the destination by rename. Failure preserves an existing file;
// cleanup of owned temporary files is attempted. This is atomic publication,
// not a promise of crash-durable storage or preserving filesystem metadata.
MeshArchiveOutcome writeMeshArchive(
    const mesh::Mesh& value, const std::string& path) noexcept;
MeshArchiveOutcome readMeshArchive(
    const std::string& path, mesh::Mesh& value) noexcept;
MeshArchiveOutcome readMeshArchive(
    const std::string& path, mesh::Mesh& value,
    const MeshArchiveLimits& limits) noexcept;
} // namespace paacs::sdk
