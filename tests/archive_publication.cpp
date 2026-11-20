// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshArchive.hpp"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

namespace {
long allocationRemaining = -1;
}
void *operator new(std::size_t size) {
    if (allocationRemaining == 0)
        throw std::bad_alloc();
    if (allocationRemaining > 0)
        --allocationRemaining;
    if (void *value = std::malloc(size ? size : 1))
        return value;
    throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *value) noexcept { std::free(value); }
void operator delete[](void *value) noexcept { std::free(value); }
void operator delete(void *value, std::size_t) noexcept { std::free(value); }
void operator delete[](void *value, std::size_t) noexcept { std::free(value); }

namespace {
void require(bool good, const char *message) {
    if (!good)
        throw std::runtime_error(message);
}
std::vector<char> read(const std::filesystem::path &path) {
    std::ifstream stream(path, std::ios::binary);
    require(bool(stream), "read held archive");
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}
void noTemporaryFiles(const std::filesystem::path &directory) {
    for (const auto &entry : std::filesystem::recursive_directory_iterator(directory))
        require(entry.path().filename().string().find(".paacs-tmp-") == std::string::npos,
                "temporary archive leaked");
}
struct Directory {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("paacs-archive-publication-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    Directory() { require(std::filesystem::create_directory(path), "exclusive test directory"); }
    ~Directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};
#if defined(_WIN32)
struct FileLock {
    HANDLE value = INVALID_HANDLE_VALUE;
    explicit FileLock(const std::filesystem::path &path)
        : value(CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr)) {
        require(value != INVALID_HANDLE_VALUE, "open locked archive");
    }
    ~FileLock() { CloseHandle(value); }
};
#endif
} // namespace

int main() {
    try {
        Directory directory;
        const auto destination = directory.path / "held.paacsmesh";
        paacs::mesh::Mesh first;
        first.addNode({1, 2, 3}, {0, 1});
        require(paacs::sdk::writeMeshArchive(first, destination.string()).completed,
                "initial archive publication");
        const auto initialBytes = read(destination);
        paacs::mesh::Mesh second;
        second.addNode({4, 5, 6}, {0, 2});
        second.addNode({7, 8, 9}, {0, 3});
        require(paacs::sdk::writeMeshArchive(second, destination.string()).completed,
                "replace existing archive");
        const auto replacedBytes = read(destination);
        require(replacedBytes != initialBytes, "archive not replaced");
        paacs::mesh::Mesh loaded;
        require(paacs::sdk::readMeshArchive(destination.string(), loaded).completed &&
                    loaded.nodeCount() == 2 && loaded.node(0).x == 4 &&
                    loaded.nodeClass(1).tag == 3,
                "published replacement readback");
        noTemporaryFiles(directory.path);

        paacs::mesh::Mesh invalid;
        invalid.addNode({std::numeric_limits<double>::quiet_NaN(), 0, 0}, {0, 1});
        const auto invalidResult = paacs::sdk::writeMeshArchive(invalid, destination.string());
        require(!invalidResult.completed &&
                    invalidResult.reason == paacs::sdk::MeshArchiveReason::InvalidMesh &&
                    read(destination) == replacedBytes,
                "invalid mesh destroyed existing archive");

        const auto missing = directory.path / "missing" / "archive.paacsmesh";
        const auto missingResult = paacs::sdk::writeMeshArchive(first, missing.string());
        require(!missingResult.completed &&
                    missingResult.reason == paacs::sdk::MeshArchiveReason::IoError &&
                    !std::filesystem::exists(missing.parent_path()),
                "missing-parent failure");
        const auto folder = directory.path / "folder.paacsmesh";
        std::filesystem::create_directory(folder);
        std::ofstream(folder / "sentinel") << "original directory contents";
        const auto folderResult = paacs::sdk::writeMeshArchive(first, folder.string());
        require(!folderResult.completed &&
                    folderResult.reason == paacs::sdk::MeshArchiveReason::IoError &&
                    std::filesystem::is_directory(folder) &&
                    std::filesystem::exists(folder / "sentinel"),
                "directory destination changed on rename failure");
        noTemporaryFiles(directory.path);
#if defined(_WIN32)
        {
            FileLock lock(destination);
            const auto locked = paacs::sdk::writeMeshArchive(first, destination.string());
            require(!locked.completed && locked.reason == paacs::sdk::MeshArchiveReason::IoError,
                    "locked archive replacement should fail");
            require(read(destination) == replacedBytes, "locked original archive changed");
            noTemporaryFiles(directory.path);
        }
        require(paacs::sdk::writeMeshArchive(first, destination.string()).completed,
                "retry after lock release");
        require(read(destination) == initialBytes, "successful retry has incorrect bytes");
#endif
        noTemporaryFiles(directory.path);
        // Sweep the actual writer until the first non-failing allocation
        // ordinal. The file and temporary-directory invariants are observed
        // independently after every injected failure, with injection disabled.
        require(paacs::sdk::writeMeshArchive(second, destination.string()).completed,
                "prepare allocation-fault destination");
        const auto destinationText = destination.string();
        unsigned allocationFailures = 0;
        bool reachedSuccess = false;
        for (long ordinal = 0; ordinal < 256; ++ordinal) {
            allocationRemaining = ordinal;
            const auto outcome = paacs::sdk::writeMeshArchive(first, destinationText);
            allocationRemaining = -1;
            noTemporaryFiles(directory.path);
            if (outcome.completed) {
                require(read(destination) == initialBytes, "fault-sweep success has wrong bytes");
                reachedSuccess = true;
                break;
            }
            require(outcome.reason == paacs::sdk::MeshArchiveReason::AllocationFailed ||
                        outcome.reason == paacs::sdk::MeshArchiveReason::IoError,
                    "allocation-fault refusal lost failure taxonomy");
            require(read(destination) == replacedBytes, "allocation failure changed destination");
            ++allocationFailures;
        }
        require(reachedSuccess && allocationFailures > 0, "writer allocation sweep did not complete");
        std::cout << "Archive publication: replacement/readback and failure preservation passed\n";
        std::cout << "Writer allocation failures preserved destination and cleaned temporaries: "
                  << allocationFailures << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
