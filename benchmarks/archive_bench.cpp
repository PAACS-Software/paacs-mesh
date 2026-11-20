// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/sdk/MeshArchive.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(int argc, char **argv) {
    using namespace paacs;
    const int count = argc > 1 ? std::stoi(argv[1]) : 200000;
    if (count < 1 || count > 10000000)
        return 2;
    mesh::Mesh m;
    (void)m.addNode({0, 0, 0}, {0, 1});
    (void)m.addNode({1, 0, 0}, {0, 2});
    (void)m.addNode({0, 1, 0}, {0, 3});
    (void)m.addNode({0, 0, 1}, {0, 4});
    m.reserveAppend(0, count, static_cast<std::size_t>(count) * 4);
    for (int i = 0; i < count; ++i)
        (void)m.addCell(mesh::CellType::Tet4, {0, 1, 2, 3}, {3, 1});
    std::vector<double> enc, dec;
    std::vector<std::uint8_t> bytes;
    for (int i = -1; i < 7; ++i) {
        const auto a = std::chrono::steady_clock::now();
        if (!sdk::serializeMeshArchive(m, bytes).completed)
            return 1;
        const auto b = std::chrono::steady_clock::now();
        mesh::Mesh decoded;
        if (!sdk::deserializeMeshArchive(bytes, decoded).completed || decoded.cellCount() != count)
            return 1;
        const auto c = std::chrono::steady_clock::now();
        if (i >= 0) {
            enc.push_back(std::chrono::duration<double>(b - a).count());
            dec.push_back(std::chrono::duration<double>(c - b).count());
        }
    }
    std::sort(enc.begin(), enc.end());
    std::sort(dec.begin(), dec.end());
    std::uint64_t hash = 1469598103934665603ULL;
    for (auto b : bytes) {
        hash ^= b;
        hash *= 1099511628211ULL;
    }
    std::cout << "case=synthetic_repeated_tet_connectivity\ncells=" << count
              << "\nbytes=" << bytes.size() << "\nencode_median_s=" << enc[enc.size() / 2]
              << "\ndecode_median_s=" << dec[dec.size() / 2] << "\narchive_hash=" << hash << '\n';
}
