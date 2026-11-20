// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "audit_reference.hpp"
#include "paacs/algo/MappedHexMesher.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
using namespace paacs;
int main(int argc, char **argv) {
    try {
        const unsigned n = argc > 1 ? static_cast<unsigned>(std::stoul(argv[1])) : 63;
        const std::string method = argc > 2 ? argv[2] : "bucket";
        const bool shuffle = argc > 3 && std::string(argv[3]) == "shuffle";
        if (!n || n > 200 || (method != "bucket" && method != "reference"))
            return 2;
        algo::MappedHexBlockRequest request;
        request.geometry.cellsU = request.geometry.cellsV = request.geometry.cellsW = n;
        request.executionThreads = 1;
        mesh::Mesh m;
        if (!algo::generateMappedHexBlock(request, m).published())
            return 3;
        if (shuffle) {
            mesh::Mesh other;
            other.reserveAppend(static_cast<std::size_t>(m.nodeCount()),
                                static_cast<std::size_t>(m.cellCount()),
                                static_cast<std::size_t>(m.connectivityCount()));
            std::vector<Index> ids(static_cast<std::size_t>(m.nodeCount()));
            std::iota(ids.begin(), ids.end(), 0);
            std::mt19937_64 rng(42);
            std::shuffle(ids.begin(), ids.end(), rng);
            std::vector<Index> inverse(ids.size());
            for (Index i = 0; i < m.nodeCount(); ++i)
                inverse[static_cast<std::size_t>(ids[static_cast<std::size_t>(i)])] = i;
            for (Index i = 0; i < m.nodeCount(); ++i) {
                const auto old = inverse[static_cast<std::size_t>(i)];
                (void)other.addNode(m.node(old), m.nodeClass(old));
            }
            for (Index i = 0; i < m.cellCount(); ++i) {
                const auto &c = m.cell(i);
                std::vector<Index> conn;
                for (int k = 0; k < mesh::nodesPerCell(c.type); ++k)
                    conn.push_back(ids[static_cast<std::size_t>(m.cellNodes(i)[k])]);
                (void)other.addCell(c.type, conn, c.onEntity);
            }
            m.swap(other);
        }
        const auto begin = std::chrono::steady_clock::now();
        const auto a = method == "reference" ? mesh::referenceLinearVolumeBoundary(m)
                                             : mesh::auditLinearVolumeBoundary(m);
        const double time =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        if (!a.valid)
            return 4;
        std::cout << std::setprecision(17) << "{\"method\":\"" << method << "\",\"n\":" << n
                  << ",\"shuffled\":" << (shuffle ? "true" : "false")
                  << ",\"cells\":" << a.volumeCells << ",\"face_uses\":" << a.faceUses
                  << ",\"seconds\":" << time << ",\"temporary_bytes\":"
                  << a.temporaryBytes
                  << ",\"maximum_bucket\":" << a.maximumBucketSize
                  << ",\"exterior_faces\":" << a.exteriorFaces
                  << ",\"interior_faces\":" << a.interiorFaces << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
