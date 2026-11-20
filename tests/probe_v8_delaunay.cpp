// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay3D.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
using paacs::Index;
using paacs::Vec3;
using paacs::algo::Tet;

double canonU01(std::mt19937& generator) {
    const std::uint64_t hi = generator();
    const std::uint64_t lo = generator();
    return static_cast<double>((hi << 21) | (lo >> 11)) * 0x1.0p-53;
}

std::uint64_t canonicalTetHash(const std::vector<Tet>& tets) {
    std::vector<std::array<Index,4>> keys;
    keys.reserve(tets.size());
    for (const auto& tet : tets) {
        auto key = tet.v;
        std::sort(key.begin(), key.end());
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    std::uint64_t hash = 1469598103934665603ULL;
    for (const auto& key : keys)
        for (Index value : key)
            for (int byte = 0; byte < 8; ++byte) {
                hash ^= (static_cast<std::uint64_t>(value) >> (8 * byte)) & 0xFF;
                hash *= 1099511628211ULL;
            }
    return hash;
}

void requireGolden() {
    {
        std::vector<Vec3> points;
        for (int k=0;k<5;++k)
            for (int j=0;j<5;++j)
                for (int i=0;i<5;++i)
                    points.push_back({double(i),double(j),double(k)});
        const auto tets=paacs::algo::delaunay3D(points);
        if(tets.size()!=384||canonicalTetHash(tets)!=0xea7550e52137ee03ULL)
            throw std::runtime_error("structured golden mismatch");
    }
    {
        std::mt19937 generator(777);
        std::vector<Vec3> points;
        for(int i=0;i<300;++i)
            points.push_back({canonU01(generator),canonU01(generator),canonU01(generator)});
        const auto tets=paacs::algo::delaunay3D(points);
        if(tets.size()!=1777||canonicalTetHash(tets)!=0x4720fb7e5d0738ccULL)
            throw std::runtime_error("random golden mismatch");
    }
    {
        std::mt19937 generator(778);
        std::vector<Vec3> points;
        for(int i=0;i<200;++i)
            points.push_back({canonU01(generator),canonU01(generator),canonU01(generator)*1e-9});
        const auto tets=paacs::algo::delaunay3D(points);
        if(tets.size()!=1109||canonicalTetHash(tets)!=0xccf396a3b736f2daULL)
            throw std::runtime_error("near-coplanar golden mismatch");
    }
}

}  // namespace

int main(int argc,char**argv){
    try {
        requireGolden();
        const int pointCount=argc>1?std::stoi(argv[1]):50000;
        const int runs=argc>2?std::stoi(argv[2]):7;
        std::mt19937 generator(20260818);
        std::vector<Vec3> points;
        points.reserve(static_cast<std::size_t>(pointCount));
        for(int i=0;i<pointCount;++i)
            points.push_back({canonU01(generator),canonU01(generator),canonU01(generator)});
        std::vector<double> times;
        times.reserve(static_cast<std::size_t>(runs));
        std::uint64_t hash=0,cells=0;
        paacs::algo::Delaunay3DStats retained{};
        for(int run=-1;run<runs;++run){
            paacs::algo::Delaunay3DStats stats{};
            const auto start=std::chrono::steady_clock::now();
            const auto tets=paacs::algo::delaunay3D(points,&stats);
            const auto stop=std::chrono::steady_clock::now();
            if(!stats.runComplete||!stats.auditPassed||stats.pointsSkipped!=0)
                throw std::runtime_error("Delaunay audit incomplete");
            const auto currentHash=canonicalTetHash(tets);
            if(run>=0){
                times.push_back(std::chrono::duration<double>(stop-start).count());
                if(hash&&hash!=currentHash)throw std::runtime_error("nondeterministic benchmark hash");
                hash=currentHash;cells=tets.size();retained=stats;
            }
        }
        std::sort(times.begin(),times.end());
        const double median=times[times.size()/2];
        std::cout<<"goldens=pass\n"
                 <<"points="<<pointCount<<"\n"
                 <<"cells="<<cells<<"\n"
                 <<"median_s="<<std::setprecision(12)<<median<<"\n"
                 <<"million_tets_per_s="<<cells/median/1.0e6<<"\n"
                 <<"canonical_hash=0x"<<std::hex<<hash<<std::dec<<"\n"
                 <<"peak_live_cells="<<retained.arenaPeakLiveCells<<"\n"
                 <<"arena_bytes="<<retained.arenaBytes<<"\n";
#ifdef PAACS_V8_PAIR_STATS
        std::cout<<"shell_pair_insertions="<<retained.shellEdgePairInsertions<<"\n"
                 <<"shell_pair_probes="<<retained.shellEdgePairProbes<<"\n"
                 <<"peak_pair_slots="<<retained.peakShellPairTableSlots<<"\n";
#endif
        return 0;
    }catch(const std::exception&error){
        std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;
    }
}
