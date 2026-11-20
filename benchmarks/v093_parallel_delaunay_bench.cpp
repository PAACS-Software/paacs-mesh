// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay3D.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using paacs::Index; using paacs::Vec3; using paacs::algo::Tet;
static double u01(std::mt19937& g){std::uint64_t h=g(),l=g();return double((h<<21)|(l>>11))*0x1.0p-53;}
static std::uint64_t canonicalHash(const std::vector<Tet>& tets){std::vector<std::array<Index,4>> keys;keys.reserve(tets.size());for(const auto&t:tets){auto k=t.v;std::sort(k.begin(),k.end());keys.push_back(k);}std::sort(keys.begin(),keys.end());std::uint64_t h=1469598103934665603ULL;for(const auto&k:keys)for(Index v:k)for(int b=0;b<8;++b){h^=(std::uint64_t(v)>>(8*b))&255U;h*=1099511628211ULL;}return h;}
int main(int ac,char**av){if(ac<5){std::cerr<<"usage: bench points serial|strict|unique threads packet [runs]\n";return 2;}const int n=std::stoi(av[1]);const std::string mode=av[2];const int threads=std::stoi(av[3]);const unsigned packet=unsigned(std::stoul(av[4]));const int runs=ac>5?std::stoi(av[5]):5;std::mt19937 g(20260818);std::vector<Vec3> p;p.reserve(n);for(int i=0;i<n;++i)p.push_back({u01(g),u01(g),u01(g)});paacs::algo::Delaunay3DOptions o;o.requestedThreads=threads;o.packetSize=packet;if(mode=="strict")o.parallelPolicy=paacs::algo::DelaunayParallelPolicy::StrictSerialEquivalent;else if(mode=="unique")o.parallelPolicy=paacs::algo::DelaunayParallelPolicy::UniqueDelaunayParallel;else o.parallelPolicy=paacs::algo::DelaunayParallelPolicy::SerialExact;paacs::algo::DelaunayWorkspace ws;std::vector<double> times;std::uint64_t h=0;std::size_t cells=0;paacs::algo::Delaunay3DStats last{};for(int r=-1;r<runs;++r){paacs::algo::Delaunay3DStats st;auto t=paacs::algo::delaunay3D(p,ws,o,&st);auto x=canonicalHash(t);if(r>=0){if(h&&h!=x){std::cerr<<"hash changed\n";return 3;}h=x;cells=t.size();times.push_back(st.totalSeconds);last=st;}}std::sort(times.begin(),times.end());const double med=times[times.size()/2];std::cout<<std::setprecision(12)<<"mode="<<mode<<"\npoints="<<n<<"\nworkers="<<threads<<"\npacket="<<packet<<"\ncells="<<cells<<"\nmedian_s="<<med<<"\ncertified_mtet_s="<<cells/med/1e6<<"\nhash=0x"<<std::hex<<h<<std::dec<<"\norder_s="<<last.normalizationOrderSeconds<<"\ninsertion_s="<<last.insertionSeconds<<"\naudit_s="<<last.adjacencyAuditSeconds<<"\nfinalization_s="<<last.finalizationSeconds<<"\nparallel_discovery_s="<<last.parallelDiscoverySeconds<<"\nparallel_commit_s="<<last.canonicalCommitSeconds<<"\nhull_prepass_s="<<last.parallelHullPrepassSeconds<<"\nrounds="<<last.parallelRounds<<"\ndiscoveries="<<last.speculativeDiscoveries<<"\nreused="<<last.speculativeCandidatesReused<<"\nrecomputed="<<last.speculativeCandidatesRecomputed<<"\nhazard_exact="<<last.parallelHazardPredicateTests<<"\nhazard_culls="<<last.parallelHazardBroadPhaseCulls<<"\npair_checks="<<last.parallelConflictPairChecks<<"\nwrite_conflicts="<<last.parallelWriteMarkConflicts<<"\ndegeneracy_fallback="<<(last.parallelDegeneracyFallback?1:0)<<"\narena_bytes="<<last.arenaBytes<<"\nscratch_bytes="<<last.scratchBytes<<"\nworkspace_bytes="<<last.workspaceRetainedBytes<<"\n";}
