// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/algo/Delaunay3D.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
using namespace paacs; using namespace paacs::algo;
static double u01(std::mt19937& g){std::uint64_t h=g(),l=g();return double((h<<21)|(l>>11))*0x1.0p-53;}
static std::vector<std::array<Index,4>> canon(const std::vector<Tet>& ts){std::vector<std::array<Index,4>> x;x.reserve(ts.size());for(const auto&t:ts){auto a=t.v;std::sort(a.begin(),a.end());x.push_back(a);}std::sort(x.begin(),x.end());return x;}
static bool sameVector(const std::vector<Tet>&a,const std::vector<Tet>&b){return a.size()==b.size()&&std::equal(a.begin(),a.end(),b.begin(),[](const Tet&x,const Tet&y){return x.v==y.v;});}
int main(){std::mt19937 g(777);std::vector<Vec3> p;for(int i=0;i<300;++i)p.push_back({u01(g),u01(g),u01(g)});const auto serial=delaunay3D(p);const auto reference=canon(serial);for(int th:{1,2,4,5})for(unsigned packet:{32u,64u,128u,256u,512u}){Delaunay3DOptions o;o.requestedThreads=th;o.packetSize=packet;o.parallelPolicy=DelaunayParallelPolicy::UniqueDelaunayParallel;Delaunay3DStats st;const auto t=delaunay3D(p,o,&st);if(canon(t)!=reference||st.parallelDegeneracyFallback){std::cerr<<"unique failure "<<th<<"/"<<packet<<"\n";return 2;}}for(int th:{1,5})for(unsigned packet:{32u,128u,512u}){Delaunay3DOptions o;o.requestedThreads=th;o.packetSize=packet;o.parallelPolicy=DelaunayParallelPolicy::StrictSerialEquivalent;const auto t=delaunay3D(p,o);if(!sameVector(t,serial)){std::cerr<<"strict failure "<<th<<"/"<<packet<<"\n";return 3;}}std::vector<Vec3> grid;for(int k=0;k<5;++k)for(int j=0;j<5;++j)for(int i=0;i<5;++i)grid.push_back({double(i),double(j),double(k)});const auto gs=delaunay3D(grid);Delaunay3DOptions o;o.requestedThreads=5;o.packetSize=128;o.parallelPolicy=DelaunayParallelPolicy::UniqueDelaunayParallel;Delaunay3DStats st;const auto gp=delaunay3D(grid,o,&st);if(!sameVector(gp,gs)||!st.parallelDegeneracyFallback){std::cerr<<"degenerate fallback failure\n";return 4;}std::cout<<"v093_parallel_probe=pass\n";}
