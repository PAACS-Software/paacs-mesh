// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#include "paacs/core/Parallel.hpp"
#include <atomic>
#include <iostream>
#include <stdexcept>
#include <thread>
int main() {
    using namespace paacs::core;
    std::atomic<unsigned long long> sum{0};
    std::jthread settings([] {for(int i=0;i<20000;++i)setThreadCount(i%5);});
    std::jthread reader([] {for(int i=0;i<20000;++i)(void)effectiveThreads();});
    for(int round=0;round<50;++round) {
        parallelForWithThreads(0,4096,4,[&](std::size_t i) {
            parallelForWithThreads(0,2,4,[&](auto){sum.fetch_add(i,std::memory_order_relaxed);});
        });
        bool caught=false;
        try {parallelForWithThreads(0,4096,4,[](auto){throw std::runtime_error("expected");});}
        catch(const std::runtime_error&){caught=true;}
        if(!caught || inPoolContext())return 1;
    }
    if(sum.load()!=50ULL*4096*4095)return 1;
    std::cout<<"parallel stress: PASS\n";
}
