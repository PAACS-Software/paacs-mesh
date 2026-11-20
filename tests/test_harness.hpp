// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

// PAACS Mesh test harness - shared by all test translation units.
#pragma once
#include <cstdio>
#include <random>

extern int gFailures;

// Project-defined canonical U(0,1): raw mt19937 draws combined into a
// 53-bit mantissa (frozen; identical across implementations).
inline double canonU01(std::mt19937& g) {
    const std::uint64_t hi = g(), lo = g();
    return double((hi << 21) | (lo >> 11)) * 0x1.0p-53;
}

#define CHECK(cond)                                                        \
    do {                                                                   \
        if (!(cond)) {                                                     \
            std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);  \
            ++gFailures;                                                   \
        }                                                                  \
    } while (0)

#define TEST(name) std::printf("[test] %s\n", name);
