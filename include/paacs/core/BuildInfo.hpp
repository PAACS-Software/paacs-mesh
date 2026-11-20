// PAACS Mesh
// Project author: Paul LeDewhurst
// Copyright PAACS Software Ltd.
// SPDX-License-Identifier: MPL-2.0

#pragma once
//
// paacs::core::BuildInfo - reproducibility metadata GENERATED from the
// actual build (directive 5.3), never written as static prose. Review
// finding: the baseline recorder hardcoded "g++-13.x -O2 Linux x86_64"
// and kept printing it under GCC 14 - recreating the record-vs-reality
// problem this project exists to prevent.
//
// CMake supplies PAACS_BUILD_TYPE / PAACS_BUILD_FLAGS / PAACS_BUILD_ID
// definitions when available; compiler/platform/arch always come from
// compiler macros so even a bare compiler invocation reports truthfully.
//
#include <cstdint>
#include <cstdio>
#include <string>

namespace paacs::core {

inline std::string buildCompiler() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#elif defined(_MSC_VER)
    return "msvc " + std::to_string(_MSC_VER);
#else
    return "unknown-compiler";
#endif
}

inline const char* buildPlatform() {
#if defined(__linux__)
    return "linux";
#elif defined(__APPLE__)
    return "macos";
#elif defined(_WIN32)
    return "windows";
#else
    return "unknown-platform";
#endif
}

inline const char* buildArch() {
#if defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    return "aarch64";
#else
    return "unknown-arch";
#endif
}

inline constexpr bool buildIsOptimized() {
#if defined(__OPTIMIZE__)
    return true;
#elif defined(_MSC_VER) && defined(NDEBUG)
    // MSVC has no __OPTIMIZE__; Release-family CMake configurations define
    // NDEBUG and carry their exact flags separately in evidence manifests.
    return true;
#else
    return false;
#endif
}

inline const char* buildType() {
#if defined(PAACS_BUILD_TYPE)
    if (PAACS_BUILD_TYPE[0] != '\0') return PAACS_BUILD_TYPE;
#endif
    return buildIsOptimized() ? "optimized (no CMake build type)"
                              : "unoptimized (no CMake build type)";
}

inline const char* buildFlags() {
#if defined(PAACS_BUILD_FLAGS)
    return PAACS_BUILD_FLAGS;
#else
    return "(flags not recorded: non-CMake build)";
#endif
}

inline const char* buildId() {
#if defined(PAACS_BUILD_ID)
    return PAACS_BUILD_ID;
#else
    return "unversioned-workspace";
#endif
}

inline const char* buildStdlib() {
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#elif defined(_MSC_VER)
    return "msvc-stl";
#else
    return "unknown-stdlib";
#endif
}

// RECORDED EFFECTIVE BUILD FINGERPRINT (round-3 hardening; renamed per
// round-4 review - this is the recorded effective configuration, NOT a
// mathematically complete toolchain description: compiler defaults,
// external toolchain files, and target options added outside the
// captured flag variables remain outside the contract): compiler
// + version, stdlib + release, platform, arch, and the FP-relevant
// effective configuration - optimization state, fast-math, salient FP
// codegen features (FMA/AVX2, which -march=native would enable),
// sanitizers - plus an FNV hash of the ACTUAL build flags when the
// build system records them (CMake always does; the canonical no-CMake
// script passes -DPAACS_BUILD_FLAGS). A gcc-13.3 build with altered
// floating-point or architecture flags therefore no longer classifies
// as the strict reference. Residual honesty note: without a recorded
// flag string, -O2 vs -O3 are not macro-distinguishable; such builds
// carry flagshash=none and are NEVER the reference environment.
inline std::string environmentFingerprint() {
    std::string fp = buildCompiler();
    fp += "|";
    fp += buildStdlib();
#if defined(_GLIBCXX_RELEASE)
    fp += "-" + std::to_string(_GLIBCXX_RELEASE);
#endif
    fp += "|";
    fp += buildPlatform();
    fp += "|";
    fp += buildArch();
    fp += "|";
#if defined(__OPTIMIZE__) || (defined(_MSC_VER) && defined(NDEBUG))
    fp += "optimized";
#else
    fp += "unoptimized";
#endif
#if defined(__FAST_MATH__)
    fp += "+fastmath";
#endif
#if defined(__FMA__)
    fp += "+fma";
#endif
#if defined(__AVX2__)
    fp += "+avx2";
#endif
#if defined(__SANITIZE_ADDRESS__)
    fp += "+asan";
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
    fp += "+asan";
#endif
#endif
#if defined(__SANITIZE_THREAD__)
    fp += "+tsan";
#elif defined(__has_feature)
#if __has_feature(thread_sanitizer)
    fp += "+tsan";
#endif
#endif
    fp += "|flags:";
#if defined(PAACS_BUILD_FLAGS)
    {
        std::uint64_t h = 1469598103934665603ULL;
        for (const char* p = PAACS_BUILD_FLAGS; *p; ++p) {
            h ^= static_cast<unsigned char>(*p);
            h *= 1099511628211ULL;
        }
        char buf[20];
        std::snprintf(buf, sizeof buf, "%016llx",
                      static_cast<unsigned long long>(h));
        fp += buf;
    }
#else
    fp += "none";
#endif
    return fp;
}

// Frozen with the 2026-07-04 baselines; a re-freeze on a new reference
// environment updates this literal in the SAME commit as the new strict
// values (docs/BASELINE_FREEZE.md rule 2).
inline constexpr const char* kStrictBaselineFingerprint =
    "gcc 13.3.0|libstdc++-13|linux|x86_64|optimized|flags:672a9d24f745392a";

inline bool isStrictBaselineReferenceEnv() {
    return environmentFingerprint() == kStrictBaselineFingerprint;
}

inline std::string buildSummary() {
    return environmentFingerprint() + " | " + buildType() +
           " | flags: " + buildFlags() + " | build-id: " + buildId();
}

} // namespace paacs::core
