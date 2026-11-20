# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30...4.4)
project(PAACSMesh VERSION 1.2.0 LANGUAGES C CXX)

set(PAACS_RELEASE_VERSION "1.2.0-rc.1")
set(PAACS_CXX_STANDARD 26 CACHE STRING "C++ language level: 26 or compatibility mode 23")
set_property(CACHE PAACS_CXX_STANDARD PROPERTY STRINGS 26 23)
if(NOT PAACS_CXX_STANDARD MATCHES "^(23|26)$")
    message(FATAL_ERROR "PAACS_CXX_STANDARD must be 23 or 26")
endif()
set(CMAKE_CXX_STANDARD ${PAACS_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# This distribution uses headers; do not pay for experimental module scans.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(PAACS_REQUIRE_GCC16 "Require GCC 16 or newer for the reference build" OFF)
if(PAACS_REQUIRE_GCC16 AND (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16))
    message(FATAL_ERROR "This preset requires GCC >= 16; select the compatibility preset for older compilers")
endif()
option(PAACS_BUILD_CLI "Build paacs-mesh" ON)
option(PAACS_BUILD_UI_BRIDGE "Build optional C ABI shared library for UI bindings" OFF)
option(PAACS_BUILD_TESTS "Build all unit tests and retained probes" ${PROJECT_IS_TOP_LEVEL})
option(PAACS_BUILD_EXTENDED_PROBES "Build retained topology allocation/rollback probes" OFF)
option(PAACS_BUILD_BENCHMARKS "Build reproducible benchmarks" OFF)
option(PAACS_WITH_OCCT "Build optional Open CASCADE 8+ native and planar STEP/BREP adapters" OFF)
option(PAACS_WARNINGS_AS_ERRORS "Treat library warnings as errors" OFF)
option(PAACS_ENABLE_IPO "Enable checked interprocedural optimization" OFF)
option(PAACS_NATIVE_ARCH "Optimize for this host (nonportable binaries)" OFF)
set(PAACS_SANITIZER "none" CACHE STRING "none, address-undefined, or thread")
set_property(CACHE PAACS_SANITIZER PROPERTY STRINGS none address-undefined thread)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
find_package(Threads REQUIRED)
include(cmake/Sources.cmake)
add_library(paacs_mesh STATIC ${PAACS_SOURCES})
add_library(PAACS::mesh ALIAS paacs_mesh)
target_compile_features(paacs_mesh PUBLIC cxx_std_${PAACS_CXX_STANDARD})
target_include_directories(paacs_mesh PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(paacs_mesh PUBLIC Threads::Threads)
# Some Apple SDKs gate the thread cancellation types behind a compiler flag.
# Carry the flag through the SDK target because public headers use stop_token.
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    include(CheckCXXSourceCompiles)
    include(CMakePushCheckState)
    set(paacs_thread_probe [=[
        #include <stop_token>
        #include <thread>
        int main() {
            std::stop_source stop;
            std::stop_callback callback(stop.get_token(), [] {});
            std::jthread worker([] {});
            return stop.request_stop() ? 0 : 1;
        }
    ]=])
    cmake_push_check_state()
    list(APPEND CMAKE_REQUIRED_LIBRARIES Threads::Threads)
    check_cxx_source_compiles("${paacs_thread_probe}" PAACS_APPLE_STD_STOP_TOKEN)
    if(NOT PAACS_APPLE_STD_STOP_TOKEN)
        string(APPEND CMAKE_REQUIRED_FLAGS " -fexperimental-library")
        list(APPEND CMAKE_REQUIRED_LINK_OPTIONS -fexperimental-library)
        check_cxx_source_compiles("${paacs_thread_probe}" PAACS_APPLE_EXPERIMENTAL_STOP_TOKEN)
        if(NOT PAACS_APPLE_EXPERIMENTAL_STOP_TOKEN)
            message(FATAL_ERROR "This Apple Clang toolchain cannot use std::stop_token and std::jthread. Use a newer Xcode SDK.")
        endif()
        target_compile_options(paacs_mesh PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-fexperimental-library>")
        target_link_options(paacs_mesh PUBLIC "$<$<LINK_LANG_AND_ID:CXX,AppleClang>:-fexperimental-library>")
    endif()
    cmake_pop_check_state()
endif()
set_target_properties(paacs_mesh PROPERTIES EXPORT_NAME mesh POSITION_INDEPENDENT_CODE ON)
# Exact floating-point filters assume IEEE semantics and unfused operations.
# Apply to consumers too: geometry and topology routines are also inline.
if(MSVC)
    target_compile_options(paacs_mesh PUBLIC /fp:strict PRIVATE /W4)
else()
    target_compile_options(paacs_mesh PUBLIC -fno-fast-math -ffp-contract=off
        PRIVATE -Wall -Wextra -Wpedantic)
endif()
set_property(TARGET paacs_mesh PROPERTY COMPILE_WARNING_AS_ERROR ${PAACS_WARNINGS_AS_ERRORS})
if(PAACS_NATIVE_ARCH)
    if(MSVC)
        message(FATAL_ERROR "PAACS_NATIVE_ARCH is only supported with GNU/Clang")
    endif()
    add_compile_options(-march=native)
    target_compile_options(paacs_mesh PRIVATE -march=native)
endif()
if(PAACS_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT PAACS_IPO_SUPPORTED OUTPUT PAACS_IPO_ERROR)
    if(NOT PAACS_IPO_SUPPORTED)
        message(FATAL_ERROR "IPO unavailable: ${PAACS_IPO_ERROR}")
    endif()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    set_property(TARGET paacs_mesh PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
endif()
if(NOT PAACS_SANITIZER STREQUAL "none")
    if(MSVC)
        message(FATAL_ERROR "Use GNU or Clang for sanitizer presets")
    endif()
    if(PAACS_SANITIZER STREQUAL "address-undefined")
        set(PAACS_SANITIZER_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer)
    elseif(PAACS_SANITIZER STREQUAL "thread")
        set(PAACS_SANITIZER_FLAGS -fsanitize=thread -fno-omit-frame-pointer)
    else()
        message(FATAL_ERROR "Unknown PAACS_SANITIZER=${PAACS_SANITIZER}")
    endif()
    target_compile_options(paacs_mesh PUBLIC ${PAACS_SANITIZER_FLAGS})
    target_link_options(paacs_mesh PUBLIC ${PAACS_SANITIZER_FLAGS})
endif()
# Metadata is private to this build; installed consumers retain their own flags.
target_compile_definitions(paacs_mesh PRIVATE
    PAACS_BUILD_TYPE="$<CONFIG>" PAACS_BUILD_ID="${PAACS_RELEASE_VERSION}")
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build-info-$<CONFIG>.txt" CONTENT
    "release=${PAACS_RELEASE_VERSION}\ncompiler=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\ncmake=${CMAKE_VERSION}\nstandard=${PAACS_CXX_STANDARD}\nconfiguration=$<CONFIG>\nflags=${CMAKE_CXX_FLAGS}\nstrict_fp=ON\nipo=${PAACS_ENABLE_IPO}\nnative=${PAACS_NATIVE_ARCH}\nsanitizer=${PAACS_SANITIZER}\nSee compile_commands.json for exact translation-unit flags.\n")
if(PAACS_WITH_OCCT)
    # OCCT's shipped version file uses exact-version compatibility even for
    # patch releases, so check the supported major after loading its config.
    find_package(OpenCASCADE REQUIRED CONFIG)
    if(OpenCASCADE_MAJOR_VERSION LESS 8)
        message(FATAL_ERROR "PAACS_WITH_OCCT requires Open CASCADE 8 or newer")
    endif()
    add_library(paacs_occt STATIC src/geometry/OcctImport.cpp
        src/geometry/OcctModel.cpp src/geometry/OcctMesher.cpp src/geometry/OcctBoundary.cpp src/geometry/OcctStarMesher.cpp
        src/geometry/OcctPartition.cpp src/geometry/OcctTrimTolerance.cpp src/geometry/OcctPartitionMesher.cpp src/geometry/OcctContactCylinder.cpp
        src/geometry/OcctBoxContact.cpp src/geometry/OcctNativeValidation.cpp src/geometry/OcctQueryContext.cpp
        src/geometry/OcctContactAssembly.cpp src/geometry/OcctContactComponents.cpp src/geometry/OcctNativeSweep.cpp
        src/geometry/OcctPairedShell.cpp src/geometry/OcctRecoveryMesher.cpp src/geometry/OcctPlcFrontRecovery.cpp src/geometry/OcctFrontMesher.cpp
        src/geometry/OcctWireRefinement.cpp src/geometry/OcctFacetRefinement.cpp)
    add_library(PAACS::occt ALIAS paacs_occt)
    # Target-name compatibility only; the older V1 CAD APIs are not restored.
    add_library(PAACS::MeshOCCT ALIAS paacs_occt)
    target_link_libraries(paacs_occt PUBLIC PAACS::mesh PRIVATE
        TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKPrim TKBO TKBool TKMesh TKXSBase TKDESTEP)
    target_include_directories(paacs_occt SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
    # Some OCCT binary SDKs export their own PCHs. Compile the adapter's explicit
    # includes with the same floating-point policy as PAACS.
    set_target_properties(paacs_occt PROPERTIES EXPORT_NAME occt POSITION_INDEPENDENT_CODE ON
        DISABLE_PRECOMPILE_HEADERS ON)
    if(MSVC)
        target_compile_options(paacs_occt PRIVATE /W4)
    else()
        target_compile_options(paacs_occt PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    set_property(TARGET paacs_occt PROPERTY COMPILE_WARNING_AS_ERROR ${PAACS_WARNINGS_AS_ERRORS})
    # Keep optional toolkit references out of the portable core export.
    install(TARGETS paacs_occt EXPORT PAACSOCCTTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
if(PAACS_BUILD_CLI)
    add_executable(paacs-mesh apps/paacs_mesh.cpp)
    target_link_libraries(paacs-mesh PRIVATE PAACS::mesh)
    add_executable(paacs-primitive-mesh apps/paacs_primitive_mesh.cpp)
    target_link_libraries(paacs-primitive-mesh PRIVATE PAACS::mesh)
    install(TARGETS paacs-primitive-mesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    if(PAACS_WITH_OCCT)
        add_executable(paacs-cad-mesh apps/paacs_cad_mesh.cpp apps/native_cad_cli.cpp)
        target_link_libraries(paacs-cad-mesh PRIVATE PAACS::occt)
        install(TARGETS paacs-cad-mesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
        add_executable(paacs-native-cad-preview apps/native_cad_preview.cpp)
        target_include_directories(paacs-native-cad-preview SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
        target_include_directories(paacs-native-cad-preview PRIVATE src/geometry)
        target_link_libraries(paacs-native-cad-preview PRIVATE PAACS::occt TKMesh)
        set_target_properties(paacs-native-cad-preview PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)
        install(TARGETS paacs-native-cad-preview RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
endif()
if(PAACS_BUILD_UI_BRIDGE)
    add_library(paacs_ui SHARED src/sdk/UIBridge.cpp)
    add_library(PAACS::ui ALIAS paacs_ui)
    target_compile_definitions(paacs_ui PRIVATE PAACS_UI_BUILD=1)
    target_link_libraries(paacs_ui PRIVATE PAACS::mesh)
    # C consumers also need the sanitizer runtime before this shared adapter.
    target_compile_options(paacs_ui INTERFACE ${PAACS_SANITIZER_FLAGS})
    target_link_options(paacs_ui INTERFACE ${PAACS_SANITIZER_FLAGS})
    target_include_directories(paacs_ui PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    set_target_properties(paacs_ui PROPERTIES EXPORT_NAME ui VERSION ${PROJECT_VERSION} SOVERSION 1
        CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
    if(NOT MSVC)
        target_compile_options(paacs_ui PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    set_property(TARGET paacs_ui PROPERTY COMPILE_WARNING_AS_ERROR ${PAACS_WARNINGS_AS_ERRORS})
    install(TARGETS paacs_ui EXPORT PAACSMeshTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
if(PAACS_BUILD_TESTS)
    enable_testing()
    # Fault-injection hooks must not enter the shipped library. The inherited
    # white-box suite needs a consistently instrumented copy of every TU.
    add_library(paacs_mesh_testing STATIC EXCLUDE_FROM_ALL ${PAACS_SOURCES})
    target_compile_features(paacs_mesh_testing PUBLIC cxx_std_${PAACS_CXX_STANDARD})
    target_include_directories(paacs_mesh_testing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(paacs_mesh_testing PUBLIC Threads::Threads)
    target_compile_definitions(paacs_mesh_testing PUBLIC PAACS_TOPO_TESTING=1
        PRIVATE $<TARGET_PROPERTY:paacs_mesh,COMPILE_DEFINITIONS>)
    target_compile_options(paacs_mesh_testing PUBLIC $<TARGET_PROPERTY:paacs_mesh,INTERFACE_COMPILE_OPTIONS>)
    target_link_options(paacs_mesh_testing PUBLIC $<TARGET_PROPERTY:paacs_mesh,INTERFACE_LINK_OPTIONS>)
    add_subdirectory(tests)
endif()
if(PAACS_BUILD_BENCHMARKS)
    add_executable(paacs-bench-ui-jobs benchmarks/ui_jobs.cpp)
    target_link_libraries(paacs-bench-ui-jobs PRIVATE PAACS::mesh)
    foreach(pair IN ITEMS "delaunay;tests/probe_v8_delaunay.cpp" "mapped-hex;benchmarks/v09_50m_bench.cpp" "pre-occt;benchmarks/v010_pre_occt_bench.cpp" "archive;benchmarks/archive_bench.cpp")
        list(GET pair 0 name)
        list(GET pair 1 source)
        add_executable(paacs-bench-${name} ${source})
        target_link_libraries(paacs-bench-${name} PRIVATE PAACS::mesh)
    endforeach()
endif()
install(TARGETS paacs_mesh EXPORT PAACSMeshTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(PAACS_BUILD_CLI)
    install(TARGETS paacs-mesh RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES LICENSE NOTICE LICENSING_STATUS.md CITATION.cff
    DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/${PROJECT_NAME})
configure_package_config_file(cmake/PAACSMeshConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PAACSMeshConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PAACSMesh)
# Pin the source SDK version; installed C++ clients use a matching build/toolchain.
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/PAACSMeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(EXPORT PAACSMeshTargets NAMESPACE PAACS:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PAACSMesh)
if(PAACS_WITH_OCCT)
    install(EXPORT PAACSOCCTTargets NAMESPACE PAACS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PAACSMesh)
endif()
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PAACSMeshConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PAACSMeshConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PAACSMesh)
set(CPACK_PACKAGE_VERSION ${PAACS_RELEASE_VERSION})
set(CPACK_SOURCE_GENERATOR ZIP)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "paacs-mesh-${PAACS_RELEASE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/out/;/evidence/;/\\.git/;/__pycache__/;CMakeUserPresets.json;\\.zip$;\\.paacsmesh$")
include(CPack)
