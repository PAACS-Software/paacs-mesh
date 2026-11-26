# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30...4.4)
project(PAACSThinHollowFixtures LANGUAGES CXX)
find_package(OpenCASCADE CONFIG REQUIRED)
if(OpenCASCADE_MAJOR_VERSION LESS 8)
    message(FATAL_ERROR "Thin hollow fixtures require Open CASCADE 8 or newer")
endif()
add_executable(paacs-make-thin-hollow make_fixtures.cpp)
target_compile_features(paacs-make-thin-hollow PRIVATE cxx_std_20)
set_target_properties(paacs-make-thin-hollow PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(paacs-make-thin-hollow SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(paacs-make-thin-hollow PRIVATE
    TKPrim TKBO TKBool TKBRep TKTopAlgo TKGeomAlgo TKGeomBase TKG3d
    TKG2d TKMath TKernel TKDESTEP TKXSBase)
if(MSVC)
    # The selected OCCT SDK supplies this tool's imports and runtime libraries.
    # Global vcpkg integration must not copy a different OCCT beside the exe.
    set_target_properties(paacs-make-thin-hollow PROPERTIES
        VS_GLOBAL_VcpkgEnabled false
        VS_GLOBAL_VcpkgApplocalDeps false)
    target_compile_options(paacs-make-thin-hollow PRIVATE /W4 /fp:strict)
else()
    target_compile_options(paacs-make-thin-hollow PRIVATE
        -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
endif()

# An explicit benchmark prototype, separate from ordinary import and meshing.
add_executable(paacs-prepare-thin-hollow prepare_fixture.cpp)
target_compile_features(paacs-prepare-thin-hollow PRIVATE cxx_std_20)
set_target_properties(paacs-prepare-thin-hollow PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(paacs-prepare-thin-hollow SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(paacs-prepare-thin-hollow PRIVATE
    TKPrim TKBO TKBool TKBRep TKTopAlgo TKGeomAlgo TKGeomBase TKG3d
    TKG2d TKMath TKernel TKDESTEP TKXSBase TKShHealing)
if(MSVC)
    set_target_properties(paacs-prepare-thin-hollow PROPERTIES
        VS_GLOBAL_VcpkgEnabled false VS_GLOBAL_VcpkgApplocalDeps false)
    target_compile_options(paacs-prepare-thin-hollow PRIVATE /W4 /fp:strict)
else()
    target_compile_options(paacs-prepare-thin-hollow PRIVATE
        -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
endif()
