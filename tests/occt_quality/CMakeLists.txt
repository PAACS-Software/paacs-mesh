# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30)
# Match the isolated native STEP build: system MSBuild package integration
# must not add unrelated dependencies or replace its app-local runtime.
set(CMAKE_VS_GLOBALS "VcpkgEnabled=false;VcpkgApplocalDeps=false")
project(PAACSNativeQualityQualification LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
get_property(paacs_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT paacs_multi_config AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()

set(PAACS_NATIVE_BUILD "" CACHE PATH "Existing Release PAACS native build directory")
if(NOT PAACS_NATIVE_BUILD)
    message(FATAL_ERROR "Set PAACS_NATIVE_BUILD to the existing PAACS native build directory")
endif()
get_filename_component(PAACS_NATIVE_BUILD "${PAACS_NATIVE_BUILD}" ABSOLUTE)
get_filename_component(paacs_source "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(PAACS_NATIVE_RUNTIME_DIR "${PAACS_NATIVE_BUILD}/../bin" CACHE PATH
    "Optional native tool directory supplying matching Windows CRT DLLs")

find_package(Threads REQUIRED)
find_package(OpenCASCADE REQUIRED CONFIG)

# Import the already-built engine rather than reconfiguring it or requiring
# optional dependencies used only by the full project test suite.
foreach(component IN ITEMS mesh occt)
    set(paacs_library_name "${CMAKE_STATIC_LIBRARY_PREFIX}paacs_${component}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(paacs_library "")
    foreach(directory IN ITEMS "${PAACS_NATIVE_BUILD}/Release" "${PAACS_NATIVE_BUILD}")
        if(EXISTS "${directory}/${paacs_library_name}")
            set(paacs_library "${directory}/${paacs_library_name}")
            break()
        endif()
    endforeach()
    if(NOT paacs_library)
        message(FATAL_ERROR "Missing ${paacs_library_name} in ${PAACS_NATIVE_BUILD}/Release or ${PAACS_NATIVE_BUILD}")
    endif()
    add_library(qualified_${component} STATIC IMPORTED)
    set_target_properties(qualified_${component} PROPERTIES IMPORTED_LOCATION "${paacs_library}")
endforeach()

# The existing Apple build records whether its public stop_token/jthread
# headers require libc++'s experimental-library flag. Carry that choice into
# these consumers of its static libraries.
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    load_cache("${PAACS_NATIVE_BUILD}" READ_WITH_PREFIX PAACS_NATIVE_
        PAACS_APPLE_EXPERIMENTAL_STOP_TOKEN)
endif()
if(WIN32)
    get_filename_component(paacs_occt_prefix "${OpenCASCADE_DIR}/../../.." ABSOLUTE)
    file(GLOB paacs_crt_dlls "${PAACS_NATIVE_RUNTIME_DIR}/*140*.dll")
endif()

enable_testing()
foreach(name IN ITEMS delaunay_cavity delaunay_quality interior_sites)
    add_executable(${name} "${paacs_source}/tests/occt_${name}.cpp")
    target_include_directories(${name} PRIVATE
        "${paacs_source}/include" "${paacs_source}/src/geometry")
    target_include_directories(${name} SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
    if(MSVC)
        target_compile_options(${name} PRIVATE /fp:strict)
    else()
        target_compile_options(${name} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
    if(PAACS_NATIVE_PAACS_APPLE_EXPERIMENTAL_STOP_TOKEN)
        target_compile_options(${name} PRIVATE -fexperimental-library)
        target_link_options(${name} PRIVATE -fexperimental-library)
    endif()
    target_link_libraries(${name} PRIVATE qualified_occt qualified_mesh Threads::Threads
        TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKPrim TKBO TKBool TKMesh TKXSBase TKDESTEP)
    add_test(NAME paacs.occt.${name} COMMAND ${name})
    set_tests_properties(paacs.occt.${name} PROPERTIES TIMEOUT 120 LABELS "occt;quality;regression")
    if(WIN32)
        set_tests_properties(paacs.occt.${name} PROPERTIES
            ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${paacs_occt_prefix}/bin")
        if(paacs_crt_dlls)
            add_custom_command(TARGET ${name} POST_BUILD
                COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                    ${paacs_crt_dlls} "$<TARGET_FILE_DIR:${name}>"
                VERBATIM)
        endif()
    endif()
endforeach()
