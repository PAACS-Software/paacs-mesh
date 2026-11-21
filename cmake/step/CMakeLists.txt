# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30)
project(PAACSStepBuild LANGUAGES C CXX)

include(ExternalProject)
enable_testing()

set(PAACS_STEP_JOBS 4 CACHE STRING "Parallel jobs for the OCCT and PAACS builds")
option(PAACS_STEP_BUILD_TOOLS "Build and test the PAACS STEP tools after OCCT" ON)
get_filename_component(paacs_source "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
include("${paacs_source}/cmake/OCCTVersion.cmake")
set(occt_install "${CMAKE_CURRENT_BINARY_DIR}/occt")
set(bin_dir "${CMAKE_CURRENT_BINARY_DIR}/bin")

set(compiler_args
    "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
)
foreach(paacs_setting CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
        CMAKE_OSX_SYSROOT CMAKE_TOOLCHAIN_FILE)
    if(DEFINED ${paacs_setting} AND NOT "${${paacs_setting}}" STREQUAL "")
        list(APPEND compiler_args "-D${paacs_setting}:STRING=${${paacs_setting}}")
    endif()
endforeach()
get_property(multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT multi_config)
    list(APPEND compiler_args "-DCMAKE_BUILD_TYPE:STRING=Release")
endif()
if(MSVC)
    # Ignore system-wide MSBuild package integration in this isolated build.
    list(APPEND compiler_args "-DCMAKE_VS_GLOBALS:STRING=VcpkgEnabled=false|VcpkgApplocalDeps=false")
endif()

set(occt_runtime_args)
if(APPLE)
    list(APPEND occt_runtime_args "-DCMAKE_INSTALL_RPATH:STRING=@loader_path")
elseif(UNIX)
    list(APPEND occt_runtime_args "-DCMAKE_INSTALL_RPATH:STRING=$ORIGIN")
endif()

# The toolkits below and their OCCT-declared dependencies form the STEP SDK.
# OCCT builds in its own CMake project and is installed as shared libraries.
ExternalProject_Add(occt
    URL "${PAACS_OCCT_URL}"
    URL_HASH "SHA512=${PAACS_OCCT_SHA512}"
    TLS_VERIFY ON
    DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/downloads"
    DOWNLOAD_NAME "occt-${PAACS_OCCT_VERSION}.tar.gz"
    INSTALL_DIR "${occt_install}"
    LIST_SEPARATOR "|"
    CMAKE_ARGS
        ${compiler_args}
        ${occt_runtime_args}
        "-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>"
        "-DCMAKE_EXPORT_NO_PACKAGE_REGISTRY:BOOL=ON"
        "-DINSTALL_DIR_LAYOUT:STRING=Unix"
        "-DINSTALL_DIR_BIN:STRING=bin"
        "-DINSTALL_DIR_LIB:STRING=lib"
        "-DINSTALL_DIR_CMAKE:STRING=lib/cmake/opencascade"
        "-DINSTALL_DIR_DOC:STRING=share/doc/opencascade"
        "-DBUILD_LIBRARY_TYPE:STRING=Shared"
        "-DBUILD_MODULE_FoundationClasses:BOOL=OFF"
        "-DBUILD_MODULE_ModelingData:BOOL=OFF"
        "-DBUILD_MODULE_ModelingAlgorithms:BOOL=OFF"
        "-DBUILD_MODULE_Visualization:BOOL=OFF"
        "-DBUILD_MODULE_ApplicationFramework:BOOL=OFF"
        "-DBUILD_MODULE_DataExchange:BOOL=OFF"
        "-DBUILD_MODULE_Draw:BOOL=OFF"
        "-DBUILD_ADDITIONAL_TOOLKITS:STRING=TKPrim TKBO TKBool TKMesh TKDESTEP"
        "-DUSE_TK:BOOL=OFF"
        "-DUSE_FREETYPE:BOOL=OFF"
        "-DUSE_XLIB:BOOL=OFF"
        "-DUSE_TBB:BOOL=OFF"
    BUILD_COMMAND "${CMAKE_COMMAND}" --build <BINARY_DIR> --config Release --parallel "${PAACS_STEP_JOBS}"
    INSTALL_COMMAND "${CMAKE_COMMAND}" --install <BINARY_DIR> --config Release
    USES_TERMINAL_DOWNLOAD ON
    USES_TERMINAL_BUILD ON
)

if(NOT PAACS_STEP_BUILD_TOOLS)
    return()
endif()

ExternalProject_Add(paacs
    DEPENDS occt
    SOURCE_DIR "${paacs_source}"
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/paacs"
    LIST_SEPARATOR "|"
    CMAKE_ARGS
        ${compiler_args}
        "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${bin_dir}"
        "-DOpenCASCADE_DIR:PATH=${occt_install}/lib/cmake/opencascade"
        "-DPAACS_CXX_STANDARD:STRING=23"
        "-DPAACS_WITH_OCCT:BOOL=ON"
        "-DPAACS_BUILD_CLI:BOOL=ON"
        "-DPAACS_BUILD_TESTS:BOOL=OFF"
        "-DPAACS_BUILD_UI_BRIDGE:BOOL=OFF"
    BUILD_COMMAND "${CMAKE_COMMAND}" --build <BINARY_DIR> --config Release
        --parallel "${PAACS_STEP_JOBS}"
        --target paacs-mesh paacs-cad-mesh paacs-native-cad-preview
    BUILD_ALWAYS ON
    INSTALL_COMMAND ""
    USES_TERMINAL_BUILD ON
)

add_test(NAME paacs.step.cylinder
    COMMAND "${bin_dir}/paacs-cad-mesh${CMAKE_EXECUTABLE_SUFFIX}"
        "${paacs_source}/tests/fixtures/cad/cylinder.step" .7
        "${CMAKE_CURRENT_BINARY_DIR}/cylinder.paacsmesh"
        --native-cad --face-size 1 .38 --size-box -.3 -.3 .7 .3 .3 1.3 .31
        --min-size .34 --growth 1.55 --report-json
)
set_tests_properties(paacs.step.cylinder PROPERTIES TIMEOUT 90)
if(WIN32)
    set_tests_properties(paacs.step.cylinder PROPERTIES
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${occt_install}/bin")
endif()
