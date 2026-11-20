# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30...4.4)
project(PAACSInstalledComponents LANGUAGES CXX)
if(PAACS_COMPONENT_MODE STREQUAL "core")
    find_package(PAACSMesh 1.2.0 EXACT CONFIG REQUIRED)
elseif(PAACS_COMPONENT_MODE STREQUAL "optional")
    find_package(PAACSMesh 1.2.0 EXACT CONFIG REQUIRED OPTIONAL_COMPONENTS OCCT)
elseif(PAACS_COMPONENT_MODE STREQUAL "required")
    find_package(PAACSMesh 1.2.0 EXACT CONFIG REQUIRED COMPONENTS OCCT)
else()
    message(FATAL_ERROR "Unknown PAACS_COMPONENT_MODE")
endif()
if(NOT TARGET PAACS::mesh OR NOT TARGET PAACS::Mesh)
    message(FATAL_ERROR "Missing portable core or compatibility target")
endif()
if(PAACS_COMPONENT_EXPECT_OCCT)
    if(NOT PAACSMesh_OCCT_FOUND OR NOT TARGET PAACS::occt OR
       NOT TARGET PAACS::MeshOCCT)
        message(FATAL_ERROR "Requested OCCT component was not loaded")
    endif()
elseif(PAACSMesh_OCCT_FOUND OR TARGET PAACS::occt)
    message(FATAL_ERROR "An unrequested or unavailable OCCT component was loaded")
endif()
if(PAACS_COMPONENT_MODE STREQUAL "core" AND
   DEFINED OpenCASCADE_CONFIG)
    message(FATAL_ERROR "Core-only lookup discovered OpenCASCADE")
endif()
add_executable(installed-core main.cpp)
target_link_libraries(installed-core PRIVATE PAACS::Mesh)
if(PAACS_COMPONENT_EXPECT_OCCT)
    # Force resolution of the optional export's real toolkit link closure.
    # CAD runtime/round-trip coverage lives in the installed consumer example.
    add_executable(installed-cad cad.cpp)
    target_link_libraries(installed-cad PRIVATE PAACS::MeshOCCT)
endif()
enable_testing()
add_test(NAME installed.core COMMAND installed-core)
