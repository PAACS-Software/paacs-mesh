# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.20)
project(paacs_plate_hole_fixture LANGUAGES CXX)

find_package(OpenCASCADE CONFIG REQUIRED)
add_executable(make-plate-hole-step make_step.cpp)
target_compile_features(make-plate-hole-step PRIVATE cxx_std_17)
target_include_directories(make-plate-hole-step PRIVATE "${OpenCASCADE_INCLUDE_DIR}")
target_link_libraries(make-plate-hole-step PRIVATE TKPrim TKBO TKDESTEP)

option(PAACS_PLATE_BUILD_MESHER "Build the specialized PAACS sweep example" ON)
if(PAACS_PLATE_BUILD_MESHER)
    find_package(PAACSMesh CONFIG REQUIRED COMPONENTS OCCT)
    add_executable(mesh-plate-hole mesh_plate_hole.cpp)
    target_link_libraries(mesh-plate-hole PRIVATE PAACS::occt)
endif()
