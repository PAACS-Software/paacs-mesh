# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.20)
project(paacs_straight_duct LANGUAGES CXX)

find_package(PAACSMesh CONFIG REQUIRED COMPONENTS OCCT)

add_executable(make-straight-duct-step make_step.cpp)
target_compile_features(make-straight-duct-step PRIVATE cxx_std_17)
target_include_directories(make-straight-duct-step PRIVATE "${OpenCASCADE_INCLUDE_DIR}")
target_link_libraries(make-straight-duct-step PRIVATE TKPrim TKDESTEP)

add_executable(mesh-straight-duct duct_mesh.cpp)
target_link_libraries(mesh-straight-duct PRIVATE PAACS::occt)

option(PAACS_DUCT_BUILD_OPPOSED_WALL_LAYERS "Build the wall-normal prism comparison" ON)
if(PAACS_DUCT_BUILD_OPPOSED_WALL_LAYERS)
    add_executable(duct-opposed-wall-layers opposed_wall_layers.cpp)
    target_link_libraries(duct-opposed-wall-layers PRIVATE PAACS::mesh)
endif()
