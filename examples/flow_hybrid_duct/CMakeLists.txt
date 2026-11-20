# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.20)
project(paacs_hybrid_duct LANGUAGES CXX)

find_package(PAACSMesh CONFIG REQUIRED COMPONENTS OCCT)
add_executable(mesh-hybrid-duct mesh_hybrid.cpp)
target_link_libraries(mesh-hybrid-duct PRIVATE PAACS::occt)
