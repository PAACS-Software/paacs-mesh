# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.20)
project(paacs_cantilever_fixture LANGUAGES CXX)

find_package(OpenCASCADE CONFIG REQUIRED)
add_executable(make-cantilever-step make_step.cpp)
target_compile_features(make-cantilever-step PRIVATE cxx_std_17)
target_include_directories(make-cantilever-step PRIVATE "${OpenCASCADE_INCLUDE_DIR}")
target_link_libraries(make-cantilever-step PRIVATE TKPrim TKDESTEP)
