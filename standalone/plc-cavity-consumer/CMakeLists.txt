# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30)
project(PaacsPlcCavityConsumer LANGUAGES CXX)
find_package(PAACSMesh CONFIG REQUIRED)
add_executable(paacs-cavity-consumer main.cpp)
target_link_libraries(paacs-cavity-consumer PRIVATE PAACS::mesh)
enable_testing()
add_test(NAME paacs.installed.cavity_recovery COMMAND paacs-cavity-consumer)
