# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

# Compatibility entry point. Prefer cmake -S . from the project root.
cmake_minimum_required(VERSION 3.30...4.4)
project(PAACSMeshStandalone LANGUAGES C CXX)
enable_testing()
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/.." "${CMAKE_CURRENT_BINARY_DIR}/paacs")
