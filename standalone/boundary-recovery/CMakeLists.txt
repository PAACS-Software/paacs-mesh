# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30)
project(PaacsRecoveryExperiments LANGUAGES CXX)
find_package(PAACSMesh CONFIG REQUIRED)
add_executable(paacs-recovery-checks recovery_checks.cpp)
target_link_libraries(paacs-recovery-checks PRIVATE PAACS::mesh)
target_include_directories(paacs-recovery-checks BEFORE PRIVATE ../../include ../../src/bnd)
enable_testing()
add_test(NAME paacs.experimental.recovery COMMAND paacs-recovery-checks)
set_tests_properties(paacs.experimental.recovery PROPERTIES TIMEOUT 30)
