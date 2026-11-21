# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

# A wheel source build uses the same pinned, shared OCCT project as the STEP SDK.
# An explicit OpenCASCADE_DIR and PAACS_PYTHON_BUILD_OCCT=OFF allow local SDK reuse.
set(PAACS_PYTHON_OCCT_JOBS 4 CACHE STRING "Parallel jobs for the wheel's OCCT build")
if(NOT PAACS_PYTHON_OCCT_JOBS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "PAACS_PYTHON_OCCT_JOBS must be a positive integer")
endif()
set(paacs_python_occt_build "${CMAKE_BINARY_DIR}/_occt-build")
set(paacs_python_occt_args
    -S "${PROJECT_SOURCE_DIR}/cmake/step"
    -B "${paacs_python_occt_build}"
    -G "${CMAKE_GENERATOR}"
    "-DCMAKE_BUILD_TYPE:STRING=Release"
    "-DPAACS_STEP_JOBS:STRING=${PAACS_PYTHON_OCCT_JOBS}"
    "-DPAACS_STEP_BUILD_TOOLS:BOOL=OFF")
foreach(paacs_setting CMAKE_GENERATOR_PLATFORM CMAKE_GENERATOR_TOOLSET
        CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM CMAKE_OSX_ARCHITECTURES
        CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_TOOLCHAIN_FILE)
    if(DEFINED ${paacs_setting} AND NOT "${${paacs_setting}}" STREQUAL "")
        list(APPEND paacs_python_occt_args "-D${paacs_setting}:STRING=${${paacs_setting}}")
    endif()
endforeach()
execute_process(COMMAND "${CMAKE_COMMAND}" ${paacs_python_occt_args}
    COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${paacs_python_occt_build}"
    --config Release --target occt --parallel "${PAACS_PYTHON_OCCT_JOBS}"
    COMMAND_ERROR_IS_FATAL ANY)
set(OpenCASCADE_DIR "${paacs_python_occt_build}/occt/lib/cmake/opencascade"
    CACHE PATH "Pinned shared OCCT SDK for Python" FORCE)
