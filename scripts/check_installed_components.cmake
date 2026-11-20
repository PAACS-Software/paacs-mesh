# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.30)
if(NOT PAACS_INSTALL_PREFIX OR NOT PAACS_COMPONENT_TEST_BINARY_DIR)
    message(FATAL_ERROR "Set PAACS_INSTALL_PREFIX and PAACS_COMPONENT_TEST_BINARY_DIR")
endif()
get_filename_component(_paacs_source "${CMAKE_CURRENT_LIST_DIR}/../tests/package_components" ABSOLUTE)
get_filename_component(PAACS_INSTALL_PREFIX "${PAACS_INSTALL_PREFIX}" ABSOLUTE)
get_filename_component(PAACS_COMPONENT_TEST_BINARY_DIR "${PAACS_COMPONENT_TEST_BINARY_DIR}" ABSOLUTE)
set(_paacs_common "-DCMAKE_PREFIX_PATH=${PAACS_INSTALL_PREFIX}")
if(PAACS_COMPONENT_TEST_GENERATOR)
    list(APPEND _paacs_common -G "${PAACS_COMPONENT_TEST_GENERATOR}")
endif()
if(PAACS_COMPONENT_TEST_GENERATOR_PLATFORM)
    list(APPEND _paacs_common -A "${PAACS_COMPONENT_TEST_GENERATOR_PLATFORM}")
endif()
if(PAACS_COMPONENT_TEST_GENERATOR_TOOLSET)
    list(APPEND _paacs_common -T "${PAACS_COMPONENT_TEST_GENERATOR_TOOLSET}")
endif()
if(PAACS_COMPONENT_TEST_CXX_COMPILER)
    list(APPEND _paacs_common "-DCMAKE_CXX_COMPILER=${PAACS_COMPONENT_TEST_CXX_COMPILER}")
endif()
if(PAACS_COMPONENT_TEST_MAKE_PROGRAM)
    list(APPEND _paacs_common "-DCMAKE_MAKE_PROGRAM=${PAACS_COMPONENT_TEST_MAKE_PROGRAM}")
endif()
if(NOT PAACS_COMPONENT_TEST_CONFIG)
    set(PAACS_COMPONENT_TEST_CONFIG Release)
endif()

function(_paacs_check name mode disable expect_occt expect_success)
    set(_paacs_build "${PAACS_COMPONENT_TEST_BINARY_DIR}/${name}")
    set(_paacs_args ${_paacs_common}
        "-DCMAKE_BUILD_TYPE=${PAACS_COMPONENT_TEST_CONFIG}"
        "-DPAACS_COMPONENT_MODE=${mode}"
        "-DPAACS_COMPONENT_EXPECT_OCCT=${expect_occt}"
        "-DCMAKE_DISABLE_FIND_PACKAGE_OpenCASCADE=${disable}")
    if(NOT disable AND PAACS_COMPONENT_TEST_OCCT_DIR)
        list(APPEND _paacs_args "-DOpenCASCADE_DIR=${PAACS_COMPONENT_TEST_OCCT_DIR}")
    endif()
    if(WIN32)
        # Global MSBuild vcpkg integration must not substitute app-local OCCT
        # DLLs from an unrelated installation for the explicitly supplied SDK.
        execute_process(COMMAND "${CMAKE_COMMAND}" -S "${_paacs_source}" -B "${_paacs_build}" ${_paacs_args}
            "-DCMAKE_VS_GLOBALS=VcpkgEnabled=false;VcpkgApplocalDeps=false"
            RESULT_VARIABLE _paacs_status OUTPUT_VARIABLE _paacs_output ERROR_VARIABLE _paacs_error)
    else()
        execute_process(COMMAND "${CMAKE_COMMAND}" -S "${_paacs_source}" -B "${_paacs_build}" ${_paacs_args}
            RESULT_VARIABLE _paacs_status OUTPUT_VARIABLE _paacs_output ERROR_VARIABLE _paacs_error)
    endif()
    file(MAKE_DIRECTORY "${_paacs_build}")
    file(WRITE "${_paacs_build}/configure-check.log" "${_paacs_output}${_paacs_error}")
    if(NOT expect_success)
        if(_paacs_status EQUAL 0)
            message(FATAL_ERROR "${name}: required unavailable OCCT unexpectedly configured")
        endif()
        if(NOT "${_paacs_output}${_paacs_error}" MATCHES "OpenCASCADE|OCCT")
            message(FATAL_ERROR "${name}: failed for an unrelated reason; see configure-check.log")
        endif()
        message(STATUS "${name}: unavailable required component rejected")
        return()
    endif()
    if(NOT _paacs_status EQUAL 0)
        message(FATAL_ERROR "${name}: configure failed\n${_paacs_output}${_paacs_error}")
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}" --build "${_paacs_build}" --config "${PAACS_COMPONENT_TEST_CONFIG}"
        RESULT_VARIABLE _paacs_status OUTPUT_VARIABLE _paacs_output ERROR_VARIABLE _paacs_error)
    file(WRITE "${_paacs_build}/build-check.log" "${_paacs_output}${_paacs_error}")
    if(NOT _paacs_status EQUAL 0)
        message(FATAL_ERROR "${name}: core consumer build failed\n${_paacs_output}${_paacs_error}")
    endif()
    get_filename_component(_paacs_cmake_bin "${CMAKE_COMMAND}" DIRECTORY)
    execute_process(COMMAND "${_paacs_cmake_bin}/ctest" --test-dir "${_paacs_build}"
        -C "${PAACS_COMPONENT_TEST_CONFIG}" --output-on-failure
        RESULT_VARIABLE _paacs_status OUTPUT_VARIABLE _paacs_output ERROR_VARIABLE _paacs_error)
    file(WRITE "${_paacs_build}/test-check.log" "${_paacs_output}${_paacs_error}")
    if(NOT _paacs_status EQUAL 0)
        message(FATAL_ERROR "${name}: core consumer test failed\n${_paacs_output}${_paacs_error}")
    endif()
    message(STATUS "${name}: core consumer compiled and ran")
endfunction()

_paacs_check(core core TRUE FALSE TRUE)
_paacs_check(optional-unavailable optional TRUE FALSE TRUE)
_paacs_check(required-unavailable required TRUE FALSE FALSE)
if(PAACS_COMPONENT_TEST_OCCT_DIR)
    _paacs_check(optional-available optional FALSE TRUE TRUE)
    _paacs_check(required-available required FALSE TRUE TRUE)
endif()
