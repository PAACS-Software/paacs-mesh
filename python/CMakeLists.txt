# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 3.1 CONFIG REQUIRED)
pybind11_add_module(_native MODULE NO_EXTRAS bindings.cpp)
# Keep experimental artifacts identifiable independently of package versions.
set(paacs_python_source_revision "unversioned")
set(paacs_python_source_dirty 1)
find_package(Git QUIET)
if(Git_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" rev-parse HEAD
        RESULT_VARIABLE paacs_python_git_result OUTPUT_VARIABLE paacs_python_git_revision
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(paacs_python_git_result EQUAL 0 AND paacs_python_git_revision MATCHES "^[0-9a-f]+$")
        set(paacs_python_source_revision "${paacs_python_git_revision}")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" status --porcelain --untracked-files=normal
            RESULT_VARIABLE paacs_python_status_result OUTPUT_VARIABLE paacs_python_git_status
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(paacs_python_status_result EQUAL 0 AND paacs_python_git_status STREQUAL "")
            set(paacs_python_source_dirty 0)
        endif()
    endif()
endif()
target_compile_definitions(_native PRIVATE
    PAACS_PYTHON_SOURCE_REVISION="${paacs_python_source_revision}"
    PAACS_PYTHON_SOURCE_DIRTY=${paacs_python_source_dirty})
target_link_libraries(_native PRIVATE PAACS::occt)
install(TARGETS _native
    LIBRARY DESTINATION paacs_mesh COMPONENT Python
    RUNTIME DESTINATION paacs_mesh COMPONENT Python)

# pybind11 headers contribute code to the extension, so retain their BSD notice.
execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "from importlib.metadata import distribution; d=distribution('pybind11'); print(next(str(d.locate_file(p)) for p in d.files if p.name == 'LICENSE'))"
    OUTPUT_VARIABLE paacs_pybind_license OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
file(TO_CMAKE_PATH "${paacs_pybind_license}" paacs_pybind_license)
install(FILES "${paacs_pybind_license}" DESTINATION paacs_mesh/_occt/licenses
    RENAME pybind11-LICENSE.txt COMPONENT Python)
if(WIN32)
    install(FILES "${PROJECT_SOURCE_DIR}/cmake/python/delvewheel-LICENSE.txt"
        DESTINATION paacs_mesh/_occt/licenses COMPONENT Python)
endif()

include("${PROJECT_SOURCE_DIR}/cmake/OCCTVersion.cmake")
if(NOT "${OpenCASCADE_MAJOR_VERSION}.${OpenCASCADE_MINOR_VERSION}.${OpenCASCADE_MAINTENANCE_VERSION}"
        VERSION_EQUAL "${PAACS_OCCT_VERSION}")
    message(FATAL_ERROR "Python wheels require the pinned OCCT ${PAACS_OCCT_VERSION} SDK")
endif()

# The SDK remains separate; only notices and runtime resources enter the wheel.
set(paacs_occt_doc "${OpenCASCADE_INSTALL_PREFIX}/share/doc/opencascade")
foreach(paacs_license LICENSE_LGPL_21.txt OCCT_LGPL_EXCEPTION.txt)
    if(NOT EXISTS "${paacs_occt_doc}/${paacs_license}")
        message(FATAL_ERROR "Missing OCCT redistribution notice: ${paacs_occt_doc}/${paacs_license}")
    endif()
    install(FILES "${paacs_occt_doc}/${paacs_license}"
        DESTINATION paacs_mesh/_occt/licenses COMPONENT Python)
endforeach()
configure_file("${PROJECT_SOURCE_DIR}/cmake/python/OCCT_NOTICE.txt.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OCCT_NOTICE.txt" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/OCCT_NOTICE.txt"
    DESTINATION paacs_mesh/_occt COMPONENT Python)
foreach(paacs_resource BOPAlgo SHMessage XSMessage XSTEPResource UnitsAPI StdResource)
    if(EXISTS "${OpenCASCADE_RESOURCE_DIR}/${paacs_resource}")
        install(DIRECTORY "${OpenCASCADE_RESOURCE_DIR}/${paacs_resource}"
            DESTINATION paacs_mesh/_occt/resources COMPONENT Python)
    endif()
endforeach()
