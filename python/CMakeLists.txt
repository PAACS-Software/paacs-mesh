# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 3.1 CONFIG REQUIRED)
pybind11_add_module(_native MODULE NO_EXTRAS bindings.cpp)
target_link_libraries(_native PRIVATE PAACS::occt)
install(TARGETS _native
    LIBRARY DESTINATION paacs_mesh COMPONENT Python
    RUNTIME DESTINATION paacs_mesh COMPONENT Python)

# pybind11 headers contribute code to the extension, so retain their BSD notice.
execute_process(COMMAND "${Python_EXECUTABLE}" -c
    "from importlib.metadata import distribution; d=distribution('pybind11'); print(next(str(d.locate_file(p)) for p in d.files if p.name == 'LICENSE'))"
    OUTPUT_VARIABLE paacs_pybind_license OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
install(FILES "${paacs_pybind_license}" DESTINATION paacs_mesh/_occt/licenses
    RENAME pybind11-LICENSE.txt COMPONENT Python)
install(FILES "${PROJECT_SOURCE_DIR}/cmake/python/delvewheel-LICENSE.txt"
    DESTINATION paacs_mesh/_occt/licenses COMPONENT Python)

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
