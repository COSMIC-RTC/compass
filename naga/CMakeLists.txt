cmake_minimum_required (VERSION 3.8)

project(naga
    VERSION 1.0
    DESCRIPTION "COMPASS NAGA lirary"
    LANGUAGES CXX CUDA
)

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

link_directories( ${MAGMA_LIBRARY_DIRS} )

add_subdirectory(../pybind11 pybind11)

add_library(naga MODULE
                    src/naga_context.cpp
                    src/naga.cpp
                    src/naga_obj_double.cpp
                    src/naga_obj_float.cpp
                    src/naga_obj.hpp
                    src/naga_obj_int.cpp
                    src/naga_obj_uint.cpp)

target_link_libraries(naga PRIVATE pybind11::module carma)

set_target_properties(naga PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
