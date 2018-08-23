cmake_minimum_required (VERSION 3.8)

project(naga
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS NAGA lirary"
    LANGUAGES CXX CUDA
)

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

link_directories( ${MAGMA_LIBRARY_DIRS} )
if(NOT TARGET wyrm)
    add_subdirectory(../wyrm wyrm)
endif()
add_library(naga MODULE
                    src/naga_context.cpp
                    src/naga.cpp
                    src/naga_obj_double.cpp
                    src/naga_obj_float.cpp
                    src/naga_sparse_obj_float.cpp
                    src/naga_sparse_obj_double.cpp
                    src/naga_obj_double2.cpp
                    src/naga_obj_float2.cpp
                    src/naga_obj_double_complex.cpp
                    src/naga_obj_float_complex.cpp
                    src/naga_obj.hpp
                    src/naga_sparse_obj.hpp
                    src/naga_obj_int.cpp
                    src/naga_obj_uint.cpp)

target_link_libraries(naga PRIVATE wyrm carma)

set_target_properties(naga PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
target_compile_features(naga PUBLIC cxx_std_11)
