cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "carmaWrap library"
    LANGUAGES CXX CUDA
)

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)
link_directories( ${MAGMA_LIBRARY_DIRS} )

add_library(carmaWrap MODULE
                    src/context.cpp
                    src/timer.cpp
                    src/carmaWrap.cpp
                    src/obj_double.cpp
                    src/obj_float.cpp
                    src/obj_half.cpp
                    src/sparse_obj_float.cpp
                    src/sparse_obj_double.cpp
                    src/obj_double2.cpp
                    src/obj_float2.cpp
                    src/obj_double_complex.cpp
                    src/obj_float_complex.cpp
                    src/obj.hpp
                    src/sparse_obj.hpp
                    src/obj_int.cpp
                    src/obj_uint.cpp)

target_link_libraries(carmaWrap PRIVATE wyrm carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
target_compile_features(carmaWrap PUBLIC cxx_std_11)
