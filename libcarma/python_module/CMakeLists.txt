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
src/obj.cpp
src/sparse_obj.cpp
src/carmaWrap.cpp
)

target_compile_features(carmaWrap PUBLIC cxx_std_14)
target_link_libraries(carmaWrap PRIVATE wyrm PUBLIC carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS carmaWrap
    PUBLIC_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib}
)
