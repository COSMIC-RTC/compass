cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

if(NOT TARGET wyrm)
    add_subdirectory($ENV{COMPASS_ROOT}/tplib/wyrm wyrm)
endif()

add_library(carmaWrap MODULE
src/context.cpp
src/timer.cpp
src/obj.cpp
src/sparse_obj.cpp
src/carmaWrap.cpp
)

target_link_libraries(carmaWrap PRIVATE wyrm PUBLIC carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS carmaWrap
    PUBLIC_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib}
)
