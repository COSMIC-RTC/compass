cmake_minimum_required (VERSION 3.18)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

find_package(pybind11 REQUIRED)
find_package(wyrm REQUIRED)

pybind11_add_module(carmaWrap MODULE
    src/context.cpp
    src/timer.cpp
    src/host_obj.cpp
    src/obj.cpp
    src/sparse_obj.cpp
    src/carmaWrap.cpp
)

target_link_libraries(carmaWrap PRIVATE wyrm::wyrm carma)

install(TARGETS carmaWrap EXPORT CarmaWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python
)
