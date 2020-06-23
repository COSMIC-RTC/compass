cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

# WYRM IMPORT #
################

find_package(pybind11 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA magma)
if(${MAGMA_FOUND})
  link_directories(${MAGMA_LIBRARY_DIRS})
else()
  message(STATUS "MAGMA not found, MAGMA features disabled")
endif()

pybind11_add_module(carmaWrap MODULE
    src/context.cpp
    src/timer.cpp
    src/host_obj.cpp
    src/obj.cpp
    src/sparse_obj.cpp
    src/carmaWrap.cpp
)

target_link_libraries(carmaWrap PRIVATE CONAN_PKG::wyrm carma)

install(TARGETS carmaWrap EXPORT CarmaWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python)
