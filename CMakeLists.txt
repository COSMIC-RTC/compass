cmake_minimum_required (VERSION 3.8)

if(NOT DEFINED VERSION_INFO)
    set(VERSION_INFO "4.0")
endif()

message(STATUS "COMPASS will be installed in: ${CMAKE_INSTALL_PREFIX}")

option(python_build "Compile Python Wrappers" ON)

project(compass
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS library"
    LANGUAGES CUDA CXX
)

add_subdirectory(libcarma)
add_subdirectory(libsutra)
# add_subdirectory(shesha)
