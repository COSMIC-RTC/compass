cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED VERSION_INFO)
  set(VERSION_INFO 5.1.0)
endif()

message(STATUS "COMPASS will be installed in: ${CMAKE_INSTALL_PREFIX}")

option(python_build "Compile Python Wrappers" ON)
option(do_half "Compile half codes" ON)

project(compass
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS library"
        LANGUAGES CUDA CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS NO_OUTPUT_DIRS)

add_subdirectory(libcarma)
add_subdirectory(libsutra)
# add_subdirectory(shesha)
