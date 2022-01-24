cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED VERSION_INFO)
    set(VERSION_INFO 5.2.1)
endif()

option(python_build "Compile Python Wrappers" ON)
option(do_half "Compile half codes" ON)

project(compass
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS library"
        LANGUAGES CUDA CXX)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(GCOV gcov)
endif()

add_subdirectory(libcarma)
add_subdirectory(libsutra)
