cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED VERSION_INFO)
    set(VERSION_INFO 5.5.0)
endif()

option(libs_build "Compile Python Wrappers" ON)
option(python_build "Compile Python Wrappers" ON)
option(do_half "Compile half codes" OFF)

project(compass
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS library"
        LANGUAGES CUDA CXX)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(GCOV gcov)
endif()

if(libs_build)
    add_subdirectory(libcarma)
    add_subdirectory(libsutra)
endif()

if(python_build)
    add_subdirectory(python_module)
endif()
