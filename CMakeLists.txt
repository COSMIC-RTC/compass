cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED VERSION_INFO)
    set(VERSION_INFO 6.0.0)
endif()


project(compass
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS library"
        LANGUAGES CUDA CXX)

option(libs_build "Compile Python Wrappers" ON)
option(python_build "Compile Python Wrappers" ON)

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
