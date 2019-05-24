cmake_minimum_required(VERSION 3.8)

if(NOT DEFINED VERSION_INFO)
  set(VERSION_INFO "4.2")
endif()

message(STATUS "COMPASS will be installed in: ${CMAKE_INSTALL_PREFIX}")

option(python_build "Compile Python Wrappers" ON)
option(hunter_build "Compile using hunter" OFF)

if(hunter_build)
    include("cmake/HunterGate.cmake")

    HunterGate(
      URL "https://github.com/ruslo/hunter/archive/v0.23.141.tar.gz"
      SHA1 "d81118e8ceb47b829dfb4bd63f5048cacb315780"
      )
endif()

project(compass
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS library"
        LANGUAGES CUDA CXX)

add_subdirectory(libcarma)
add_subdirectory(libsutra)
# add_subdirectory(shesha)
