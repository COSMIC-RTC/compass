cmake_minimum_required (VERSION 3.8)

project(compass
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS library"
    LANGUAGES CUDA CXX
)

add_subdirectory(libcarma)
add_subdirectory(libsutra)
# add_subdirectory(shesha)
