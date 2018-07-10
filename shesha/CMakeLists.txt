cmake_minimum_required (VERSION 3.8)

project(shesha_bind
    VERSION 1.0
    DESCRIPTION "COMPASS shesha_atmos lirary"
    LANGUAGES CXX CUDA
)

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

link_directories( ${MAGMA_LIBRARY_DIRS} )
if(NOT TARGET wyrm)
    add_subdirectory(../wyrm wyrm)
endif()
add_library(shesha_bind MODULE
                    shesha/sutra_pybind/shesha_atmos.cpp
                    shesha/sutra_pybind/shesha_bind.cpp)
#link_directories($ENV{COMPASS_ROOT}/libsutra)

target_link_libraries(shesha_bind PRIVATE wyrm carma sutra)

set_target_properties(shesha_bind PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
