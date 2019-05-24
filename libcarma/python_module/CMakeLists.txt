cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

# WYRM IMPORT #
################

if(NOT TARGET wyrm::wyrm)
    if(hunter_build)
        hunter_add_package(wyrm)
        find_package(wyrm CONFIG REQUIRED)
    else()
        if(IS_DIRECTORY $ENV{WYRM_ROOT})
            add_subdirectory($ENV{WYRM_ROOT} wyrm_dir)
            add_library(wyrm::wyrm ALIAS wyrm) # to match exported target
        else()
            message(FATAL_ERROR "Wyrm not found please define WYRM_ROOT or use hunter")
        endif()
    endif()
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)
link_directories( ${MAGMA_LIBRARY_DIRS} )

add_library(carmaWrap MODULE
src/context.cpp
src/timer.cpp
src/host_obj.cpp
src/obj.cpp
src/sparse_obj.cpp
src/carmaWrap.cpp
)

target_link_libraries(carmaWrap PRIVATE wyrm::wyrm carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS carmaWrap EXPORT CarmaWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python)
