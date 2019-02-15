cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

# WYRM IMPORT #
################

if(NOT TARGET wyrm)
    if(NOT DEFINED WYRM_ROOT)
        set(WYRM_ROOT $ENV{WYRM_ROOT})
    endif()
    if(NOT DEFINED WYRM_ROOT)
        message(STATUS "wyrm path set WYRM_ROOT to $COMPASS_ROOT/tplib/wyrm" )
        set(WYRM_ROOT $ENV{COMPASS_ROOT}/tplib/wyrm)
    endif()
    if(IS_DIRECTORY ${WYRM_ROOT})
		add_subdirectory(${WYRM_ROOT} wyrm_dir)
    else()
        message(ERROR "wyrm path not found, specify it by set WYRM_ROOT to the right path" )
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

target_link_libraries(carmaWrap PRIVATE wyrm PUBLIC carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS carmaWrap EXPORT CarmaWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python)
