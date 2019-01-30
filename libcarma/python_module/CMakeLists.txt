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
        message(WARNING "set WYRM_ROOT to the default path : $COMPASS_ROOT/tplib/wyrm" )
        set(WYRM_ROOT $ENV{COMPASS_ROOT}/tplib/wyrm)
    endif()
	if(IS_DIRECTORY ${WYRM_ROOT})
		add_subdirectory(${WYRM_ROOT} wyrm_dir)
    else()
        message(FATAL_ERROR "wyrm path not found, specify it by setting WYRM_ROOT to the right path" )
    endif()
endif()

add_library(carmaWrap MODULE
src/context.cpp
src/timer.cpp
src/obj.cpp
src/sparse_obj.cpp
src/carmaWrap.cpp
)

target_link_libraries(carmaWrap PRIVATE wyrm PUBLIC carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS carmaWrap
    PUBLIC_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib}
)
