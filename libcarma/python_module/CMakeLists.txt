cmake_minimum_required (VERSION 3.8)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "carmaWrap library"
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

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

add_library(carmaWrap MODULE
src/context.cpp
src/timer.cpp
src/carmaWrap.cpp
src/obj_double.cpp
src/obj_float.cpp
src/obj_half.cpp
src/sparse_obj_float.cpp
src/sparse_obj_double.cpp
src/obj_double2.cpp
src/obj_float2.cpp
src/obj_double_complex.cpp
src/obj_float_complex.cpp
src/obj.hpp
src/sparse_obj.hpp
src/obj_int.cpp
src/obj_uint.cpp)

target_compile_features(carmaWrap PUBLIC cxx_std_14)
target_link_libraries(carmaWrap PRIVATE wyrm carma)

set_target_properties(carmaWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
