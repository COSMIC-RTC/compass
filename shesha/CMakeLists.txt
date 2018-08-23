cmake_minimum_required (VERSION 3.8)

project(shesha_bind
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS shesha library"
    LANGUAGES CXX CUDA
)

if (DEFINED ENV{BRAHMA_ROOT})
    SET(BRAHMA_DIRS $ENV{ACE_ROOT}/lib $ENV{TAO_ROOT}/tao $ENV{DDS_ROOT}/lib $ENV{BRAHMA_ROOT})
    link_directories( ${BRAHMA_DIRS} )
    SET(BRAHMA_CPP shesha/sutra_pybind/shesha_rtc_brahma.cpp
                shesha/sutra_pybind/shesha_target_brahma.cpp)
endif()

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

link_directories( ${MAGMA_LIBRARY_DIRS} )
if(NOT TARGET wyrm)
    add_subdirectory(../wyrm wyrm)
endif()

add_library(shesha_bind MODULE
                    shesha/sutra_pybind/shesha_telescope.cpp
                    shesha/sutra_pybind/shesha_tscreen.cpp
                    shesha/sutra_pybind/shesha_atmos.cpp
                    shesha/sutra_pybind/shesha_dms.cpp
                    shesha/sutra_pybind/shesha_source.cpp
                    shesha/sutra_pybind/shesha_lgs.cpp
                    shesha/sutra_pybind/shesha_target.cpp
                    shesha/sutra_pybind/shesha_sensors.cpp
                    shesha/sutra_pybind/shesha_wfs.cpp
                    shesha/sutra_pybind/shesha_wfs_sh.cpp
                    shesha/sutra_pybind/shesha_wfs_pyrhr.cpp
                    shesha/sutra_pybind/shesha_centroider.cpp
                    shesha/sutra_pybind/shesha_centroider_tcog.cpp
                    shesha/sutra_pybind/shesha_centroider_wcog.cpp
                    shesha/sutra_pybind/shesha_centroider_bpcog.cpp
                    shesha/sutra_pybind/shesha_centroider_corr.cpp
                    shesha/sutra_pybind/shesha_centroider_pyr.cpp
                    shesha/sutra_pybind/shesha_controller.cpp
                    shesha/sutra_pybind/shesha_controller_ls.cpp
                    shesha/sutra_pybind/shesha_controller_mv.cpp
                    shesha/sutra_pybind/shesha_controller_geo.cpp
                    shesha/sutra_pybind/shesha_controller_generic.cpp
                    shesha/sutra_pybind/shesha_controller_cured.cpp
                    shesha/sutra_pybind/shesha_rtc.cpp
                    shesha/sutra_pybind/shesha_gamora.cpp
                    shesha/sutra_pybind/shesha_groot.cpp
                    shesha/sutra_pybind/shesha_bind.cpp
                    ${BRAHMA_CPP}
                    )
#link_directories($ENV{COMPASS_ROOT}/libsutra)

target_link_libraries(shesha_bind PRIVATE wyrm sutra)

set_target_properties(shesha_bind PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")
target_compile_features(shesha_bind PUBLIC cxx_std_11)
