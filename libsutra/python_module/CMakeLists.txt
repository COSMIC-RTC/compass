cmake_minimum_required (VERSION 3.8)

project(sutraWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS sutraWrap library"
    LANGUAGES CXX CUDA
)

# set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# add_definitions(${sloth_DEFINITIONS})

if (DEFINED ENV{BRAHMA_ROOT})
    SET(BRAHMA_DIRS $ENV{ACE_ROOT}/lib $ENV{TAO_ROOT}/tao $ENV{DDS_ROOT}/lib $ENV{BRAHMA_ROOT})
    link_directories( ${BRAHMA_DIRS} )
    SET(BRAHMA_CPP src/rtc_brahma.cpp
                src/target_brahma.cpp)
endif()

if (DEFINED ENV{OCTOPUS_ROOT})
    SET(OCTOPUS_CPP src/rtc_cacao.cpp)
endif()

add_library(sutraWrap MODULE
                    src/telescope.cpp
                    src/tscreen.cpp
                    src/atmos.cpp
                    src/dms.cpp
                    src/source.cpp
                    src/lgs.cpp
                    src/target.cpp
                    src/sensors.cpp
                    src/wfs.cpp
                    src/wfs_sh.cpp
                    src/wfs_pyrhr.cpp
                    src/centroider.cpp
                    src/centroider_tcog.cpp
                    src/centroider_wcog.cpp
                    src/centroider_bpcog.cpp
                    src/centroider_corr.cpp
                    src/centroider_pyr.cpp
                    src/controller.cpp
                    src/controller_ls.cpp
                    src/controller_mv.cpp
                    src/controller_geo.cpp
                    src/controller_generic.cpp
                    src/controller_cured.cpp
                    src/rtc.cpp
                    src/gamora.cpp
                    src/groot.cpp
                    src/sutraWrap.cpp
                    ${BRAHMA_CPP}
                    ${OCTOPUS_CPP}
                    )
#link_directories($ENV{COMPASS_ROOT}/libsutra)

target_link_libraries(sutraWrap PRIVATE wyrm carma sutra)
set_target_properties(sutraWrap PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                      SUFFIX "${PYTHON_MODULE_EXTENSION}")

install(TARGETS sutraWrap EXPORT SutraWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python)
