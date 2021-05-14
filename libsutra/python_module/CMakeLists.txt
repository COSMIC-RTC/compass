cmake_minimum_required (VERSION 3.8)

project(sutraWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS sutraWrap library"
    LANGUAGES CXX CUDA
)

find_package(pybind11 REQUIRED)

# cmake-format: off
# if (DEFINED ENV{BRAHMA_ROOT})
#     SET(BRAHMA_DIRS $ENV{ACE_ROOT}/lib $ENV{TAO_ROOT}/tao $ENV{DDS_ROOT}/lib $ENV{BRAHMA_ROOT})
#     link_directories( ${BRAHMA_DIRS} )
#     SET(BRAHMA_CPP src/rtc_brahma.cpp
#                 src/target_brahma.cpp)
# endif()
# cmake-format: on

if (${OCTOPUS_FOUND})
    SET(OCTOPUS_CPP src/rtc_cacao.cpp)
endif()

pybind11_add_module(sutraWrap MODULE
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
                    src/centroider_maskedPix.cpp
                    src/controller.cpp
                    src/controller_ls.cpp
                    src/controller_mv.cpp
                    src/controller_geo.cpp
                    src/controller_generic.cpp
                    src/controller_generic_linear.cpp
                    src/controller_cured.cpp
                    src/rtc.cpp
                    src/gamora.cpp
                    src/groot.cpp
                    src/sutraWrap.cpp
                    ${BRAHMA_CPP}
                    ${OCTOPUS_CPP}
                    )
#link_directories($ENV{COMPASS_ROOT}/libsutra)

target_link_libraries(sutraWrap PRIVATE CONAN_PKG::wyrm carma sutra)

install(TARGETS sutraWrap EXPORT SutraWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python)
