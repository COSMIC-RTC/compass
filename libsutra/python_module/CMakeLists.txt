cmake_minimum_required (VERSION 3.18)

project(sutraWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS sutraWrap library"
    LANGUAGES CXX CUDA
)

find_package(pybind11 REQUIRED)
find_package(wyrm REQUIRED)

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
    src/coronagraph.cpp
    src/perfect_coronagraph.cpp
    src/stellar_coronagraph.cpp
    src/sutraWrap.cpp
    ${BRAHMA_CPP}
)

if (${OCTOPUS_FOUND})
    target_sources(sutraWrap PRIVATE src/rtc_cacao.cpp)
endif()

target_link_libraries(sutraWrap PRIVATE wyrm::wyrm sutra)

install(TARGETS sutraWrap EXPORT SutraWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python
)
