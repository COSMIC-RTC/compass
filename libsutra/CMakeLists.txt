cmake_minimum_required(VERSION 3.18)

project(libsutra VERSION ${VERSION_INFO} DESCRIPTION "COMPASS SUTRA lirary" LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(GCOV gcov)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")

add_library(libsutra
    src.cpp/sutra_acquisim.cpp
    src.cpp/sutra_atmos.cpp
    src.cpp/sutra_centroider_bpcog.cpp
    src.cpp/sutra_centroider_cog.cpp
    src.cpp/sutra_centroider_corr.cpp
    src.cpp/sutra_centroider.cpp
    src.cpp/sutra_centroider_pyr.cpp
    src.cpp/sutra_centroider_maskedPix.cpp
    src.cpp/sutra_centroider_tcog.cpp
    src.cpp/sutra_centroider_wcog.cpp
    src.cpp/sutra_controller.cpp
    src.cpp/sutra_controller_generic.cpp
    src.cpp/sutra_controller_generic_linear.cpp
    src.cpp/sutra_controller_geo.cpp
    src.cpp/sutra_controller_ls.cpp
    src.cpp/sutra_controller_mv.cpp
    src.cpp/sutra_dm.cpp
    src.cpp/sutra_gamora.cpp
    src.cpp/sutra_groot.cpp
    src.cpp/sutra_kl.cpp
    src.cpp/sutra_lgs.cpp
    src.cpp/sutra_phase.cpp
    # src.cpp/sutra_roket.cpp
    src.cpp/sutra_rtc.cpp
    src.cpp/sutra_source.cpp
    src.cpp/sutra_target.cpp
    src.cpp/sutra_telemetry.cpp
    src.cpp/sutra_telescope.cpp
    src.cpp/sutra_tscreen.cpp
    src.cpp/sutra_sensors.cpp
    src.cpp/sutra_wfs.cpp
    # src.cpp/sutra_wfs_geom.cpp
    src.cpp/sutra_wfs_pyr_pyrhr.cpp
    src.cpp/sutra_wfs_sh.cpp
    src.cpp/sutra_coronagraph.cpp
    src.cpp/sutra_perfectCoronagraph.cpp
    src.cpp/sutra_stellarCoronagraph.cpp
    src.cu/sutra_acquisim.cu
    src.cu/sutra_utils.cu
    src.cu/sutra_centroider.cu
    src.cu/sutra_centroider_cog.cu
    src.cu/sutra_centroider_tcog.cu
    src.cu/sutra_centroider_wcog.cu
    src.cu/sutra_centroider_bpcog.cu
    src.cu/sutra_centroider_corr.cu
    src.cu/sutra_centroider_pyr.cu
    src.cu/sutra_centroider_maskedPix.cu
    src.cu/sutra_controller.cu
    src.cu/sutra_controller_utils.cu
    src.cu/sutra_dm.cu
    src.cu/sutra_gamora.cu
    src.cu/sutra_groot.cu
    src.cu/sutra_kl.cu
    src.cu/sutra_lgs.cu
    # src.cu/sutra_brightest_pix.cu
    # src.cu/sutra_roket.cu
    src.cu/sutra_target.cu
    src.cu/sutra_tscreen.cu
    src.cu/sutra_wfs.cu
    src.cu/sutra_coronagraph.cu
)

target_include_directories(libsutra
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>
        $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/src.cu>
)


set_target_properties(libsutra PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_features(libsutra PUBLIC cxx_std_14)

target_link_libraries(libsutra PUBLIC
    libcarma
    CUDA::cudart
    ${GCOV}
)

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig"
    CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libsutra.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/libsutra.pc @ONLY)

install(TARGETS libsutra
        EXPORT SutraConfig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(DIRECTORY include DESTINATION include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libsutra.pc
        DESTINATION "${INSTALL_PKGCONFIG_DIR}")
