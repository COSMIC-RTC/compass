cmake_minimum_required (VERSION 3.8)

project(sutra
    VERSION 1.0
    DESCRIPTION "COMPASS SUTRA lirary"
    LANGUAGES CXX CUDA
)

include(GNUInstallDirs)

# EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
# message( STATUS "Architecture: ${ARCHITECTURE}" )
# ifeq ($(OS_ARCH),armv7l)
#   ifeq ($(abi),gnueabi)
#     CCFLAGS += -mfloat-abi=softfp
#   else
#     # default to gnueabihf
#     override abi := gnueabihf
#     #      LDFLAGSS += --dynamic-linker=/lib/ld-linux-armhf.so.3
#     CCFLAGS   += -marm -mfpu=vfpv3  -mfloat-abi=hard -march=armv7-a
#     NVCCFLAGS := -target-cpu-arch ARM -Xlinker --dynamic-linker=/lib/ld-linux-armhf.so.3 -Xcompiler "$(CCFLAGS)"
#   endif
# else
#   ifeq ($(OS_ARCH),aarch64)
#     # CCFLAGS += -march=armv8-a+crc -mtune=cortex-a50  # Armageddon
#     CCFLAGS += -march=armv8-a+crc+crypto+fp+simd -mtune=cortex-a57  # TX2
#   else
#     CCFLAGS += -m64
#   endif
# endif

find_package(CUDA)

add_library(sutra SHARED
    src.cpp/sutra_acquisim.cpp
    src.cpp/sutra_aotemplate.cpp
    src.cpp/sutra_atmos.cpp
    src.cpp/sutra_centroider_bpcog.cpp
    src.cpp/sutra_centroider_cog.cpp
    src.cpp/sutra_centroider_corr.cpp
    src.cpp/sutra_centroider.cpp
    src.cpp/sutra_centroider_pyr.cpp
    src.cpp/sutra_centroider_roof.cpp
    src.cpp/sutra_centroider_tcog.cpp
    src.cpp/sutra_centroider_wcog.cpp
    src.cpp/sutra_controller.cpp
    src.cpp/sutra_controller_cured.cpp
    src.cpp/sutra_controller_generic.cpp
    src.cpp/sutra_controller_geo.cpp
    src.cpp/sutra_controller_kalman.cpp
    src.cpp/sutra_controller_ls.cpp
    src.cpp/sutra_controller_mv.cpp
    src.cpp/sutra_dm.cpp
    src.cpp/sutra_gamora.cpp
    src.cpp/sutra_groot.cpp
    src.cpp/sutra_kl.cpp
    src.cpp/sutra_lgs.cpp
    src.cpp/sutra_phase.cpp
    src.cpp/sutra_roket.cpp
    # src.cpp/sutra_rtc_brama.cpp
    # src.cpp/sutra_rtc_bramaListenerImpl.cpp
    src.cpp/sutra_rtc.cpp
    # src.cpp/sutra_target_brama.cpp
    # src.cpp/sutra_target_bramaListenerImpl.cpp
    src.cpp/sutra_target.cpp
    src.cpp/sutra_telemetry.cpp
    src.cpp/sutra_telescope.cpp
    src.cpp/sutra_tscreen.cpp
    src.cpp/sutra_wfs.cpp
    src.cpp/sutra_wfs_geom.cpp
    src.cpp/sutra_wfs_pyr.cpp
    src.cpp/sutra_wfs_pyr_pyr4.cpp
    src.cpp/sutra_wfs_pyr_pyrhr.cpp
    src.cpp/sutra_wfs_pyr_roof.cpp
    src.cpp/sutra_wfs_sh.cpp
    src.cu/sutra_acquisim.cu
    src.cu/sutra_aotemplate.cu
    src.cu/sutra_ao_utils.cu
    src.cu/sutra_centroider.cu
    src.cu/sutra_controller.cu
    src.cu/sutra_controller_utils.cu
    src.cu/sutra_dm.cu
    src.cu/sutra_gamora.cu
    src.cu/sutra_groot.cu
    src.cu/sutra_kl.cu
    src.cu/sutra_lgs.cu
    src.cu/sutra_roket.cu
    src.cu/sutra_target.cu
    src.cu/sutra_tscreen.cu
    src.cu/sutra_wfs.cu
    aao/cured/cured.cpp
)

set_target_properties(sutra
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    # CUDA_SEPARABLE_COMPILATION ON
)

target_include_directories(sutra
    PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${CMAKE_CURRENT_SOURCE_DIR}/aao/cured>
    $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>
)

target_compile_options(sutra
    PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fPIC -use_fast_math -std=c++11 -gencode $ENV{GENCODE}>
)

target_link_libraries(sutra PRIVATE carma ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_curand_LIBRARY} ${CUDA_cusparse_LIBRARY} )

configure_file(sutra.pc.in sutra.pc @ONLY)

install(TARGETS sutra
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CMAKE_BINARY_DIR}/sutra.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
