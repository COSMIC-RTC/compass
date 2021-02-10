cmake_minimum_required(VERSION 3.8)

project(sutra
        VERSION ${VERSION_INFO} DESCRIPTION "COMPASS SUTRA lirary"
        LANGUAGES CXX CUDA)

# cmake-format: off
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
# cmake-format: on

find_package(PkgConfig REQUIRED)

#pkg_check_modules(MAGMA magma)
if(${MAGMA_FOUND})
  link_directories(${MAGMA_LIBRARY_DIRS})
else()
message(STATUS "MAGMA not found, MAGMA features disabled")
endif()

#pkg_check_modules(OCTOPUS octopus)
if(${OCTOPUS_FOUND})
  set(OCTOPUS_CPP src.cpp/sutra_rtc_cacao.cpp)
  set(OCTOPUS_H include.h/sutra_rtc_cacao.h)
  link_directories(${OCTOPUS_LIBRARY_DIRS})
else()
  message(STATUS "OCTOPUS_ROOT not find, CACAO RTC will not be available")
endif()

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(GCOV gcov)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")

# cmake-format: off
# if(DEFINED ENV{BRAHMA_ROOT})
#   set(BRAHMA_FLAGS
#       -fvisibility-inlines-hidden
#       -fno-strict-aliasing
#       -pthread
#       -pipe
#       -Wno-deprecated-declarations
#       -D_GNU_SOURCE
#       -D__ACE_INLINE__
#       -DUSE_BRAHMA)

#   set(BRAHMA_DIRS
#       $ENV{ACE_ROOT}/lib
#       $ENV{TAO_ROOT}/tao
#       $ENV{DDS_ROOT}/lib
#       $ENV{BRAHMA_ROOT})
#   link_directories(${BRAHMA_DIRS})

#   set(BRAHMA_INC
#       $ENV{BRAHMA_ROOT}
#       $ENV{ACE_ROOT}
#       $ENV{TAO_ROOT}
#       $ENV{TAO_ROOT}/orbsvcs
#       $ENV{DDS_ROOT})

#   set(BRAHMA_CPP src.cpp/SutraRtc_brahma.cpp
#       # src.cpp/SutraRtcBrahmaListenerImpl.cpp
#       src.cpp/sutra_target_brahma.cpp
#       # src.cpp/sutra_target_brahmaListenerImpl.cpp
#       )

#   set(BRAHMA_H include.h/SutraRtc_brahma.h
#       # include.h/SutraRtcBrahmaListenerImpl.h
#       include.h/sutra_target_brahma.h
#       # include.h/sutra_target_brahmaListenerImpl.h
#       )

#   set(BRAHMA_LIBS
#       BRAHMACommon
#       OpenDDS_InfoRepoDiscovery
#       OpenDDS_Dcps
#       TAO_PortableServer
#       TAO_AnyTypeCode
#       TAO
#       ACE
#       dl
#       rt)
# else()
#   message(STATUS "BRAHMA_ROOT not find, BRAHMA RTC will not be available")
# endif()
# cmake-format: on

add_library(sutra
            SHARED
            src.cpp/sutra_acquisim.cpp
            src.cpp/sutra_template.cpp
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
            src.cpp/sutra_controller_cured.cpp
            src.cpp/sutra_controller_generic.cpp
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
            src.cu/sutra_acquisim.cu
            src.cu/sutra_template.cu
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
            # src.cu/sutra_roket.cu
            src.cu/sutra_target.cu
            src.cu/sutra_tscreen.cu
            src.cu/sutra_wfs.cu
            aao/cured/cured.cpp
            ${BRAHMA_CPP}
            ${OCTOPUS_CPP})

if(${OCTOPUS_FOUND})
  target_compile_options(sutra PUBLIC -DUSE_CACAO ${OCTOPUS_CFLAGS})
  target_include_directories(
    sutra
    PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${OCTOPUS_INCLUDE_DIRS}>)
  target_link_libraries(sutra PUBLIC ${OCTOPUS_LIBRARIES})
else()
  message(STATUS "OCTOPUS_ROOT not find, CACAO RTC will not be available")
endif()

set_target_properties(sutra
                      PROPERTIES VERSION
                                 ${PROJECT_VERSION}
                                 SOVERSION
                                 ${PROJECT_VERSION_MAJOR}
                                 # CUDA_SEPARABLE_COMPILATION OFF
                      )

target_include_directories(
  sutra
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include.h ${BRAHMA_INC}
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/aao/cured>
          $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.h
          ${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>)

target_compile_options(sutra
                       PUBLIC $<$<COMPILE_LANGUAGE:CXX>: ${CACAO_FLAGS}
                              ${BRAHMA_FLAGS}>)

target_link_libraries(sutra PUBLIC carma ${BRAHMA_LIBS} ${GCOV})

if(python_build)
  add_subdirectory(python_module)
endif()

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig"
    CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sutra.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/sutra.pc @ONLY)

install(TARGETS sutra
        EXPORT SutraConfig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES include.h/sutra_acquisim.h
              include.h/sutra_template.h
              include.h/sutra_atmos.h
              include.h/sutra_centroider_bpcog.h
              include.h/sutra_centroider_cog.h
              include.h/sutra_centroider_corr.h
              include.h/sutra_centroider.h
              include.h/sutra_centroider_maskedPix.h
              include.h/sutra_centroider_pyr.h
              include.h/sutra_centroider_tcog.h
              include.h/sutra_centroider_wcog.h
              include.h/sutra_controller_cured.h
              include.h/sutra_controller_generic.h
              include.h/sutra_controller_geo.h
              include.h/sutra_controller.h
              include.h/sutra_controller_ls.h
              include.h/sutra_controller_mv.h
              include.h/sutra_controller_utils.h
              include.h/sutra_dm.h
              include.h/sutra_gamora.h
              include.h/sutra_groot.h
              include.h/sutra_kl.h
              include.h/sutra_lgs.h
              include.h/sutra_phase.h
              include.h/sutra_roket.h
              include.h/sutra_rtc.h
              include.h/sutra_sensors.h
              include.h/sutra_source.h
              include.h/sutra_target.h
              include.h/sutra_telemetry.h
              include.h/sutra_telescope.h
              include.h/sutra_tscreen.h
              include.h/sutra_utils.h
              include.h/sutra_wfs_geom.h
              include.h/sutra_wfs.h
              include.h/sutra_wfs_pyr_pyrhr.h
              include.h/sutra_wfs_sh.h
              ${BRAHMA_H}
              ${OCTOPUS_H}
        DESTINATION include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sutra.pc
        DESTINATION "${INSTALL_PKGCONFIG_DIR}")
