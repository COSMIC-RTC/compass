cmake_minimum_required(VERSION 3.18)

project(carma VERSION ${VERSION_INFO} DESCRIPTION "COMPASS CARMA lirary" LANGUAGES CUDA CXX)

find_package(PkgConfig REQUIRED)
find_package(CUDAToolkit REQUIRED)

pkg_check_modules(MAGMA   IMPORTED_TARGET magma)
pkg_check_modules(OCTOPUS IMPORTED_TARGET octopus)

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(GCOV gcov)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")

add_library(carma
    src.cpp/carma_context.cpp
    src.cpp/carma_cublas.cpp
    src.cpp/carma_cusolver.cpp
    src.cpp/carma_cusparse.cpp
    src.cpp/carma_fft.cpp
    src.cpp/carma_fft_conv.cpp
    src.cpp/carma_host_obj.cpp
    src.cpp/carma_ipcs.cpp
    src.cpp/carma_ksparse.cpp
    src.cpp/carma_magma.cpp
    src.cpp/carma_multithread.cpp
    src.cpp/carma_sparse_host_obj.cpp
    src.cpp/carma_sparse_obj.cpp
    src.cpp/carma_streams.cpp
    src.cpp/carma_utils.cpp
    src.cpp/carma_obj.cpp
    src.cu/carma_fft_conv.cu
    src.cu/carma_timer.cu
    src.cu/carma_rng.cu
    src.cu/carma_svd.cu
    src.cu/carma_utils.cu
    src.cu/convolutionFFT2D.cu
    src.cu/carma_obj.cu
    src.cu/carma_sum.cu
    src.cu/carma_transpose.cu
)

target_include_directories(carma PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include.h
    $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>
)

set_target_properties(carma PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_features(carma PUBLIC cxx_std_14)

target_compile_options(carma PUBLIC -Wno-deprecated-declarations
    $<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>
)

target_compile_definitions(carma PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:VERSION_INFO=${VERSION_INFO} USE_UVA>
)
target_link_libraries(carma PUBLIC
    CUDA::cudart
    CUDA::cufft
    CUDA::cublas
    CUDA::curand
    CUDA::cusparse
    CUDA::cusolver
    CONAN_PKG::emu
    ${GCOV}
)

if (${CUDAToolkit_VERSION_MAJOR} LESS 11)
    target_link_libraries(carma PUBLIC CONAN_PKG::cub)
endif()

if(${OCTOPUS_FOUND})
  target_link_libraries(carma PUBLIC PkgConfig::OCTOPUS)
  target_compile_definitions(carma PUBLIC USE_OCTOPUS=1)
else()
  message(STATUS "OCTOPUS_ROOT not found, CACAO features disabled")
endif()

if(${MAGMA_FOUND})
    target_link_libraries(carma PUBLIC PkgConfig::MAGMA)
    target_compile_definitions(carma PUBLIC $<$<COMPILE_LANGUAGE:CXX>:USE_MAGMA=1>)
else()
    message(STATUS "MAGMA not found, MAGMA features disabled")
endif()


if(do_half)
  if(${CUDA_LOWEREST_SM} GREATER_EQUAL 60)
    target_compile_definitions(carma PUBLIC CAN_DO_HALF=1)
  else()
    message(FATAL_ERROR "(at least one of) Your NVIDIA card can not support half precision (SM should be >=60)")
  endif()
endif()

if(python_build)
  add_subdirectory(python_module)
endif()

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig"
    CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/carma.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/carma.pc @ONLY)

install(TARGETS carma
        EXPORT CarmaConfig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(FILES include.h/carma_context.h
              include.h/carma_cublas.h
              include.h/carma_cusolver.h
              include.h/carma_cusparse.h
              include.h/carma_exception.h
              include.h/carma_fft.h
              include.h/carma.h
              include.h/carma_host_obj.h
              include.h/carma_ipcs.h
              include.h/carma_magma.h
              include.h/carma_multithread.h
              include.h/carma_obj.h
              include.h/carma_sparse_host_obj.h
              include.h/carma_sparse_obj.h
              include.h/carma_streams.h
              include.h/carma_timer.h
              include.h/carma_utils.h
              include.h/type_list.hpp
        DESTINATION include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/carma.pc
        DESTINATION "${INSTALL_PKGCONFIG_DIR}")
