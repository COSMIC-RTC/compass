cmake_minimum_required(VERSION 3.18)

project(libcarma VERSION ${VERSION_INFO} DESCRIPTION "COMPASS CARMA lirary" LANGUAGES CUDA CXX)

find_package(PkgConfig REQUIRED)
find_package(CUDAToolkit REQUIRED)

add_library(libcarma
    src.cpp/carma_context.cpp
    src.cpp/carma_cublas.cpp
    src.cpp/carma_cusolver.cpp
    src.cpp/carma_cusparse.cpp
    src.cpp/carma_fft.cpp
    # src.cpp/carma_fft_conv.cpp
    src.cpp/carma_host_obj.cpp
    src.cpp/carma_sparse_host_obj.cpp
    src.cpp/carma_sparse_obj.cpp
    src.cpp/carma_streams.cpp
    src.cpp/carma_utils.cpp
    src.cpp/carma_obj.cpp
    # src.cu/carma_fft_conv.cu
    src.cu/carma_timer.cu
    src.cu/carma_rng.cu
    src.cu/carma_svd.cu
    src.cu/carma_utils.cu
    # src.cu/convolutionFFT2D.cu
    src.cu/carma_obj.cu
    src.cu/carma_sum.cu
    src.cu/carma_transpose.cu
)

target_include_directories(libcarma PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>
)

set_target_properties(libcarma PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_features(libcarma PUBLIC cxx_std_14)

target_compile_options(libcarma PUBLIC -Wno-deprecated-declarations
    $<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>
)

target_compile_definitions(libcarma PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:VERSION_INFO=${VERSION_INFO} USE_P2P_ACCESS>
)

target_link_libraries(libcarma PUBLIC
    CUDA::cudart
    CUDA::cufft
    CUDA::cublas
    CUDA::curand
    CUDA::cusparse
    CUDA::cusolver
    ${GCOV}
)

if (${CUDAToolkit_VERSION_MAJOR} LESS 11)
    find_package(cub REQUIRED)
    target_link_libraries(libcarma PUBLIC cub::cub)
endif()

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig"
    CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libcarma.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/libcarma.pc @ONLY)

install(TARGETS libcarma
        EXPORT CarmaConfig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(DIRECTORY include DESTINATION include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcarma.pc
        DESTINATION "${INSTALL_PKGCONFIG_DIR}")
