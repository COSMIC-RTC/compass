cmake_minimum_required(VERSION 3.18)

project(carma VERSION ${VERSION_INFO} DESCRIPTION "COMPASS CARMA lirary" LANGUAGES CUDA CXX)

find_package(PkgConfig REQUIRED)
find_package(CUDAToolkit REQUIRED)

pkg_check_modules(OCTOPUS IMPORTED_TARGET octopus)
link_directories(${OCTOPUS_LIBRARY_DIRS})

add_library(carma
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

target_include_directories(carma PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
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
    $<$<COMPILE_LANGUAGE:CXX>:VERSION_INFO=${VERSION_INFO} USE_P2P_ACCESS>
)

target_link_libraries(carma PUBLIC
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
    target_link_libraries(carma PUBLIC cub::cub)
endif()

if(${OCTOPUS_FOUND})
  target_compile_definitions(carma PUBLIC USE_OCTOPUS=1)
  target_compile_options(carma PUBLIC ${OCTOPUS_CFLAGS})
  target_include_directories(carma PUBLIC ${OCTOPUS_INCLUDE_DIRS})
  target_link_libraries(carma PUBLIC ${OCTOPUS_LIBRARIES})
else()
  message(STATUS "OCTOPUS_ROOT not found, CACAO features disabled")
endif()

if(do_half)
    message(STATUS "COMPASS half precision enabled")
    target_compile_definitions(carma PUBLIC CAN_DO_HALF=1)
else()
    message(STATUS "COMPASS half precision disabled")
endif()

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig"
    CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/carma.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/carma.pc @ONLY)

install(TARGETS carma
        EXPORT CarmaConfig
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(DIRECTORY include DESTINATION include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/carma.pc
        DESTINATION "${INSTALL_PKGCONFIG_DIR}")
