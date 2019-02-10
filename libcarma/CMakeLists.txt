cmake_minimum_required (VERSION 3.8)

if(DEFINED ENV{CUDA_SM})
    set(CUDA_SM $ENV{CUDA_SM})
else()
    message("Probing GPUs for their compute capabilities...")
    execute_process(COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_cuda_sm.sh"  OUTPUT_VARIABLE CUDA_SM)
endif()
message("Building for Compute Capability ${CUDA_SM}")

if( ${CUDA_SM} STREQUAL "ALL")
unset(CUDA_SM)
list(APPEND CUDA_SM "30" "35" "37" "50" "52" "60" "61" "70")
endif()

foreach(SM ${CUDA_SM})
list( APPEND CUDA_ARCH -gencode arch=compute_${SM},code=sm_${SM})
endforeach(SM)

#string(APPEND CMAKE_CUDA_FLAGS ${CUDA_ARCH})

project(carma
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS CARMA lirary"
    LANGUAGES CUDA CXX
)

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

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

link_directories( ${MAGMA_LIBRARY_DIRS} )

add_library(carma SHARED
    ./src.cpp/carma_context.cpp
    ./src.cpp/carma_cublas.cpp
    ./src.cpp/carma_cula.cpp
    ./src.cpp/carma_cusparse.cpp
    ./src.cpp/carma_fft.cpp
    ./src.cpp/carma_fft_conv.cpp
    ./src.cpp/carma_host_obj.cpp
    ./src.cpp/carma_ipcs.cpp
    ./src.cpp/carma_ksparse.cpp
    ./src.cpp/carma_magma.cpp
    ./src.cpp/carma_multithread.cpp
    ./src.cpp/carma_sparse_host_obj.cpp
    ./src.cpp/carma_sparse_obj.cpp
    ./src.cpp/carma_streams.cpp
    ./src.cpp/carma_utils.cpp
    ./src.cpp/carma_obj.cpp
    ./src.cu/carma_fft_conv.cu
    ./src.cu/carma_timer.cu
    ./src.cu/carma_rng.cu
    ./src.cu/carma_svd.cu
    ./src.cu/carma_utils.cu
    ./src.cu/convolutionFFT2D.cu
    ./src.cu/carma_obj.cu
    ./src.cu/carma_sum.cu
    ./src.cu/carma_transpose.cu
)

set_target_properties(carma
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    # CUDA_SEPARABLE_COMPILATION OFF
)
list(LENGTH CUDA_SM CUDA_SM_LEN)
math(EXPR CUDA_SM_LEN ${CUDA_SM_LEN}-1)
list(GET CUDA_SM ${CUDA_SM_LEN} CUDA_SM_LAST)
target_compile_definitions(carma PUBLIC CUDA_HIGHEST_SM=${CUDA_SM_LAST})
target_compile_features(carma PUBLIC cxx_std_14)

# string (REPLACE ";" " " CMAKE_CUDA_FLAGS "${CUDA_NVCC_FLAGS};${CUDA_ARCH}")

target_include_directories(carma
    PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${MAGMA_INCLUDE_DIRS}>
    $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${CMAKE_CURRENT_SOURCE_DIR}/include.cuh $ENV{COMPASS_ROOT}/tplib/cub>
)

target_compile_options(carma
    PUBLIC
    $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_ARCH} -use_fast_math >
    PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${MAGMA_CFLAGS_OTHER}>
)

target_compile_definitions(carma
    PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:VERSION_INFO=${VERSION_INFO}> #USE_UVA>
    PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:USE_MAGMA>
)

find_package(CUDA REQUIRED)
target_link_libraries(carma
    PUBLIC ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_curand_LIBRARY} ${CUDA_cusparse_LIBRARY}
    PRIVATE ${MAGMA_LIBRARIES}
)

add_subdirectory(python_module)

set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig" CACHE PATH "Installation directory for pkgconfig (.pc) files")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/carma.pc.in ${CMAKE_CURRENT_BINARY_DIR}/carma.pc @ONLY)

install(TARGETS carma
    PUBLIC_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib}
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/carma.pc DESTINATION "${INSTALL_PKGCONFIG_DIR}")
