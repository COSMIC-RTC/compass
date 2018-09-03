cmake_minimum_required (VERSION 3.8)

project(carma
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS CARMA lirary"
    LANGUAGES CXX CUDA
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

find_package(CUDA)

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
    ./src.cpp/carma_rng.cpp
    ./src.cpp/carma_sparse_host_obj.cpp
    ./src.cpp/carma_sparse_obj.cpp
    ./src.cpp/carma_streams.cpp
    ./src.cpp/carma_utils.cpp
    ./src.cpp/carma_obj.cpp
    ./src.cu/carma_fft_conv.cu
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
    CUDA_SEPARABLE_COMPILATION ON
)
target_compile_features(carma PUBLIC cxx_std_11)

target_include_directories(carma
    PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${MAGMA_INCLUDE_DIRS}>
    $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CURRENT_SOURCE_DIR}/include.h ${CMAKE_CURRENT_SOURCE_DIR}/include.cuh>
)


target_compile_options(carma
    PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:-DVERSION_INFO=${VERSION_INFO}>
    PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-DUSE_MAGMA ${MAGMA_CFLAGS_OTHER}>
    $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_SM_FLAG} -Xcompiler -fPIC -use_fast_math >
)

target_link_libraries(carma PRIVATE ${MAGMA_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_curand_LIBRARY} ${CUDA_cusparse_LIBRARY} )

include(GNUInstallDirs)
configure_file(carma.pc.in carma.pc @ONLY)

install(TARGETS carma
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CMAKE_BINARY_DIR}/carma.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
