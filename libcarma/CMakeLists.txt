cmake_minimum_required (VERSION 3.5)
project(carma 
    VERSION 1.0 
    DESCRIPTION "COMPASS CARMA lirary"
    LANGUAGES CXX CUDA)
include(GNUInstallDirs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/include.cuh)

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

# add_compile_options("-fPIC -Wall -fno-common -std=c++11 ")

find_package( CUDA REQUIRED)

# nvcc flags                                                                                                                                                                                               
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -Xcompiler -fPIC -use_fast_math -std=c++11 -gencode $ENV{GENCODE})

find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma)

message(STATUS ${MAGMA_LIBRARY_DIRS})
link_directories( ${MAGMA_LIBRARY_DIRS} )

cuda_compile(CU_O 
    ./src.cu/carma_fft_conv.cu
    ./src.cu/carma_rng.cu
    ./src.cu/carma_svd.cu
    ./src.cu/carma_utils.cu
    ./src.cu/convolutionFFT2D.cu
    ./src.cu/carma_obj.cu
    ./src.cu/carma_sum.cu
    ./src.cu/carma_transpose.cu
)

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
    ${CU_O}
    )    

set_target_properties(carma PROPERTIES 
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER include.h/carma.h
    PUBLIC_HEADER include.h/carma_context.h
    PUBLIC_HEADER include.h/carma_cublas.h
    PUBLIC_HEADER include.h/carma_cusparse.h
    PUBLIC_HEADER include.h/carma_exception.h
    PUBLIC_HEADER include.h/carma_fft.h
    PUBLIC_HEADER include.h/carma_host_obj.h
    PUBLIC_HEADER include.h/carma_multithread.h
    PUBLIC_HEADER include.h/carma_sparse_host_obj.h
    PUBLIC_HEADER include.h/carma_sparse_obj.h
    PUBLIC_HEADER include.h/carma_streams.h
    PUBLIC_HEADER include.h/carma_timer.h
    PUBLIC_HEADER include.h/carma_utils.h
    PUBLIC_HEADER include.h/convolutionFFT2D_common.h
    PUBLIC_HEADER include.h/carma_ipcs.h
    PUBLIC_HEADER include.h/carma_obj.h
    PUBLIC_HEADER include.cuh/carma_utils.cuh
    PUBLIC_HEADER include.cuh/convolutionFFT2D.cuh
    PUBLIC_HEADER include.cuh/carma_utils.cuh
    PUBLIC_HEADER include.cuh/convolutionFFT2D.cuh 
)    

target_include_directories(carma PRIVATE .)
target_link_libraries(carma PRIVATE ${MAGMA_LIBRARIES})
target_include_directories(carma PUBLIC ${MAGMA_INCLUDE_DIRS})
target_compile_options(carma PUBLIC -DUSE_MAGMA ${MAGMA_CFLAGS_OTHER})

configure_file(carma.pc.in carma.pc @ONLY)

install(TARGETS carma
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)    

install(FILES ${CMAKE_BINARY_DIR}/carma.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)