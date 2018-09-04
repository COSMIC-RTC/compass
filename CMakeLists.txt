cmake_minimum_required (VERSION 3.8)

project(compass
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS lirary"
    LANGUAGES CUDA CXX
)


if(DEFINED ENV{CUDA_SM})
    set(CUDA_SM $ENV{CUDA_SM})
else()
    message("Probing the (first) GPU for its compute capability...")
    execute_process(COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_cuda_sm.sh"  RESULT_VARIABLE CUDA_SM)
endif()
message("Building for Compute Capability ${CUDA_SM}")

if( ${CUDA_SM} STREQUAL "ALL")
    set(CUDA_SM "35" "37" "50" "52" "60" "61" "70")
endif()

foreach(SM ${CUDA_SM})
list( APPEND CUDA_ARCH "-gencode" "arch=compute_${SM},code=sm_${SM}")
endforeach(SM)

string (REPLACE ";" " " CMAKE_CUDA_FLAGS "${CUDA_NVCC_FLAGS};${CUDA_ARCH}")
# set(CMAKE_CUDA_FLAGS ${CUDA_ARCH} VERBATIM)

# string (REPLACE ";" " " CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS};${CUDA_ARCH}")
# set(CMAKE_SHARED_LINKER_FLAGS ${CUDA_ARCH} VERBATIM)

add_subdirectory(libcarma)
add_subdirectory(naga)
add_subdirectory(libsutra)
add_subdirectory(shesha)
