cmake_minimum_required (VERSION 3.8)

set(ENV{CUDA_INCLUDE_DIRS} "${CUDA_INCLUDE_DIRS}")
set(ENV{CUDA_CUDART_LIBRARY} "${CUDA_CUDART_LIBRARY}")

if(DEFINED ENV{CUDA_SM})
    set(CUDA_SM $ENV{CUDA_SM})
else()
    message("Probing the (first) GPU for its compute capability...")
    execute_process(COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_cuda_sm.sh"  RESULT_VARIABLE CUDA_SM)
endif()
message("Building for Compute Capability ${CUDA_SM}")

add_subdirectory(libcarma)
add_subdirectory(naga)
# add_subdirectory(libsutra)
