cmake_minimum_required (VERSION 3.8)

project(compass
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS lirary"
    LANGUAGES CUDA CXX
)

if( ${CUDA_SM} STREQUAL "ALL")
    set(CUDA_SM "30" "35" "37" "50" "52" "60" "61" "70")
endif()

# set(CMAKE_CUDA_FLAGS ${CUDA_SM_FLAG} )

foreach(SM ${CUDA_SM})
    list( APPEND CUDA_SM_FLAG "-gencode=arch=compute_${SM},code=sm_${SM}")
endforeach(SM)
string (REPLACE ";" " " CMAKE_CUDA_FLAGS "${CUDA_SM_FLAG}")

add_subdirectory(libcarma)
add_subdirectory(naga)
add_subdirectory(libsutra)
add_subdirectory(shesha)
