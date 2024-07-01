cmake_minimum_required (VERSION 3.18)

project(carmaWrap
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carmaWrap library"
    LANGUAGES CXX CUDA
)

find_package(PkgConfig REQUIRED)
find_package(pybind11 REQUIRED)
find_package(CUDAToolkit REQUIRED)
pkg_check_modules(OCTOPUS IMPORTED_TARGET octopus)
link_directories(${OCTOPUS_LIBRARY_DIRS})

pybind11_add_module(carmaWrap MODULE
    src/context.cpp
    src/timer.cpp
    src/host_obj.cpp
    src/obj.cpp
    src/sparse_obj.cpp
    src/carmaWrap.cpp
)

target_compile_features(carmaWrap PUBLIC cxx_std_14)

if (NOT TARGET carma)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CARMA REQUIRED IMPORTED_TARGET carma)
    message(STATUS "Found carma: ${CARMA_VERSION}")
    target_include_directories(carmaWrap PRIVATE ${CARMA_INCLUDE_DIRS})
    target_link_libraries(carmaWrap PUBLIC 
        PkgConfig::CARMA
        CUDA::cudart
        CUDA::cufft
        CUDA::cublas
        CUDA::curand
        CUDA::cusparse
        CUDA::cusolver
    )
else()
    message(STATUS "Found carma target in the project")
    target_link_libraries(carmaWrap PRIVATE carma)
endif()



install(TARGETS carmaWrap EXPORT CarmaWrapConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python
)
