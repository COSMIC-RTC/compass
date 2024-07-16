cmake_minimum_required (VERSION 3.18)

project(pycarma
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS carma library"
    LANGUAGES CXX CUDA
)

find_package(PkgConfig REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG)
find_package(CUDAToolkit REQUIRED)

pybind11_add_module(carma MODULE
    src/context.cpp
    src/timer.cpp
    src/host_obj.cpp
    src/obj.cpp
    src/sparse_obj.cpp
    src/carma.cpp
)

target_compile_features(carma PUBLIC cxx_std_14)

if (NOT TARGET libcarma)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CARMA REQUIRED IMPORTED_TARGET libcarma)
    message(STATUS "Found carma: ${CARMA_VERSION}")
    target_include_directories(carma PRIVATE ${CARMA_INCLUDE_DIRS})
    target_link_libraries(carma PUBLIC 
        PkgConfig::CARMA
        CUDA::cudart
        CUDA::cufft
        CUDA::cublas
        CUDA::curand
        CUDA::cusparse
        CUDA::cusolver
    )
else()
    message(STATUS "Found libcarma target in the project")
    target_link_libraries(carma PRIVATE libcarma)
endif()



install(TARGETS carma EXPORT carmaConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python
)
