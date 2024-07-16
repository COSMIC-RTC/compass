cmake_minimum_required (VERSION 3.18)

project(pysutra
    VERSION ${VERSION_INFO}
    DESCRIPTION "COMPASS sutra library"
    LANGUAGES CXX CUDA
)

find_package(PkgConfig REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG)
find_package(CUDAToolkit REQUIRED)

pybind11_add_module(sutra MODULE
    src/telescope.cpp
    src/tscreen.cpp
    src/atmos.cpp
    src/dms.cpp
    src/source.cpp
    src/lgs.cpp
    src/target.cpp
    src/sensors.cpp
    src/wfs.cpp
    src/wfs_sh.cpp
    src/wfs_pyrhr.cpp
    src/centroider.cpp
    src/centroider_tcog.cpp
    src/centroider_wcog.cpp
    src/centroider_bpcog.cpp
    src/centroider_corr.cpp
    src/centroider_pyr.cpp
    src/centroider_maskedPix.cpp
    src/controller.cpp
    src/controller_ls.cpp
    src/controller_mv.cpp
    src/controller_geo.cpp
    src/controller_generic.cpp
    src/controller_generic_linear.cpp
    src/rtc.cpp
    src/gamora.cpp
    src/groot.cpp
    src/coronagraph.cpp
    src/perfect_coronagraph.cpp
    src/stellar_coronagraph.cpp
    src/sutra.cpp
)

target_compile_features(sutra PUBLIC cxx_std_14)

if (NOT TARGET libsutra)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(SUTRA REQUIRED IMPORTED_TARGET libsutra)
    message(STATUS "Found libsutra: ${SUTRA_VERSION}")
    target_include_directories(sutra PRIVATE ${SUTRA_INCLUDE_DIRS})
    target_link_libraries(sutra PRIVATE 
        PkgConfig::SUTRA
        CUDA::cudart
        CUDA::cufft
        CUDA::cublas
        CUDA::curand
        CUDA::cusparse
        CUDA::cusolver
    )
else()
    message(STATUS "Found libsutra target in the project")
    target_link_libraries(sutra PRIVATE libsutra)
endif()



install(TARGETS sutra EXPORT sutraConfig
    ARCHIVE  DESTINATION python
    LIBRARY  DESTINATION python
)
