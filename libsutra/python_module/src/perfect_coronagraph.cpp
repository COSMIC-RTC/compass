// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      perfect_coronagraph.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraPerfectCoronagraph
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <wyrm>

#include <sutra_perfectCoronagraph.h>

namespace py = pybind11;

std::unique_ptr<SutraPerfectCoronagraph> perfect_coronagraph_init(CarmaContext &context, SutraSource *d_source,int im_dimx, 
                                int im_dimy, float *wavelength, int nWavelength, int device) {
  return std::unique_ptr<SutraPerfectCoronagraph>(new SutraPerfectCoronagraph(&context, d_source, 
                                                    im_dimx, im_dimy, wavelength, 
                                                    nWavelength, device));
};

void declare_perfect_coronagraph(py::module &mod) {
  py::class_<SutraPerfectCoronagraph, SutraCoronagraph>(mod, "PerfectCoronagraph")

      .def(py::init(wy::colCast(perfect_coronagraph_init)), R"pbdoc(
        Instantiates a PerfectCoronagraph object

        Args:
            context: (CarmaContext): context

            d_source: (SutraSource): Coronagraph source input

            im_dimx: (int): Coronagraphic image dimension along x axis

            im_dimy: (int): Coronagraphic image dimension along y axis

            wavelength: (np.ndarray[ndim=1, dtype=np.float32]): vector of wavelengths

            nWavelength: (int): number of wavelength

            device: (int): GPU device index
      )pbdoc",
      py::arg("context"), py::arg("d_source"), py::arg("im_dimx"), py::arg("im_dimy"), 
      py::arg("wavelength"), py::arg("nWavelength"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "AA", [](SutraPerfectCoronagraph &sc) { return sc.AA; }, "A MFT matrix")
      .def_property_readonly(
          "BB", [](SutraPerfectCoronagraph &sc) { return sc.BB; }, "B MFT matrix")
      .def_property_readonly(
          "norm", [](SutraPerfectCoronagraph &sc) { return sc.norm; }, "MFT normalization")
      .def_property_readonly(
          "AA_psf", [](SutraPerfectCoronagraph &sc) { return sc.AA_psf; }, "A MFT matrix")
      .def_property_readonly(
          "BB_psf", [](SutraPerfectCoronagraph &sc) { return sc.BB_psf; }, "B MFT matrix")
      .def_property_readonly(
          "norm_psf", [](SutraPerfectCoronagraph &sc) { return sc.norm_psf; }, "MFT normalization")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("set_mft", wy::colCast(&SutraPerfectCoronagraph::set_mft),
           R"pbdoc(
    Set MFT matrices for coronagraphic image computation

    Args:
        A : (np.ndarray[dtype=np.complex32, ndims=3]): A MFT matrix for each wavelength
        B : (np.ndarray[dtype=np.complex32, ndims=3]): B MFT matrix for each wavelength
        norm : (np.ndarray[dtype=np.complex32, ndims=3]): MFT normalization for each wavelength
    )pbdoc",
           py::arg("A"), py::arg("B"), py::arg("norm"))

      .def("set_mft_psf", wy::colCast(&SutraPerfectCoronagraph::set_mft_psf),
           R"pbdoc(
    Set MFT matrices for PSF computation

    Args:
        A : (np.ndarray[dtype=np.complex32, ndims=3]): A MFT matrix for each wavelength
        B : (np.ndarray[dtype=np.complex32, ndims=3]): B MFT matrix for each wavelength
        norm : (np.ndarray[dtype=np.complex32, ndims=3]): MFT normalization for each wavelength
    )pbdoc",
           py::arg("A"), py::arg("B"), py::arg("norm"));

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
};
