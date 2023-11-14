// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      stellar_coronagraph.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraStellarCoronagraph
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <wyrm>

#include <sutra_stellarCoronagraph.h>

namespace py = pybind11;

std::unique_ptr<SutraStellarCoronagraph> stellar_coronagraph_init(CarmaContext &context,
                                SutraSource *d_source,int im_dimx, int im_dimy,
                                int fpm_dimx, int fpm_dimy,
                                float *wavelength, int nWavelength, int babinet, int device) {
  return std::unique_ptr<SutraStellarCoronagraph>(new SutraStellarCoronagraph(&context, d_source,
                                                    im_dimx, im_dimy, fpm_dimx, fpm_dimy,
                                                    wavelength, nWavelength, babinet, device));
};

void declare_stellar_coronagraph(py::module &mod) {
  py::class_<SutraStellarCoronagraph, SutraCoronagraph>(mod, "StellarCoronagraph")

      .def(py::init(wy::colCast(stellar_coronagraph_init)), R"pbdoc(
        Instantiates a StellarCoronagraph object

        Args:
            context: (CarmaContext): context

            d_source: (SutraSource): Coronagraph source input

            im_dimx: (int): Coronagraphic image dimension along x axis

            im_dimy: (int): Coronagraphic image dimension along y axis

            fpm_dimx: (int): Focal plane dimension along x axis

            fpm_dimy: (int): Focal plane dimension along y axis

            wavelength: (np.ndarray[ndim=1, dtype=np.float32]): vector of wavelengths

            nWavelength: (int): number of wavelength

            babinet: (bool): Flag to enable Babinet trick

            device: (int): GPU device index
      )pbdoc",
      py::arg("context"), py::arg("d_source"), py::arg("im_dimx"), py::arg("im_dimy"),
      py::arg("fpm_dimx"), py::arg("fpm_dimy"),py::arg("wavelength"),
      py::arg("nWavelength"), py::arg("babinet"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "babinet", [](SutraStellarCoronagraph &sc) { return sc.babinet; }, "Babinet trick flag")
      .def_property_readonly(
          "fpmDimx", [](SutraStellarCoronagraph &sc) { return sc.fpmDimx; }, "Focal plane dimension")
      .def_property_readonly(
          "fpmDimy", [](SutraStellarCoronagraph &sc) { return sc.fpmDimy; }, "Focal plane dimension")
      .def_property_readonly(
          "AA", [](SutraStellarCoronagraph &sc) { return sc.AA; }, "A MFT matrices")
      .def_property_readonly(
          "BB", [](SutraStellarCoronagraph &sc) { return sc.BB; }, "B MFT matrices")
      .def_property_readonly(
          "norm", [](SutraStellarCoronagraph &sc) { return sc.norm; }, "MFT normalization")
      .def_property_readonly(
          "d_lyot_stop", [](SutraStellarCoronagraph &sc) { return sc.d_lyot_stop; }, "Lyot stop")
      .def_property_readonly(
          "d_apodizer", [](SutraStellarCoronagraph &sc) { return sc.d_apodizer; }, "Apodizer")
      .def_property_readonly(
          "focal_plane_mask", [](SutraStellarCoronagraph &sc) { return sc.focal_plane_mask; }, "Focal plane mask")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("set_mft", wy::colCast(&SutraStellarCoronagraph::set_mft),
           R"pbdoc(
    Set MFT matrices for coronagraphic image computation

    Args:
        A : (np.ndarray[dtype=np.complex32, ndims=3]): A MFT matrix for each wavelength

        B : (np.ndarray[dtype=np.complex32, ndims=3]): B MFT matrix for each wavelength

        norm : (np.ndarray[dtype=np.complex32, ndims=3]): MFT normalization for each wavelength

        mft_type : (str): MFT matrices to set, i.e. "img" or "psf"

    )pbdoc",
           py::arg("A"), py::arg("B"), py::arg("norm"), py::arg("mft_type"))

      .def("set_focal_plane_mask", wy::colCast(&SutraStellarCoronagraph::set_focal_plane_mask),
           R"pbdoc(
    Set focal plane mask for coronagraphic image computation

    Args:
        mask: (np.ndarray[ndim=3, dtype=np.float32]): Focal plane mask for each wavelength
    )pbdoc",
           py::arg("mask"))

      .def("set_apodizer", wy::colCast(&SutraStellarCoronagraph::set_apodizer),
           R"pbdoc(
    Set apodizer for coronagraphic image computation

    Args:
        mask: (np.ndarray[ndim=2, dtype=np.float32]): apodizer
    )pbdoc",
           py::arg("mask"))

      .def("set_lyot_stop", wy::colCast(&SutraStellarCoronagraph::set_lyot_stop),
           R"pbdoc(
    Set lyot_stop for coronagraphic image computation

    Args:
        mask: (np.ndarray[ndim=2, dtype=np.float32]): lyot_stop
    )pbdoc",
           py::arg("mask"))

      .def("compute_image_normalization", &SutraStellarCoronagraph::compute_image_normalization,
           R"pbdoc(
    Compute image for normalization
        )pbdoc");

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
};
