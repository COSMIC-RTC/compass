// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      centroider.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroider
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_centroider.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_impl(py::module &mod, const char *name) {
  using centroider = SutraCentroider<Tin, Tcomp>;

  py::class_<centroider>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("context",
                             [](centroider &sc) { return sc.current_context; },
                             "GPU context")

      .def_property_readonly("device", [](centroider &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](centroider &sc) { return sc.get_type(); },
                             "Centroider type")

      .def_property_readonly("nslopes",
                             [](centroider &sc) { return sc.nslopes; },
                             "Number of slopes")

      .def_property_readonly("wfs", [](centroider &sc) { return sc.wfs; },
                             "SutraWfs handled by this centroider")

      .def_property_readonly("nvalid", [](centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("nxsub", [](centroider &sc) { return sc.nxsub; },
                             "Number of ssp across the pupil diameter")

      .def_property_readonly("npix", [](centroider &sc) { return sc.npix; },
                             "Number of pixels along a side of WFS subap.")

      .def_property_readonly("offset", [](centroider &sc) { return sc.offset; },
                             "Offset for centroiding computation")

      .def_property_readonly("scale", [](centroider &sc) { return sc.scale; },
                             "Scale factor to get slopes in arcsec")

      .def_property_readonly("nvalid", [](centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("d_bincube",
                             [](centroider &sc) { return sc.d_bincube; },
                             "Bincube of the WFS image")

      .def_property_readonly("d_intensities",
                             [](centroider &sc) { return sc.d_intensities; },
                             "intensities of the WFS image")

      .def_property_readonly("d_centroids_ref",
                             [](centroider &sc) { return sc.d_centroids_ref; },
                             "Reference slopes vector")

      .def_property_readonly("d_img", [](centroider &sc) { return sc.d_img; },
                             "Calibrated WFS image")

      .def_property_readonly("d_img_raw",
                             [](centroider &sc) { return sc.d_img_raw; },
                             "Raw WFS image")

      .def_property_readonly("d_validx",
                             [](centroider &sc) { return sc.d_validx; },
                             "X positions of the valid ssp")

      .def_property_readonly("d_validy",
                             [](centroider &sc) { return sc.d_validy; },
                             "Y positions of the valid ssp")

      .def_property_readonly("d_dark", [](centroider &sc) { return sc.d_dark; },
                             "Dark frame for calibration")

      .def_property_readonly("d_flat", [](centroider &sc) { return sc.d_flat; },
                             "Flat frame for calibration")

      .def_property_readonly("d_lutPix",
                             [](centroider &sc) { return sc.d_lutPix; },
                             "Lookup Table of pixels for calibration")

      .def_property_readonly("d_validMask",
                             [](centroider &sc) {
                               sc.get_validMask();
                               return sc.d_validMask;
                             },
                             "Flat frame for calibration")

      .def_property_readonly("filter_TT",
                             [](centroider &sc) { return sc.filter_TT; },
                             "Tip/tilt filtering flag")

      .def_property_readonly("d_ref_Tip",
                             [](centroider &sc) { return sc.d_ref_Tip; },
                             "Tip mode reference for filtering")

      .def_property_readonly("d_ref_Tilt",
                             [](centroider &sc) { return sc.d_ref_Tilt; },
                             "Tilt mode reference for filtering")

      .def_property_readonly("d_TT_slopes",
                             [](centroider &sc) { return sc.d_TT_slopes; },
                             "Tip/tilt slopes removed after filtering")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("get_cog",
           wy::colCast((int (centroider::*)(void)) & centroider::get_cog),
           "Computes centroids and stores it in d_slopes of the WFS")

      .def("load_validpos", wy::colCast(&centroider::load_validpos),
           R"pbdoc(
     Load the validx and validy arrays

     Args:
        validx: (np.array[ndim=1,dtype=np.float32]): X positions of the valid ssp

        validy: (np.array[ndim=1,dtype=np.float32]): Y positions of the valid ssp

        N: (int): arrays size
    )pbdoc",
           py::arg("validx"), py::arg("validy"), py::arg("N"))

      .def("set_npix", wy::colCast(&centroider::set_npix),
           R"pbdoc(
     Set the number of pixels per subap for a RTC standalone

     Args:
        npix: (int): number of pixel along a subap. side
    )pbdoc",
           py::arg("npix"))

      .def("set_nxsub", wy::colCast(&centroider::set_nxsub),
           R"pbdoc(
     Set the number of ssp across the pupil diameter for a RTC standalone

     Args:
        nxsub: (int): number of ssp across the pupil diameter
    )pbdoc",
           py::arg("nxsub"))

      .def("load_img",
           wy::colCast((int (centroider::*)(Tin *, int, int, int)) &
                       centroider::load_img),
           R"pbdoc(
     Load an image in a RTC standalone (host to device)

     Args:
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image

        m: (int): Image support size X

        n: (int): Image support size Y

        location: (int): If -1, image is located on the CPU (hostToDevice). Else, it is the GPU index where the image is located
    )pbdoc",
           py::arg("img"), py::arg("m"), py::arg("n"), py::arg("location"))

      .def("load_img",
           wy::colCast((int (centroider::*)(Tin *, int, int)) &
                       centroider::load_img),
           R"pbdoc(
     Load a square image (n, n) in a RTC standalone

     Args:
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image

        n: (int): Image support size along one axis

     Kwargs:
        location: (int): (optionnal) If -1 (default), image is located on the CPU (hostToDevice). Else, it is the GPU index where the image is located
    )pbdoc",
           py::arg("img"), py::arg("n"), py::arg("location") = -1)

      .def("load_img",
           wy::colCast((int (centroider::*)(CarmaObj<Tin> *)) &
                       centroider::load_img),
           R"pbdoc(
     Load an image in a RTC standalone from a CarmaObj

     Args:
        img: (CarmaObj): SH image
    )pbdoc",
           py::arg("img"))

      .def("calibrate_img", wy::colCast((int (centroider::*)(void)) &
                       centroider::calibrate_img), R"pbdoc(
           Performs the raw WFS frame calibration
           )pbdoc")

      .def("calibrate_img_validPix", wy::colCast((int (centroider::*)(void)) &
                       centroider::calibrate_img_validPix), R"pbdoc(
           Performs the raw WFS frame calibration only on useful pixels
           )pbdoc")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_centroids_ref", wy::colCast(&centroider::set_centroids_ref),
           R"pbdoc(
     Set the references slopes

     Args:
          refslopes: (np.array[ndim1,dtype=np.float32]): reference slopes to set
    )pbdoc",
           py::arg("refslopes"))

      .def("set_scale", &centroider::set_scale, R"pbdoc(
     Set the controider scale factor

     Args:
        scale: (float): new scale factor
     )pbdoc",
           py::arg("scale"))

      .def("set_offset", &centroider::set_offset, R"pbdoc(
     Set the controider offset [pixels]

     Args:
        offset: (float): new offset [pixels]
    )pbdoc",
           py::arg("offset"))

      .def("set_dark", wy::colCast(&centroider::set_dark), R"pbdoc(
     Set the dark frame for calibration

     Args:
        dark: (np.ndarray[ndim=2, dtype=np.float32_t): dark frame (size n by n)

        n: (int): image support size
    )pbdoc",
           py::arg("dark"), py::arg("n"))

      .def("set_flat", wy::colCast(&centroider::set_flat), R"pbdoc(
     Set the flat frame for calibration

     Args:
        flat: (np.ndarray[ndim=2, dtype=np.float32_t): flat frame (size n by n)

        n: (int): image support size
    )pbdoc",
           py::arg("flat"), py::arg("n"))

      .def("init_calib", wy::colCast(&centroider::init_calib), R"pbdoc(
     Initialize data used for calibration

     Args:
        n: (int): image support height

        m: (int): image support width
    )pbdoc",
           py::arg("n"), py::arg("m"))

      .def("init_img_raw", wy::colCast(&centroider::init_img_raw), R"pbdoc(
     Initialize array to store raw WFS image in RTC standalone mode

     Args:
        n: (int): image support height

        m: (int): image support width
    )pbdoc",
           py::arg("n"), py::arg("m"))

      .def("set_lutPix", wy::colCast(&centroider::set_lutPix), R"pbdoc(
     Set the lookup Table Pixel vector for calibration

     Args:
        lutPix: (np.ndarray[ndim=1, dtype=np.float32_t): lutPix vector

        n: (int): image pixel size
    )pbdoc",
           py::arg("lutPix"), py::arg("n"));
};

void declare_centroider(py::module &mod) {
  centroider_impl<float, float>(mod, "Centroider_FF");
  centroider_impl<uint16_t, float>(mod, "Centroider_UF");
#ifdef CAN_DO_HALF
  centroider_impl<float, half>(mod, "Centroider_FH");
  centroider_impl<uint16_t, half>(mod, "Centroider_UH");

#endif
}
