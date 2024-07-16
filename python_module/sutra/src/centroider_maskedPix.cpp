// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      centroider_maskedpix.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_centroider_maskedpix
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_centroider_maskedPix.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_maskedPix_impl(py::module &mod, const char *name) {
  using centroider_maskedPix = SutraCentroiderMaskedPix<Tin, Tcomp>;

  py::class_<centroider_maskedPix, SutraCentroider<Tin, Tcomp>>(mod, name)

      .def_property_readonly(
          "d_selected_pix",
          [](centroider_maskedPix &sc) { return sc.d_selected_pix; },
          "Selected pixels as an image")

      .def_property_readonly(
          "d_mask", [](centroider_maskedPix &sc) { return sc.d_mask; },
          "Mask of valid pixels")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def(
          "fill_selected_pix",
          [](centroider_maskedPix &sc, CarmaObj<Tcomp> &centroids) {
            sc.fill_selected_pix(centroids);
          },
          R"pbdoc(
    Return the given pixels vector in an image format

    Args:
        centroids: (CarmaObj): Pixels to map on the image
    )pbdoc",
          py::arg("centroids"))

      .def("fill_mask", &centroider_maskedPix::fill_mask,
           "Fill the mask of valid pixels")

      ;
};
void declare_centroider_maskedPix(py::module &mod) {
  centroider_maskedPix_impl<float, float>(mod, "CentroiderMASKEDPIX_FF");
  centroider_maskedPix_impl<uint16_t, float>(mod, "CentroiderMASKEDPIX_UF");
}
