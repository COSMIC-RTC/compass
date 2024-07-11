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

//! \file      centroider_corr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroiderCorr
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_centroider_corr.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
int32_t init_corr(SutraCentroiderCorr<Tin, Tcomp> &sc, int32_t isizex, int32_t isizey, ArrayFStyle<float> &interpmat) {
    return sc.init_corr(isizex, isizey, interpmat.mutable_data());
}

template <typename Tin, typename Tcomp>
int32_t load_corr(SutraCentroiderCorr<Tin, Tcomp> &sc, ArrayFStyle<Tcomp> &corr, ArrayFStyle<Tcomp> &corr_norm, int32_t ndim) {
    return sc.load_corr(corr.mutable_data(), corr_norm.mutable_data(), ndim);
}

template <typename Tin, typename Tcomp>
void centroider_corr_impl(py::module &mod, const char *name) {
  using centroider_corr = SutraCentroiderCorr<Tin, Tcomp>;

  py::class_<centroider_corr, SutraCentroider<Tin, Tcomp>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "npix", [](centroider_corr &sc) { return sc.npix; },
          "TODO: docstring")

      .def_property_readonly(
          "interp_sizex", [](centroider_corr &sc) { return sc.interp_sizex; },
          "TODO: docstring")

      .def_property_readonly(
          "interp_sizey", [](centroider_corr &sc) { return sc.interp_sizey; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrfnct", [](centroider_corr &sc) { return sc.d_corrfnct; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrspot", [](centroider_corr &sc) { return sc.d_corrspot; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrnorm", [](centroider_corr &sc) { return sc.d_corrnorm; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrmax", [](centroider_corr &sc) { return sc.d_corrmax; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corr", [](centroider_corr &sc) { return sc.d_corr; },
          "TODO: docstring")

      .def_property_readonly(
          "d_interpmat", [](centroider_corr &sc) { return sc.d_interpmat; },
          "TODO: docstring")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def("init_corr", &init_corr<Tin, Tcomp>, R"pbdoc(
    Initializes corr computation

    Args:
        isizex: (int): TODO: docstring

        isizey: (int):

        interpmat: (np.array[ndim=,dtype=np.float32]):
        )pbdoc",
           py::arg("isizex"), py::arg("isizey"), py::arg("interpmat"))

      .def("load_corr", &load_corr<Tin, Tcomp>, R"pbdoc(
    Load arrays for correlation computation

    Args:
        corr: (np.array[ndim=,dtype=np.float32): TODO: docstring

        corr_norm: (np.array[ndim=,dtype=np.float32):

        ndim: (int):
        )pbdoc",
           py::arg("corr"), py::arg("corr_norm"), py::arg("ndim"))


      ;
};

void declare_centroider_corr(py::module &mod) {
  centroider_corr_impl<float, float>(mod, "CentroiderCORR_FF");
  centroider_corr_impl<uint16_t, float>(mod, "CentroiderCORR_UF");
}
