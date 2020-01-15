// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      centroider_corr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_centroider_corr
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_centroider_corr.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_corr_impl(py::module &mod, const char *name) {
  using centroider_corr = sutra_centroider_corr<Tin, Tcomp>;

  py::class_<centroider_corr, sutra_centroider<Tin, Tcomp>>(mod, name)
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
      //
      .def("init_corr", wy::colCast(&centroider_corr::init_corr), R"pbdoc(
            Initializes corr computation

            Parameters
            ------------
            isizex: (int): TODO: docstring
            isizey: (int):
            interpmat: (np.array[ndim=,dtype=np.float32]):
            )pbdoc",
           py::arg("isizex"), py::arg("isizey"), py::arg("interpmat"))

      .def("load_corr", wy::colCast(&centroider_corr::load_corr), R"pbdoc(
            Load arrays for correlation computation

            Parameters
            ------------
            corr: (np.array[ndim=,dtype=np.float32): TODO: docstring
            corr_norm: (np.array[ndim=,dtype=np.float32):
            ndim: (int):
            )pbdoc",
           py::arg("corr"), py::arg("corr_norm"), py::arg("ndim"))

      .def("set_npix", wy::colCast(&centroider_corr::set_npix),
           R"pbdoc(
               Set the number of pixels per subap.
            Parameters
            ------------
            npix: (int): number of pixels per subap
            )pbdoc",
           py::arg("npix"))

      ;
};

void declare_centroider_corr(py::module &mod) {
  centroider_corr_impl<float, float>(mod, "CentroiderCORR_FF");
  centroider_corr_impl<uint16_t, float>(mod, "CentroiderCORR_UF");
}
