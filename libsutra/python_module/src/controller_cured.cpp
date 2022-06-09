// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      controller_cured.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraControllerCured
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <sutra_controller_cured.h>

#include <wyrm>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_cured_impl(py::module &mod, const char *name) {
  using controller_cured = SutraControllerCured<Tcomp, Tout>;

  py::class_<controller_cured, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "ndivs", [](controller_cured &sc) { return sc.ndivs; },
          "Number of subdivision levels")

      .def_property_readonly(
          "tt_flag", [](controller_cured &sc) { return sc.tt_flag; },
          "Flag to separate TT")

      .def_property_readonly(
          "h_centroids", [](controller_cured &sc) { return sc.h_centroids; },
          "Centroids")

      .def_property_readonly(
          "h_err", [](controller_cured &sc) { return sc.h_err; },
          "Increment error")

      .def_property_readonly(
          "d_err", [](controller_cured &sc) { return sc.d_err; },
          "Increment error")

      .def_property_readonly(
          "d_cenbuff", [](controller_cured &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly(
          "d_imat", [](controller_cured &sc) { return sc.d_imat; },
          "Interaction matrix")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("init_cured", wy::colCast(&controller_cured::init_cured),
           R"pbdoc(
    Initialize CURED

    Args:
        nxsub: (int): TODO: docstring

        isvalid:

        ndivs: (int):

        tt: (int):
    )pbdoc",
           py::arg("nxsub"), py::arg("isvalid"), py::arg("ndivs"),
           py::arg("tt"));
};

void declare_controller_cured(py::module &mod) {
  controller_cured_impl<float, float>(mod, "ControllerCURED_FF");
  controller_cured_impl<float, uint16_t>(mod, "ControllerCURED_FU");
}
