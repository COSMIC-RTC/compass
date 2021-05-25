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

//! \file      controller_generic_linear.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_controller_generic_linear
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2021/05/12
//! \copyright GNU Lesser General Public License

#include <sutra_controller_generic_linear.h>

#include <wyrm>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_generic_linear_impl(py::module &mod, const char *name) {
  using controller_generic_linear = sutra_controller_generic_linear<Tcomp, Tout>;

  py::class_<controller_generic_linear, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "polc", [](controller_generic_linear &sc) { return sc.polc(); },
          "polc flag (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "modal", [](controller_generic_linear &sc) { return sc.modal(); },
          "modal flag (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_slope_buffers", [](controller_generic_linear &sc) { return sc.n_slope_buffers(); },
          "number of slope vectors to store in buffer\n(see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_state_buffers", [](controller_generic_linear &sc) { return sc.n_state_buffers(); },
          "number of state vectors to store in buffer\n(see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_mode_buffers", [](controller_generic_linear &sc) { return sc.n_mode_buffers(); },
          "number of mode vectors to store in buffer\n(see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_states", [](controller_generic_linear &sc) { return sc.n_states(); },
          "number of states in state vector (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_modes", [](controller_generic_linear &sc) { return sc.n_modes(); },
          "number of modes in mode vector (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_iir_in", [](controller_generic_linear &sc) { return sc.n_iir_in(); },
          "number of iir inputs to use (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_iir_out", [](controller_generic_linear &sc) { return sc.n_iir_out(); },
          "number of iir outputs to use (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_com", [](controller_generic_linear &sc) { return sc.d_com; },
          "command vector (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_circular_coms", [](controller_generic_linear &sc) { return sc.d_circular_coms; },
          "circular command buffer (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_centroids", [](controller_generic_linear &sc) { return sc.d_centroids; },
          "centroid/slope vector (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_circular_x", [](controller_generic_linear &sc) { return sc.d_circular_x; },
          "circular state buffer (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_circular_s", [](controller_generic_linear &sc) { return sc.d_circular_s; },
          "circular slope buffer (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_circular_u_in", [](controller_generic_linear &sc) { return sc.d_circular_u_in; },
          "circular buffer of iir inputs (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_circular_u_out", [](controller_generic_linear &sc) { return sc.d_circular_u_out; },
          "circular buffer of iir outputs (see compass.io details on the generic linear controller)")
      
      
      
      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗    ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝    ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗  ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝  ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗    ██║     ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝    ██║     ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗  ██║     ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝  ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_polc", wy::colCast(&controller_generic_linear::set_polc),
           R"pbdoc(
    Set the polc flag
	TODO ALL DOCS
    Args:
    )pbdoc",
           py::arg("polc"))

      .def("set_A", wy::colCast(&controller_generic_linear::set_A),
           R"pbdoc(
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_L", wy::colCast(&controller_generic_linear::set_L),
           R"pbdoc(
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_K", wy::colCast(&controller_generic_linear::set_K),
           R"pbdoc(
    )pbdoc",
           py::arg("M"))

      .def("set_D", wy::colCast(&controller_generic_linear::set_D),
           R"pbdoc(
    )pbdoc",
           py::arg("M"))

      .def("set_F", wy::colCast(&controller_generic_linear::set_F),
           R"pbdoc(
    )pbdoc",
           py::arg("M"))

      .def("set_iir_a", wy::colCast(&controller_generic_linear::set_iir_a),
           R"pbdoc(
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_iir_b", wy::colCast(&controller_generic_linear::set_iir_b),
           R"pbdoc(
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      ;
};

void declare_controller_generic_linear(py::module &mod) {
  controller_generic_linear_impl<float, float>(mod, "ControllerGENERICLINEAR_FF");
  controller_generic_linear_impl<float, uint16_t>(mod, "ControllerGENERICLINEAR_FU");
/*#ifdef CAN_DO_HALF
  controller_generic_linear_impl<half, float>(mod, "ControllerGENERICLINEAR_HF");
  controller_generic_linear_impl<half, uint16_t>(mod, "ControllerGENERICLINEAR_HU");
#endif*/
}
