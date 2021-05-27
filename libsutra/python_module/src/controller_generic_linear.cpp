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
          "polc", [](controller_generic_linear &sc) { 
              return sc.polc(); }, 
          "polc flag\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "modal", [](controller_generic_linear &sc) { 
              return sc.modal(); }, 
          "modal flag\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_slope_buffers", [](controller_generic_linear &sc) { 
              return sc.n_slope_buffers(); }, 
          "number of slope vectors to store in buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_state_buffers", [](controller_generic_linear &sc) { 
              return sc.n_state_buffers(); },
          "number of state vectors to store in buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_mode_buffers", [](controller_generic_linear &sc) { 
              return sc.n_mode_buffers(); },
          "number of mode vectors to store in buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_states", [](controller_generic_linear &sc) { 
              return sc.n_states(); },
          "number of states in state vector\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_modes", [](controller_generic_linear &sc) { 
              return sc.n_modes(); },
          "number of modes in mode vector\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_iir_in", [](controller_generic_linear &sc) { 
              return sc.n_iir_in(); },
          "number of iir inputs to use\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "n_iir_out", [](controller_generic_linear &sc) { 
              return sc.n_iir_out(); },
          "number of iir outputs to use\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_com", [](controller_generic_linear &sc) { 
              return sc.d_com; },
          "command vector\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_circular_coms", [](controller_generic_linear &sc) { 
              return sc.d_circular_coms; },
          "circular command buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_centroids", [](controller_generic_linear &sc) { 
              return sc.d_centroids; },
          "centroid/slope vector\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_circular_x", [](controller_generic_linear &sc) { 
              return sc.d_circular_x; },
          "circular state buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_circular_s", [](controller_generic_linear &sc) { 
              return sc.d_circular_s; },
          "circular slope buffer\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_circular_u_in", [](controller_generic_linear &sc) { 
              return sc.d_circular_u_in; },
          "circular buffer of iir inputs\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_circular_u_out", [](controller_generic_linear &sc) { 
              return sc.d_circular_u_out; },
          "circular buffer of iir outputs\n\
          (see compass.io details on the generic linear controller)")


      .def_property_readonly(
          "d_A", [](controller_generic_linear &sc) { 
              return sc.d_A; },
          "List of A matrices used for state recursion\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_L", [](controller_generic_linear &sc) { 
              return sc.d_L; },
          "List of L matrices used for innovation from measurements\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_K", [](controller_generic_linear &sc) { 
              return sc.d_K.get(); },
          "K matrix, used for projection from state to modes\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_D", [](controller_generic_linear &sc) { 
              return sc.d_D.get(); },
          "D matrix, typically set to the interaction matrix\n\
          (see compass.io details on the generic linear controller)")
      
      .def_property_readonly(
          "d_imat", [](controller_generic_linear &sc) { 
              return sc.d_D.get(); },
          "D matrix, typically set to the interaction matrix\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_F", [](controller_generic_linear &sc) { 
              return sc.d_F.get(); },
          "F matrix, used to project from modes to actuator voltages\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_iir_a", [](controller_generic_linear &sc) { 
              return sc.d_iir_a; },
          "List of iir coefficient vectors for outputs\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_iir_b", [](controller_generic_linear &sc) { 
              return sc.d_iir_b; },
          "List of iir coefficient vectors for inputs\n\
          (see compass.io details on the generic linear controller)")



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
    Set the POLC flag

    Args:
      polc: (bool): POLC flag
    )pbdoc",
           py::arg("polc"))

      .def("set_A", wy::colCast(&controller_generic_linear::set_A),
           R"pbdoc(
    Set a single A matrix within the list of A matrices (state recursions)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be sliced
      i:  (int): index of list to slice matrix 
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_L", wy::colCast(&controller_generic_linear::set_L),
           R"pbdoc(
    Set a single L matrix within the list of L matrices (innovations)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be sliced
      i:  (int): index of list to slice matrix 
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_K", wy::colCast(&controller_generic_linear::set_K),
           R"pbdoc(
    Set K matrix (state to mode projection)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_D", wy::colCast(&controller_generic_linear::set_D),
           R"pbdoc(
    Set D matrix (interaction matrix)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_F", wy::colCast(&controller_generic_linear::set_F),
           R"pbdoc(
    Set F matrix (mode to actuator voltage projection)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_iir_a", wy::colCast(&controller_generic_linear::set_iir_a),
           R"pbdoc(
    Set a single iir 'a' vector within list (combined with iir outputs)

    Args:
      M:  (np.array[ndim=1,dtype=np.float32]): Vector to be set
      i:  (int): index of list to slice vector 
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_iir_b", wy::colCast(&controller_generic_linear::set_iir_b),
           R"pbdoc(
    Set a single iir 'b' vector within list (combined with iir inputs)

    Args:
      M:  (np.array[ndim=1,dtype=np.float32]): Vector to be set
      i:  (int): index of list to slice vector 
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      ;
};

void declare_controller_generic_linear(py::module &mod) {
  controller_generic_linear_impl<float, float>(mod,
   "ControllerGENERICLINEAR_FF");
  controller_generic_linear_impl<float, uint16_t>(mod,
   "ControllerGENERICLINEAR_FU");
/*#ifdef CAN_DO_HALF
  controller_generic_linear_impl<half, float>(mod, 
   "ControllerGENERICLINEAR_HF");
  controller_generic_linear_impl<half, uint16_t>(mod,
   "ControllerGENERICLINEAR_HU");
#endif*/
}
