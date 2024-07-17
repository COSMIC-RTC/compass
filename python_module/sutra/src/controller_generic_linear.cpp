// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      controller_generic_linear.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraControllerGenericLinear
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2021/05/12

#include "sutraUtils.hpp"

#include <sutra_controller_generic_linear.hpp>

namespace py = pybind11;


template <typename Tcomp, typename Tout>
int32_t set_matA(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &A, int32_t i){
    return scgl.set_matA(A.mutable_data(), i);
}

template <typename Tcomp, typename Tout>
int32_t set_matL(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &L, int32_t i){
    return scgl.set_matL(L.mutable_data(), i);
}

template <typename Tcomp, typename Tout>
int32_t set_matK(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &K){
    return scgl.set_matK(K.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_matD(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &D){
    return scgl.set_matD(D.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_matF(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &F){
    return scgl.set_matF(F.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_iir_a(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &iir_a, int32_t i){
    return scgl.set_iir_a(iir_a.mutable_data(), i);
}

template <typename Tcomp, typename Tout>
int32_t set_iir_b(SutraControllerGenericLinear<Tcomp, Tout> &scgl, ArrayFStyle<float> &iir_b, int32_t i){
    return scgl.set_iir_b(iir_b.mutable_data(), i);
}

template <typename Tcomp, typename Tout>
void controller_generic_linear_impl(py::module &mod, const char *name) {
  using controller_generic_linear = SutraControllerGenericLinear<Tcomp, Tout>;

  py::class_<controller_generic_linear, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

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
          "d_x_now", [](controller_generic_linear &sc) {
              return sc.d_x_now.get(); },
          "temporary state vector used for control calcs\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_s_now", [](controller_generic_linear &sc) {
              return sc.d_s_now.get(); },
          "temporary slope vector used for control calcs\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_u_now", [](controller_generic_linear &sc) {
              return sc.d_u_now.get(); },
          "temporary command vector used for control calcs\n\
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
          "d_matA", [](controller_generic_linear &sc) {
              return sc.d_matA; },
          "List of A matrices used for state recursion\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_matL", [](controller_generic_linear &sc) {
              return sc.d_matL; },
          "List of L matrices used for innovation from measurements\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_matK", [](controller_generic_linear &sc) {
              return sc.d_matK.get(); },
          "K matrix, used for projection from state to modes\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_matD", [](controller_generic_linear &sc) {
              return sc.d_matD.get(); },
          "D matrix, typically set to the interaction matrix\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_imat", [](controller_generic_linear &sc) {
              return sc.d_matD.get(); },
          "D matrix, typically set to the interaction matrix\n\
          (see compass.io details on the generic linear controller)")

      .def_property_readonly(
          "d_matF", [](controller_generic_linear &sc) {
              return sc.d_matF.get(); },
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

      .def("set_polc", &controller_generic_linear::set_polc,
           R"pbdoc(
    Set the POLC flag

    Args:
      polc: (bool): POLC flag
    )pbdoc",
           py::arg("polc"))

      .def("set_matA", &set_matA< Tcomp, Tout >,
           R"pbdoc(
    Set a single A matrix within the list of A matrices (state recursions)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be sliced
      i:  (int): index of list to slice matrix
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_matL", &set_matL< Tcomp, Tout >,
           R"pbdoc(
    Set a single L matrix within the list of L matrices (innovations)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be sliced
      i:  (int): index of list to slice matrix
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_matK", &set_matK< Tcomp, Tout >,
           R"pbdoc(
    Set K matrix (state to mode projection)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_matD", &set_matD< Tcomp, Tout >,
           R"pbdoc(
    Set D matrix (interaction matrix)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_matF", &set_matF< Tcomp, Tout >,
           R"pbdoc(
    Set F matrix (mode to actuator voltage projection)

    Args:
      M:  (np.array[ndim=2,dtype=np.float32]): Matrix to be set
    )pbdoc",
           py::arg("M"))

      .def("set_iir_a", &set_iir_a< Tcomp, Tout >,
           R"pbdoc(
    Set a single iir 'a' vector within list (combined with iir outputs)

    Args:
      M:  (np.array[ndim=1,dtype=np.float32]): Vector to be set
      i:  (int): index of list to slice vector
    )pbdoc",
           py::arg("M"),
	   py::arg("i"))

      .def("set_iir_b", &set_iir_b< Tcomp, Tout >,
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
}
