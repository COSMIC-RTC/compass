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

//! \file      controller_generic.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraControllerGeneric
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_controller_generic.hpp>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
int32_t set_decayFactor(SutraControllerGeneric<Tcomp, Tout> &sgc, ArrayFStyle<float> &decayFactor) {
    return sgc.set_decayFactor(decayFactor.mutable_data());
}


template <typename Tcomp, typename Tout>
int32_t set_modal_gains(SutraControllerGeneric<Tcomp, Tout> &sgc, ArrayFStyle<float> &gain) {
    return sgc.set_modal_gains(gain.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_cmat(SutraControllerGeneric<Tcomp, Tout> &sgc, ArrayFStyle<float> &cmat) {
    return sgc.set_cmat(cmat.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_matE(SutraControllerGeneric<Tcomp, Tout> &sgc, ArrayFStyle<float> &matE) {
    return sgc.set_matE(matE.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_imat(SutraControllerGeneric<Tcomp, Tout> &sgc, ArrayFStyle<float> &imat) {
    return sgc.set_imat(imat.mutable_data());
}

template <typename Tcomp, typename Tout>
void controller_generic_impl(py::module &mod, const char *name) {
  using controller_generic = SutraControllerGeneric<Tcomp, Tout>;

  py::class_<controller_generic, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "d_matE", [](controller_generic &sc) { return sc.d_matE; },
          "E matrix (see compass.io details on generic controller)")

      .def_property_readonly(
          "d_decayFactor",
          [](controller_generic &sc) { return sc.d_decayFactor; },
          "decayFactor vector (see compass.io details on generic controller)")

      .def_property_readonly(
          "d_cmat", [](controller_generic &sc) { return sc.d_cmat; },
          "Control matrix")

      .def_property_readonly(
          "d_imat", [](controller_generic &sc) { return sc.d_imat; },
          "Control matrix")

      .def_property_readonly(
          "d_gain", [](controller_generic &sc) { return sc.d_gain; },
          "vector of modal gains")

      .def_property_readonly(
          "polc", [](controller_generic &sc) { return sc.polc; }, "POLC flag")

      .def_property_readonly(
          "d_err_ngpu", [](controller_generic &sc) { return sc.d_err_ngpu; },
          "")

      .def_property_readonly(
          "d_compbuff", [](controller_generic &sc) { return sc.d_compbuff; },
          "Computation buffer buffer")

      .def_property_readonly(
          "command_law", [](controller_generic &sc) { return sc.command_law; },
          "Command law currently used")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def("set_decayFactor", &set_decayFactor<Tcomp, Tout>,
           R"pbdoc(
    Set the decay factor vector

    Args:
      decayFactor: (np.array[ndim=1,dtype=np.float32]): decay factor
    )pbdoc",
           py::arg("decayFactor"))

      .def("set_polc", &controller_generic::set_polc,
           R"pbdoc(
    Set the polc flag

    Args:
      polc: (bool): polc flag
    )pbdoc",
           py::arg("polc"))

      .def("set_matE", &set_matE<Tcomp, Tout>,
           R"pbdoc(
    Set the E matrix

    Args:
      E: (np.array[ndim=2,dtype=np.float32]): E matrix to set
    )pbdoc",
           py::arg("E"))

      .def("set_commandlaw", &controller_generic::set_commandlaw,
           R"pbdoc(
    Set the command law to use

    Args:
      commandlaw: (str): command law "integrator", "modal_integrator" or "2matrices"
    )pbdoc",
           py::arg("commandlaw"))

      .def("set_modal_gains", set_modal_gains<Tcomp, Tout>,
           R"pbdoc(
    Set the controller modal gains

    Args:
      mgain: (np.array[ndim=1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_imat", set_imat<Tcomp, Tout>, R"pbdoc(
    Set the interaction matrix

    Args:
      imat: (np.array[ndim=2,dtype=np.float32]): interaction matrix to set
    )pbdoc",
           py::arg("imat"))

      .def("set_cmat", set_cmat<Tcomp, Tout>,
           R"pbdoc(
    Set the command matrix

    Args:
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      ;
};

void declare_controller_generic(py::module &mod) {
  controller_generic_impl<float, float>(mod, "ControllerGENERIC_FF");
  controller_generic_impl<float, uint16_t>(mod, "ControllerGENERIC_FU");
}
