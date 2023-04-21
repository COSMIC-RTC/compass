// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      controller_generic.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_controller_generic
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include <sutra_controller_generic.h>

#include <wyrm>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_generic_impl(py::module &mod, const char *name) {
  using controller_generic = sutra_controller_generic<Tcomp, Tout>;

  py::class_<controller_generic, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
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
      //
      .def("set_decayFactor", wy::colCast(&controller_generic::set_decayFactor),
           R"pbdoc(
    Set the decay factor vector

    Args:
      decayFactor: (np.array[ndim=1,dtype=np.float32]): decay factor
    )pbdoc",
           py::arg("decayFactor"))

      .def("set_polc", wy::colCast(&controller_generic::set_polc),
           R"pbdoc(
    Set the polc flag

    Args:
      polc: (bool): polc flag
    )pbdoc",
           py::arg("polc"))

      .def("set_matE", wy::colCast(&controller_generic::set_matE),
           R"pbdoc(
    Set the E matrix

    Args:
      E: (np.array[ndim=2,dtype=np.float32]): E matrix to set
    )pbdoc",
           py::arg("E"))

      .def("set_commandlaw", wy::colCast(&controller_generic::set_commandlaw),
           R"pbdoc(
    Set the command law to use

    Args:
      commandlaw: (str): command law "integrator", "modal_integrator" or "2matrices"
    )pbdoc",
           py::arg("commandlaw"))

      .def("set_modal_gains", wy::colCast(&controller_generic::set_modal_gains),
           R"pbdoc(
    Set the controller modal gains

    Args:
      mgain: (np.array[ndim=1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_imat", wy::colCast(&controller_generic::set_imat), R"pbdoc(
    Set the interaction matrix

    Args:
      imat: (np.array[ndim=2,dtype=np.float32]): interaction matrix to set
    )pbdoc",
           py::arg("imat"))

      .def("set_cmat", wy::colCast(&controller_generic::set_cmat),
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
#ifdef CAN_DO_HALF
  controller_generic_impl<half, float>(mod, "ControllerGENERIC_HF");
  controller_generic_impl<half, uint16_t>(mod, "ControllerGENERIC_HU");
#endif
}
