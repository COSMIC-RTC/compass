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

//! \file      groot.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraGroot
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_groot.hpp>

namespace py = pybind11;

std::unique_ptr<SutraGroot> groot_init(CarmaContext &context, int32_t device,
                                        int32_t nactus, int32_t nlayers, float gsangle,
                                        ArrayFStyle<float> &vdt, ArrayFStyle<float> &Htheta, ArrayFStyle<float> &L0,
                                        ArrayFStyle<float> &winddir, ArrayFStyle<float> &scale,
                                        ArrayFStyle<float> &pzt2tt, ArrayFStyle<float> &TTPfilter,
                                        ArrayFStyle<float> &Nact, ArrayFStyle<float> &xpos, ArrayFStyle<float> &ypos,
                                        float fc) {
  return std::unique_ptr<SutraGroot>(new SutraGroot(
      &context, device, nactus, nlayers, gsangle, vdt.mutable_data(), Htheta.mutable_data(), L0.mutable_data(),
      winddir.mutable_data(), scale.mutable_data(), pzt2tt.mutable_data(), TTPfilter.mutable_data(),
      Nact.mutable_data(), xpos.mutable_data(), ypos.mutable_data(), fc));
};

std::unique_ptr<SutraGroot> groot_init_alias(CarmaContext &context,
                                              int32_t device, int32_t nssp,
                                              ArrayFStyle<float> &weights, float scale,
                                              ArrayFStyle<float> &xpos, ArrayFStyle<float> &ypos,
                                              float fc, float d, int32_t npts) {
  return std::unique_ptr<SutraGroot>(new SutraGroot(
      &context, device, nssp, weights.mutable_data(), scale, xpos.mutable_data(), ypos.mutable_data(), fc, d, npts));
};

void declare_groot(py::module &mod) {
  py::class_<SutraGroot>(mod, "Groot")
      .def(py::init(&groot_init), R"pbdoc(
    Initializes Groot to compute aniso and bandwidth model

    Args:
          context: (CarmaContext): context

          device: (int): context active device

          nssp : (str) :

          nlayers:

          gsangle:

          vdt:

          Htheta:

          L0:

          winddir:

          scale:

          pzt2tt:

          TTPfilter:

          Nact:

          xpos:

          ypos:

          fc:
           )pbdoc",
           py::arg("context"), py::arg("device"), py::arg("nssp"),
           py::arg("nlayers"), py::arg("gsangle"), py::arg("vdt"),
           py::arg("Htheta"), py::arg("L0"), py::arg("winddir"),
           py::arg("scale"), py::arg("pzt2tt"), py::arg("TTPfilter"),
           py::arg("Nact"), py::arg("xpos"), py::arg("ypos"), py::arg("fc"))

      .def(py::init(&groot_init_alias), R"pbdoc(
    Initializes Groot to compute aliasing model

    Args:
          context: (CarmaContext): context

          device: (int): context active device

          nssp : (str) :

          weights:

          scale:

          xpos:

          ypos:

          fc:

          d:

          npts:
           )pbdoc",
           py::arg("context"), py::arg("device"), py::arg("nssp"),
           py::arg("weights"), py::arg("scale"), py::arg("xpos"),
           py::arg("ypos"), py::arg("fc"), py::arg("d"), py::arg("npts"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "device", [](SutraGroot &sg) { return sg.device; }, "GPU device")

      .def_property_readonly(
          "nactus", [](SutraGroot &sg) { return sg.nactus; },
          "Number of actuators")

      .def_property_readonly(
          "nssp", [](SutraGroot &sg) { return sg.nssp; }, "number of subap")

      .def_property_readonly(
          "nlayers", [](SutraGroot &sg) { return sg.nlayers; },
          "number of turbulent layers")

      .def_property_readonly(
          "npts", [](SutraGroot &sg) { return sg.npts; },
          "number of samples for aliasig computation")

      .def_property_readonly(
          "gsangle", [](SutraGroot &sg) { return sg.gsangle; },
          "Guide star angle [rad]")

      .def_property_readonly(
          "fc", [](SutraGroot &sg) { return sg.fc; },
          "DM cut-off frequency [m]")

      .def_property_readonly(
          "d", [](SutraGroot &sg) { return sg.d; }, "DM pitch")

      .def_property_readonly(
          "d_Cerr", [](SutraGroot &sg) { return sg.d_Cerr; },
          "Model of aniso and bandwidth covariance error matrix")

      .def_property_readonly(
          "d_CaXX", [](SutraGroot &sg) { return sg.d_CaXX; },
          "XX component of the aliasing model")

      .def_property_readonly(
          "d_CaYY", [](SutraGroot &sg) { return sg.d_CaYY; },
          "YY component of the aliasing model")

      .def_property_readonly(
          "d_TT", [](SutraGroot &sg) { return sg.d_TT; }, "tip-tilt IF matrix")

      .def_property_readonly(
          "scale", [](SutraGroot &sg) { return sg.scale; }, "Scale factor")

      .def_property_readonly(
          "d_TTPfilter", [](SutraGroot &sg) { return sg.d_TTPfilter; },
          "Tip-tilt and piston filter matrix (= Btt.dot(P))")

      .def_property_readonly(
          "d_pzt2tt", [](SutraGroot &sg) { return sg.d_pzt2tt; },
          "pzt to TT matrix")

      .def_property_readonly(
          "d_Nact", [](SutraGroot &sg) { return sg.d_Nact; },
          "Coupling matrix")

      .def_property_readonly(
          "d_xpos", [](SutraGroot &sg) { return sg.d_xpos; },
          "X-positions of DM actuators or ssp [m]")

      .def_property_readonly(
          "d_ypos", [](SutraGroot &sg) { return sg.d_ypos; },
          "Y-positions of DM actuators or ssp [m]")

      .def_property_readonly(
          "d_tab_int_x", [](SutraGroot &sg) { return sg.d_tab_int_x; },
          "Tabulated integral")

      .def_property_readonly(
          "d_tab_int_y", [](SutraGroot &sg) { return sg.d_tab_int_y; },
          "Tabulated integral")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("compute_Cerr", &SutraGroot::compute_Cerr,
           "Computes the aniso and bandwidth error covariance matrix")

      .def("compute_Calias", &SutraGroot::compute_Calias,
           "Computes the aliasing error covariance matrix");
};
