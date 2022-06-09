// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      groot.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraGroot
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_groot.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

std::unique_ptr<SutraGroot> groot_init(CarmaContext &context, int device,
                                        int nactus, int nlayers, float gsangle,
                                        float *vdt, float *Htheta, float *L0,
                                        float *winddir, float *scale,
                                        float *pzt2tt, float *TTPfilter,
                                        float *Nact, float *xpos, float *ypos,
                                        float fc) {
  return std::unique_ptr<SutraGroot>(new SutraGroot(
      &context, device, nactus, nlayers, gsangle, vdt, Htheta, L0, winddir,
      scale, pzt2tt, TTPfilter, Nact, xpos, ypos, fc));
};

std::unique_ptr<SutraGroot> groot_init_alias(CarmaContext &context,
                                              int device, int nssp,
                                              float *weights, float scale,
                                              float *xpos, float *ypos,
                                              float fc, float d, int npts) {
  return std::unique_ptr<SutraGroot>(new SutraGroot(
      &context, device, nssp, weights, scale, xpos, ypos, fc, d, npts));
};

void declare_groot(py::module &mod) {
  py::class_<SutraGroot>(mod, "Groot")
      .def(py::init(wy::colCast(groot_init)), R"pbdoc(
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

      .def(py::init(wy::colCast(groot_init_alias)), R"pbdoc(
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
      //

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
      .def("compute_Cerr", wy::colCast(&SutraGroot::compute_Cerr),
           "Computes the aniso and bandwidth error covariance matrix")

      .def("compute_Calias", wy::colCast(&SutraGroot::compute_Calias),
           "Computes the aliasing error covariance matrix");
};
