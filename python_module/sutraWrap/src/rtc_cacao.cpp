// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      rtc_cacao.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraRtcCacao
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_rtc_cacao.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp, typename Tout>
std::unique_ptr<SutraRtcCacao<Tin, Tcomp, Tout>> rtc_cacao_init(
    std::string interface_cal_frame_name, std::string interface_loop_frame_name) {
  return std::unique_ptr<SutraRtcCacao<Tin, Tcomp, Tout>>(
      new SutraRtcCacao<Tin, Tcomp, Tout>(interface_cal_frame_name, interface_loop_frame_name));
}

template <typename Tin, typename Tcomp, typename Tout>
void rtc_cacao_impl(py::module &mod, const char *name) {
  using rtc = SutraRtc<Tin, Tcomp, Tout>;
  using rtc_cacao = SutraRtcCacao<Tin, Tcomp, Tout>;

  py::class_<rtc_cacao, rtc>(mod, name)
      .def(py::init(&rtc_cacao_init<Tin, Tcomp, Tout>), R"pbdoc(
    Create and initialise a cacao rtc object

    Args:
        interface_cal_frame_name:

        interface_loop_frame_name:
            )pbdoc",
           py::arg("interface_cal_frame_name"), py::arg("interface_loop_frame_name"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      // .def_property_readonly("framecounter",
      //                        [](rtc_cacao &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](rtc_cacao &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](rtc_cacao &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](rtc_cacao &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](rtc_cacao &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](rtc_cacao &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](rtc_cacao &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](rtc_cacao &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("publish", &rtc_cacao::publish)

      ;
}

void declare_rtc_cacao(py::module &mod) {
#ifdef CAN_DO_HALF
  rtc_cacao_impl<float, half, float>(mod, "Rtc_cacao_FHF");
  rtc_cacao_impl<uint16_t, half, float>(mod, "Rtc_cacao_UHF");
  rtc_cacao_impl<float, half, uint16_t>(mod, "Rtc_cacao_FHU");
  rtc_cacao_impl<uint16_t, half, uint16_t>(mod, "Rtc_cacao_UHU");
#endif
  rtc_cacao_impl<float, float, float>(mod, "Rtc_cacao_FFF");
  rtc_cacao_impl<uint16_t, float, float>(mod, "Rtc_cacao_UFF");
  rtc_cacao_impl<float, float, uint16_t>(mod, "Rtc_cacao_FFU");
  rtc_cacao_impl<uint16_t, float, uint16_t>(mod, "Rtc_cacao_UFU");
};
