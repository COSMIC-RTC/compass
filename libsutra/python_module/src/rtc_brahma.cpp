// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      rtc_brahma.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraRtcBrahma
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <wyrm>

#include <SutraRtc_brahma.h>

namespace py = pybind11;

template <typename T>
std::unique_ptr<SutraRtcBrahma<T>> rtc_brahma_init(CarmaContext *context,
                                                     SutraSensors *wfs,
                                                     SutraTarget *target,
                                                     ACE_TCHAR *name) {
  return std::unique_ptr<SutraRtcBrahma<T>>(
      new SutraRtcBrahma<T>(context, wfs, target, name));
}

template <typename T>
void rtc_brahma_impl(py::module &mod, const char *name) {
  using rtc = SutraRtc<float, T, float>;
  using rtc_brahma = SutraRtcBrahma<T>;

  py::class_<rtc_brahma, rtc>(mod, name)
      .def(py::init(wy::colCast(rtc_brahma_init<T>)), R"pbdoc(
    Create and initialise a brahma rtc object

    Args:
        context: (CarmaContext) : current carma context

        wfs:

        target:

        name:
        )pbdoc",
           py::arg("context"), py::arg("wfs"), py::arg("target"),
           py::arg("name"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      // .def_property_readonly("framecounter",
      //                        [](rtc_brahma &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](rtc_brahma &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](rtc_brahma &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](rtc_brahma &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](rtc_brahma &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](rtc_brahma &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](rtc_brahma &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](rtc_brahma &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &rtc_brahma::publish)

      ;
}

void declare_rtc_brahma(py::module &mod) {
#ifdef CAN_DO_HALF
  rtc_brahma_impl<half>(mod, "Rtc_brahmaH");
#endif
  rtc_brahma_impl<float>(mod, "Rtc_brahma");
};
