// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      timer.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaTimer
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

namespace py = pybind11;

void declare_carmaWrap_timer(py::module &mod) {
  py::class_<CarmaTimer>(mod, "timer")
      .def(py::init())
      .def_property_readonly("total_time",
                             [](CarmaTimer &ct) { return ct.elapsed(); })
      .def("reset", &CarmaTimer::reset)
      .def("start", &CarmaTimer::start)
      .def("stop", &CarmaTimer::stop)
      .def("set_stream", [](CarmaTimer &ct, CarmaDevice &cd) {
        ct.set_stream(cd.get_stream());
      });
  //  .def("stop",
  //       [](CarmaTimer &ct, CarmaDevice &cd) { ct.stop(cd.get_stream()); });

  py::class_<CarmaClock>(mod, "clock")
      .def(py::init([](CarmaContext &context, int i) {
             return std::unique_ptr<CarmaClock>(new CarmaClock(&context, i));
           }),
           R"pbdoc(
    Create a CarmaClock object which provides timing based on GPU clock count

    Args:
        context: (CarmaContext): carma context

        i: (int): time buffer size
      )pbdoc",
           py::arg("context"), py::arg("i"))

      .def_property_readonly("time_buffer",
                             [](CarmaClock &clk) { return clk.time_buffer; })
      .def_property_readonly("cc", [](CarmaClock &clk) { return clk.cc; })
      .def_property_readonly("gpu_freq",
                             [](CarmaClock &clk) { return clk.gpu_freq; })

      .def("tic", &CarmaClock::tic)
      .def("toc", &CarmaClock::toc);
}
