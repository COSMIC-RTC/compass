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

//! \file      timer.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaTimer
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.hpp>

namespace py = pybind11;

void declare_carma_timer(py::module &mod) {
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
      .def(py::init([](CarmaContext &context, int32_t i) {
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
