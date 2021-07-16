// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      timer.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaTimer
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
