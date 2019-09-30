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
//! \brief     this file provides pybind wrapper for carma_timer
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

namespace py = pybind11;

void declare_carmaWrap_timer(py::module &mod) {
  py::class_<carma_timer>(mod, "timer")
      .def(py::init())
      .def_property_readonly("total_time",
                             [](carma_timer &ct) { return ct.elapsed(); })
      .def("reset", &carma_timer::reset)
      .def("start", &carma_timer::start)
      .def("stop", &carma_timer::stop)
      .def("setStream", [](carma_timer &ct, carma_device &cd) {
        ct.setStream(cd.get_stream());
      });
  //  .def("stop",
  //       [](carma_timer &ct, carma_device &cd) { ct.stop(cd.get_stream()); });

  py::class_<carma_clock>(mod, "clock")
      .def(py::init([](carma_context &context, int i) {
             return std::unique_ptr<carma_clock>(new carma_clock(&context, i));
           }),
           R"pbdoc(
        Create a carma_clock object which provides timing based on GPU clock count

        Parameters
        ------------
        context: (carma_context): carma context
        i: (int): time buffer size
      )pbdoc",
           py::arg("context"), py::arg("i"))

      .def_property_readonly("timeBuffer",
                             [](carma_clock &clk) { return clk.timeBuffer; })
      .def_property_readonly("cc", [](carma_clock &clk) { return clk.cc; })
      .def_property_readonly("GPUfreq",
                             [](carma_clock &clk) { return clk.GPUfreq; })

      .def("tic", &carma_clock::tic)
      .def("toc", &carma_clock::toc);
}
