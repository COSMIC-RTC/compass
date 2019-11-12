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

//! \file      tscreen.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_tscreen
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_tscreen.h>

namespace py = pybind11;

void declare_tscreen(py::module &mod) {
  py::class_<sutra_tscreen>(mod, "Tscreen")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](sutra_tscreen &st) { return st.device; }, "Device index")

      .def_property_readonly(
          "d_A", [](sutra_tscreen &st) { return st.d_A; },
          "A matrix for extrusion")

      .def_property_readonly(
          "d_B", [](sutra_tscreen &st) { return st.d_B; },
          "B matrix for extrusion")

      .def_property_readonly(
          "d_istencilx", [](sutra_tscreen &st) { return st.d_istencilx; },
          "stencil for row extrusion")

      .def_property_readonly(
          "d_istencily", [](sutra_tscreen &st) { return st.d_istencily; },
          "stencil for column extrusion")

      .def_property_readonly(
          "d_z", [](sutra_tscreen &st) { return st.d_z; },
          "tmp array for extrusion process")

      .def_property_readonly(
          "d_noise", [](sutra_tscreen &st) { return st.d_noise; },
          "random numbers for extrusion")

      .def_property_readonly(
          "d_ytmp", [](sutra_tscreen &st) { return st.d_ytmp; },
          "contains the extrude update")

      .def_property_readonly(
          "screen_size", [](sutra_tscreen &st) { return st.screen_size; },
          "size of phase screen")

      .def_property_readonly(
          "r0", [](sutra_tscreen &st) { return st.r0; }, "layer r0 in pixels")

      .def_property_readonly(
          "amplitude", [](sutra_tscreen &st) { return st.amplitude; },
          "amplitude for extrusion (r0**(-5/6)")

      .def_property_readonly(
          "altitude", [](sutra_tscreen &st) { return st.altitude; },
          "altitude of the phase screen")

      .def_property_readonly(
          "windspeed", [](sutra_tscreen &st) { return st.windspeed; },
          "wind speed of phase screen")

      .def_property_readonly(
          "winddir", [](sutra_tscreen &st) { return st.winddir; },
          "wind direction of phase screen")

      .def_property_readonly(
          "deltax", [](sutra_tscreen &st) { return st.deltax; },
          "number of columns to extrude per iteration")

      .def_property_readonly(
          "deltay", [](sutra_tscreen &st) { return st.deltay; },
          "number of rows to extrude per iteration")

      .def_property_readonly(
          "accumx", [](sutra_tscreen &st) { return st.accumx; },
          "accumulate columns to extrude")

      .def_property_readonly(
          "accumy", [](sutra_tscreen &st) { return st.accumy; },
          "accumulate rows to extrude")

      .def_property_readonly(
          "d_screen", [](sutra_tscreen &st) { return st.d_tscreen->d_screen; },
          "Turbulent phase screen");
};
