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

//! \file      atmos.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_atmos
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_atmos.h>

namespace py = pybind11;

std::unique_ptr<sutra_atmos> atmos_init(carma_context &context, int nscreens,
                                        float global_r0, float *r0_per_layers,
                                        long *dim_screens, long *stencil_size,
                                        float *altitude, float *windspeed,
                                        float *winddir, float *deltax,
                                        float *deltay, int device) {
  return std::unique_ptr<sutra_atmos>(new sutra_atmos(
      &context, nscreens, global_r0, r0_per_layers, dim_screens, stencil_size,
      altitude, windspeed, winddir, deltax, deltay, device));
}

void declare_atmos(py::module &mod) {
  py::class_<sutra_atmos>(mod, "Atmos")
      .def(py::init(wy::colCast(atmos_init)), R"pbdoc(
        Create and initialise an atmos object
        Parameters
        ------------
        context: (carma_context) : current carma context
        nscreens: (float) : number of turbulent layers
        global_r0: (float): global r0
        r0_per_layer: (float) : r0 per layer
        size1: (np.ndarray[ndim=1, dtype=np.int64_t]) : First size of screens
        size2: (np.ndarray[ndim=1, dtype=np.int64_t]) : Second size of screens
        altitude: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes [m]
        windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind speed [m/s]
        winddir: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind direction [deg]
        deltax: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltax pixels in the x-direction at each iteration
        deltay: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltay pixels in the y-direction at each iteration
        device: (int): GPU device index
        )pbdoc",
           py::arg("context"), py::arg("nscreens"), py::arg("global_r0"),
           py::arg("r0_per_layer"), py::arg("size1"), py::arg("size2"),
           py::arg("altitude"), py::arg("windspeed"), py::arg("winddir"),
           py::arg("deltax"), py::arg("deltay"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "nscreens", [](sutra_atmos &sa) { return sa.nscreens; },
          "Number of turbulent screens")

      .def_property_readonly(
          "r0", [](sutra_atmos &sa) { return sa.r0; }, "Global r0")

      .def_property_readonly(
          "d_screens",
          [](sutra_atmos &sa) -> vector<sutra_tscreen *> & {
            return sa.d_screens;
          },
          "Vector of tscreens")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("move_atmos", &sutra_atmos::move_atmos, R"pbdoc(
        Move the turbulence in the atmos screen following loaded
        parameters such as windspeed and wind direction
        )pbdoc")

      .def("init_screen", wy::colCast(&sutra_atmos::init_screen), R"pbdoc(
        Initialize an newly allocated screen

        Parameters
        ------------
        idx: (int): index of the screen
        A: (np.ndarray[ndim=2, dtype=np.float32]): A matrix (cf. Assemat)
        B: (np.ndarray[ndim=2, dtype=np.float32]): B matrix (cf. Assemat)
        istencilx: (np.ndarray[ndim=2, dtype=int32]): X stencil index
        istencily: (np.ndarray[ndim=2, dtype=int32]): Y stencil index
        seed: (int): seed for RNG
        )pbdoc",
           py::arg("idx"), py::arg("A"), py::arg("B"), py::arg("istencilx"),
           py::arg("istencily"), py::arg("seed"))

      .def("add_screen", wy::colCast(&sutra_atmos::add_screen), R"pbdoc(
        Add a screen to the atmos object.

        Parameters
        ------------
        altitude: (float) : altitude of the screen in meters
        size: (long) : dimension of the screen (size x size)
        stencil_size: (long): dimension of the stencil
        r0: (float) : frac of r0**(5/3)
        windspeed: (float) : windspeed of the screen [m/s]
        winddir: (float) : wind direction (deg)
        deltax: (float) : extrude deltax pixels in the x-direction at each iteration
        deltay: (float) : extrude deltay pixels in the y-direction at each iteration
        device: (int) : device number
        )pbdoc",
           py::arg("altitude"), py::arg("size"), py::arg("stencil_size"),
           py::arg("r0"), py::arg("windspeed"), py::arg("winddir"),
           py::arg("deltax"), py::arg("deltay"), py::arg("device"))

      .def("refresh_screen", wy::colCast(&sutra_atmos::refresh_screen),
           R"pbdoc(
        Refresh the selected screen by extrusion
        Parameters
        ------------
        idx: (int): index of the screen
        )pbdoc",
           py::arg("idx"))

      .def("del_screen", wy::colCast(&sutra_atmos::del_screen), R"pbdoc(
        Delete the selected screen
        Parameters
        ------------
        idx: (int): index of the screen
        )pbdoc",
           py::arg("idx"))

      .def("__str__",
           [](sutra_atmos &sa) {
             std::cout << "Screen # | alt.(m) | speed (m/s) | dir.(deg) | r0 "
                          "(pix) | deltax | deltay"
                       << std::endl;
             vector<sutra_tscreen *>::iterator it = sa.d_screens.begin();
             int i = 0;
             while (it != sa.d_screens.end()) {
               std::cout << i << " | " << (*it)->altitude << " | "
                         << (*it)->windspeed << " | " << (*it)->winddir << " | "
                         << powf((*it)->amplitude, -6 / 5) << " | "
                         << (*it)->deltax << " | " << (*it)->deltay
                         << std::endl;
               i++;
               it++;
             }
             return "";
           })

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_global_r0", wy::colCast(&sutra_atmos::set_global_r0),
           R"pbdoc(
        Change the current global r0 of all layers

        Parameters
        ------------
        r0: (float): r0 @ 0.5 microns
        )pbdoc",
           py::arg("r0"))

      .def("set_frac", wy::colCast(&sutra_atmos::set_frac),
           R"pbdoc(
        Change the current fraction of r0 for all layers

        Parameters
        ------------
        r0: np.ndarray[ndim=1,dtype=np.float32_t]): fraction of r0 for all layers (%)
        )pbdoc",
           py::arg("frac"))

      .def("set_seed", wy::colCast(&sutra_atmos::set_seed), R"pbdoc(
        Set the seed of the selected screen RNG

        Parameters
        ------------
        idx: (int): index of the screen
        seed: (int) :new seed
        )pbdoc",
           py::arg("idx"), py::arg("seed"));
};
