// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      atmos.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraAtmos
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_atmos.hpp>

namespace py = pybind11;

std::unique_ptr<SutraAtmos> atmos_init(CarmaContext &context, int32_t nscreens,
                                        float global_r0, ArrayFStyle<float> &r0_per_layers,
                                        ArrayFStyle<int64_t> &dim_screens, ArrayFStyle<int64_t> &stencil_size,
                                        ArrayFStyle<float> &altitude, ArrayFStyle<float> &windspeed,
                                        ArrayFStyle<float> &winddir, ArrayFStyle<float> &deltax,
                                        ArrayFStyle<float> &deltay, int32_t device) {
  return std::unique_ptr<SutraAtmos>(new SutraAtmos(
      &context, nscreens, global_r0, r0_per_layers.mutable_data(), dim_screens.mutable_data(), stencil_size.mutable_data(),
      altitude.mutable_data(), windspeed.mutable_data(), winddir.mutable_data(), deltax.mutable_data(), deltay.mutable_data(), device));
}

int32_t init_screen(SutraAtmos &sa, int32_t idx, ArrayFStyle<float> &h_A, ArrayFStyle<float> &h_B, ArrayFStyle<uint32_t> &h_istencilx,
                ArrayFStyle<uint32_t> &h_istencily, int32_t seed) {
  return sa.init_screen(idx, h_A.mutable_data(), h_B.mutable_data(), h_istencilx.mutable_data(), h_istencily.mutable_data(), seed);
}


void declare_atmos(py::module &mod) {
  py::class_<SutraAtmos>(mod, "Atmos")
      .def(py::init(&atmos_init),
    R"pbdoc(
    Create and initialise an atmos object

    Attributes:
        context: (CarmaContext) : current carma context

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

      .def_property_readonly(
          "nscreens", [](SutraAtmos &sa) { return sa.nscreens; },
          "Number of turbulent screens")

      .def_property_readonly(
          "r0", [](SutraAtmos &sa) { return sa.r0; }, "Global r0")

      .def_property_readonly(
          "d_screens",
          [](SutraAtmos &sa) -> vector<SutraTurbuScreen *> & {
            return sa.d_screens;
          },
          "Vector of tscreens")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("move_atmos", &SutraAtmos::move_atmos, R"pbdoc(
        Move the turbulence in the atmos screen following loaded
        parameters such as windspeed and wind direction
        )pbdoc")

      .def("init_screen", &init_screen, R"pbdoc(
    Initialize an newly allocated screen

    Args:
        idx: (int): index of the screen

        A: (np.ndarray[ndim=2, dtype=np.float32]): A matrix (cf. Assemat)

        B: (np.ndarray[ndim=2, dtype=np.float32]): B matrix (cf. Assemat)

        istencilx: (np.ndarray[ndim=2, dtype=int32]): X stencil index

        istencily: (np.ndarray[ndim=2, dtype=int32]): Y stencil index

        seed: (int): seed for RNG
        )pbdoc",
           py::arg("idx"), py::arg("A"), py::arg("B"), py::arg("istencilx"),
           py::arg("istencily"), py::arg("seed"))

      .def("add_screen", &SutraAtmos::add_screen, R"pbdoc(
    Add a screen to the atmos object.

    Args:
        altitude: (float) : altitude of the screen in meters

        size: (int64_t) : dimension of the screen (size x size)

        stencil_size: (int64_t): dimension of the stencil

        r0: (float) : frac of r0**(5/3)

        windspeed: (float) : windspeed of the screen [m/s]

        winddir: (float) : wind direction (deg)

        deltax: (float) : extrude deltax pixels in the x-direction at each iteration

        deltay: (float) : extrude deltay pixels in the y-direction at each iteration

        device: (int32_t) : device number
        )pbdoc",
           py::arg("altitude"), py::arg("size"), py::arg("stencil_size"),
           py::arg("r0"), py::arg("windspeed"), py::arg("winddir"),
           py::arg("deltax"), py::arg("deltay"), py::arg("device"))

      .def("refresh_screen", &SutraAtmos::refresh_screen,
           R"pbdoc(
    Refresh the selected screen by extrusion

    Args:
        idx: (int): index of the screen
        )pbdoc",
           py::arg("idx"))

      .def("del_screen", &SutraAtmos::del_screen, R"pbdoc(
    Delete the selected screen

    Args:
        idx: (int): index of the screen
        )pbdoc",
           py::arg("idx"))

      .def("__str__",
           [](SutraAtmos &sa) {
             std::cout << "Screen # | alt.(m) | speed (m/s) | dir.(deg) | r0 "
                          "(pix) | deltax | deltay"
                       << std::endl;
             vector<SutraTurbuScreen *>::iterator it = sa.d_screens.begin();
             int32_t i = 0;
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

      .def("set_r0", &SutraAtmos::set_r0,
           R"pbdoc(
    Change the current global r0 of all layers

    Args:
        r0: (float): r0 @ 0.5 microns
        )pbdoc",
           py::arg("r0"))

      .def("set_seed", &SutraAtmos::set_seed, R"pbdoc(
    Set the seed of the selected screen RNG

    Args:
        idx: (int): index of the screen

        seed: (int32_t) :new seed
        )pbdoc",
           py::arg("idx"), py::arg("seed"));
};
