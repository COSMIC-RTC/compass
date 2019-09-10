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

//! \file      centroider_pyr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_centroider_pyr
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_centroider_pyr.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_pyr_impl(py::module &mod, const char *name) {
  using centroider_pyr = sutra_centroider_pyr<Tin, Tcomp>;

  py::class_<centroider_pyr, sutra_centroider<Tin, Tcomp>>(mod, name)

      .def_property_readonly(
          "pyr_method", [](centroider_pyr &sc) { return sc.get_method_str(); },
          "Method used for pyramid slopes compuation")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def(
          "set_pyr_method",
          [](centroider_pyr &sc, uint8_t method) {
            return sc.set_method(method);
          },
          R"pbdoc(
            Set the pyramid method for slopes computation

            Parameters
            ------------
            method : (int) : new centroiding method (0: nosinus global
                                                    1: sinus global
                                                    2: nosinus local
                                                    3: sinus local)
                            favor use of shesha_constant.PyrCentroiderMethod

        )pbdoc",
          py::arg("method"))

      .def(
          "set_pyr_thresh",
          [](centroider_pyr &sc, float thresh) {
            return sc.set_valid_thresh(thresh);
          },
          R"pbdoc(
            Set the pyramid threshold value

            Parameters
            ------------
            thresh : (float) : threshold value
        )pbdoc",
          py::arg("thresh"));
};
void declare_centroider_pyr(py::module &mod) {
  centroider_pyr_impl<float, float>(mod, "CentroiderPYR_FF");
  centroider_pyr_impl<uint16_t, float>(mod, "CentroiderPYR_UF");
}
