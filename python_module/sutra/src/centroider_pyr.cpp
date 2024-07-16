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

//! \file      centroider_pyr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroiderPyr
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_centroider_pyr.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_pyr_impl(py::module &mod, const char *name) {
  using centroider_pyr = SutraCentroiderPyr<Tin, Tcomp>;

  py::class_<centroider_pyr, SutraCentroider<Tin, Tcomp>>(mod, name)

      .def_property_readonly(
          "pyr_method", [](centroider_pyr &sc) { return sc.get_method_str(); },
          "Method used for pyramid slopes compuation")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def(
          "set_pyr_method",
          [](centroider_pyr &sc, uint8_t method) {
            return sc.set_method(method);
          },
          R"pbdoc(
    Set the pyramid method for slopes computation

    Args:
            method : (int32_t) : new centroiding method (0: nosinus global
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

    Args:
            thresh : (float) : threshold value
        )pbdoc",
          py::arg("thresh"))

      ;
};
void declare_centroider_pyr(py::module &mod) {
  centroider_pyr_impl<float, float>(mod, "CentroiderPYR_FF");
  centroider_pyr_impl<uint16_t, float>(mod, "CentroiderPYR_UF");
}
