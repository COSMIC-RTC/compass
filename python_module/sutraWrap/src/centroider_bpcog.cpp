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

//! \file      centroider_bpcog.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroiderBpcog
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_centroider_bpcog.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_bpcog_impl(py::module &mod, const char *name) {
  using centroider_bpcog = SutraCentroiderBpcog<Tin, Tcomp>;

  py::class_<centroider_bpcog, SutraCentroider<Tin, Tcomp>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "nmax", [](centroider_bpcog &sc) { return sc.nmax; },
          "Number of brightest pixels")

      .def_property_readonly(
          "d_bpix", [](centroider_bpcog &sc) { return sc.d_bpix; },
          "Brightest pixels")

      .def_property_readonly(
          "d_bpind", [](centroider_bpcog &sc) { return sc.d_bpind; },
          "Brightest pixels indices")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def("set_nmax", &centroider_bpcog::set_nmax,
           R"pbdoc(
    Set the number of brightest pixels considered for COG computation

    Args:
        nmax : (float) : nmax value to set
        )pbdoc",
           py::arg("nmax"))

      ;
};

void declare_centroider_bpcog(py::module &mod) {
  centroider_bpcog_impl<float, float>(mod, "CentroiderBPCOG_FF");
  centroider_bpcog_impl<uint16_t, float>(mod, "CentroiderBPCOG_UF");
#ifdef CAN_DO_HALF
  centroider_bpcog_impl<float, half>(mod, "CentroiderBPCOG_FH");
  centroider_bpcog_impl<uint16_t, half>(mod, "CentroiderBPCOG_UH");

#endif
}
