// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      centroider_tcog.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroiderTcog
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_centroider_tcog.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_tcog_impl(py::module &mod, const char *name) {
  using centroider_tcog = SutraCentroiderTcog<Tin, Tcomp>;

  py::class_<centroider_tcog, SutraCentroider<Tin, Tcomp>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "threshold", [](centroider_tcog &sc) { return sc.threshold; },
          "Threshold value")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

      .def("set_threshold", &centroider_tcog::set_threshold,
           R"pbdoc(
    Set the threshold value of a TCOG centroider

    Args:
      thresh : (float) : threshold value to set

        )pbdoc",
           py::arg("thresh"))

      ;
};
void declare_centroider_tcog(py::module &mod) {
  centroider_tcog_impl<float, float>(mod, "CentroiderTCOG_FF");
  centroider_tcog_impl<uint16_t, float>(mod, "CentroiderTCOG_UF");
#ifdef CAN_DO_HALF
  centroider_tcog_impl<float, half>(mod, "CentroiderTCOG_FH");
  centroider_tcog_impl<uint16_t, half>(mod, "CentroiderTCOG_UH");

#endif
}
