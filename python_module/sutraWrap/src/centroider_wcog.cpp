// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      centroider_wcog.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCentroiderWcog
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_centroider_wcog.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
int32_t load_weights(SutraCentroiderWcog<Tin, Tcomp> &scw, ArrayFStyle<float> &weights, int32_t ndim) {
    return scw.load_weights(weights.mutable_data(), ndim);
}

template <typename Tin, typename Tcomp>
void centroider_wcog_impl(py::module &mod, const char *name) {
  using centroider_wcog = SutraCentroiderWcog<Tin, Tcomp>;

  py::class_<centroider_wcog, SutraCentroider<Tin, Tcomp>>(mod, name)
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "d_weights", [](centroider_wcog &sc) { return sc.d_weights; },
          "Weights applied")

      .def_property_readonly(
          "threshold", [](centroider_wcog &sc) { return sc.threshold; },
          "Threshold value")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("init_weights", &centroider_wcog::init_weights,
           "Initializes WCOG computation")

      .def("load_weights", &load_weights<Tin, Tcomp>,
           R"pbdoc(
    Load weights on WCOG

    Args:
            weight: (np.array[ndim=, dtype=np.float32]: weights

            ndim: (int): TODO: docstring
        )pbdoc",
           py::arg("weights"), py::arg("ndim"))

      .def("set_threshold", &centroider_wcog::set_threshold,
           R"pbdoc(
    Set the threshold value of a TCOG centroider

    Args:
      thresh : (float) : threshold value to set

        )pbdoc",
           py::arg("thresh"))

      ;
};
void declare_centroider_wcog(py::module &mod) {
  centroider_wcog_impl<float, float>(mod, "CentroiderWCOG_FF");
  centroider_wcog_impl<uint16_t, float>(mod, "CentroiderWCOG_UF");
}
