// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      target.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraTarget
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_target.hpp>

namespace py = pybind11;

std::unique_ptr<SutraTarget> target_init(CarmaContext &context,
                                          SutraTelescope *d_tel, int32_t ntargets,
                                          ArrayFStyle<float> &xpos, ArrayFStyle<float> &ypos,
                                          ArrayFStyle<float> &lambda, ArrayFStyle<float> &mag,
                                          float zerop, ArrayFStyle<int64_t> &sizes, int32_t Npts,
                                          int32_t device) {
  return std::unique_ptr<SutraTarget>(
      new SutraTarget(&context, d_tel, ntargets, xpos.mutable_data(), ypos.mutable_data(),
                      lambda.mutable_data(), mag.mutable_data(), zerop, sizes.mutable_data(), Npts,
                      device));
}

void declare_target(py::module &mod) {
  py::class_<SutraTarget>(mod, "Target")
      .def(py::init(&target_init), R"pbdoc(
    Create and initialise an target object

    Args:
        context: (CarmaContext) : current carma context

        d_tel: (SutraTelescope) : SutraTelescope object

        ntargets: (int): number of targets

        xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : X positions of each target in arcsec

        ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : Y positions of each target in arcsec

        lambda_um: (np.ndarray[ndim=1,dtype=np.float32_t]) : Wavelength of each target in µm

        mag: (np.ndarray[ndim=1,dtype=np.float32_t]) : magnitude of each target

        zerop: (float) : Flux at magnitude 0 in photons/m²/s

        sizes: (np.ndarray[ndim=1,dtype=np.int64_t]) : Support size of each target

        Npts : (int): number of points in the pupil

        device: (int): GPU device index
        )pbdoc",
           py::arg("context"), py::arg("d_tel"), py::arg("ntargets"),
           py::arg("xpos"), py::arg("ypos"), py::arg("lambda_um"), py::arg("mag"),
           py::arg("zerop"), py::arg("sizes"), py::arg("Npts"),
           py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "ntargets", [](SutraTarget &st) { return st.ntargets; },
          "Number of targets")

      .def_property_readonly(
          "d_targets",
          [](SutraTarget &st) -> vector<SutraSource *> & {
            return st.d_targets;
          },
          "Vector of targets")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("__str__",
           [](SutraTarget &st) {
             std::cout << "Source # | position(\") |  Mag | Lambda (mic.)"
                       << std::endl;
             vector<SutraSource *>::iterator it = st.d_targets.begin();
             int32_t i = 0;
             while (it != st.d_targets.end()) {
               std::cout << i << " | "
                         << "(" << (*it)->posx << "," << (*it)->posy << ") | "
                         << (*it)->mag << " | " << (*it)->lambda << std::endl;
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
      ;
}
