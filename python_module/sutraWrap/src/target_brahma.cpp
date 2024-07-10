// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      target_brahma.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraTargetBrahma
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_target_brahma.hpp>

namespace py = pybind11;

std::unique_ptr<SutraTargetBrahma> target_brahma_init(
    CarmaContext *context, ACE_TCHAR *name, SutraTelescope *d_tel,
    int32_t subsample, int32_t ntargets, ArrayFStyle<float> &xpos, ArrayFStyle<float> &ypos, ArrayFStyle<float> &lambda,
    ArrayFStyle<float> &mag, float zerop, int64_t *sizes, int32_t Npts, int32_t device) {
  return std::unique_ptr<SutraTargetBrahma>(
      new SutraTargetBrahma(context, name, d_tel, subsample, ntargets, xpos.mutable_data(), ypos.mutable_data(),
                      lambda.mutable_data(), mag.mutable_data(), zerop, sizes.mutable_data(), Npts,
                      device));
}

void declare_target_brahma(py::module &mod) {
  py::class_<SutraTargetBrahma, SutraTarget>(mod, "Target_brahma")
      .def(py::init(&target_brahma_init), R"pbdoc(
    Create and initialise a brahma target object

    Args:
        context: (CarmaContext) : current carma context

        name:

        subsample:

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
           py::arg("context"), py::arg("name"), py::arg("subsample"),
           py::arg("d_tel"), py::arg("ntargets"), py::arg("xpos"),
           py::arg("ypos"), py::arg("lambda_um"), py::arg("mag"), py::arg("zerop"),
           py::arg("sizes"), py::arg("Npts"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      // .def_property_readonly("framecounter",
      //                        [](SutraTargetBrahma &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("samplecounter",
      //                        [](SutraTargetBrahma &st) {
      //                          return st.samplecounter;
      //                        },
      //                        "Sample counter")

      // .def_property_readonly("subsample",
      //                        [](SutraTargetBrahma &st) {
      //                          return st.subsample;
      //                        },
      //                        "Subsample")

      // .def_property_readonly("is_initialised",
      //                        [](SutraTargetBrahma &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "is_initialised flag")
      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &SutraTargetBrahma::publish)

      ;
};
