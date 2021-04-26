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

//! \file      target_brahma.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraTargetBrahma
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_target_brahma.h>

namespace py = pybind11;

std::unique_ptr<SutraTargetBrahma> target_brahma_init(
    CarmaContext *context, ACE_TCHAR *name, SutraTelescope *d_tel,
    int subsample, int ntargets, float *xpos, float *ypos, float *lambda,
    float *mag, float zerop, long *sizes, int Npts, int device) {
  return std::unique_ptr<SutraTargetBrahma>(
      new SutraTargetBrahma(context, name, d_tel, subsample, ntargets, xpos,
                              ypos, lambda, mag, zerop, sizes, Npts, device));
}

void declare_target_brahma(py::module &mod) {
  py::class_<SutraTargetBrahma, SutraTarget>(mod, "Target_brahma")
      .def(py::init(wy::colCast(target_brahma_init)), R"pbdoc(
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
