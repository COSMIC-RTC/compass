// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      centroider_maskedpix.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_centroider_maskedpix
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_centroider_maskedPix.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_maskedPix_impl(py::module &mod, const char *name) {
  using centroider_maskedPix = sutra_centroider_maskedPix<Tin, Tcomp>;

  py::class_<centroider_maskedPix, sutra_centroider<Tin, Tcomp>>(mod, name)

      .def_property_readonly(
          "d_selected_pix",
          [](centroider_maskedPix &sc) { return sc.d_selected_pix; },
          "Selected pixels as an image")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def("get_selected_pix",
           wy::colCast(&centroider_maskedPix::get_selected_pix),
           R"pbdoc(
        Return the given pixels vector in an image format

        Parameters
        ------------
        pix: (np.array[ndim=1,dtype=np.float32]): Pixels to map on the image
    )pbdoc",
           py::arg("pix"));
};
void declare_centroider_maskedPix(py::module &mod) {
  centroider_maskedPix_impl<float, float>(mod, "CentroiderMASKEDPIX_FF");
  centroider_maskedPix_impl<uint16_t, float>(mod, "CentroiderMASKEDPIX_UF");
}
