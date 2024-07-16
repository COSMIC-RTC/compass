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

//! \file      coronagraph.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraCoronagraph
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_coronagraph.hpp>

namespace py = pybind11;

int32_t set_amplitude(SutraCoronagraph &sc, ArrayFStyle<float> &amplitude) {
    return sc.set_amplitude(amplitude.mutable_data());
}

void declare_coronagraph(py::module &mod) {
  py::class_<SutraCoronagraph>(mod, "Coronagraph")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "device", [](SutraCoronagraph &sc) { return sc.device; }, "GPU device index")

      .def_property_readonly(
          "type", [](SutraCoronagraph &sc) { return sc.type; }, "Coronagraph type")

      .def_property_readonly(
          "cntPsf", [](SutraCoronagraph &sc) { return sc.cntPsf; },
          "Accumulation counter of psf")

      .def_property_readonly(
          "cntImg", [](SutraCoronagraph &sc) { return sc.cntImg; },
          "Accumulation counter of coronagraphic images")

      .def_property_readonly(
          "imageDimx", [](SutraCoronagraph &sc) { return sc.imageDimx; },
          "Coronagraphic image dimension along X-axis")

      .def_property_readonly(
          "imageDimy", [](SutraCoronagraph &sc) { return sc.imageDimy; },
          "Coronagraphic image dimension along Y-axis")

      .def_property_readonly(
          "pupDimx", [](SutraCoronagraph &sc) { return sc.pupDimx; },
          "Pupil dimension along X-axis")

      .def_property_readonly(
          "pupDimy", [](SutraCoronagraph &sc) { return sc.pupDimy; },
          "Pupil dimension along Y-axis")

      .def_property_readonly(
          "wavelength", [](SutraCoronagraph &sc) { return sc.wavelength; },
          "Vector of wavelength used to compute coronagraphic image")

      .def_property_readonly(
          "d_image_se", [](SutraCoronagraph &sc) { return sc.d_image_se; }, "Short exposure coronagraphic image")

      .def_property_readonly(
          "d_image_le", [](SutraCoronagraph &sc) { return sc.d_image_le; }, "Long exposure coronagraphic image")

      .def_property_readonly(
          "d_psf_se", [](SutraCoronagraph &sc) { return sc.d_psf_se; }, "Short exposure coronagraphic psf")

      .def_property_readonly(
          "d_psf_le", [](SutraCoronagraph &sc) { return sc.d_psf_le; }, "Long exposure coronagraphic psf")

      .def_property_readonly(
          "amplitude", [](SutraCoronagraph &sc) { return sc.amplitude; }, "Electric field amplitude in the pupil")

      .def_property_readonly(
          "d_electric_field", [](SutraCoronagraph &sc) { return sc.d_electric_field; }, "Electric field in the pupil")

      .def_property_readonly(
          "d_complex_image", [](SutraCoronagraph &sc) { return sc.d_complex_image; }, "Complex coronagraphic image")

      .def_property_readonly(
          "d_pupil", [](SutraCoronagraph &sc) { return sc.d_pupil; }, "Telescope pupil")

      .def_property_readonly(
          "d_source", [](SutraCoronagraph &sc) { return sc.d_source; },
          "SutraSource used as OPD input")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("compute_image", &SutraCoronagraph::compute_image,
           R"pbdoc(
    Computes the coronagraphic image from the source phase screen

    Args:
        accumulate: (bool, optionnal): If True (default), accumulate the short exposure image in the int64_t exposure one
    )pbdoc",
           py::arg("accumulate") = true)

      .def("compute_psf", &SutraCoronagraph::compute_psf,
           R"pbdoc(
    Computes the psf from the source phase screen

    Args:
        accumulate: (bool, optionnal): If True (default), accumulate the short exposure psf in the int64_t exposure one
    )pbdoc",
           py::arg("accumulate") = true)

      .def("reset", &SutraCoronagraph::reset,
           R"pbdoc(
    Reset int64_t exposure image and counter
    )pbdoc")

      .def("set_amplitude", &set_amplitude,
           R"pbdoc(
    Set the electric field amplitude

    Args:
        amplitude: (np.ndarray[ndim=3, dtype=np.float32]): electric field amplitude
    )pbdoc")

      .def("compute_electric_field", &SutraCoronagraph::compute_electric_field,
           R"pbdoc(
    Computes the electric field from the specified wavelength

    Args:
        wavelengthIndex: (int): Index of the wavelength to use
    )pbdoc");

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
};
