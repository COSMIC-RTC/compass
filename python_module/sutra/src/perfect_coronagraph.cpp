// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      perfect_coronagraph.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraPerfectCoronagraph
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_perfectCoronagraph.hpp>

namespace py = pybind11;

std::unique_ptr<SutraPerfectCoronagraph> perfect_coronagraph_init(CarmaContext &context, SutraSource *d_source,int32_t im_dimx,
                                int32_t im_dimy, ArrayFStyle<float> &wavelength, int32_t nWavelength, int32_t device) {
  return std::unique_ptr<SutraPerfectCoronagraph>(new SutraPerfectCoronagraph(&context, d_source,
                                                    im_dimx, im_dimy, wavelength.mutable_data(),
                                                    nWavelength, device));
};

int32_t set_mft(SutraPerfectCoronagraph &spc, ArrayFStyle<std::complex<float>> &A, ArrayFStyle<std::complex<float>> &B, ArrayFStyle<float> &norm, std::string &mftType) {
    return spc.set_mft(reinterpret_cast<cuFloatComplex*>(A.mutable_data()), reinterpret_cast<cuFloatComplex*>(B.mutable_data()), norm.mutable_data(), mftType);
}


void declare_perfect_coronagraph(py::module &mod) {
  py::class_<SutraPerfectCoronagraph, SutraCoronagraph>(mod, "PerfectCoronagraph")

      .def(py::init(&perfect_coronagraph_init), R"pbdoc(
        Instantiates a PerfectCoronagraph object

        Args:
            context: (CarmaContext): context

            d_source: (SutraSource): Coronagraph source input

            im_dimx: (int): Coronagraphic image dimension along x axis

            im_dimy: (int): Coronagraphic image dimension along y axis

            wavelength: (np.ndarray[ndim=1, dtype=np.float32]): vector of wavelengths

            nWavelength: (int): number of wavelength

            device: (int): GPU device index
      )pbdoc",
      py::arg("context"), py::arg("d_source"), py::arg("im_dimx"), py::arg("im_dimy"),
      py::arg("wavelength"), py::arg("nWavelength"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "AA", [](SutraPerfectCoronagraph &sc) { return sc.AA; }, "A MFT matrix")
      .def_property_readonly(
          "BB", [](SutraPerfectCoronagraph &sc) { return sc.BB; }, "B MFT matrix")
      .def_property_readonly(
          "norm", [](SutraPerfectCoronagraph &sc) { return sc.norm; }, "MFT normalization")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("set_mft", &set_mft,
           R"pbdoc(
    Set MFT matrices for coronagraphic image computation

    Args:
        A : (np.ndarray[dtype=np.complex32, ndims=3]): A MFT matrix for each wavelength

        B : (np.ndarray[dtype=np.complex32, ndims=3]): B MFT matrix for each wavelength

        norm : (np.ndarray[dtype=np.complex32, ndims=3]): MFT normalization for each wavelength

        mft_type : (str): MFT matrices to set, i.e. "img" or "psf"

    )pbdoc",
           py::arg("A"), py::arg("B"), py::arg("norm"), py::arg("mft_type"));

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
};
