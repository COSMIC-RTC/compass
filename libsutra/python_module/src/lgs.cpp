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

//! \file      lgs.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraLGS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_lgs.h>

namespace py = pybind11;

void declare_lgs(py::module &mod) {
  py::class_<SutraLGS>(mod, "LGS")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](SutraLGS &sl) { return sl.device; }, "GPU device index")

      .def_property_readonly(
          "nvalid", [](SutraLGS &sl) { return sl.nvalid; }, "TODO: docstring")

      .def_property_readonly(
          "npix", [](SutraLGS &sl) { return sl.npix; }, "TODO: docstring")

      .def_property_readonly(
          "nmaxhr", [](SutraLGS &sl) { return sl.nmaxhr; },
          "Size of HR support")

      .def_property_readonly(
          "hg", [](SutraLGS &sl) { return sl.hg; }, "TODO: docstring")

      .def_property_readonly(
          "h0", [](SutraLGS &sl) { return sl.h0; }, "TODO: docstring")

      .def_property_readonly(
          "deltah", [](SutraLGS &sl) { return sl.deltah; }, "TODO: docstring")

      .def_property_readonly(
          "pixsize", [](SutraLGS &sl) { return sl.pixsize; },
          "Pixel size on sky[arcsec]")

      .def_property_readonly(
          "nprof", [](SutraLGS &sl) { return sl.nprof; }, "TODO: docstring")

      .def_property_readonly(
          "d_doffaxis", [](SutraLGS &sl) { return sl.d_doffaxis; },
          "TODO: docstring")

      .def_property_readonly(
          "d_azimuth", [](SutraLGS &sl) { return sl.d_azimuth; },
          "TODO: docstring")

      .def_property_readonly(
          "d_prof1d", [](SutraLGS &sl) { return sl.d_prof1d; },
          "TODO: docstring")

      .def_property_readonly(
          "d_profcum", [](SutraLGS &sl) { return sl.d_profcum; },
          "TODO: docstring")

      .def_property_readonly(
          "d_prof2d", [](SutraLGS &sl) { return sl.d_prof2d; },
          "TODO: docstring")

      .def_property_readonly(
          "d_beam", [](SutraLGS &sl) { return sl.d_beam; }, "TODO: docstring")

      .def_property_readonly(
          "d_ftbeam", [](SutraLGS &sl) { return sl.d_ftbeam; },
          "TODO: docstring")

      .def_property_readonly(
          "d_lgskern", [](SutraLGS &sl) { return sl.d_lgskern; },
          "TODO: docstring")

      .def_property_readonly(
          "d_ftlgskern", [](SutraLGS &sl) { return sl.d_ftlgskern; },
          "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("lgs_init", wy::colCast(&SutraLGS::lgs_init), R"pbdoc(
    Initialize LGS object

    Args:
        nprof: (int): TODO: docstring

        hg: (float):

        h0: (float):

        deltah: (float):

        pixsize: (float):

        doffaxis:(np.array[ndim= , dtype=np.float32]):

        prof1d:(np.array[ndim= , dtype=np.float32]):

        profcum:(np.array[ndim= , dtype=np.float32]):

        beam:(np.array[ndim= , dtype=np.float32]):

        ftbeam:(np.array[ndim= , dtype=np.complex64]):

        azimuth:(np.array[ndim= , dtype=np.float32]):
    )pbdoc",
           py::arg("nprof"), py::arg("hg"), py::arg("h0"), py::arg("deltah"),
           py::arg("pixsize"), py::arg("doffaxis"), py::arg("prof1d"),
           py::arg("profcum"), py::arg("beam"), py::arg("ftbeam"),
           py::arg("azimuth"))
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      ;
};
