// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      telescope.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraTelescope
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_telescope.h>

namespace py = pybind11;

std::unique_ptr<SutraTelescope> telescope_init(CarmaContext &context,
                                                long n_pup, long npos,
                                                float *pupil, long n_pup_m,
                                                float *pupil_m) {
  return std::unique_ptr<SutraTelescope>(
      new SutraTelescope(&context, n_pup, npos, pupil, n_pup_m, pupil_m));
}

void declare_telescope(py::module &mod) {
  py::class_<SutraTelescope>(mod, "Telescope")
      .def(py::init(wy::colCast(telescope_init)), R"pbdoc(
    Create and initialise a Telescope object

    Args:
        context: (CarmaContext) : current carma context

        n_pup: (long) : spupil size

        npos : (long): number of points in the pupil

        pupil: (np.ndarray[ndim=2, dtype=np.float32_t]) : spupil

        n_pup_m: (long) : mpupil size

        pupil_m: (np.ndarray[ndim=2, dtype=np.float32_t]) : mpupil
        )pbdoc",
           py::arg("context"), py::arg("n_pup"), py::arg("npos"),
           py::arg("pupil"), py::arg("n_pup_m"), py::arg("pupil_m"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](SutraTelescope &sp) { return sp.device; },
          "Device number")

      .def_property_readonly(
          "pup_size", [](SutraTelescope &sp) { return sp.pup_size; },
          "Small Pupil size")

      .def_property_readonly(
          "pup_size_m", [](SutraTelescope &sp) { return sp.pup_size_m; },
          "Medium Pupil size")

      .def_property_readonly(
          "num_eleme_pup", [](SutraTelescope &sp) { return sp.num_eleme_pup; },
          "number of points in the pupil")

      .def_property_readonly(
          "d_pupil", [](SutraTelescope &sp) { return sp.d_pupil; },
          "Small pupil of the Telescope")

      .def_property_readonly(
          "d_pupil_m", [](SutraTelescope &sp) { return sp.d_pupil_m; },
          "Medium pupil of the Telescope")

      .def_property_readonly(
          "d_phase_ab_M1", [](SutraTelescope &sp) { return sp.d_phase_ab_M1; },
          "M1 aberrations on the small pupil")

      .def_property_readonly(
          "d_phase_ab_M1_m",
          [](SutraTelescope &sp) { return sp.d_phase_ab_M1_m; },
          "M1 aberrations on the medium pupil")

      .def_property_readonly(
          "d_input_phase",
          [](SutraTelescope &sp) { return sp.d_input_phase; },
          "Cube of user-defined input phase screens")

      .def_property_readonly(
          "input_phase_counter",
          [](SutraTelescope &sp) { return sp.input_phase_counter; },
          "Index of the current phase screen in the cube d_input_phase")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      //
      .def(
          "update_input_phase", &SutraTelescope::update_input_phase,
          R"pbdoc(
              Update input_phase_counter to take the next phase screen in the circular buffer d_input_phase
          )pbdoc"
      )
      
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def(
          "set_pupil",
          [](SutraTelescope &sp,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (sp.d_pupil->get_nb_elements() == data.size())
              sp.d_pupil->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the small pupil

    Args:
        pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  small pupil
      )pbdoc",
          py::arg("pup"))

      .def(
          "set_pupil_m",
          [](SutraTelescope &sp,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (sp.d_pupil_m->get_nb_elements() == data.size())
              sp.d_pupil_m->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the medium pupil

    Args:
        pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  medium pupil
      )pbdoc",
          py::arg("pup"))

      .def(
          "set_phase_ab_M1",
          [](SutraTelescope &sp,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            return sp.set_phase_ab_M1(data.mutable_data(), data.size());
          },
          R"pbdoc(
    Set the M1 phase aberration in the small pupil

    Args:
        phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the small pupil
      )pbdoc",
          py::arg("phase_ab"))

      .def(
          "set_phase_ab_M1_m",
          [](SutraTelescope &sp,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            return sp.set_phase_ab_M1_m(data.mutable_data(), data.size());
          },
          R"pbdoc(
    Set the M1 phase aberration in the medium pupil

    Args:
        phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the medium pupil
      )pbdoc",
          py::arg("phase_ab"))

      .def(
          "set_input_phase",
          [](SutraTelescope &sp,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            return sp.set_input_phase(data.mutable_data(), data.shape(0) * data.shape(1), data.shape(2));
          },
          R"pbdoc(
        Set a 3D cube of phase screens to be played. Each phase screen is shown to sources as an additional layer to be raytraced. Each phase screen must have the same dimensions as m_pupil
    
    Args:
        input_hase: (np.ndarray[ndim=3,dtype=np.float32_t]) : Cube of input phase screens
      )pbdoc",
          py::arg("phase_ab"))

      ;
};
