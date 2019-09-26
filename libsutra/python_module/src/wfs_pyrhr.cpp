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

//! \file      wfs_pyrhr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_wfs_pyrhr
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_wfs_pyr_pyrhr.h>

namespace py = pybind11;

typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;
typedef py::array_t<int, py::array::f_style | py::array::forcecast> F_arrayI;

void declare_wfs_pyrhr(py::module &mod) {
  py::class_<sutra_wfs_pyr_pyrhr, sutra_wfs>(mod, "PYRWFS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "npupils", [](sutra_wfs_pyr_pyrhr &sp) { return sp.npupils; },
          "Number of pupil images")

      .def_property_readonly(
          "d_hrimg", [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_hrimg; },
          "TODO: docstring")

      .def_property_readonly(
          "d_submask", [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_submask; },
          "TODO: docstring")

      .def_readwrite("compute_pyrfocalplane",
                     &sutra_wfs_pyr_pyrhr::compute_pyrfocalplane,
                     "TODO: docstring")

      .def_property_readonly(
          "d_pyrfocalplane",
          [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_pyrfocalplane; },
          "TODO: docstring")

      .def_property_readonly(
          "d_psum", [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_psum; },
          "TODO: docstring")

      .def_property_readonly(
          "d_phalfxy", [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_phalfxy; },
          "TODO: docstring")

      .def_property_readonly(
          "d_poffsets", [](sutra_wfs_pyr_pyrhr &sp) { return sp.d_poffsets; },
          "TODO: docstring")

      .def_property_readonly(
          "pyr_cx", [](sutra_wfs_pyr_pyrhr &sp) { return sp.pyr_cx; },
          "TODO: docstring")

      .def_property_readonly(
          "pyr_cy", [](sutra_wfs_pyr_pyrhr &sp) { return sp.pyr_cy; },
          "Modulation points X-positions")

      .def_property_readonly(
          "pyr_mod_weights",
          [](sutra_wfs_pyr_pyrhr &sp) { return sp.pyr_mod_weights; },
          "Ponderation weights for each modulation points")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("loadarrays", wy::colCast(&sutra_wfs_pyr_pyrhr::loadarrays), R"pbdoc(
      Load PYRHR WFS arrays

      Parameters
      ------------
      halfxy:
      cx:
      cy:
      weights:
      sincar:
      submask:
      validsubsx:
      validsubsy:
      phasemap:
      fluxPerSub:
    )pbdoc",
           py::arg("halfxy"), py::arg("cx"), py::arg("cy"), py::arg("weights"),
           py::arg("sincar"), py::arg("submask"), py::arg("validsubsx"),
           py::arg("validsubsy"), py::arg("phasemap"), py::arg("fluxPerSub"))

      .def("comp_nphot", &sutra_wfs_pyr_pyrhr::comp_nphot, R"pbdoc(
      Compute the currect number of photons for a given system

      Parameters
      ------------
      ittime: (float): 1/loop frequency [s].
      optthroughput: (float): wfs global throughput.
      diam: (float): telescope diameter.
      cobs: (float): telescope central obstruction.
      zerop: (float): (optional for LGS)  detector zero point expressed in ph/m**2/s in the bandwidth of the WFS.
      gsmag: (float): (optional for LGS)  magnitude of guide star.
    )pbdoc",
           py::arg("ittime"), py::arg("optthroughput"), py::arg("diam"),
           py::arg("cobs"), py::arg("zerop") = 0, py::arg("gsmag") = 0)

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def(
          "set_pyrimg",
          [](sutra_wfs_pyr_pyrhr &sp, F_arrayS data) {
            if (data.size() == sp.d_binimg->getNbElem())
              sp.d_binimg->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
        Set the image of the PWFS

        Parameters
        ------------
        img: (np.array[ndim=2,dtype=np.float32]): new image to set
      )pbdoc",
          py::arg("img"))

      .def("set_pyr_modulation",
           wy::colCast((int (sutra_wfs_pyr_pyrhr::*)(float *, float *, int)) &
                       sutra_wfs_pyr_pyrhr::set_pyr_modulation),
           R"pbdoc(
        Set the modulation points of a PWFS

        Parameters
        ------------
        cx: (np.ndarray[ndim=1, dtype=np.float32_t]): X position of modulation points
        cy: (np.ndarray[ndim=1, dtype=np.float32_t]): Y position of modulation points
        npts: (int): number of modulation points
      )pbdoc",
           py::arg("cx"), py::arg("cy"), py::arg("npts"))

      .def("set_pyr_modulation",
           wy::colCast(
               (int (sutra_wfs_pyr_pyrhr::*)(float *, float *, float *, int)) &
               sutra_wfs_pyr_pyrhr::set_pyr_modulation),
           R"pbdoc(
        Set the modulation points and weights of a PWFS

        Parameters
        ------------
        cx: (np.ndarray[ndim=1, dtype=np.float32_t]): X position of modulation points
        cy: (np.ndarray[ndim=1, dtype=np.float32_t]): Y position of modulation points
        weights: (np.ndarray[ndim=1, dtype=np.float32_t]): modulation points weights ponderation
        npts: (int): number of modulation points
      )pbdoc",
           py::arg("cx"), py::arg("cy"), py::arg("weights"), py::arg("npts"))

      .def("set_pyr_mod_weights",
           wy::colCast(&sutra_wfs_pyr_pyrhr::set_pyr_mod_weights), R"pbdoc(
        Set the modulation points weights of a PWFS

        Parameters
        ------------
        weights: (np.ndarray[ndim=1, dtype=np.float32_t]): modulation points weights ponderation
        npts: (int): number of modulation points
      )pbdoc",
           py::arg("weights"), py::arg("npts"))

      .def(
          "set_submask",
          [](sutra_wfs_pyr_pyrhr &sp, F_arrayS data) {
            if (data.size() == sp.d_submask->getNbElem())
              sp.d_submask->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
        Set the field stop of the PWFS

        Parameters
        ------------
        mask: (np.array[ndim=2,dtype=np.float32]): new field stop to set
      )pbdoc",
          py::arg("mask"))

      .def(
          "set_validpix",
          [](sutra_wfs_pyr_pyrhr &sp, F_arrayI datax, F_arrayI datay) {
            if (datax.size() == sp.d_validsubsx->getNbElem() &&
                datay.size() == sp.d_validsubsy->getNbElem()) {
              sp.d_validsubsx->host2device(datax.mutable_data());
              sp.d_validsubsy->host2device(datay.mutable_data());
            } else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
        Set the valid pixels of the PWFS

        Parameters
        ------------
        datax: (np.array[ndim=2,dtype=np.float32]): new X positions of valid pixels
        datay: (np.array[ndim=2,dtype=np.float32]): new Y positions of valid pixels
      )pbdoc",
          py::arg("datax"), py::arg("datay"))

      .def("copyValidPix", wy::colCast(&sutra_wfs_pyr_pyrhr::copyValidPix),
           R"pbdoc(
        Copy the given pixels on the right place in the binimg of PWFS

        Parameters
        ------------
        data:
        validx:
        validy:
        dim:
      )pbdoc",
           py::arg("data"), py::arg("validx"), py::arg("validy"),
           py::arg("dim"))

      ;
};
