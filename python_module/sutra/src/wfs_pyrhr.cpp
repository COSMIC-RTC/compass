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

//! \file      wfs_pyrhr.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_wfs_pyrhr
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_wfs_pyr_pyrhr.hpp>

namespace py = pybind11;

int32_t set_pyr_modulation_points(SutraWfs_PyrHR &sp, ArrayFStyle<float> &cx, ArrayFStyle<float> &cy, int32_t npts)
{
  return sp.set_pyr_modulation_points(cx.mutable_data(), cy.mutable_data(), npts);
}

int32_t set_pyr_modulation_points_with_weights(SutraWfs_PyrHR &sp, ArrayFStyle<float> &cx, ArrayFStyle<float> &cy, ArrayFStyle<float> &weights, int32_t npts)
{
  return sp.set_pyr_modulation_points(cx.mutable_data(), cy.mutable_data(), weights.mutable_data(), npts);
}

int32_t load_arrays(SutraWfs_PyrHR &sp, ArrayFStyle<std::complex<float>> &halfxy, 
                ArrayFStyle<float> &cx, ArrayFStyle<float> &cy, ArrayFStyle<float> &weights,
                ArrayFStyle<float> &sincar, ArrayFStyle<float> &submask, ArrayFStyle<int32_t> &validsubsx,
                ArrayFStyle<int32_t> &validsubsy, ArrayFStyle<int32_t> &phasemap, ArrayFStyle<float> &fluxPerSub,
                ArrayFStyle<float> &ttprojmat)
{
  return sp.load_arrays(reinterpret_cast<cuFloatComplex *>(halfxy.mutable_data()), cx.mutable_data(), cy.mutable_data(), weights.mutable_data(), sincar.mutable_data(), submask.mutable_data(), validsubsx.mutable_data(), validsubsy.mutable_data(), phasemap.mutable_data(), fluxPerSub.mutable_data(), ttprojmat.mutable_data());
}

int32_t set_phalfxy(SutraWfs_PyrHR &sp, ArrayFStyle<std::complex<float>> &phalfxy)
{
  return sp.set_phalfxy(reinterpret_cast<cuFloatComplex *>(phalfxy.mutable_data()));
}

int32_t set_pyr_mod_weights(SutraWfs_PyrHR &sp, ArrayFStyle<float> &weights, int32_t npts)
{
  return sp.set_pyr_mod_weights(weights.mutable_data(), npts);
}

int32_t copy_valid_pix(SutraWfs_PyrHR &sp, ArrayFStyle<float> &img, ArrayFStyle<int32_t> &validx, ArrayFStyle<int32_t> &validy, int32_t im_dim)
{
  return sp.copy_valid_pix(img.mutable_data(), validx.mutable_data(), validy.mutable_data(), im_dim);
}

void declare_wfs_pyrhr(py::module &mod)
{
  py::class_<SutraWfs_PyrHR, SutraWfs>(mod, "PYRWFS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "npupils", [](SutraWfs_PyrHR &sp)
          { return sp.npupils; },
          "Number of pupil images")

      .def_property_readonly(
          "d_hrimg", [](SutraWfs_PyrHR &sp)
          { return sp.d_hrimg; },
          "TODO: docstring")

      .def_readwrite("compute_pyrfocalplane",
                     &SutraWfs_PyrHR::compute_pyrfocalplane,
                     "TODO: docstring")

      .def_property_readonly(
          "d_pyrfocalplane",
          [](SutraWfs_PyrHR &sp)
          { return sp.d_pyrfocalplane; },
          "TODO: docstring")

      .def_property_readonly(
          "d_psum", [](SutraWfs_PyrHR &sp)
          { return sp.d_psum; },
          "TODO: docstring")

      .def_property_readonly(
          "d_phalfxy", [](SutraWfs_PyrHR &sp)
          { return sp.d_phalfxy; },
          "TODO: docstring")

      .def_property_readonly(
          "d_poffsets", [](SutraWfs_PyrHR &sp)
          { return sp.d_poffsets; },
          "TODO: docstring")

      .def_property_readonly(
          "pyr_cx", [](SutraWfs_PyrHR &sp)
          { return sp.pyr_cx; },
          "TODO: docstring")

      .def_property_readonly(
          "pyr_cy", [](SutraWfs_PyrHR &sp)
          { return sp.pyr_cy; },
          "Modulation points X-positions")

      .def_property_readonly(
          "pyr_mod_weights",
          [](SutraWfs_PyrHR &sp)
          { return sp.pyr_mod_weights; },
          "Ponderation weights for each modulation points")

      .def_property_readonly(
          "d_hrimg_modu",
          [](SutraWfs_PyrHR &sp)
          { return sp.d_hrimg_modu; },
          "Pyramid HR img on each modulation point")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("load_arrays", &load_arrays, R"pbdoc(
  Load PYRHR WFS arrays

  Args:
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

      ttprojmat: (np.array[ndim=2, dtype=np.float32]): slope projection matrix
                 for geom wfs.

    )pbdoc",
           py::arg("halfxy"), py::arg("cx"), py::arg("cy"), py::arg("weights"),
           py::arg("sincar"), py::arg("submask"), py::arg("validsubsx"),
           py::arg("validsubsy"), py::arg("phasemap"), py::arg("fluxPerSub"),
           py::arg("ttprojmat"))

      .def("comp_nphot", &SutraWfs_PyrHR::comp_nphot, R"pbdoc(
    Compute the currect number of photons for a given system

    Args:
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

      .def(
          "set_pyrimg",
          [](SutraWfs_PyrHR &sp, ArrayFStyle<float> data)
          {
            if (data.size() == sp.d_binimg->get_nb_elements())
              sp.d_binimg->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the image of the PWFS

    Args:
        img: (np.array[ndim=2,dtype=np.float32]): new image to set
      )pbdoc",
          py::arg("img"))

      .def("set_pyr_modulation_points", &set_pyr_modulation_points,
           R"pbdoc(
    Set the modulation points of a PWFS

    Args:
        cx: (np.ndarray[ndim=1, dtype=np.float32_t]): X position of modulation points

        cy: (np.ndarray[ndim=1, dtype=np.float32_t]): Y position of modulation points

        npts: (int): number of modulation points
      )pbdoc",
           py::arg("cx"), py::arg("cy"), py::arg("npts"))

      .def("set_pyr_modulation_points", &set_pyr_modulation_points_with_weights,
           R"pbdoc(
    Set the modulation points and weights of a PWFS

    Args:
        cx: (np.ndarray[ndim=1, dtype=np.float32_t]): X position of modulation points

        cy: (np.ndarray[ndim=1, dtype=np.float32_t]): Y position of modulation points

        weights: (np.ndarray[ndim=1, dtype=np.float32_t]): modulation points weights ponderation

        npts: (int): number of modulation points
      )pbdoc",
           py::arg("cx"), py::arg("cy"), py::arg("weights"), py::arg("npts"))

      .def("set_phalfxy",
           &set_phalfxy, R"pbdoc(
    Set the pyramid mask for each modulation point

    Args:
        phalfxy: (np.ndarray[ndim=2, dtype=np.complex64]): pyramid mask for each modulation point
      )pbdoc",
           py::arg("phalfxy"))

      .def("set_pyr_mod_weights",
           &set_pyr_mod_weights, R"pbdoc(
    Set the modulation points weights of a PWFS

    Args:
        weights: (np.ndarray[ndim=1, dtype=np.float32_t]): modulation points weights ponderation

        npts: (int): number of modulation points
      )pbdoc",
           py::arg("weights"), py::arg("npts"))

      .def(
          "set_submask",
          [](SutraWfs_PyrHR &sp, ArrayFStyle<float> data)
          {
            if (data.size() == sp.d_submask->get_nb_elements())
              sp.d_submask->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the field stop of the PWFS

    Args:
        mask: (np.array[ndim=2,dtype=np.float32]): new field stop to set
      )pbdoc",
          py::arg("mask"))

      .def(
          "set_validpix",
          [](SutraWfs_PyrHR &sp, ArrayFStyle<int32_t> &datax, ArrayFStyle<int32_t> &datay)
          {
            if (datax.size() == sp.d_validsubsx->get_nb_elements() &&
                datay.size() == sp.d_validsubsy->get_nb_elements())
            {
              sp.d_validsubsx->host2device(datax.mutable_data());
              sp.d_validsubsy->host2device(datay.mutable_data());
            }
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the valid pixels of the PWFS

    Args:
        datax: (np.array[ndim=2,dtype=np.float32]): new X positions of valid pixels

        datay: (np.array[ndim=2,dtype=np.float32]): new Y positions of valid pixels
      )pbdoc",
          py::arg("datax"), py::arg("datay"))

      .def("copy_valid_pix", &copy_valid_pix,
           R"pbdoc(
    Copy the given pixels on the right place in the binimg of PWFS

    Args:
        data:

        validx:

        validy:

        dim:
      )pbdoc",
           py::arg("data"), py::arg("validx"), py::arg("validy"),
           py::arg("dim"))

      ;
}
