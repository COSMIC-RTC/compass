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

//! \file      wfs.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_wfs
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_wfs.h>

namespace py = pybind11;

void declare_wfs(py::module &mod) {
  py::class_<sutra_wfs>(mod, "Wfs")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](sutra_wfs &sw) { return sw.device; }, "GPU device index")

      .def_property_readonly(
          "type", [](sutra_wfs &sw) { return sw.type; }, "WFS type")

      .def_property_readonly(
          "is_low_order", [](sutra_wfs &sw) { return sw.is_low_order; },
          "WFS for low order")

      .def_property_readonly(
          "nxsub", [](sutra_wfs &sw) { return sw.nxsub; },
          "Number of ssp in the diameter")

      .def_property_readonly(
          "nvalid", [](sutra_wfs &sw) { return sw.nvalid; },
          "Number of valid ssp")

      .def_property_readonly(
          "npix", [](sutra_wfs &sw) { return sw.npix; }, "Pixels per ssp")

      .def_property_readonly(
          "nrebin", [](sutra_wfs &sw) { return sw.nrebin; }, "Rebin factor")

      .def_property_readonly(
          "nfft", [](sutra_wfs &sw) { return sw.nfft; }, "FFT support size")

      .def_property_readonly(
          "ntot", [](sutra_wfs &sw) { return sw.ntot; }, "HR support size")

      .def_property_readonly(
          "npup", [](sutra_wfs &sw) { return sw.npup; }, "Pupil support size")

      .def_property_readonly(
          "nphase", [](sutra_wfs &sw) { return sw.nphase; },
          "Number of phase point per ssp")

      .def_property_readonly(
          "nmaxhr", [](sutra_wfs &sw) { return sw.nmaxhr; }, "TODO: docstring")

      .def_property_readonly(
          "nffthr", [](sutra_wfs &sw) { return sw.nffthr; }, "TODO: docstring")

      .def_property_readonly(
          "subapd", [](sutra_wfs &sw) { return sw.subapd; },
          "ssp diameter in pixels")

      .def_property_readonly(
          "nphot", [](sutra_wfs &sw) { return sw.nphot; },
          "Number of photons/ssp/iter")

      .def_property_readonly(
          "nphot4imat", [](sutra_wfs &sw) { return sw.nphot4imat; },
          "Number of photons/ssp/iter used for imat computation")

      .def_property_readonly(
          "noise", [](sutra_wfs &sw) { return sw.noise; }, "WFS noise [e-]")

      .def_property_readonly(
          "lgs", [](sutra_wfs &sw) { return sw.lgs; }, "Is the WFS a LGS one ?")

      .def_property_readonly(
          "kernconv", [](sutra_wfs &sw) { return sw.kernconv; },
          "Convolution kernel for spot computation")

      .def_property_readonly(
          "is_low_order", [](sutra_wfs &sw) { return sw.is_low_order; },
          "Flag for low order WFS")

      .def_property_readonly(
          "fakecam", [](sutra_wfs &sw) { return sw.fakecam; },
          "Flag for uint16 image")

      .def_property_readonly(
          "maxFluxPerPix", [](sutra_wfs &sw) { return sw.maxFluxPerPix; },
          "Maximum number of photons allowed before pixel saturation")

      .def_property_readonly(
          "maxPixValue", [](sutra_wfs &sw) { return sw.maxPixValue; },
          "Maximum number of ADU allowed in the uint16 image")

      .def_property_readonly(
          "roket", [](sutra_wfs &sw) { return sw.roket; },
          "Is the WFS a LGS one ?")

      .def_property_readonly(
          "d_camplipup", [](sutra_wfs &sw) { return sw.d_camplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_camplifoc", [](sutra_wfs &sw) { return sw.d_camplifoc; },
          "Complex amplitude in the focal plane")

      .def_property_readonly(
          "d_fttotim", [](sutra_wfs &sw) { return sw.d_fttotim; },
          "Buffer for FFT computation")

      .def_property_readonly(
          "d_pupil", [](sutra_wfs &sw) { return sw.d_pupil; }, "Pupil")

      .def_property_readonly(
          "d_bincube", [](sutra_wfs &sw) { return sw.d_bincube; },
          "WFS spots as a 3D array")

      .def_property_readonly(
          "d_binimg", [](sutra_wfs &sw) { return sw.d_binimg; }, "WFS image")

      .def_property_readonly(
          "d_binimg_notnoisy",
          [](sutra_wfs &sw) { return sw.d_binimg_notnoisy; },
          "WFS image without noise (ROKET only)")

      .def_property_readonly(
          "d_intensities", [](sutra_wfs &sw) { return sw.d_intensities; },
          "Sum of intensities in each ssp")

      .def_property_readonly(
          "d_offsets", [](sutra_wfs &sw) { return sw.d_offsets; },
          "TODO: docstring")

      .def_property_readonly(
          "d_fluxPerSub", [](sutra_wfs &sw) { return sw.d_fluxPerSub; },
          "Normalized flux per ssp")

      .def_property_readonly(
          "d_sincar", [](sutra_wfs &sw) { return sw.d_sincar; },
          "TODO: docstring")

      .def_property_readonly(
          "d_hrmap", [](sutra_wfs &sw) { return sw.d_hrmap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_camimg", [](sutra_wfs &sw) { return sw.d_camimg; },
          "uint16 WFS image")

      .def_property_readonly(
          "d_dark", [](sutra_wfs &sw) { return sw.d_dark; }, "Dark WFS frame")

      .def_property_readonly(
          "d_flat", [](sutra_wfs &sw) { return sw.d_flat; }, "Flat WFS frame")

      .def_property_readonly(
          "d_slopes", [](sutra_wfs &sw) { return sw.d_slopes; },
          "Slopes vector")

      .def_property_readonly(
          "d_gs", [](sutra_wfs &sw) { return sw.d_gs; },
          "WGS GS (sutra_source object)")

      .def_property_readonly(
          "d_phasemap", [](sutra_wfs &sw) { return sw.d_phasemap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_validsubsx", [](sutra_wfs &sw) { return sw.d_validsubsx; },
          "X-position of valid ssp")

      .def_property_readonly(
          "d_validsubsy", [](sutra_wfs &sw) { return sw.d_validsubsy; },
          "Y-position of valid ssp")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("comp_image", &sutra_wfs::comp_image,
           R"pbdoc(
               Computes the WFS image from the WFS phase

               Parameters
               ------------
               noise: (bool): take noise into account or not
               )pbdoc",
           py::arg("noise") = true)
      .def("slopes_geom",
           wy::colCast((int (sutra_wfs::*)(int)) & sutra_wfs::slopes_geom),
           R"pbdoc(
          Computes theoretical slopes in wfs.d_slopes

          Parameters
          ------------
          type: (int): method to use (0: reduce, 1: derive)
        )pbdoc",
           py::arg("type") = 0)

      .def("slopes_geom",
           wy::colCast((int (sutra_wfs::*)(float *, int)) &
                       sutra_wfs::slopes_geom),
           R"pbdoc(
          Computes theoretical slopes in given array

          Parameters
          ------------
          slopes: (np.array(ndim=1, dtype=np.float32)):
          type: (int): method to use (0: reduce, 1: derive)
        )pbdoc",
           py::arg("slopes"), py::arg("type") = 0)

      .def("fill_binimage", &sutra_wfs::fill_binimage,
           "Fill d_binimg from d_bincube")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def(
          "set_pupil",
          [](sutra_wfs &sw,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (data.size() == sw.d_pupil->getNbElem())
              sw.set_pupil(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
                      Set the pupil seen by the WFS

                      Parameters
                      ------------
                      pupil: (np.array(ndim=2,dtype=np.float32)): pupil to set
                  )pbdoc",
          py::arg("pupil"))

      .def("set_noise", &sutra_wfs::set_noise, R"pbdoc(
            Set the noise of the WFS

            Parameters
            ------------
            noise: (float): desired noise (< 0 = no noise
                                           0 = photon only
                                           > 0 = photon + ron in e-)
            seed: (int): seed for the RNG
        )pbdoc",
           py::arg("noise"), py::arg("seed"))

      .def("set_fakecam", &sutra_wfs::set_fakecam, R"pbdoc(
            Enable or disable uint16 computation for the WFS

            Parameters
            ------------
            fakecam: (bool): fakecam flag
        )pbdoc",
           py::arg("fakecam"))

      .def("set_maxFluxPerPix", &sutra_wfs::set_maxFluxPerPix, R"pbdoc(
            Set the maximum number of photons allowed before pixel saturation

            Parameters
            ------------
            maxFluxPerPix: (int): maximum number of photons allowed before pixel saturation
        )pbdoc",
           py::arg("maxFluxPerPix"))

      .def("set_maxPixValue", &sutra_wfs::set_maxPixValue, R"pbdoc(
            Set the maximum number of ADU allowed in the uint16 image

            Parameters
            ------------
            maxPixValue: (int): maximum number of ADU allowed in the uint16 image
        )pbdoc",
           py::arg("maxPixValue"))

      .def("set_binimg", wy::colCast(&sutra_wfs::set_binimg), R"pbdoc(
        Set the binimg of the SH WFS

        Parameters
        ------------
        binimg: (np.array[ndim=3, dtype=np.float32]) : cube of subap. images
        nElem: (int): Number of elements in binimg
      )pbdoc",
           py::arg("binimg"), py::arg("nElem"))

      .def("set_dark", wy::colCast(&sutra_wfs::set_dark), R"pbdoc(
        Set the dark of the SH WFS

        Parameters
        ------------
        dark: (np.array[ndim=2, dtype=np.float32]) : dark image
        nElem: (int): Number of elements in dark
      )pbdoc",
           py::arg("dark"), py::arg("nElem"))

      .def("set_flat", wy::colCast(&sutra_wfs::set_flat), R"pbdoc(
        Set the flat of the SH WFS

        Parameters
        ------------
        flat: (np.array[ndim=2, dtype=np.float32]) : flat image
        nElem: (int): Number of elements in flat
      )pbdoc",
           py::arg("flat"), py::arg("nElem"));
};
