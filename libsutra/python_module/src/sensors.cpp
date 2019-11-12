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

//! \file      sensors.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_sensors
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_sensors.h>

namespace py = pybind11;

std::unique_ptr<sutra_sensors> sensors_init(
    carma_context &context, sutra_telescope *d_tel, vector<string> type,
    int nwfs, long *nxsub, long *nvalid, long *npupils, long *npix,
    long *nphase, long *nrebin, long *nfft, long *ntot, long *npup,
    float *pdiam, float *nphot, float *nphot4imat, int *lgs, bool *fakecam,
    int *maxFluxPerPix, int *maxPixValue, int device, bool roket) {
  return std::unique_ptr<sutra_sensors>(new sutra_sensors(
      &context, d_tel, type, nwfs, nxsub, nvalid, npupils, npix, nphase, nrebin,
      nfft, ntot, npup, pdiam, nphot, nphot4imat, lgs, fakecam, maxFluxPerPix,
      maxPixValue, device, roket));
}

void declare_sensors(py::module &mod) {
  py::class_<sutra_sensors>(mod, "Sensors")
      .def(py::init(wy::colCast(sensors_init)), R"pbdoc(
        Create and initialise a sensors object
        Parameters
        ------------
        context: (carma_context) : current carma context
        d_tel: (sutra_telescope) : sutra_telescope object
        type: (list of string): WFS types
        nwfs: (int) : number of WFS
        nxsub: (np.ndarray[ndim=1, dtype=np.int64]) : number of ssp in the diameter for each WFS
        nvalid: (np.ndarray[ndim=1, dtype=np.int64]) : number of valid ssp for each WFS
        npupils: (np.ndarray[ndim=1, dtype=np.int64]) : number of pupil images for each WFS
        npix: (np.ndarray[ndim=1,dtype=np.int64]) : number of pix per ssp for each WFS
        nphase: (np.ndarray[ndim=1,dtype=np.int64]) : number of phase points per ssp for each WFS
        nrebin: (np.ndarray[ndim=1,dtype=np.int64]) : rebin factor for each WFS
        nfft: (np.ndarray[ndim=1,dtype=np.int64]) : FFT support size for each WFS
        ntot: (np.ndarray[ndim=1,dtype=np.int64]) : HR support size for each WFS
        npup: (np.ndarray[ndim=1,dtype=np.int64]) : Pupil support size for each WFS
        pdiam: (np.ndarray[ndim=1,dtype=np.float32]) : ssp diameter in pixels for each WFS
        nphot: (np.ndarray[ndim=1,dtype=np.float32]) : photons per subap per iter for each WFS
        nphot4imat: (np.ndarray[ndim=1,dtype=np.float32]) : photons per subap per iter for each WFS (for imat computation only)
        lgs: (np.ndarray[ndim=1,dtype=np.int64]) : LGS flag for each WFS
        fakecam: (bool): if True, image is computed in uint16
        maxFluxPerPix: (np.ndarray[ndim=1, dtype=np.int32]): maximum number of photons a pixel can handle before saturation
        maxPixValue: (np.ndarray[ndim=1, dtype=np.int32]): maximum number of ADU possible in the uint16 image
        device: (int): GPU device index
        roket : (bool): flag for enabling ROKET

        )pbdoc",
           py::arg("context"), py::arg("d_tel"), py::arg("type"),
           py::arg("nwfs"), py::arg("nxsub"), py::arg("nvalid"),
           py::arg("npupils"), py::arg("npix"), py::arg("nphase"),
           py::arg("nrebin"), py::arg("nfft"), py::arg("ntot"), py::arg("npup"),
           py::arg("pdiam"), py::arg("nphot"), py::arg("nphot4imat"),
           py::arg("lgs"), py::arg("fakecam"), py::arg("maxFluxPerPix"),
           py::arg("maxPixValue"), py::arg("device"), py::arg("roket"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](sutra_sensors &ss) { return ss.device; },
          "GPU device index")

      .def_property_readonly(
          "roket", [](sutra_sensors &ss) { return ss.roket; }, "ROKET flag")

      .def_property_readonly(
          "nsensors", [](sutra_sensors &ss) { return ss.nsensors(); },
          "Number of WFS")

      .def_property_readonly(
          "d_wfs",
          [](sutra_sensors &ss) -> vector<sutra_wfs *> & { return ss.d_wfs; },
          "Vector of WFS")

      .def_property_readonly(
          "d_camplipup", [](sutra_sensors &ss) { return ss.d_camplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_camplifoc", [](sutra_sensors &ss) { return ss.d_camplifoc; },
          "Complex amplitude in the focal plane")

      .def_property_readonly(
          "d_fttotim", [](sutra_sensors &ss) { return ss.d_fttotim; },
          "Buffer for FFT computation")

      .def_property_readonly(
          "d_ftlgskern", [](sutra_sensors &ss) { return ss.d_ftlgskern; },
          "Convolution kernel for LGS spot")

      .def_property_readonly(
          "d_lgskern", [](sutra_sensors &ss) { return ss.d_lgskern; },
          "LGS spot")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("initgs",
           wy::colCast(
               (int (sutra_sensors::*)(float *, float *, float *, float *,
                                       float, long *, float *, long *, float *,
                                       float *, float *, float *)) &
               sutra_sensors::initgs),
           R"pbdoc(
                         Initializes the guide stars of all WFS

                         Parameters
                         ------------
                         xpos: (np.array(ndim=1,dtype=np.float32)): X position of the GSs [arcsec]
                         ypos: (np.array(ndim=1,dtype=np.float32)): Y position of the GSs [arcsec]
                         lambda: (np.array(ndim=1,dtype=np.float32)): Wavelength of the GSs [µm]
                         mag: (np.array(ndim=1,dtype=np.float32)): Magnitude of the GSs
                         zerop: (float): Flux at magnitude 0
                         sizes: (np.array(ndim=1,dtype=np.int64)): Support size of the GSs
                         noise: (np.array(ndim=1,dtype=np.float32)): Noise of the WFS [e-]
                         seeds: (np.array(ndim=1,dtype=np.int64)): seeds for noise generation
                         G: (np.array(ndim=1,dtype=np.float32)): Magnification factors for WFS misalignment
                         thetaML: (np.array(ndim=1,dtype=np.float32)): Pupil rotation angle for WFS misalignment
                         dx: (np.array(ndim=1,dtype=np.float32)): X axis misalignment for WFS
                         dy: (np.array(ndim=1,dtype=np.float32)): Y axis misalignment for WFS
                     )pbdoc",
           py::arg("xpos"), py::arg("ypos"), py::arg("lambda"), py::arg("mag"),
           py::arg("zerop"), py::arg("sizes"), py::arg("noise"),
           py::arg("seeds"), py::arg("G"), py::arg("thetaML"), py::arg("dx"),
           py::arg("dy"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      ;
};
