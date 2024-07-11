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

//! \file      sensors.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraSensors
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_sensors.hpp>

namespace py = pybind11;

std::unique_ptr<SutraSensors> sensors_init(
    CarmaContext &context, SutraTelescope *d_tel, vector<string> type, int32_t nwfs,
    ArrayFStyle<int64_t> &nxsub, ArrayFStyle<int64_t> &nvalid, ArrayFStyle<int64_t> &npupils, ArrayFStyle<int64_t> &npix, ArrayFStyle<int64_t> &nphase,
    ArrayFStyle<int64_t> &nrebin, ArrayFStyle<int64_t> &nfft, ArrayFStyle<int64_t> &ntot, ArrayFStyle<int64_t> &npup, ArrayFStyle<float> &pdiam,
    ArrayFStyle<float> &nphot, ArrayFStyle<float> &nphot4imat, ArrayFStyle<int32_t> &lgs, ArrayFStyle<bool> &fakecam,
    ArrayFStyle<int32_t> &max_flux_per_pix, ArrayFStyle<int32_t> &max_pix_value, int32_t device, bool roket) {
  return std::unique_ptr<SutraSensors>(new SutraSensors(
      &context, d_tel, type, nwfs, nxsub.mutable_data(), nvalid.mutable_data(), npupils.mutable_data(), npix.mutable_data(), nphase.mutable_data(),
      nrebin.mutable_data(), nfft.mutable_data(), ntot.mutable_data(), npup.mutable_data(), pdiam.mutable_data(), nphot.mutable_data(), nphot4imat.mutable_data(),
      lgs.mutable_data(), fakecam.mutable_data(), max_flux_per_pix.mutable_data(), max_pix_value.mutable_data(), device, roket));
}

int32_t initgs(SutraSensors &ss, ArrayFStyle<float> &xpos, ArrayFStyle<float> &ypos, ArrayFStyle<float> &lambda, ArrayFStyle<float> &mag, float zerop,
             ArrayFStyle<int64_t> &size, ArrayFStyle<float> &noise, ArrayFStyle<int64_t> &seed, ArrayFStyle<float> &G, ArrayFStyle<float> &thetaML,
             ArrayFStyle<float> &dx, ArrayFStyle<float> &dy) {
    return ss.initgs(xpos.mutable_data(), ypos.mutable_data(), lambda.mutable_data(), mag.mutable_data(), zerop, size.mutable_data(), noise.mutable_data(),
                     seed.mutable_data(), G.mutable_data(), thetaML.mutable_data(), dx.mutable_data(), dy.mutable_data());
}

int32_t set_field_stop(SutraSensors &ss, int32_t nwfs, ArrayFStyle<float> &field_stop, int32_t N) {
    return ss.set_field_stop(nwfs, field_stop.mutable_data(), N);
}

void declare_sensors(py::module &mod) {
  py::class_<SutraSensors>(mod, "Sensors")
      .def(py::init(&sensors_init), R"pbdoc(
        Create and initialise a sensors object

    Args:
        context: (CarmaContext) : current carma context

        d_tel: (SutraTelescope) : SutraTelescope object

        type: (list of string): WFS types

        nwfs: (int32_t) : number of WFS

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

        max_flux_per_pix: (np.ndarray[ndim=1, dtype=np.int32]): maximum number of photons a pixel can handle before saturation

        max_pix_value: (np.ndarray[ndim=1, dtype=np.int32]): maximum number of ADU possible in the uint16 image

        device: (int): GPU device index

        roket : (bool): flag for enabling ROKET
        )pbdoc",
           py::arg("context"), py::arg("d_tel"), py::arg("type"),
           py::arg("nwfs"), py::arg("nxsub"), py::arg("nvalid"),
           py::arg("npupils"), py::arg("npix"), py::arg("nphase"),
           py::arg("nrebin"), py::arg("nfft"), py::arg("ntot"), py::arg("npup"),
           py::arg("pdiam"), py::arg("nphot"), py::arg("nphot4imat"),
           py::arg("lgs"), py::arg("fakecam"), py::arg("max_flux_per_pix"),
           py::arg("max_pix_value"), py::arg("device"), py::arg("roket"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "device", [](SutraSensors &ss) { return ss.device; },
          "GPU device index")

      .def_property_readonly(
          "roket", [](SutraSensors &ss) { return ss.roket; }, "ROKET flag")

      .def_property_readonly(
          "nsensors", [](SutraSensors &ss) { return ss.nsensors(); },
          "Number of WFS")

      .def_property_readonly(
          "d_wfs",
          [](SutraSensors &ss) -> vector<SutraWfs *> & { return ss.d_wfs; },
          "Vector of WFS")

      .def_property_readonly(
          "d_camplipup", [](SutraSensors &ss) { return ss.d_camplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_camplifoc", [](SutraSensors &ss) { return ss.d_camplifoc; },
          "Complex amplitude in the focal plane")

      .def_property_readonly(
          "d_fttotim", [](SutraSensors &ss) { return ss.d_fttotim; },
          "Buffer for FFT computation")

      .def_property_readonly(
          "d_ftlgskern", [](SutraSensors &ss) { return ss.d_ftlgskern; },
          "Convolution kernel for LGS spot")

      .def_property_readonly(
          "d_lgskern", [](SutraSensors &ss) { return ss.d_lgskern; },
          "LGS spot")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("initgs",&initgs,
           R"pbdoc(
    Initializes the guide stars of all WFS

    Args:
        xpos: (np.array(ndim=1,dtype=np.float32)): X position of the GSs [arcsec]

        ypos: (np.array(ndim=1,dtype=np.float32)): Y position of the GSs [arcsec]

        lambda_um: (np.array(ndim=1,dtype=np.float32)): Wavelength of the GSs [µm]

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
           py::arg("xpos"), py::arg("ypos"), py::arg("lambda_um"),
           py::arg("mag"), py::arg("zerop"), py::arg("sizes"), py::arg("noise"),
           py::arg("seeds"), py::arg("G"), py::arg("thetaML"), py::arg("dx"),
           py::arg("dy"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

        .def("set_field_stop", &set_field_stop, R"pbdoc(
    Set new field stop for the specified SH WFS

    Args:
        nwfs: (int): WFS index

        field_stop: (np.array(ndim=2,dtype=np.float32)): Field stop to use

        N : (int): Size of the field stop array along one axis
    )pbdoc",
            py::arg("nwfs"), py::arg("field_stop"), py::arg("N"))

      ;
};
