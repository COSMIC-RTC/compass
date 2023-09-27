// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      wfs.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraWfs
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <wyrm>

#include <sutra_wfs.h>

namespace py = pybind11;

void declare_wfs(py::module &mod) {
  py::class_<SutraWfs>(mod, "Wfs")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](SutraWfs &sw) { return sw.device; }, "GPU device index")

      .def_property_readonly(
          "type", [](SutraWfs &sw) { return sw.type; }, "WFS type")

      .def_property_readonly(
          "is_low_order", [](SutraWfs &sw) { return sw.is_low_order; },
          "WFS for low order")

      .def_property_readonly(
          "nxsub", [](SutraWfs &sw) { return sw.nxsub; },
          "Number of ssp in the diameter")

      .def_property_readonly(
          "nvalid", [](SutraWfs &sw) { return sw.nvalid; },
          "Number of valid ssp")

      .def_property_readonly(
          "npix", [](SutraWfs &sw) { return sw.npix; }, "Pixels per ssp")

      .def_property_readonly(
          "nrebin", [](SutraWfs &sw) { return sw.nrebin; }, "Rebin factor")

      .def_property_readonly(
          "nfft", [](SutraWfs &sw) { return sw.nfft; }, "FFT support size")

      .def_property_readonly(
          "ntot", [](SutraWfs &sw) { return sw.ntot; }, "HR support size")

      .def_property_readonly(
          "npup", [](SutraWfs &sw) { return sw.npup; }, "Pupil support size")

      .def_property_readonly(
          "nphase", [](SutraWfs &sw) { return sw.nphase; },
          "Number of phase point per ssp")

      .def_property_readonly(
          "nmaxhr", [](SutraWfs &sw) { return sw.nmaxhr; }, "TODO: docstring")

      .def_property_readonly(
          "nffthr", [](SutraWfs &sw) { return sw.nffthr; }, "TODO: docstring")

      .def_property_readonly(
          "subapd", [](SutraWfs &sw) { return sw.subapd; },
          "ssp diameter in pixels")

      .def_property_readonly(
          "nphot", [](SutraWfs &sw) { return sw.nphot; },
          "Number of photons/ssp/iter")

      .def_property_readonly(
          "nphot4imat", [](SutraWfs &sw) { return sw.nphot4imat; },
          "Number of photons/ssp/iter used for imat computation")

      .def_property_readonly(
          "noise", [](SutraWfs &sw) { return sw.noise; }, "WFS noise [e-]")

      .def_property_readonly(
          "lgs", [](SutraWfs &sw) { return sw.lgs; }, "Is the WFS a LGS one ?")

      .def_property_readonly(
          "kernconv", [](SutraWfs &sw) { return sw.kernconv; },
          "Convolution kernel for spot computation")

      .def_property_readonly(
          "is_low_order", [](SutraWfs &sw) { return sw.is_low_order; },
          "Flag for low order WFS")

      .def_property_readonly(
          "fakecam", [](SutraWfs &sw) { return sw.fakecam; },
          "Flag for uint16 image")

      .def_property_readonly(
          "max_flux_per_pix", [](SutraWfs &sw) { return sw.max_flux_per_pix; },
          "Maximum number of photons allowed before pixel saturation")

      .def_property_readonly(
          "max_pix_value", [](SutraWfs &sw) { return sw.max_pix_value; },
          "Maximum number of ADU allowed in the uint16 image")

      .def_property_readonly(
          "roket", [](SutraWfs &sw) { return sw.roket; },
          "Is the WFS a LGS one ?")

      .def_property_readonly(
          "d_camplipup", [](SutraWfs &sw) { return sw.d_camplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_camplifoc", [](SutraWfs &sw) { return sw.d_camplifoc; },
          "Complex amplitude in the focal plane")

      .def_property_readonly(
          "d_fttotim", [](SutraWfs &sw) { return sw.d_fttotim; },
          "Buffer for FFT computation")

      .def_property_readonly(
          "d_pupil", [](SutraWfs &sw) { return sw.d_pupil; }, "Pupil")

      .def_property_readonly(
          "d_bincube", [](SutraWfs &sw) { return sw.d_bincube; },
          "WFS spots as a 3D array")

      .def_property_readonly(
          "d_binimg", [](SutraWfs &sw) { return sw.d_binimg; }, "WFS image")

      .def_property_readonly(
          "d_binimg_notnoisy",
          [](SutraWfs &sw) { return sw.d_binimg_notnoisy; },
          "WFS image without noise (ROKET only)")

      .def_property_readonly(
          "d_intensities", [](SutraWfs &sw) { return sw.d_intensities; },
          "Sum of intensities in each ssp")

      .def_property_readonly(
          "d_offsets", [](SutraWfs &sw) { return sw.d_offsets; },
          "TODO: docstring")

      .def_property_readonly(
          "d_fluxPerSub", [](SutraWfs &sw) { return sw.d_fluxPerSub; },
          "Normalized flux per ssp")

      .def_property_readonly(
          "d_sincar", [](SutraWfs &sw) { return sw.d_sincar; },
          "TODO: docstring")

      .def_property_readonly(
          "d_hrmap", [](SutraWfs &sw) { return sw.d_hrmap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_camimg", [](SutraWfs &sw) { return sw.d_camimg; },
          "uint16 WFS image")

      .def_property_readonly(
          "d_dark", [](SutraWfs &sw) { return sw.d_dark; }, "Dark WFS frame")

      .def_property_readonly(
          "d_flat", [](SutraWfs &sw) { return sw.d_flat; }, "Flat WFS frame")

      .def_property_readonly(
          "d_slopes", [](SutraWfs &sw) { return sw.d_slopes; },
          "Slopes vector")

      .def_property_readonly(
          "d_gs", [](SutraWfs &sw) { return sw.d_gs; },
          "WGS GS (SutraSource object)")

      .def_property_readonly(
          "d_phasemap", [](SutraWfs &sw) { return sw.d_phasemap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_ttprojmat", [](SutraWfs &sw) { return sw.d_ttprojmat; },
          "TT projection matrix from subap phase to slopes (geom wfs type 2)")

      .def_property_readonly(
          "d_ttprojvec", [](SutraWfs &sw) { return sw.d_ttprojvec; },
          "Input vector for TT projection from subap phase to slopes (geom wfs type 2)")

      .def_property_readonly(
          "d_validsubsx", [](SutraWfs &sw) { return sw.d_validsubsx; },
          "X-position of valid ssp")

      .def_property_readonly(
          "d_validsubsy", [](SutraWfs &sw) { return sw.d_validsubsy; },
          "Y-position of valid ssp")
      .def_property_readonly(
          "d_submask", [](SutraWfs &sw) { return sw.d_submask; },
          "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("comp_image", &SutraWfs::comp_image,
           R"pbdoc(
    Computes the WFS image from the WFS phase

    Args:
        noise: (bool): take noise into account or not
    )pbdoc",
           py::arg("noise") = true)
      .def("slopes_geom",
           wy::colCast((int (SutraWfs::*)(int)) & SutraWfs::slopes_geom),
           R"pbdoc(
    Computes theoretical slopes in wfs.d_slopes

    Args:
        type: (int): method to use (0: reduce, 1: derive)
    )pbdoc",
           py::arg("type") = 0)

      .def("slopes_geom",
           wy::colCast((int (SutraWfs::*)(float *, int)) &
                       SutraWfs::slopes_geom),
           R"pbdoc(
    Computes theoretical slopes in given array

    Args:
        slopes: (np.array(ndim=1, dtype=np.float32)):

        type: (int): method to use (0: reduce, 1: derive)
    )pbdoc",
           py::arg("slopes"), py::arg("type") = 0)

      .def("fill_binimage", &SutraWfs::fill_binimage,
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
          [](SutraWfs &sw,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (data.size() == sw.d_pupil->get_nb_elements())
              sw.set_pupil(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the pupil seen by the WFS

    Args:
        pupil: (np.array(ndim=2,dtype=np.float32)): pupil to set
)pbdoc",
          py::arg("pupil"))

      .def("set_noise", &SutraWfs::set_noise, R"pbdoc(
    Set the noise of the WFS

    Args:
        noise: (float): desired noise (< 0 = no noise
                                        0 = photon only
                                        > 0 = photon + ron in e-)
        seed: (int): seed for the RNG
        )pbdoc",
           py::arg("noise"), py::arg("seed"))

      .def("set_fakecam", &SutraWfs::set_fakecam, R"pbdoc(
    Enable or disable uint16 computation for the WFS

    Args:
        fakecam: (bool): fakecam flag
        )pbdoc",
           py::arg("fakecam"))

      .def("set_max_flux_per_pix", &SutraWfs::set_max_flux_per_pix, R"pbdoc(
    Set the maximum number of photons allowed before pixel saturation

    Args:
        max_flux_per_pix: (int): maximum number of photons allowed before pixel saturation
        )pbdoc",
           py::arg("max_flux_per_pix"))

      .def("set_max_pix_value", &SutraWfs::set_max_pix_value, R"pbdoc(
    Set the maximum number of ADU allowed in the uint16 image

    Args:
        max_pix_value: (int): maximum number of ADU allowed in the uint16 image
        )pbdoc",
           py::arg("max_pix_value"))

      .def("set_binimg", wy::colCast(&SutraWfs::set_binimg), R"pbdoc(
    Set the binimg of the SH WFS

    Args:
        binimg: (np.array[ndim=3, dtype=np.float32]) : cube of subap. images

        nElem: (int): Number of elements in binimg
      )pbdoc",
           py::arg("binimg"), py::arg("nElem"))

      .def("set_dark", wy::colCast(&SutraWfs::set_dark), R"pbdoc(
    Set the dark of the SH WFS

    Args:
        dark: (np.array[ndim=2, dtype=np.float32]) : dark image

        nElem: (int): Number of elements in dark
      )pbdoc",
           py::arg("dark"), py::arg("nElem"))

      .def("set_flat", wy::colCast(&SutraWfs::set_flat), R"pbdoc(
    Set the flat of the SH WFS

    Args:
        flat: (np.array[ndim=2, dtype=np.float32]) : flat image

        nElem: (int): Number of elements in flat
      )pbdoc",
           py::arg("flat"), py::arg("nElem"));
};
