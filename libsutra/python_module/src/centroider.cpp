#include <wyrm>

#include <sutra_centroider.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_centroider(py::module &mod) {
  py::class_<sutra_centroider<float>>(mod, "Centroider")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "context",
          [](sutra_centroider<float> &sc) { return sc.current_context; },
          "GPU context")

      .def_property_readonly(
          "device", [](sutra_centroider<float> &sc) { return sc.device; },
          "GPU device index")

      .def_property_readonly(
          "type", [](sutra_centroider<float> &sc) { return sc.get_type(); },
          "Centroider type")

      .def_property_readonly(
          "nslopes", [](sutra_centroider<float> &sc) { return sc.nslopes; },
          "Number of slopes")

      .def_property_readonly("wfs",
                             [](sutra_centroider<float> &sc) { return sc.wfs; },
                             "sutra_wfs handled by this centroider")

      .def_property_readonly(
          "nvalid", [](sutra_centroider<float> &sc) { return sc.nvalid; },
          "Number of valid ssp of the WFS")

      .def_property_readonly(
          "nxsub", [](sutra_centroider<float> &sc) { return sc.nxsub; },
          "Number of ssp across the pupil diameter")

      .def_property_readonly(
          "npix", [](sutra_centroider<float> &sc) { return sc.npix; },
          "Number of pixels along a side of WFS subap.")

      .def_property_readonly(
          "offset", [](sutra_centroider<float> &sc) { return sc.offset; },
          "Offset for centroiding computation")

      .def_property_readonly(
          "scale", [](sutra_centroider<float> &sc) { return sc.scale; },
          "Scale factor to get slopes in arcsec")

      .def_property_readonly(
          "nvalid", [](sutra_centroider<float> &sc) { return sc.nvalid; },
          "Number of valid ssp of the WFS")

      .def_property_readonly(
          "d_bincube", [](sutra_centroider<float> &sc) { return sc.d_bincube; },
          "Bincube of the WFS image")

      .def_property_readonly(
          "d_intensities",
          [](sutra_centroider<float> &sc) { return sc.d_intensities; },
          "intensities of the WFS image")

      .def_property_readonly(
          "d_centroids_ref",
          [](sutra_centroider<float> &sc) { return sc.d_centroids_ref; },
          "Reference slopes vector")

      .def_property_readonly(
          "d_img", [](sutra_centroider<float> &sc) { return sc.d_img; },
          "Calibrated WFS image")

      .def_property_readonly(
          "d_img_raw", [](sutra_centroider<float> &sc) { return sc.d_img_raw; },
          "Raw WFS image")

      .def_property_readonly(
          "d_validx", [](sutra_centroider<float> &sc) { return sc.d_validx; },
          "X positions of the valid ssp")

      .def_property_readonly(
          "d_validy", [](sutra_centroider<float> &sc) { return sc.d_validy; },
          "Y positions of the valid ssp")

      .def_property_readonly(
          "d_dark", [](sutra_centroider<float> &sc) { return sc.d_dark; },
          "Dark frame for calibration")

      .def_property_readonly(
          "d_flat", [](sutra_centroider<float> &sc) { return sc.d_flat; },
          "Flat frame for calibration")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("get_cog",
           wy::colCast((int (sutra_centroider<float>::*)(void)) &
                       sutra_centroider<float>::get_cog),
           "Computes centroids and stores it in d_slopes of the WFS")

      .def("get_cog",
           wy::colCast(
               (int (sutra_centroider<float>::*)(float *, float *, bool)) &
               sutra_centroider<float>::get_cog),
           R"pbdoc(
        Computes centroids and stores it in the given array

        Parameters
        ------------
        intensities: (np.array[ndim=1,dtype=np.float32]): array where to compute the sum of each ssp
        slopes: (np.array[ndim=1,dtype=np.float32]): array where to store centroids
        noise: (bool): Computing the centroids of the image with or without noise
        )pbdoc",
           py::arg("intensities"), py::arg("slopes"), py::arg("noise"))

      .def("load_validpos",
           wy::colCast(&sutra_centroider<float>::load_validpos),
           R"pbdoc(
        Load the validx and validy arrays

        Parameters
        ------------
        validx: (np.array[ndim=1,dtype=np.float32]): X positions of the valid ssp
        validy: (np.array[ndim=1,dtype=np.float32]): Y positions of the valid ssp
        N: (int): arrays size
    )pbdoc",
           py::arg("validx"), py::arg("validy"), py::arg("N"))

      .def("set_npix", wy::colCast(&sutra_centroider<float>::set_npix),
           R"pbdoc(
            Set the number of pixels per subap for a RTC standalone

        Parameters
        ------------
        npix: (int): number of pixel along a subap. side
    )pbdoc",
           py::arg("npix"))

      .def("set_nxsub", wy::colCast(&sutra_centroider<float>::set_nxsub),
           R"pbdoc(
            Set the number of ssp across the pupil diameter for a RTC standalone

        Parameters
        ------------
        nxsub: (int): number of ssp across the pupil diameter
    )pbdoc",
           py::arg("nxsub"))

      .def("load_img", wy::colCast(&sutra_centroider<float>::load_img), R"pbdoc(
            Load a SH image in a RTC standalone (host to device)

        Parameters
        ------------
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image
        n: (int): Image support size
    )pbdoc",
           py::arg("img"), py::arg("n"))

      .def("calibrate_img", &sutra_centroider<float>::calibrate_img, R"pbdoc(
           Performs the raw WFS frame calibration

           Parameters
           ------------
           save_raw: (bool): (default=false) If True, saves the raw WFS image in d_img_raw before calibration
           )pbdoc",
           py::arg("save_raw"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_centroids_ref",
           wy::colCast(&sutra_centroider<float>::set_centroids_ref), R"pbdoc(
      Set the references slopes

      Parameters
      ------------
      refslopes: (np.array[ndim1,dtype=np.float32]): reference slopes to set
    )pbdoc",
           py::arg("refslopes"))

      .def("set_scale", &sutra_centroider<float>::set_scale, R"pbdoc(
        Set the controider scale factor

        Parameters
        ------------
        scale: (float): new scale factor
    )pbdoc",
           py::arg("scale"))

      .def("set_dark", wy::colCast(&sutra_centroider<float>::set_dark), R"pbdoc(
        Set the dark frame for calibration

        Parameters
        ------------
        dark: (np.ndarray[ndim=2, dtype=np.float32_t): dark frame
        n: (int): image support size
    )pbdoc",
           py::arg("dark"), py::arg("n"))

      .def("set_flat", wy::colCast(&sutra_centroider<float>::set_flat), R"pbdoc(
        Set the flat frame for calibration

        Parameters
        ------------
        flat: (np.ndarray[ndim=2, dtype=np.float32_t): flat frame
        n: (int): image support size
    )pbdoc",
           py::arg("flat"), py::arg("n"));
};
