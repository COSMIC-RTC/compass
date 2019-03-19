#include <wyrm>

#include <sutra_centroider.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_impl(py::module &mod, const char *name) {
  using centroider = sutra_centroider<Tin, Tcomp>;

  py::class_<centroider>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("context",
                             [](centroider &sc) { return sc.current_context; },
                             "GPU context")

      .def_property_readonly("device", [](centroider &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](centroider &sc) { return sc.get_type(); },
                             "Centroider type")

      .def_property_readonly("nslopes",
                             [](centroider &sc) { return sc.nslopes; },
                             "Number of slopes")

      .def_property_readonly("wfs", [](centroider &sc) { return sc.wfs; },
                             "sutra_wfs handled by this centroider")

      .def_property_readonly("nvalid", [](centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("nxsub", [](centroider &sc) { return sc.nxsub; },
                             "Number of ssp across the pupil diameter")

      .def_property_readonly("npix", [](centroider &sc) { return sc.npix; },
                             "Number of pixels along a side of WFS subap.")

      .def_property_readonly("offset", [](centroider &sc) { return sc.offset; },
                             "Offset for centroiding computation")

      .def_property_readonly("scale", [](centroider &sc) { return sc.scale; },
                             "Scale factor to get slopes in arcsec")

      .def_property_readonly("nvalid", [](centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("d_bincube",
                             [](centroider &sc) { return sc.d_bincube; },
                             "Bincube of the WFS image")

      .def_property_readonly("d_intensities",
                             [](centroider &sc) { return sc.d_intensities; },
                             "intensities of the WFS image")

      .def_property_readonly("d_centroids_ref",
                             [](centroider &sc) { return sc.d_centroids_ref; },
                             "Reference slopes vector")

      .def_property_readonly("d_img", [](centroider &sc) { return sc.d_img; },
                             "Calibrated WFS image")

      .def_property_readonly("d_img_raw",
                             [](centroider &sc) { return sc.d_img_raw; },
                             "Raw WFS image")

      .def_property_readonly("d_validx",
                             [](centroider &sc) { return sc.d_validx; },
                             "X positions of the valid ssp")

      .def_property_readonly("d_validy",
                             [](centroider &sc) { return sc.d_validy; },
                             "Y positions of the valid ssp")

      .def_property_readonly("d_dark", [](centroider &sc) { return sc.d_dark; },
                             "Dark frame for calibration")

      .def_property_readonly("d_flat", [](centroider &sc) { return sc.d_flat; },
                             "Flat frame for calibration")

      .def_property_readonly("d_validMask",
                             [](centroider &sc) {
                               sc.get_validMask();
                               return sc.d_validMask;
                             },
                             "Flat frame for calibration")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("get_cog",
           wy::colCast((int (centroider::*)(void)) & centroider::get_cog),
           "Computes centroids and stores it in d_slopes of the WFS")

      .def("load_validpos", wy::colCast(&centroider::load_validpos),
           R"pbdoc(
        Load the validx and validy arrays

        Parameters
        ------------
        validx: (np.array[ndim=1,dtype=np.float32]): X positions of the valid ssp
        validy: (np.array[ndim=1,dtype=np.float32]): Y positions of the valid ssp
        N: (int): arrays size
    )pbdoc",
           py::arg("validx"), py::arg("validy"), py::arg("N"))

      .def("set_npix", wy::colCast(&centroider::set_npix),
           R"pbdoc(
            Set the number of pixels per subap for a RTC standalone

        Parameters
        ------------
        npix: (int): number of pixel along a subap. side
    )pbdoc",
           py::arg("npix"))

      .def("set_nxsub", wy::colCast(&centroider::set_nxsub),
           R"pbdoc(
            Set the number of ssp across the pupil diameter for a RTC standalone

        Parameters
        ------------
        nxsub: (int): number of ssp across the pupil diameter
    )pbdoc",
           py::arg("nxsub"))

      .def("load_img", wy::colCast(&centroider::load_img), R"pbdoc(
            Load a SH image in a RTC standalone (host to device)

        Parameters
        ------------
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image
        n: (int): Image support size
    )pbdoc",
           py::arg("img"), py::arg("n"))

      .def("load_img_gpu", &centroider::load_img_gpu, R"pbdoc(
            Load a SH image in a RTC standalone (device to device)

        Parameters
        ------------
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image
        n: (int): Image support size
    )pbdoc",
           py::arg("img"))

      .def("calibrate_img", &centroider::calibrate_img, R"pbdoc(
           Performs the raw WFS frame calibration
           )pbdoc")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_centroids_ref", wy::colCast(&centroider::set_centroids_ref),
           R"pbdoc(
      Set the references slopes

      Parameters
      ------------
      refslopes: (np.array[ndim1,dtype=np.float32]): reference slopes to set
    )pbdoc",
           py::arg("refslopes"))

      .def("set_scale", &centroider::set_scale, R"pbdoc(
        Set the controider scale factor

        Parameters
        ------------
        scale: (float): new scale factor
    )pbdoc",
           py::arg("scale"))

      .def("set_dark", wy::colCast(&centroider::set_dark), R"pbdoc(
        Set the dark frame for calibration

        Parameters
        ------------
        dark: (np.ndarray[ndim=2, dtype=np.float32_t): dark frame
        n: (int): image support size
    )pbdoc",
           py::arg("dark"), py::arg("n"))

      .def("set_flat", wy::colCast(&centroider::set_flat), R"pbdoc(
        Set the flat frame for calibration

        Parameters
        ------------
        flat: (np.ndarray[ndim=2, dtype=np.float32_t): flat frame
        n: (int): image support size
    )pbdoc",
           py::arg("flat"), py::arg("n"));
};

void declare_centroider(py::module &mod) {
  centroider_impl<float, float>(mod, "Centroider_FF");
  centroider_impl<uint16_t, float>(mod, "Centroider_UF");
#ifdef CAN_DO_HALF
  centroider_impl<float, half>(mod, "Centroider_FH");
  centroider_impl<uint16_t, half>(mod, "Centroider_UH");

#endif
}
