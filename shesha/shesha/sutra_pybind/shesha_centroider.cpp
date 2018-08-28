#include <wyrm>

#include <sutra_centroider.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_shesha_centroider(py::module &mod) {
  py::class_<sutra_centroider>(mod, "Centroider")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly("device",
                             [](sutra_centroider &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](sutra_centroider &sc) { return sc.get_type(); },
                             "Centroider type")

      .def_property_readonly("nslopes",
                             [](sutra_centroider &sc) { return sc.nslopes; },
                             "Number of slopes")
      .def_property_readonly("wfs", [](sutra_centroider &sc) { return sc.wfs; },
                             "sutra_wfs handled by this centroider")

      .def_property_readonly("nvalid",
                             [](sutra_centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("offset",
                             [](sutra_centroider &sc) { return sc.offset; },
                             "Offset for centroiding computation")

      .def_property_readonly("scale",
                             [](sutra_centroider &sc) { return sc.scale; },
                             "Scale factor to get slopes in arcsec")

      .def_property_readonly("nvalid",
                             [](sutra_centroider &sc) { return sc.nvalid; },
                             "Number of valid ssp of the WFS")

      .def_property_readonly("d_bincube",
                             [](sutra_centroider &sc) { return sc.d_bincube; },
                             "Bincube of the WFS image")

      .def_property_readonly("d_subsum",
                             [](sutra_centroider &sc) { return sc.d_subsum; },
                             "subsum of the WFS image")

      .def_property_readonly(
          "d_img", [](sutra_centroider &sc) { return sc.d_img; }, "WFS image")

      .def_property_readonly("d_validx",
                             [](sutra_centroider &sc) { return sc.d_validx; },
                             "X positions of the valid ssp")

      .def_property_readonly("d_validy",
                             [](sutra_centroider &sc) { return sc.d_validy; },
                             "Y positions of the valid ssp")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("get_cog",
           wy::colCast((int (sutra_centroider::*)(void)) &
                       sutra_centroider::get_cog),
           "Computes centroids and stores it in d_slopes of the WFS")

      .def("get_cog",
           wy::colCast((int (sutra_centroider::*)(float *, float *, bool)) &
                       sutra_centroider::get_cog),
           R"pbdoc(
        Computes centroids and stores it in the given array

        Parameters
        ------------
        subsum: (np.array[ndim=1,dtype=np.float32]): array where to compute the sum of each ssp
        slopes: (np.array[ndim=1,dtype=np.float32]): array where to store centroids
        noise: (bool): Computing the centroids of the image with or without noise
        )pbdoc",
           py::arg("subsum"), py::arg("slopes"), py::arg("noise"))

      .def("load_validpos", wy::colCast(&sutra_centroider::load_validpos),
           R"pbdoc(
        Load the validx and validy arrays

        Parameters
        ------------
        validx: (np.array[ndim=1,dtype=np.float32]): X positions of the valid ssp
        validy: (np.array[ndim=1,dtype=np.float32]): Y positions of the valid ssp
        N: (int): arrays size
    )pbdoc",
           py::arg("validx"), py::arg("validy"), py::arg("N"))

      .def("fill_bincube", wy::colCast(&sutra_centroider::fill_bincube),
           R"pbdoc(
            Fill the bincube from the previously loaded image.
            Only use it with a SH RTC standalone, after using load_img

        Parameters
        ------------
        npix: (int): number of pixel along a subap. side
    )pbdoc",
           py::arg("npix"))

      .def("load_img", wy::colCast(&sutra_centroider::load_img), R"pbdoc(
            Load a SH image in a RTC standalone (host to device)

        Parameters
        ------------
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image
        n: (int): Image support size
    )pbdoc",
           py::arg("img"), py::arg("n"))

      .def("load_img_gpu", wy::colCast(&sutra_centroider::load_img_gpu),
           R"pbdoc(
            Load a SH image in a RTC standalone (device to device)

        Parameters
        ------------
        img: (np.ndarray[ndim=2, dtype=np.float32_t]): SH image
        n: (int): Image support size
    )pbdoc",
           py::arg("img"), py::arg("n"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_scale", &sutra_centroider::set_scale, R"pbdoc(
        Set the controider scale factor

        Parameters
        ------------
        scale: (float): new scale factor
    )pbdoc",
           py::arg("scale"));
};
