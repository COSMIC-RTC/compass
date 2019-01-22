#include <wyrm>

#include <sutra_centroider_corr.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_centroider_corr(py::module &mod) {
  py::class_<sutra_centroider_corr, sutra_centroider<float>>(mod,
                                                             "CentroiderCORR")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly("npix",
                             [](sutra_centroider_corr &sc) { return sc.npix; },
                             "TODO: docstring")

      .def_property_readonly(
          "interp_sizex",
          [](sutra_centroider_corr &sc) { return sc.interp_sizex; },
          "TODO: docstring")

      .def_property_readonly(
          "interp_sizey",
          [](sutra_centroider_corr &sc) { return sc.interp_sizey; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrfnct", [](sutra_centroider_corr &sc) { return sc.d_corrfnct; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrspot", [](sutra_centroider_corr &sc) { return sc.d_corrspot; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrnorm", [](sutra_centroider_corr &sc) { return sc.d_corrnorm; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrmax", [](sutra_centroider_corr &sc) { return sc.d_corrmax; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corr", [](sutra_centroider_corr &sc) { return sc.d_corr; },
          "TODO: docstring")

      .def_property_readonly(
          "d_interpmat",
          [](sutra_centroider_corr &sc) { return sc.d_interpmat; },
          "TODO: docstring")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("init_corr", wy::colCast(&sutra_centroider_corr::init_corr), R"pbdoc(
            Initializes corr computation

            Parameters
            ------------
            isizex: (int): TODO: docstring
            isizey: (int):
            interpmat: (np.array[ndim=,dtype=np.float32]):
            )pbdoc",
           py::arg("isizex"), py::arg("isizey"), py::arg("interpmat"))

      .def("load_corr", wy::colCast(&sutra_centroider_corr::load_corr), R"pbdoc(
            Load arrays for correlation computation

            Parameters
            ------------
            corr: (np.array[ndim=,dtype=np.float32): TODO: docstring
            corr_norm: (np.array[ndim=,dtype=np.float32):
            ndim: (int):
            )pbdoc",
           py::arg("corr"), py::arg("corr_norm"), py::arg("ndim"))

      .def("set_npix", wy::colCast(&sutra_centroider_corr::set_npix),
           R"pbdoc(
               Set the number of pixels per subap.
            Parameters
            ------------
            npix: (int): number of pixels per subap
            )pbdoc",
           py::arg("npix"))

      ;
};
