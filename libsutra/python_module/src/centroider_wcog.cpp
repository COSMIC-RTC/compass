#include <wyrm>

#include <sutra_centroider_wcog.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_centroider_wcog(py::module &mod) {
  py::class_<sutra_centroider_wcog, sutra_centroider>(mod, "CentroiderWCOG")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly("npix",
                             [](sutra_centroider_wcog &sc) { return sc.npix; },
                             "TODO: docstring")

      .def_property_readonly(
          "d_weights", [](sutra_centroider_wcog &sc) { return sc.d_weights; },
          "Weights applied")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("init_weights", &sutra_centroider_wcog::init_weights,
           "Initializes WCOG computation")

      .def("set_npix", &sutra_centroider_wcog::set_npix, R"pbdoc(
            Set the number of pixels per subap.
            Parameters
            ------------
            npix: (int): number of pixels per subap
            )pbdoc",
           py::arg("npix"))

      .def("load_weights", wy::colCast(&sutra_centroider_wcog::load_weights),
           R"pbdoc(
            Load weights on WCOG

            Parameters
            ------------
            weight: (np.array[ndim=, dtype=np.float32]: weights
            ndim: (int): TODO: docstring
        )pbdoc",
           py::arg("weights"), py::arg("ndim"))

      ;
};
