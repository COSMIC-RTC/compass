#include <wyrm>

#include <sutra_centroider_wcog.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_wcog_impl(py::module &mod, const char *name) {
  using centroider_wcog = sutra_centroider_wcog<Tin, Tcomp>;

  py::class_<centroider_wcog, sutra_centroider<Tin, Tcomp>>(mod, name)
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "npix", [](centroider_wcog &sc) { return sc.npix; },
          "TODO: docstring")

      .def_property_readonly(
          "d_weights", [](centroider_wcog &sc) { return sc.d_weights; },
          "Weights applied")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("init_weights", &centroider_wcog::init_weights,
           "Initializes WCOG computation")

      .def("set_npix", &centroider_wcog::set_npix, R"pbdoc(
            Set the number of pixels per subap.
            Parameters
            ------------
            npix: (int): number of pixels per subap
            )pbdoc",
           py::arg("npix"))

      .def("load_weights", wy::colCast(&centroider_wcog::load_weights),
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
void declare_centroider_wcog(py::module &mod) {
  centroider_wcog_impl<float, float>(mod, "CentroiderWCOG_FF");
  centroider_wcog_impl<uint16_t, float>(mod, "CentroiderWCOG_UF");
}
