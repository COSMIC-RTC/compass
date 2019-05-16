#include <wyrm>

#include <sutra_centroider_corr.h>
#include "declare_name.hpp"

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_corr_impl(py::module &mod, const char *name) {
  using centroider_corr = sutra_centroider_corr<Tin, Tcomp>;

  py::class_<centroider_corr, sutra_centroider<Tin, Tcomp>>(mod, name)
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "npix", [](centroider_corr &sc) { return sc.npix; },
          "TODO: docstring")

      .def_property_readonly(
          "interp_sizex", [](centroider_corr &sc) { return sc.interp_sizex; },
          "TODO: docstring")

      .def_property_readonly(
          "interp_sizey", [](centroider_corr &sc) { return sc.interp_sizey; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrfnct", [](centroider_corr &sc) { return sc.d_corrfnct; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrspot", [](centroider_corr &sc) { return sc.d_corrspot; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrnorm", [](centroider_corr &sc) { return sc.d_corrnorm; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corrmax", [](centroider_corr &sc) { return sc.d_corrmax; },
          "TODO: docstring")

      .def_property_readonly(
          "d_corr", [](centroider_corr &sc) { return sc.d_corr; },
          "TODO: docstring")

      .def_property_readonly(
          "d_interpmat", [](centroider_corr &sc) { return sc.d_interpmat; },
          "TODO: docstring")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("init_corr", wy::colCast(&centroider_corr::init_corr), R"pbdoc(
            Initializes corr computation

            Parameters
            ------------
            isizex: (int): TODO: docstring
            isizey: (int):
            interpmat: (np.array[ndim=,dtype=np.float32]):
            )pbdoc",
           py::arg("isizex"), py::arg("isizey"), py::arg("interpmat"))

      .def("load_corr", wy::colCast(&centroider_corr::load_corr), R"pbdoc(
            Load arrays for correlation computation

            Parameters
            ------------
            corr: (np.array[ndim=,dtype=np.float32): TODO: docstring
            corr_norm: (np.array[ndim=,dtype=np.float32):
            ndim: (int):
            )pbdoc",
           py::arg("corr"), py::arg("corr_norm"), py::arg("ndim"))

      .def("set_npix", wy::colCast(&centroider_corr::set_npix),
           R"pbdoc(
               Set the number of pixels per subap.
            Parameters
            ------------
            npix: (int): number of pixels per subap
            )pbdoc",
           py::arg("npix"))

      ;
};

void declare_centroider_corr(py::module &mod) {
  centroider_corr_impl<float, float>(mod, "CentroiderCORR_FF");
  centroider_corr_impl<uint16_t, float>(mod, "CentroiderCORR_UF");
}
