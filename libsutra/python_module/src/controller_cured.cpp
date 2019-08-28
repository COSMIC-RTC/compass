#include <wyrm>

#include <sutra_controller_cured.h>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_cured_impl(py::module &mod, const char *name) {
  using controller_cured = sutra_controller_cured<Tcomp, Tout>;

  py::class_<controller_cured, sutra_controller<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "ndivs", [](controller_cured &sc) { return sc.ndivs; },
          "Number of subdivision levels")

      .def_property_readonly(
          "tt_flag", [](controller_cured &sc) { return sc.tt_flag; },
          "Flag to separate TT")

      .def_property_readonly(
          "h_centroids", [](controller_cured &sc) { return sc.h_centroids; },
          "Centroids")

      .def_property_readonly(
          "h_err", [](controller_cured &sc) { return sc.h_err; },
          "Increment error")

      .def_property_readonly(
          "d_err", [](controller_cured &sc) { return sc.d_err; },
          "Increment error")

      .def_property_readonly(
          "d_cenbuff", [](controller_cured &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly(
          "d_imat", [](controller_cured &sc) { return sc.d_imat; },
          "Interaction matrix")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("init_cured", wy::colCast(&controller_cured::init_cured),
           R"pbdoc(
        Initialize CURED

        Parameters
        ------------
        nxsub: (int): TODO: docstring
        isvalid:
        ndivs: (int):
        tt: (int):
    )pbdoc",
           py::arg("nxsub"), py::arg("isvalid"), py::arg("ndivs"),
           py::arg("tt"));
};

void declare_controller_cured(py::module &mod) {
  controller_cured_impl<float, float>(mod, "ControllerCURED_FF");
  controller_cured_impl<float, uint16_t>(mod, "ControllerCURED_FU");
}
