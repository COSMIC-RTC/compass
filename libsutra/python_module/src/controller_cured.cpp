#include <wyrm>

#include <sutra_controller_cured.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_controller_cured(py::module &mod) {
  py::class_<sutra_controller_cured, sutra_controller>(mod, "ControllerCURED")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("gain",
                             [](sutra_controller_cured &sc) { return sc.gain; },
                             "Controller gain")

      .def_property_readonly(
          "ndivs", [](sutra_controller_cured &sc) { return sc.ndivs; },
          "Number of subdivision levels")

      .def_property_readonly(
          "tt_flag", [](sutra_controller_cured &sc) { return sc.tt_flag; },
          "Flag to separate TT")

      .def_property_readonly(
          "h_centroids",
          [](sutra_controller_cured &sc) { return sc.h_centroids; },
          "Centroids")

      .def_property_readonly(
          "h_err", [](sutra_controller_cured &sc) { return sc.h_err; },
          "Increment error")

      .def_property_readonly(
          "d_err", [](sutra_controller_cured &sc) { return sc.d_err; },
          "Increment error")

      .def_property_readonly(
          "d_cenbuff", [](sutra_controller_cured &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly(
          "d_imat", [](sutra_controller_cured &sc) { return sc.d_imat; },
          "Interaction matrix")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("init_cured", wy::colCast(&sutra_controller_cured::init_cured),
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
           py::arg("tt"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_gain", &sutra_controller_cured::set_gain, R"pbdoc(
      Set the gain
      Parameters
      ------------
      gain: (float): gain
    )pbdoc",
           py::arg("gain"));
};
