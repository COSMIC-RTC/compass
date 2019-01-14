#include <wyrm>

#include <sutra_controller.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_controller(py::module &mod) {
  py::class_<sutra_controller>(mod, "Controller")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly("device",
                             [](sutra_controller &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](sutra_controller &sc) { return sc.get_type(); },
                             "Controller type")

      .def_property_readonly("nactu",
                             [](sutra_controller &sc) { return sc.nactu(); },
                             "Number of actuators to control")

      .def_property_readonly("nslope",
                             [](sutra_controller &sc) { return sc.nslope(); },
                             "Number of slopes")

      .def_property_readonly("open_loop",
                             [](sutra_controller &sc) { return sc.open_loop; },
                             "Open loop flag")

      .def_property_readonly(
          "delay", [](sutra_controller &sc) { return sc.delay; }, "Loop delay")

      .def_property_readonly("d_dmseen",
                             [](sutra_controller &sc) { return sc.d_dmseen; },
                             "Vector of sutra_dm commanded")

      .def_property_readonly("d_subsum",
                             [](sutra_controller &sc) { return sc.d_subsum; },
                             "Array to store ssp intensities sum")

      .def_property_readonly(
          "d_centroids", [](sutra_controller &sc) { return sc.d_centroids; },
          "Slopes vector")

      .def_property_readonly(
          "d_centroids_ref",
          [](sutra_controller &sc) { return sc.d_centroids_ref; },
          "Reference slopes vector")

      .def_property_readonly("d_com",
                             [](sutra_controller &sc) { return sc.d_com; },
                             "Current command vector")

      .def_property_readonly("d_com1",
                             [](sutra_controller &sc) { return sc.d_com1; },
                             "Command vector at iteration k-1")

      .def_property_readonly("d_com1",
                             [](sutra_controller &sc) { return sc.d_com1; },
                             "Command vector at iteration k-2")

      .def_property_readonly(
          "d_perturb_map",
          [](sutra_controller &sc)
              -> map<string, tuple<carma_obj<float> *, int, bool>> & {
            return sc.d_perturb_map;
          },
          "Perturbation voltage buffers")

      .def_property_readonly("d_voltage",
                             [](sutra_controller &sc) { return sc.d_voltage; },
                             "Total voltage to apply on the DMs")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("add_perturb", &sutra_controller::add_perturb,
           "Add the perturbation voltage to the command")

      .def("command_delay", &sutra_controller::command_delay,
           "Delay the command")

      .def("enable_perturb_voltage",
           wy::colCast(&sutra_controller::enable_perturb_voltage),
           R"pbdoc(
      Enable a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to enable
    )pbdoc",
           py::arg("name"))

      .def("disable_perturb_voltage",
           wy::colCast(&sutra_controller::disable_perturb_voltage),
           R"pbdoc(
      Disable a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to enable
    )pbdoc",
           py::arg("name"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_centroids_ref",
           wy::colCast(&sutra_controller::set_centroids_ref), R"pbdoc(
      Set the references slopes

      Parameters
      ------------
      refslopes: (np.array[ndim1,dtype=np.float32]): reference slopes to set
    )pbdoc",
           py::arg("refslopes"))

      .def("add_perturb_voltage",
           wy::colCast(&sutra_controller::add_perturb_voltage),
           R"pbdoc(
      Add a new perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the new buffer
      perturb: (np.array[ndim=2,dtype=np.float32]): perturbation voltage to set
      N: (int): Number of perturb voltage vectors in the buffer
    )pbdoc",
           py::arg("name"), py::arg("perturb"), py::arg("N"))

      .def("remove_perturb_voltage",
           wy::colCast(&sutra_controller::remove_perturb_voltage),
           R"pbdoc(
      Remove a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to remove
    )pbdoc",
           py::arg("name"))

      .def("reset_perturb_voltage", &sutra_controller::reset_perturb_voltage,
           R"pbdoc(
      Remove all perturbation voltage buffers
    )pbdoc")

      .def("set_openloop", wy::colCast(&sutra_controller::set_openloop),
           R"pbdoc(
      Open (1) or close (0) the loop

      Parameters
      ------------
      status: (int): open loop status
      rst: (bool): reset integrator if True
    )pbdoc",
           py::arg("status"), py::arg("rst") = true)

      .def("set_com", wy::colCast(&sutra_controller::set_com), R"pbdoc(
        Set the command vector of the controller
        Parameters
        ------------
        com: (np.array[ndim=3, dtype=np.float32]) : command vector
        nElem: (int): Number of elements in com
      )pbdoc",
           py::arg("com"), py::arg("nElem"))

      ;
};
