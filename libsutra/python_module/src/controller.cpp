#include <wyrm>

#include <sutra_controller.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;
using controller = sutra_controller<float>;
using controllerH = sutra_controller<half>;

void declare_controller(py::module &mod) {
  py::class_<controller>(mod, "Controller")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("context",
                             [](controller &sc) { return sc.current_context; },
                             "GPU context")

      .def_property_readonly("device", [](controller &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](controller &sc) { return sc.get_type(); },
                             "Controller type")

      .def_property_readonly("nactu", [](controller &sc) { return sc.nactu(); },
                             "Number of actuators to control")

      .def_property_readonly("nslope",
                             [](controller &sc) { return sc.nslope(); },
                             "Number of slopes")

      .def_property_readonly("open_loop",
                             [](controller &sc) { return sc.open_loop; },
                             "Open loop flag")

      .def_property_readonly("delay", [](controller &sc) { return sc.delay; },
                             "Loop delay")

      .def_property_readonly("d_dmseen",
                             [](controller &sc) { return sc.d_dmseen; },
                             "Vector of sutra_dm commanded")

      .def_property_readonly("d_centroids",
                             [](controller &sc) { return sc.d_centroids; },
                             "Slopes vector")

      .def_property_readonly("d_com", [](controller &sc) { return sc.d_com; },
                             "Current command vector")

      .def_property_readonly("d_com1", [](controller &sc) { return sc.d_com1; },
                             "Command vector at iteration k-1")

      .def_property_readonly("d_com2", [](controller &sc) { return sc.d_com2; },
                             "Command vector at iteration k-2")

      .def_property_readonly(
          "d_perturb_map",
          [](controller &sc)
              -> map<string, tuple<carma_obj<float> *, int, bool>> & {
            return sc.d_perturb_map;
          },
          "Perturbation voltage buffers")

      .def_property_readonly("d_voltage",
                             [](controller &sc) { return sc.d_voltage; },
                             "Total voltage to apply on the DMs")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("add_perturb", &controller::add_perturb,
           "Add the perturbation voltage to the command")

      .def("command_delay", &controller::command_delay, "Delay the command")

      .def("enable_perturb_voltage",
           wy::colCast(&controller::enable_perturb_voltage),
           R"pbdoc(
      Enable a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to enable
    )pbdoc",
           py::arg("name"))

      .def("disable_perturb_voltage",
           wy::colCast(&controller::disable_perturb_voltage),
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
      .def("add_perturb_voltage", wy::colCast(&controller::add_perturb_voltage),
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
           wy::colCast(&controller::remove_perturb_voltage),
           R"pbdoc(
      Remove a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to remove
    )pbdoc",
           py::arg("name"))

      .def("reset_perturb_voltage", &controller::reset_perturb_voltage,
           R"pbdoc(
      Remove all perturbation voltage buffers
    )pbdoc")

      .def("set_openloop", wy::colCast(&controller::set_openloop),
           R"pbdoc(
      Open (1) or close (0) the loop

      Parameters
      ------------
      status: (int): open loop status
      rst: (bool): reset integrator if True
    )pbdoc",
           py::arg("status"), py::arg("rst") = true)

      .def("set_com", wy::colCast(&controller::set_com), R"pbdoc(
        Set the command vector of the controller
        Parameters
        ------------
        com: (np.array[ndim=3, dtype=np.float32]) : command vector
        nElem: (int): Number of elements in com
      )pbdoc",
           py::arg("com"), py::arg("nElem"))

      ;

  py::class_<controllerH>(mod, "ControllerH")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("context",
                             [](controllerH &sc) { return sc.current_context; },
                             "GPU context")

      .def_property_readonly("device",
                             [](controllerH &sc) { return sc.device; },
                             "GPU device index")

      .def_property_readonly("type",
                             [](controllerH &sc) { return sc.get_type(); },
                             "Controller type")

      .def_property_readonly("nactu",
                             [](controllerH &sc) { return sc.nactu(); },
                             "Number of actuators to control")

      .def_property_readonly("nslope",
                             [](controllerH &sc) { return sc.nslope(); },
                             "Number of slopes")

      .def_property_readonly("open_loop",
                             [](controllerH &sc) { return sc.open_loop; },
                             "Open loop flag")

      .def_property_readonly("delay", [](controllerH &sc) { return sc.delay; },
                             "Loop delay")

      .def_property_readonly("d_dmseen",
                             [](controllerH &sc) { return sc.d_dmseen; },
                             "Vector of sutra_dm commanded")

      .def_property_readonly("d_centroids",
                             [](controllerH &sc) { return sc.d_centroids; },
                             "Slopes vector")

      .def_property_readonly("d_com", [](controllerH &sc) { return sc.d_com; },
                             "Current command vector")

      .def_property_readonly("d_com1",
                             [](controllerH &sc) { return sc.d_com1; },
                             "Command vector at iteration k-1")

      .def_property_readonly("d_com2",
                             [](controllerH &sc) { return sc.d_com2; },
                             "Command vector at iteration k-2")

      .def_property_readonly(
          "d_perturb_map",
          [](controllerH &sc)
              -> map<string, tuple<carma_obj<half> *, int, bool>> & {
            return sc.d_perturb_map;
          },
          "Perturbation voltage buffers")

      .def_property_readonly("d_voltage",
                             [](controllerH &sc) { return sc.d_voltage; },
                             "Total voltage to apply on the DMs")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("add_perturb", &controllerH::add_perturb,
           "Add the perturbation voltage to the command")

      .def("command_delay", &controllerH::command_delay, "Delay the command")

      .def("enable_perturb_voltage",
           wy::colCast(&controllerH::enable_perturb_voltage),
           R"pbdoc(
      Enable a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to enable
    )pbdoc",
           py::arg("name"))

      .def("disable_perturb_voltage",
           wy::colCast(&controllerH::disable_perturb_voltage),
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
      .def("add_perturb_voltage",
           wy::colCast(&controllerH::add_perturb_voltage),
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
           wy::colCast(&controllerH::remove_perturb_voltage),
           R"pbdoc(
      Remove a perturbation voltage buffer

      Parameters
      ------------
      name: (str): name of the buffer to remove
    )pbdoc",
           py::arg("name"))

      .def("reset_perturb_voltage", &controllerH::reset_perturb_voltage,
           R"pbdoc(
      Remove all perturbation voltage buffers
    )pbdoc")

      .def("set_openloop", wy::colCast(&controllerH::set_openloop),
           R"pbdoc(
      Open (1) or close (0) the loop

      Parameters
      ------------
      status: (int): open loop status
      rst: (bool): reset integrator if True
    )pbdoc",
           py::arg("status"), py::arg("rst") = true)

      .def("set_com", wy::colCast(&controllerH::set_com), R"pbdoc(
        Set the command vector of the controller
        Parameters
        ------------
        com: (np.array[ndim=3, dtype=np.float32]) : command vector
        nElem: (int): Number of elements in com
      )pbdoc",
           py::arg("com"), py::arg("nElem"))

      ;
};
