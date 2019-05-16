#include <wyrm>

#include <sutra_controller.h>
#include "declare_name.hpp"

namespace py = pybind11;

template <typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, float>::type
get_gain(sutra_controller<Tcomp, Tout> &sc) {
  return float(sc.gain);
}

template <typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, float>::type get_gain(
    sutra_controller<Tcomp, Tout> &sc) {
  return __half2float(sc.gain);
}

template <typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, float>::type
get_delay(sutra_controller<Tcomp, Tout> &sc) {
  return float(sc.delay);
}

template <typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, float>::type
get_delay(sutra_controller<Tcomp, Tout> &sc) {
  return __half2float(sc.delay);
}

template <typename Tcomp, typename Tout>
void controller_impl(py::module &mod, const char *name) {
  using controller = sutra_controller<Tcomp, Tout>;

  py::class_<controller>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "context", [](controller &sc) { return sc.current_context; },
          "GPU context")

      .def_property_readonly(
          "device", [](controller &sc) { return sc.device; },
          "GPU device index")

      .def_property_readonly(
          "type", [](controller &sc) { return sc.get_type(); },
          "Controller type")

      .def_property_readonly(
          "nactu", [](controller &sc) { return sc.nactu(); },
          "Number of actuators to control")

      .def_property_readonly(
          "nslope", [](controller &sc) { return sc.nslope(); },
          "Number of slopes")

      .def_property_readonly(
          "gain", [](controller &sc) { return get_gain(sc); },
          "Controller gain")
      .def_property_readonly(
          "open_loop", [](controller &sc) { return sc.open_loop; },
          "Open loop flag")

      .def_property_readonly(
          "delay", [](controller &sc) { return get_delay(sc); }, "Loop delay")

      .def_property_readonly(
          "d_dmseen", [](controller &sc) { return sc.d_dmseen; },
          "Vector of sutra_dm commanded")

      .def_property_readonly(
          "d_centroids", [](controller &sc) { return sc.d_centroids; },
          "Slopes vector")

      .def_property_readonly(
          "d_com", [](controller &sc) { return sc.d_com; },
          "Current command vector")

      .def_property_readonly(
          "d_com1", [](controller &sc) { return sc.d_com1; },
          "Command vector at iteration k-1")

      .def_property_readonly(
          "d_com2", [](controller &sc) { return sc.d_com2; },
          "Command vector at iteration k-2")
      .def_property_readonly(
          "comRange",
          [](controller &sc) { return std::make_tuple(sc.Vmin, sc.Vmax); },
          "Tuple (Vmin, Vmax) used for command clipping")
      .def_property_readonly(
          "valMax", [](controller &sc) { return sc.valMax; },
          "Maximum value for d_voltage (ADU). Only used if "
          "output is expected in uint16")

      .def_property_readonly(
          "d_perturb_map",
          [](controller &sc)
              -> map<string, tuple<carma_obj<Tcomp> *, int, bool>> & {
            return sc.d_perturb_map;
          },
          "Perturbation voltage buffers")

      .def_property_readonly(
          "d_voltage", [](controller &sc) { return sc.d_voltage; },
          "Total voltage to apply on the DMs")

      .def_property_readonly(
          "d_comClipped", [](controller &sc) { return sc.d_comClipped; },
          "Delayed commands")

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

      .def("comp_voltage", wy::colCast(&controller::comp_voltage),
           "Computes the final voltage to send to the DM")

      .def("clip_commands", wy::colCast(&controller::clip_commands),
           "Clip the commands between Vmin and Vmax (values set in the "
           "controller)")

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

      .def("set_perturb_voltage", wy::colCast(&controller::set_perturb_voltage),
           R"pbdoc(
          Set an existing perturbation voltage buffer

          Parameters
          ------------
          name: (str): name of the buffer
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

      .def("set_delay", wy::colCast(&controller::set_delay),
           R"pbdoc(
          Set the delay

          Parameters
          ------------
          delay: (float): loop delay in frames
     )pbdoc",
           py::arg("delay"))

      .def("set_gain", wy::colCast(&controller::set_gain),
           R"pbdoc(
          Set the gain

          Parameters
          ------------
          gain: (float): loop gain
     )pbdoc",
           py::arg("gain"))

      .def(
          "set_comRange",
          [](controller &sc, float Vmin, float Vmax) {
            sc.set_Vmax(Vmax);
            sc.set_Vmin(Vmin);
          },
          R"pbdoc(
          Set the Vmin and Vmax value for command clipping

          Parameters
          ------------
          Vmin: (float): Vmin value for clipping
          Vmax: (float): Vmax value for clipping
     )pbdoc",
          py::arg("Vmin"), py::arg("Vmax"))

      .def("set_Vmax", wy::colCast(&controller::set_Vmax),
           R"pbdoc(
          Set the Vmax value for command clipping

          Parameters
          ------------
          Vmax: (float): Vmax value for clipping
     )pbdoc",
           py::arg("Vmax"))

      .def("set_valMax", wy::colCast(&controller::set_valMax),
           R"pbdoc(
          Set the valMax value for command conversion

          Parameters
          ------------
          valMax: (float): valMax value for conversion
     )pbdoc",
           py::arg("valMax"))

      .def("set_com", wy::colCast(&controller::set_com), R"pbdoc(
          Set the command vector of the controller
          Parameters
          ------------
          com: (np.array[ndim=3, dtype=np.float32]) : command vector
          nElem: (int): Number of elements in com
          )pbdoc",
           py::arg("com"), py::arg("nElem"));
};

void declare_controller(py::module &mod) {
  controller_impl<float, float>(mod, "ControllerFF");
  controller_impl<float, uint16_t>(mod, "ControllerFU");
#ifdef CAN_DO_HALF
  controller_impl<half, float>(mod, "ControllerHF");
  controller_impl<half, uint16_t>(mod, "ControllerHU");
#endif
}
