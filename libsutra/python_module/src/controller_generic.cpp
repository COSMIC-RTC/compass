#include <wyrm>

#include <sutra_controller_generic.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;
using controller_generic = sutra_controller_generic<float>;
using controller_genericH = sutra_controller_generic<half>;

void declare_controller_generic(py::module &mod) {
  py::class_<controller_generic, sutra_controller<float>>(mod,
                                                          "ControllerGENERIC")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "d_matE", [](controller_generic &sc) { return sc.d_matE; },
          "E matrix (see compass.io details on generic controller)")

      .def_property_readonly(
          "d_decayFactor",
          [](controller_generic &sc) { return sc.d_decayFactor; },
          "decayFactor vector (see compass.io details on generic controller)")

      .def_property_readonly("d_cmat",
                             [](controller_generic &sc) { return sc.d_cmat; },
                             "Control matrix")

      .def_property_readonly("d_imat",
                             [](controller_generic &sc) { return sc.d_imat; },
                             "Control matrix")

      .def_property_readonly("d_gain",
                             [](controller_generic &sc) { return sc.d_gain; },
                             "vector of modal gains")

      .def_property_readonly("gain",
                             [](controller_generic &sc) { return sc.gain; },
                             "Integrator loop gain")

      .def_property_readonly("polc",
                             [](controller_generic &sc) { return sc.polc; },
                             "POLC flag")

      .def_property_readonly(
          "d_compbuff", [](controller_generic &sc) { return sc.d_compbuff; },
          "Computation buffer buffer")

      .def_property_readonly(
          "command_law", [](controller_generic &sc) { return sc.command_law; },
          "Command law currently used")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_decayFactor", wy::colCast(&controller_generic::set_decayFactor),
           R"pbdoc(
      Set the decay factor vector

      Parameters
      ------------
      decayFactor: (np.array[ndim1,dtype=np.float32]): decay factor
    )pbdoc",
           py::arg("decayFactor"))

      .def("set_polc", wy::colCast(&controller_generic::set_polc),
           R"pbdoc(
      Set the polc flag

      Parameters
      ------------
      polc: (bool): polc flag
    )pbdoc",
           py::arg("polc"))

      .def("set_matE", wy::colCast(&controller_generic::set_matE),
           R"pbdoc(
      Set the E matrix

      Parameters
      ------------
      E: (np.array[ndim=2,dtype=np.float32]): E matrix to set
    )pbdoc",
           py::arg("E"))

      .def("set_commandlaw", wy::colCast(&controller_generic::set_commandlaw),
           R"pbdoc(
      Set the command law to use

      Parameters
      ------------
      commandlaw: (str): command law "integrator" or "2matrices"
    )pbdoc",
           py::arg("commandlaw"))

      .def("set_mgain", wy::colCast(&controller_generic::set_mgain),
           R"pbdoc(
      Set the controller modal gains

      Parameters
      ------------
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_gain", &controller_generic::set_gain,
           R"pbdoc(
      Set the controller loop gain

      Parameters
      ------------
      gain: (float): loop gain to set
    )pbdoc",
           py::arg("gain"))

      .def("set_imat", wy::colCast(&controller_generic::set_imat), R"pbdoc(
      Set the interaction matrix

      Parameters
      ------------
      imat: (np.array[ndim=2,dtype=np.float32]): interaction matrix to set
    )pbdoc",
           py::arg("imat"))

      .def("set_cmat", wy::colCast(&controller_generic::set_cmat),
           R"pbdoc(
      Set the command matrix

      Parameters
      ------------
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      ;

  py::class_<controller_genericH, sutra_controller<half>>(mod,
                                                          "ControllerGENERIC_H")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "d_matE", [](controller_genericH &sc) { return sc.d_matE; },
          "E matrix (see compass.io details on generic controller)")

      .def_property_readonly(
          "d_decayFactor",
          [](controller_genericH &sc) { return sc.d_decayFactor; },
          "decayFactor vector (see compass.io details on generic controller)")

      .def_property_readonly("d_cmat",
                             [](controller_genericH &sc) { return sc.d_cmat; },
                             "Control matrix")

      .def_property_readonly("d_imat",
                             [](controller_genericH &sc) { return sc.d_imat; },
                             "Control matrix")

      .def_property_readonly("d_gain",
                             [](controller_genericH &sc) { return sc.d_gain; },
                             "vector of modal gains")

      .def_property_readonly("gain",
                             [](controller_genericH &sc) { return sc.gain; },
                             "Integrator loop gain")

      .def_property_readonly("polc",
                             [](controller_genericH &sc) { return sc.polc; },
                             "POLC flag")

      .def_property_readonly(
          "d_compbuff", [](controller_genericH &sc) { return sc.d_compbuff; },
          "Computation buffer buffer")

      .def_property_readonly(
          "command_law", [](controller_genericH &sc) { return sc.command_law; },
          "Command law currently used")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_decayFactor",
           wy::colCast(&controller_genericH::set_decayFactor),
           R"pbdoc(
      Set the decay factor vector

      Parameters
      ------------
      decayFactor: (np.array[ndim1,dtype=np.float32]): decay factor
    )pbdoc",
           py::arg("decayFactor"))

      .def("set_polc", wy::colCast(&controller_genericH::set_polc),
           R"pbdoc(
      Set the polc flag

      Parameters
      ------------
      polc: (bool): polc flag
    )pbdoc",
           py::arg("polc"))

      .def("set_matE", wy::colCast(&controller_genericH::set_matE),
           R"pbdoc(
      Set the E matrix

      Parameters
      ------------
      E: (np.array[ndim=2,dtype=np.float32]): E matrix to set
    )pbdoc",
           py::arg("E"))

      .def("set_commandlaw", wy::colCast(&controller_genericH::set_commandlaw),
           R"pbdoc(
      Set the command law to use

      Parameters
      ------------
      commandlaw: (str): command law "integrator" or "2matrices"
    )pbdoc",
           py::arg("commandlaw"))

      .def("set_mgain", wy::colCast(&controller_genericH::set_mgain),
           R"pbdoc(
      Set the controller modal gains

      Parameters
      ------------
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_gain",
           [](controller_genericH &sc, float gain) {
             sc.set_gain(__float2half(gain));
           },
           R"pbdoc(
      Set the controller loop gain

      Parameters
      ------------
      gain: (float): loop gain to set
    )pbdoc",
           py::arg("gain"))

      .def("set_imat", wy::colCast(&controller_genericH::set_imat), R"pbdoc(
      Set the interaction matrix

      Parameters
      ------------
      imat: (np.array[ndim=2,dtype=np.float32]): interaction matrix to set
    )pbdoc",
           py::arg("imat"))

      .def("set_cmat", wy::colCast(&controller_genericH::set_cmat),
           R"pbdoc(
      Set the command matrix

      Parameters
      ------------
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      ;
};
