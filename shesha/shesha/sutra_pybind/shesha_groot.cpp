#include <wyrm>

#include <sutra_groot.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

std::unique_ptr<sutra_groot> groot_init(carma_context &context, int device,
                                        int nactus, int nlayers, float gsangle,
                                        float *vdt, float *Htheta, float *L0,
                                        float *winddir, float *scale,
                                        float *pzt2tt, float *TTPfilter,
                                        float *Nact, float *xpos, float *ypos,
                                        float fc) {
  return std::unique_ptr<sutra_groot>(new sutra_groot(
      &context, device, nactus, nlayers, gsangle, vdt, Htheta, L0, winddir,
      scale, pzt2tt, TTPfilter, Nact, xpos, ypos, fc));
};

std::unique_ptr<sutra_groot> groot_init_alias(carma_context &context,
                                              int device, int nssp,
                                              float *weights, float scale,
                                              float *xpos, float *ypos,
                                              float fc, float d, int npts) {
  return std::unique_ptr<sutra_groot>(new sutra_groot(
      &context, device, nssp, weights, scale, xpos, ypos, fc, d, npts));
};

void declare_shesha_groot(py::module &mod) {
  py::class_<sutra_groot>(mod, "Groot")
      .def(py::init(wy::colCast(groot_init)), R"pbdoc(
          Initializes Groot to compute aniso and bandwidth model

          Parameters
          ------------
          context: (carma_context): context
          device: (int): context active device
          nssp : (str) :
          nlayers:
          gsangle:
          vdt:
          Htheta:
          L0:
          winddir:
          scale:
          pzt2tt:
          TTPfilter:
          Nact:
          xpos:
          ypos:
          fc:
           )pbdoc",
           py::arg("context"), py::arg("device"), py::arg("nssp"),
           py::arg("nlayers"), py::arg("gsangle"), py::arg("vdt"),
           py::arg("Htheta"), py::arg("L0"), py::arg("winddir"),
           py::arg("scale"), py::arg("pzt2tt"), py::arg("TTPfilter"),
           py::arg("Nact"), py::arg("xpos"), py::arg("ypos"), py::arg("fc"))

      .def(py::init(wy::colCast(groot_init_alias)), R"pbdoc(
          Initializes Groot to compute aliasing model

          Parameters
          ------------
          context: (carma_context): context
          device: (int): context active device
          nssp : (str) :
          weights:
          scale:
          xpos:
          ypos:
          fc:
          d:
          npts:
           )pbdoc",
           py::arg("context"), py::arg("device"), py::arg("nssp"),
           py::arg("weights"), py::arg("scale"), py::arg("xpos"),
           py::arg("ypos"), py::arg("fc"), py::arg("d"), py::arg("npts"))
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "device", [](sutra_groot &sg) { return sg.device; }, "GPU device")

      .def_property_readonly("nactus",
                             [](sutra_groot &sg) { return sg.nactus; },
                             "Number of actuators")

      .def_property_readonly("nssp", [](sutra_groot &sg) { return sg.nssp; },
                             "number of subap")

      .def_property_readonly("nlayers",
                             [](sutra_groot &sg) { return sg.nlayers; },
                             "number of turbulent layers")

      .def_property_readonly("npts", [](sutra_groot &sg) { return sg.npts; },
                             "number of samples for aliasig computation")

      .def_property_readonly("gsangle",
                             [](sutra_groot &sg) { return sg.gsangle; },
                             "Guide star angle [rad]")

      .def_property_readonly("fc", [](sutra_groot &sg) { return sg.fc; },
                             "DM cut-off frequency [m]")

      .def_property_readonly("d", [](sutra_groot &sg) { return sg.d; },
                             "DM pitch")

      .def_property_readonly(
          "d_Cerr", [](sutra_groot &sg) { return sg.d_Cerr; },
          "Model of aniso and bandwidth covariance error matrix")

      .def_property_readonly("d_CaXX",
                             [](sutra_groot &sg) { return sg.d_CaXX; },
                             "XX component of the aliasing model")

      .def_property_readonly("d_CaYY",
                             [](sutra_groot &sg) { return sg.d_CaYY; },
                             "YY component of the aliasing model")

      .def_property_readonly("d_TT", [](sutra_groot &sg) { return sg.d_TT; },
                             "tip-tilt IF matrix")

      .def_property_readonly("scale", [](sutra_groot &sg) { return sg.scale; },
                             "Scale factor")

      .def_property_readonly("d_TTPfilter",
                             [](sutra_groot &sg) { return sg.d_TTPfilter; },
                             "Tip-tilt and piston filter matrix (= Btt.dot(P))")

      .def_property_readonly("d_pzt2tt",
                             [](sutra_groot &sg) { return sg.d_pzt2tt; },
                             "pzt to TT matrix")

      .def_property_readonly("d_Nact",
                             [](sutra_groot &sg) { return sg.d_Nact; },
                             "Coupling matrix")

      .def_property_readonly("d_xpos",
                             [](sutra_groot &sg) { return sg.d_xpos; },
                             "X-positions of DM actuators or ssp [m]")

      .def_property_readonly("d_ypos",
                             [](sutra_groot &sg) { return sg.d_ypos; },
                             "Y-positions of DM actuators or ssp [m]")

      .def_property_readonly("d_tab_int_x",
                             [](sutra_groot &sg) { return sg.d_tab_int_x; },
                             "Tabulated integral")

      .def_property_readonly("d_tab_int_y",
                             [](sutra_groot &sg) { return sg.d_tab_int_y; },
                             "Tabulated integral")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("compute_Cerr", wy::colCast(&sutra_groot::compute_Cerr),
           "Computes the aniso and bandwidth error covariance matrix")

      .def("compute_Calias", wy::colCast(&sutra_groot::compute_Calias),
           "Computes the aliasing error covariance matrix");
};
