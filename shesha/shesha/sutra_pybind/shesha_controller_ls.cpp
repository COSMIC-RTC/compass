#include <wyrm>

#include <sutra_controller_ls.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_shesha_controller_ls(py::module &mod) {
  py::class_<sutra_controller_ls>(mod, "ControllerLS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("gain",
                             [](sutra_controller_ls &sc) { return sc.gain; },
                             "Controller gain")

      .def_property_readonly("d_imat",
                             [](sutra_controller_ls &sc) { return sc.d_imat; },
                             "Interaction matrix")

      .def_property_readonly("d_cmat",
                             [](sutra_controller_ls &sc) { return sc.d_cmat; },
                             "Control matrix")

      .def_property_readonly("d_gain",
                             [](sutra_controller_ls &sc) { return sc.d_gain; },
                             "vector of modal gains")

      .def_property_readonly(
          "d_eigenvals", [](sutra_controller_ls &sc) { return sc.d_eigenvals; },
          "Eigen values of the imat")

      .def_property_readonly("d_U",
                             [](sutra_controller_ls &sc) { return sc.d_U; },
                             "Eigen modes of the imat")

      .def_property_readonly(
          "d_cenbuff", [](sutra_controller_ls &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly("d_err",
                             [](sutra_controller_ls &sc) { return sc.d_err; },
                             "Current increment on the command")

      .def_property_readonly(
          "is_modopti", [](sutra_controller_ls &sc) { return sc.is_modopti; },
          "Falg for modal optimization")

      .def_property_readonly(
          "nrec", [](sutra_controller_ls &sc) { return sc.nrec; },
          "Number of open loop slopes to take for modal optimization")

      .def_property_readonly("d_M2V",
                             [](sutra_controller_ls &sc) { return sc.nmodes; },
                             "Number of modes for modal optimization")

      .def_property_readonly("gmin",
                             [](sutra_controller_ls &sc) { return sc.gmin; },
                             "Minimal gain for modal optimization")

      .def_property_readonly("gmax",
                             [](sutra_controller_ls &sc) { return sc.gmax; },
                             "Maximal gain for modal optimization")

      .def_property_readonly("ngain",
                             [](sutra_controller_ls &sc) { return sc.ngain; },
                             "Number of gain values to test between gmin and "
                             "gmax for modal optimization")

      .def_property_readonly("Fs",
                             [](sutra_controller_ls &sc) { return sc.Fs; },
                             "Sampling frequency for modal optimization")

      .def_property_readonly("cpt_rec",
                             [](sutra_controller_ls &sc) { return sc.cpt_rec; },
                             "Counter for modal gains refresh")

      .def_property_readonly("d_M2V",
                             [](sutra_controller_ls &sc) { return sc.d_M2V; },
                             "Modes to volt matrix for modal optimization")

      .def_property_readonly("d_S2M",
                             [](sutra_controller_ls &sc) { return sc.d_S2M; },
                             "Slopes to modes matrix for modal optimization")

      .def_property_readonly("d_slpol",
                             [](sutra_controller_ls &sc) { return sc.d_slpol; },
                             "Open loop slopes for modal optimization")

      .def_property_readonly("d_Hcor",
                             [](sutra_controller_ls &sc) { return sc.d_Hcor; },
                             "Transfer function for modal optimization")

      .def_property_readonly(
          "d_compbuff", [](sutra_controller_ls &sc) { return sc.d_compbuff; },
          "Buffer for POLC computation")

      .def_property_readonly(
          "d_compbuff2", [](sutra_controller_ls &sc) { return sc.d_compbuff2; },
          "Buffer for POLC computation")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("svdec_imat", wy::colCast(&sutra_controller_ls::svdec_imat),
           "Performs interaction matrix SVD")

      .def("build_cmat",
           wy::colCast((int (sutra_controller_ls::*)(int)) &
                       sutra_controller_ls::build_cmat),
           R"pbdoc(
        Computes the command matrix after imat SVD

        Parameters
        ------------
        nfilt: (int): number of modes to filter
    )pbdoc",
           py::arg("nfilt"))

      .def("init_modalOpti", wy::colCast(&sutra_controller_ls::init_modalOpti),
           R"pbdoc(
      Initialize modal optimization control

      Parameters
      ------------
      nmodes: (int): number of modes to control
      nrec: (int): number of open loop slopes to consider
      M2V: (np.array[ndim=2,dtype=np.float32]): Modes to Volt matrix
      gmin: (float): Minimal gain
      gmax: (float): Maximal gain
      ngain: (int): Number of gain values to test between gmin and gmax
      Fs: (float): Sampling frequency [Hz]
    )pbdoc",
           py::arg("nmodes"), py::arg("nrec"), py::arg("M2V"), py::arg("gmin"),
           py::arg("gmax"), py::arg("ngain"), py::arg("Fs"))

      .def("loadOpenLoopSlp",
           wy::colCast(&sutra_controller_ls::loadOpenLoopSlp), R"pbdoc(
      Load recorded open loop slopes for modal optimization initialization

      Parameters
      ------------
      slopes: (np.array[ndim=2,dtype=np.float32]): Open loop slopes
    )pbdoc",
           py::arg("slopes"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_centroids_ref",
           wy::colCast(&sutra_controller_ls::set_centroids_ref), R"pbdoc(
      Set the references slopes

      Parameters
      ------------
      refslopes: (np.array[ndim1,dtype=np.float32]): reference slopes to set
    )pbdoc",
           py::arg("refslopes"))

      .def("set_gain", wy::colCast(&sutra_controller_ls::set_gain), R"pbdoc(
      Set the controller gain

      Parameters
      ------------
      gain: (float): gain to set
    )pbdoc",
           py::arg("gain"))

      .def("set_delay", wy::colCast(&sutra_controller_ls::set_delay), R"pbdoc(
      Set the loop delay

      Parameters
      ------------
      delay: (float): delay to set
    )pbdoc",
           py::arg("delay"))

      .def("set_mgain", wy::colCast(&sutra_controller_ls::set_mgain), R"pbdoc(
      Set the controller modal gains

      Parameters
      ------------
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_cmat", wy::colCast(&sutra_controller_ls::set_cmat), R"pbdoc(
      Set the command matrix

      Parameters
      ------------
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      ;
};
