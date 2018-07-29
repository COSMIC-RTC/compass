#include <wyrm>

#include <sutra_controller_mv.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_shesha_controller_mv(py::module &mod) {
  py::class_<sutra_controller_mv, sutra_controller>(mod, "ControllerMV")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("gain",
                             [](sutra_controller_mv &sc) { return sc.gain; },
                             "Controller gain")

      .def_property_readonly("d_imat",
                             [](sutra_controller_mv &sc) { return sc.d_imat; },
                             "Interaction matrix")

      .def_property_readonly("d_cmat",
                             [](sutra_controller_mv &sc) { return sc.d_cmat; },
                             "Control matrix")

      .def_property_readonly("d_gain",
                             [](sutra_controller_mv &sc) { return sc.d_gain; },
                             "vector of modal gains")

      .def_property_readonly(
          "d_covmat", [](sutra_controller_mv &sc) { return sc.d_covmat; },
          "TODO: docstring")

      .def_property_readonly(
          "d_KLbasis", [](sutra_controller_mv &sc) { return sc.d_KLbasis; },
          "KL basis")

      .def_property_readonly(
          "d_noisemat", [](sutra_controller_mv &sc) { return sc.d_noisemat; },
          "Noise on WFS measurements matrix")

      .def_property_readonly("d_Cmm",
                             [](sutra_controller_mv &sc) { return sc.d_Cmm; },
                             "Slope covariance matrix")

      .def_property_readonly("d_Cphim",
                             [](sutra_controller_mv &sc) { return sc.d_Cphim; },
                             "Actuators-Slopes covariance marix")

      .def_property_readonly(
          "h_Cmmeigenvals",
          [](sutra_controller_mv &sc) { return sc.h_Cmmeigenvals; },
          "Cmm eigenvalues")

      .def_property_readonly(
          "h_eigenvals", [](sutra_controller_mv &sc) { return sc.h_eigenvals; },
          "Eigen values")

      .def_property_readonly(
          "d_cenbuff", [](sutra_controller_mv &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly("d_com1",
                             [](sutra_controller_mv &sc) { return sc.d_com1; },
                             "Commands at iteration k-1 (for POLC)")

      .def_property_readonly("d_com2",
                             [](sutra_controller_mv &sc) { return sc.d_com2; },
                             "Commands at iteration k-2 (for POLC)")

      .def_property_readonly(
          "d_olmeas", [](sutra_controller_mv &sc) { return sc.d_olmeas; },
          "Reconstructed open loop measurement")

      .def_property_readonly("d_err",
                             [](sutra_controller_mv &sc) { return sc.d_err; },
                             "Increment error")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("build_cmat",
           (int (sutra_controller_mv::*)(float)) &
               sutra_controller_mv::build_cmat,
           R"pbdoc(
        Computes the command matrix

        Parameters
        ------------
        cond: (float): Conditioning number for inversion
    )pbdoc",
           py::arg("cond"))

      .def("load_noisemat", wy::colCast(&sutra_controller_mv::load_noisemat),
           R"pbdoc(
        Load the noise covariance matrix

        Parameters
        ------------
        noisemat: (float): noise covariance marix
    )pbdoc",
           py::arg("noisemat"))

      .def("filter_cmat", &sutra_controller_mv::filter_cmat, R"pbdoc(
        Filter command matrix from TT

        Parameters
        ------------
        cond: (float): TODO: docstring
    )pbdoc",
           py::arg("cond"))

      .def("filter_cphim", wy::colCast(&sutra_controller_mv::filter_cphim),
           R"pbdoc(
        Filter Cphim from piston and apply coupling

        Parameters
        ------------
        F: (np.array[ndim=2,dtype=np.float32]): Piston filter matrix
        Nact: (np.array[ndim=2,dtype=np.float32]): Coupling matrix

    )pbdoc",
           py::arg("F"), py::arg("Nact"))

      .def("compute_Cmm", wy::colCast(&sutra_controller_mv::compute_Cmm),
           R"pbdoc(
        Compute the Cmm matrix

        Parameters
        ------------
        atmos : (sutra_atmos): sutra_atmos object
        sensors: (sutra_sensors): sutra_sensors object
        LO: (np.array[ndim=1,dtype=np.float64]): outer scale of each layer
        Cn2: (np.array[ndim=1,dtype=np.float64]): Cn2 profile
        alphaX: (np.array[ndim=1,dtype=np.float64]): X position of each WFS
        alphaY: (np.array[ndim=1,dtype=np.float64]): Y position of each WFS
        diamTel: (double): Telescope diameter
        cobs: (double): Central obstruction ratio

    )pbdoc",
           py::arg("atmos"), py::arg("sensors"), py::arg("LO"), py::arg("Cn2"),
           py::arg("alphaX"), py::arg("alphaY"), py::arg("diamTel"),
           py::arg("cobs"))

      .def("compute_Cphim", wy::colCast(&sutra_controller_mv::compute_Cphim),
           R"pbdoc(
        Compute the Cphim matrix

        Parameters
        ------------
        atmos : (sutra_atmos): sutra_atmos object
        sensors: (sutra_sensors): sutra_sensors object
        dms: (sutra_dms): sutra_dms object
        LO: (np.array[ndim=1,dtype=np.float64]): outer scale of each layer
        Cn2: (np.array[ndim=1,dtype=np.float64]): Cn2 profile
        alphaX: (np.array[ndim=1,dtype=np.float64]): X position of each WFS
        alphaY: (np.array[ndim=1,dtype=np.float64]): Y position of each WFS
        X: (np.array[ndim=1,dtype=np.float64]): X position of each subaperture
        Y: (np.array[ndim=1,dtype=np.float64]): Y position of each subaperture
        xactu: (np.array[ndim=1,dtype=np.float64]): X position of each actuators
        yactu: (np.array[ndim=1,dtype=np.float64]): Y position of each actuators
        diamTel: (double): Telescope diameter
        k2: (np.array[ndim=1,dtype=np.float64]): scales
        NlayerDm: (np.array[ndim=1,dtype=np.int64]): Number of layers handled by each DM
        indLayerDm: (np.array[ndim=1,dtype=np.int64]): Indices of layers handled by each DM
        FoV: (double): Field of view
        pitch: (np.array[ndim=1,dtype=np.int64]): Pitch of each DM
        alt_dm: (np.array[ndim=1,dtype=np.int64]): Altitude of each DM
    )pbdoc",
           py::arg("atmos"), py::arg("sensors"), py::arg("dms"), py::arg("LO"),
           py::arg("Cn2"), py::arg("alphaX"), py::arg("alphaY"), py::arg("X"),
           py::arg("Y"), py::arg("xactu"), py::arg("yactu"), py::arg("diamTel"),
           py::arg("k2"), py::arg("NlayerDm"), py::arg("indLayerDm"),
           py::arg("FoV"), py::arg("pitch"), py::arg("alt_dm"))
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def("set_gain", &sutra_controller_mv::set_gain, R"pbdoc(
      Set the controller gain

      Parameters
      ------------
      gain: (float): gain to set
    )pbdoc",
           py::arg("gain"))

      .def("set_delay", &sutra_controller_mv::set_delay, R"pbdoc(
      Set the loop delay

      Parameters
      ------------
      delay: (float): delay to set
    )pbdoc",
           py::arg("delay"))

      .def("set_mgain", wy::colCast(&sutra_controller_mv::set_mgain), R"pbdoc(
      Set the controller modal gains

      Parameters
      ------------
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_cmat", wy::colCast(&sutra_controller_mv::set_cmat), R"pbdoc(
      Set the command matrix

      Parameters
      ------------
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      .def("set_imat", wy::colCast(&sutra_controller_mv::set_imat), R"pbdoc(
      Set the interaction matrix

      Parameters
      ------------
      imat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("imat"))

      ;
};
