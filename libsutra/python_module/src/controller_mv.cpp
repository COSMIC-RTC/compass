// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      controller_mv.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_controller_mv
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#include <wyrm>

#include <sutra_controller_mv.h>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_mv_impl(py::module &mod, const char *name) {
  using controller_mv = sutra_controller_mv<Tcomp, Tout>;

  py::class_<controller_mv, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "d_imat", [](controller_mv &sc) { return sc.d_imat; },
          "Interaction matrix")

      .def_property_readonly(
          "d_cmat", [](controller_mv &sc) { return sc.d_cmat; },
          "Control matrix")

      .def_property_readonly(
          "d_gain", [](controller_mv &sc) { return sc.d_gain; },
          "vector of modal gains")

      .def_property_readonly(
          "d_covmat", [](controller_mv &sc) { return sc.d_covmat; },
          "TODO: docstring")

      .def_property_readonly(
          "d_KLbasis", [](controller_mv &sc) { return sc.d_KLbasis; },
          "KL basis")

      .def_property_readonly(
          "d_noisemat", [](controller_mv &sc) { return sc.d_noisemat; },
          "Noise on WFS measurements matrix")

      .def_property_readonly(
          "d_Cmm", [](controller_mv &sc) { return sc.d_Cmm; },
          "Slope covariance matrix")

      .def_property_readonly(
          "d_Cphim", [](controller_mv &sc) { return sc.d_Cphim; },
          "Actuators-Slopes covariance marix")

      .def_property_readonly(
          "h_Cmmeigenvals", [](controller_mv &sc) { return sc.h_Cmmeigenvals; },
          "Cmm eigenvalues")

      .def_property_readonly(
          "h_eigenvals", [](controller_mv &sc) { return sc.h_eigenvals; },
          "Eigen values")

      .def_property_readonly(
          "d_cenbuff", [](controller_mv &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly(
          "d_com1", [](controller_mv &sc) { return sc.d_com1; },
          "Commands at iteration k-1 (for POLC)")

      .def_property_readonly(
          "d_olmeas", [](controller_mv &sc) { return sc.d_olmeas; },
          "Reconstructed open loop measurement")

      .def_property_readonly(
          "d_err", [](controller_mv &sc) { return sc.d_err; },
          "Increment error")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("build_cmat",
           (int (controller_mv::*)(float)) & controller_mv::build_cmat,
           R"pbdoc(
    Computes the command matrix

    Args:
        cond: (float): Conditioning number for inversion
    )pbdoc",
           py::arg("cond"))

      .def("load_noisemat", wy::colCast(&controller_mv::load_noisemat),
           R"pbdoc(
    Load the noise covariance matrix

    Args:
        noisemat: (float): noise covariance marix
    )pbdoc",
           py::arg("noisemat"))

      .def("filter_cmat", &controller_mv::filter_cmat, R"pbdoc(
    Filter command matrix from TT

    Args:
        cond: (float): TODO: docstring
    )pbdoc",
           py::arg("cond"))

      .def("filter_cphim", wy::colCast(&controller_mv::filter_cphim),
           R"pbdoc(
    Filter Cphim from piston and apply coupling

    Args:
        F: (np.array[ndim=2,dtype=np.float32]): Piston filter matrix

        Nact: (np.array[ndim=2,dtype=np.float32]): Coupling matrix

    )pbdoc",
           py::arg("F"), py::arg("Nact"))

      .def("compute_Cmm", wy::colCast(&controller_mv::compute_Cmm),
           R"pbdoc(
    Compute the Cmm matrix

    Args:
        atmos : (SutraAtmos): SutraAtmos object

        sensors: (SutraSensors): SutraSensors object

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

      .def("compute_Cphim", wy::colCast(&controller_mv::compute_Cphim),
           R"pbdoc(
    Compute the Cphim matrix

    Args:
        atmos : (SutraAtmos): SutraAtmos object

        sensors: (SutraSensors): SutraSensors object

        dms: (SutraDms): SutraDms object

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

      .def("set_modal_gains", wy::colCast(&controller_mv::set_modal_gains), R"pbdoc(
    Set the controller modal gains

    Args:
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_cmat", wy::colCast(&controller_mv::set_cmat), R"pbdoc(
    Set the command matrix

    Args:
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      .def("set_imat", wy::colCast(&controller_mv::set_imat), R"pbdoc(
    Set the interaction matrix

    Args:
      imat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("imat"))

      ;
};

void declare_controller_mv(py::module &mod) {
  controller_mv_impl<float, float>(mod, "ControllerMV_FF");
  controller_mv_impl<float, uint16_t>(mod, "ControllerMV_FU");
}
