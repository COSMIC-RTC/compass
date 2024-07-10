// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      controller_mv.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraControllerMv
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_controller_mv.hpp>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
int32_t load_noisemat(SutraControllerMv<Tcomp, Tout> &scmv, ArrayFStyle<Tcomp> &noise) {
    return scmv.load_noisemat(noise.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t filter_cphim(SutraControllerMv<Tcomp, Tout> &scmv, ArrayFStyle<Tcomp> &F, ArrayFStyle<Tcomp> &Nact) {
    return scmv.filter_cphim(F.mutable_data(), Nact.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t compute_Cmm(SutraControllerMv<Tcomp, Tout> &scmv, SutraAtmos *atmos, SutraSensors *sensors, ArrayFStyle<double> &L0,
                ArrayFStyle<double> &cn2, ArrayFStyle<double> &alphaX, ArrayFStyle<double> &alphaY, double diamTel,
                double cobs){
    return scmv.compute_Cmm(atmos, sensors, L0.mutable_data(), cn2.mutable_data(), alphaX.mutable_data(), alphaY.mutable_data(), diamTel, cobs);
}

template <typename Tcomp, typename Tout>
int32_t compute_Cphim(SutraControllerMv<Tcomp, Tout> &scmv, SutraAtmos *atmos, SutraSensors *sensors, SutraDms *dms,
                ArrayFStyle<double> &L0, ArrayFStyle<double> &cn2, ArrayFStyle<double> &alphaX, ArrayFStyle<double> &alphaY,
                ArrayFStyle<double> &X, ArrayFStyle<double> &Y, ArrayFStyle<double> &xactu, ArrayFStyle<double> &yactu,
                double diamTel, ArrayFStyle<double> &k2, ArrayFStyle<int64_t> &NlayerDm,
                ArrayFStyle<int64_t> &indLayerDm, double FoV, ArrayFStyle<double> &pitch,
                ArrayFStyle<double> &alt_dm){
    return scmv.compute_Cphim(atmos, sensors, dms, L0.mutable_data(), cn2.mutable_data(), alphaX.mutable_data(), alphaY.mutable_data(), X.mutable_data(),
                            Y.mutable_data(), xactu.mutable_data(), yactu.mutable_data(), diamTel, k2.mutable_data(), NlayerDm.mutable_data(),
                            indLayerDm.mutable_data(), FoV, pitch.mutable_data(), alt_dm.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_modal_gains(SutraControllerMv<Tcomp, Tout> &scmv, ArrayFStyle<Tcomp> &mgain) {
    return scmv.set_modal_gains(mgain.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_cmat(SutraControllerMv<Tcomp, Tout> &scmv, ArrayFStyle<Tcomp> &cmat) {
    return scmv.set_cmat(cmat.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_imat(SutraControllerMv<Tcomp, Tout> &scmv, ArrayFStyle<Tcomp> &imat) {
    return scmv.set_imat(imat.mutable_data());
}


template <typename Tcomp, typename Tout>
void controller_mv_impl(py::module &mod, const char *name) {
  using controller_mv = SutraControllerMv<Tcomp, Tout>;

  py::class_<controller_mv, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

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
           [](controller_mv &sc, float cond) { return sc.build_cmat(cond); },
           R"pbdoc(
    Computes the command matrix

    Args:
        cond: (float): Conditioning number for inversion
    )pbdoc",
           py::arg("cond"))

      .def("load_noisemat", &load_noisemat<Tcomp, Tout>,
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

      .def("filter_cphim", &filter_cphim<Tcomp, Tout>,
           R"pbdoc(
    Filter Cphim from piston and apply coupling

    Args:
        F: (np.array[ndim=2,dtype=np.float32]): Piston filter matrix

        Nact: (np.array[ndim=2,dtype=np.float32]): Coupling matrix

    )pbdoc",
           py::arg("F"), py::arg("Nact"))

      .def("compute_Cmm", &compute_Cmm<Tcomp, Tout>,
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

      .def("compute_Cphim", &compute_Cphim<Tcomp, Tout>,
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

      .def("set_modal_gains", &set_modal_gains<Tcomp, Tout>, R"pbdoc(
    Set the controller modal gains

    Args:
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      .def("set_cmat", &set_cmat<Tcomp, Tout>, R"pbdoc(
    Set the command matrix

    Args:
      cmat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("cmat"))

      .def("set_imat", &set_imat<Tcomp, Tout>, R"pbdoc(
    Set the interaction matrix

    Args:
      imat: (np.array[ndim=2,dtype=np.float32]): command matrix to set
    )pbdoc",
           py::arg("imat"))

      ;
}

void declare_controller_mv(py::module &mod) {
  controller_mv_impl<float, float>(mod, "ControllerMV_FF");
  controller_mv_impl<float, uint16_t>(mod, "ControllerMV_FU");
}
