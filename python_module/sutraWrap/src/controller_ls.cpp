// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      controller_ls.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraControllerLs
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_controller_ls.hpp>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
int32_t init_modalOpti(SutraControllerLs<Tcomp, Tout> &sc, int32_t nmodes, int32_t nrec, ArrayFStyle<Tcomp> &M2V, Tcomp gmin, Tcomp gmax,
                    int32_t ngain, Tcomp Fs) {
    return sc.init_modalOpti(nmodes, nrec, M2V.mutable_data(), gmin, gmax, ngain, Fs);
}

template <typename Tcomp, typename Tout>
int32_t loadopen_loopSlp(SutraControllerLs<Tcomp, Tout> &sc, ArrayFStyle<Tcomp> &ol_slopes) {
    return sc.loadopen_loopSlp(ol_slopes.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_modal_gains(SutraControllerLs<Tcomp, Tout> &sc, ArrayFStyle<Tcomp> &mgain) {
    return sc.set_modal_gains(mgain.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_cmat(SutraControllerLs<Tcomp, Tout> &sc, ArrayFStyle<Tcomp> &cmat) {
    return sc.set_cmat(cmat.mutable_data());
}

template <typename Tcomp, typename Tout>
int32_t set_imat(SutraControllerLs<Tcomp, Tout> &sc, ArrayFStyle<Tcomp> &imat) {
    return sc.set_imat(imat.mutable_data());
}


template <typename Tcomp, typename Tout>
void controller_ls_impl(py::module &mod, const char *name) {
  using controller_ls = SutraControllerLs<Tcomp, Tout>;

  py::class_<controller_ls, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "d_imat", [](controller_ls &sc) { return sc.d_imat; },
          "Interaction matrix")

      .def_property_readonly(
          "d_cmat", [](controller_ls &sc) { return sc.d_cmat; },
          "Control matrix")

      .def_property_readonly(
          "d_gain", [](controller_ls &sc) { return sc.d_gain; },
          "vector of modal gains")

      .def_property_readonly(
          "d_eigenvals", [](controller_ls &sc) { return sc.d_eigenvals; },
          "Eigen values of the imat")

      .def_property_readonly(
          "d_U", [](controller_ls &sc) { return sc.d_U; },
          "Eigen modes of the imat")

      .def_property_readonly(
          "d_cenbuff", [](controller_ls &sc) { return sc.d_cenbuff; },
          "Centroids circular buffer")

      .def_property_readonly(
          "d_err", [](controller_ls &sc) { return sc.d_err; },
          "Current increment on the command")

      .def_property_readonly(
          "is_modopti", [](controller_ls &sc) { return sc.is_modopti; },
          "Falg for modal optimization")

      .def_property_readonly(
          "nrec", [](controller_ls &sc) { return sc.nrec; },
          "Number of open loop slopes to take for modal optimization")

      .def_property_readonly(
          "nmodes", [](controller_ls &sc) { return sc.nmodes; },
          "Number of modes for modal optimization")

      .def_property_readonly(
          "gmin", [](controller_ls &sc) { return sc.gmin; },
          "Minimal gain for modal optimization")

      .def_property_readonly(
          "gmax", [](controller_ls &sc) { return sc.gmax; },
          "Maximal gain for modal optimization")

      .def_property_readonly(
          "ngain", [](controller_ls &sc) { return sc.ngain; },
          "Number of gain values to test between gmin and "
          "gmax for modal optimization")

      .def_property_readonly(
          "Fs", [](controller_ls &sc) { return sc.Fs; },
          "Sampling frequency for modal optimization")

      .def_property_readonly(
          "cpt_rec", [](controller_ls &sc) { return sc.cpt_rec; },
          "Counter for modal gains refresh")

      .def_property_readonly(
          "d_M2V", [](controller_ls &sc) { return sc.d_M2V; },
          "Modes to volt matrix for modal optimization")

      .def_property_readonly(
          "d_S2M", [](controller_ls &sc) { return sc.d_S2M; },
          "Slopes to modes matrix for modal optimization")

      .def_property_readonly(
          "d_slpol", [](controller_ls &sc) { return sc.d_slpol; },
          "Open loop slopes for modal optimization")

      .def_property_readonly(
          "d_Hcor", [](controller_ls &sc) { return sc.d_Hcor; },
          "Transfer function for modal optimization")

      .def_property_readonly(
          "d_compbuff", [](controller_ls &sc) { return sc.d_compbuff; },
          "Buffer for POLC computation")

      .def_property_readonly(
          "d_compbuff2", [](controller_ls &sc) { return sc.d_compbuff2; },
          "Buffer for POLC computation")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("svdec_imat", &controller_ls::svdec_imat,
           "Performs interaction matrix SVD")

      .def("build_cmat",
           (int32_t (controller_ls::*)(int32_t)) &controller_ls::build_cmat,
           R"pbdoc(
    Computes the command matrix after imat SVD

    Args:
        nfilt: (int): number of modes to filter
    )pbdoc",
           py::arg("nfilt"))

      .def("init_modalOpti", &init_modalOpti<Tcomp, Tout>,
           R"pbdoc(
    Initialize modal optimization control

    Args:
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

      .def("loadopen_loopSlp", &loadopen_loopSlp<Tcomp, Tout>,
           R"pbdoc(
    Load recorded open loop slopes for modal optimization initialization

    Args:
      slopes: (np.array[ndim=2,dtype=np.float32]): Open loop slopes
    )pbdoc",
           py::arg("slopes"))

      .def("modalControlOptimization", &controller_ls::modalControlOptimization,
           "TODO: docstring")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
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
      imat: (np.array[ndim=2,dtype=np.float32]): interaction matrix to set
    )pbdoc",
           py::arg("imat"))

      ;
};

void declare_controller_ls(py::module &mod) {
  controller_ls_impl<float, float>(mod, "ControllerLS_FF");
  controller_ls_impl<float, uint16_t>(mod, "ControllerLS_FU");
}
