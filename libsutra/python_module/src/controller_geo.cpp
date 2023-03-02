// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      controller_geo.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_controller_geo
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <sutra_controller_geo.h>

#include <wyrm>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
void controller_geo_impl(py::module &mod, const char *name) {
  using controller_geo = sutra_controller_geo<Tcomp, Tout>;

  py::class_<controller_geo, SutraController<Tcomp, Tout>>(mod, name)

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "Nphi", [](controller_geo &sc) { return sc.Nphi; },
          "Number of points in the pupil")

      .def_property_readonly(
          "Ntt", [](controller_geo &sc) { return sc.Ntt; },
          "Number of tip-tilt mirror")

      .def_property_readonly(
          "d_gain", [](controller_geo &sc) { return sc.d_gain; },
          "vector of modal gains")

      .def_property_readonly(
          "d_phi", [](controller_geo &sc) { return sc.d_phi; },
          "Phase in the pupil without piston (double)")

      .def_property_readonly(
          "d_phif", [](controller_geo &sc) { return sc.d_phif; },
          "Phase in the pupil without piston (float)")

      .def_property_readonly(
          "d_indx_pup", [](controller_geo &sc) { return sc.d_indx_pup; },
          "Indices of the valid pixels in spupil")

      .def_property_readonly(
          "d_indx_mpup", [](controller_geo &sc) { return sc.d_indx_mpup; },
          "Indices of the valid pixels in mpupil")

      .def_property_readonly(
          "d_IFsparse", [](controller_geo &sc) { return sc.d_IFsparse; },
          "Influence functions in the pupil (sparse representation")

      .def_property_readonly(
          "d_geocov", [](controller_geo &sc) { return sc.d_geocov; },
          "Geometric covariance matrix")

      .def_property_readonly(
          "d_compdouble", [](controller_geo &sc) { return sc.d_compdouble; },
          "Buffer for computation (double precision)")

      .def_property_readonly(
          "d_compfloat", [](controller_geo &sc) { return sc.d_compfloat; },
          "Buffer for computation (simple precision)")

      .def_property_readonly(
          "d_TT", [](controller_geo &sc) { return sc.d_TT; },
          "Tip-tilt influence functions")

      .def_property_readonly(
          "d_geocovTT", [](controller_geo &sc) { return sc.d_geocovTT; },
          "Geometric covariance matrix for TT mirror")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("load_Btt", wy::colCast(&controller_geo::load_Btt),
           R"pbdoc(
    Load the Btt modal basis in the geo controller for ROKET

    Args:
        Btt_pzt: (np.array[ndim=2,dtype=np.float32]) : PZT DM component

        Btt_tt: (np.array[ndim=2,dtype=np.float32]) : TT mirror component
        )pbdoc",
           py::arg("Btt_pzt"), py::arg("Btt_tt"))

      .def("init_proj_sparse", wy::colCast(&controller_geo::init_proj_sparse),
           R"pbdoc(
    Initializes projection matrices

    Args:
        dms: (SutraDms): SutraDms object

        indx_dm: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the pupil in the DM support
        unitpervolt: (np.array[ndim=1,dtype=np.float32]): Unit per volt of each DM

        indx_pup: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the small pupil

        indx_mpup: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the medium pupil

        roket: (bool): ROKET flag
    )pbdoc",
           py::arg("dms"), py::arg("indx_dm"), py::arg("unitpervolt"),
           py::arg("indx_pup"), py::arg("indx_mpup"), py::arg("roket"))

      .def("comp_dphi", wy::colCast(&controller_geo::comp_dphi),
           R"pbdoc(
    Get the pupil phase and remove piston before projection

    Args:
      source: (SutraSource): Phase source

      wfs_direction: (bool): Must be True if the source is a WFS GS
    )pbdoc",
           py::arg("source"), py::arg("wfs_direction"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def("load_mgain", wy::colCast(&controller_geo::load_mgain),
           R"pbdoc(
    Set the controller modal gains

    Args:
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      ;
};

void declare_controller_geo(py::module &mod) {
  controller_geo_impl<float, float>(mod, "ControllerGEO_FF");
  controller_geo_impl<float, uint16_t>(mod, "ControllerGEO_FU");
}
