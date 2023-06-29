// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      wfs_sh.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraWfsSH
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#include <sutra_wfs_sh.h>

#include <wyrm>

namespace py = pybind11;

void declare_wfs_sh(py::module &mod) {
  py::class_<SutraWfsSH, SutraWfs>(mod, "SHWFS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "d_binmap", [](SutraWfsSH &ssh) { return ssh.d_binmap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_validpuppixx", [](SutraWfsSH &ssh) { return ssh.d_validpuppixx; },
          "X position of the bottom left corner of each ssp")

      .def_property_readonly(
          "d_validpuppixy", [](SutraWfsSH &ssh) { return ssh.d_validpuppixy; },
          "Y position of the bottom left corner of each ssp")
      .def_property_readonly(
          "d_fsamplipup", [](SutraWfsSH &ss) { return ss.d_fsamplipup; },
          "Complex amplitude in the pupil in the field stop array size")
      .def_property_readonly(
          "d_fsamplifoc", [](SutraWfsSH &ss) { return ss.d_fsamplifoc; },
          "Focal plane with field stop")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("load_arrays", wy::colCast(&SutraWfsSH::load_arrays), R"pbdoc(
    Load SH WFS arrays

    Args:
      phasemap: TODO: docstring

      hrmap:

      binmap:

      offsets:

      fluxPerSub: (np.array[ndim=2,dtype=np.float32]): Normalized flux per ssp

      validsubsx: (np.array[ndim=1, dtype=np.int64]): X position of each valid ssp

      validsubsy: (np.array[ndim=1, dtype=np.int64]): Y position of each valid ssp

      istart:

      jstart:

      ttprojmat: (np.array[ndim=2, dtype=np.float32]): slope projection matrix
                 for geom wfs.

      kernel:
    )pbdoc",
           py::arg("phasemap"), py::arg("hrmap"), py::arg("binmap"),
           py::arg("offsets"), py::arg("fluxPerSub"), py::arg("validsubsx"),
           py::arg("validsubsy"), py::arg("istart"), py::arg("jstart"),
           py::arg("ttprojmat"), py::arg("kernel"))

      .def("comp_nphot", &SutraWfsSH::comp_nphot, R"pbdoc(
    Compute the currect number of photons for a given system

    Args:
      ittime: (float): 1/loop frequency [s].

      optthroughput: (float): wfs global throughput.

      diam: (float):  telescope diameter.

      nxsub: (float): linear number of subaps.

      zerop: (float): (optional for LGS)  detector zero point expressed in ph/m**2/s in the bandwidth of the WFS.

      gsmag: (float): (optional for LGS)  magnitude of guide star.

      lgsreturnperwatt: (float): (optional for NGS) return per watt factor (high season : 10 ph/cm2/s/W).

      laserpower: (float): (optional for NGS) laser power in W.
    )pbdoc",
           py::arg("ittime"), py::arg("optthroughput"), py::arg("diam"),
           py::arg("nxsub"), py::arg("zerop") = 0, py::arg("gsmag") = 0,
           py::arg("lgsreturnperwatt") = 0, py::arg("laserpower") = 0)

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_bincube", wy::colCast(&SutraWfsSH::set_bincube), R"pbdoc(
    Set the bincube of the SH WFS

    Args:
        bincube: (np.array[ndim=3, dtype=np.float32]) : cube of subap. images

        nElem: (int): Number of elements in bincube
      )pbdoc",
           py::arg("bincube"), py::arg("nElem"));
};
