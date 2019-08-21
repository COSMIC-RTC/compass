#include <wyrm>

#include <sutra_wfs_sh.h>

namespace py = pybind11;

void declare_wfs_sh(py::module &mod) {
  py::class_<sutra_wfs_sh, sutra_wfs>(mod, "SHWFS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "d_binmap", [](sutra_wfs_sh &ssh) { return ssh.d_binmap; },
          "TODO: docstring")

      .def_property_readonly(
          "d_validpuppixx",
          [](sutra_wfs_sh &ssh) { return ssh.d_validpuppixx; },
          "X position of the bottom left corner of each ssp")

      .def_property_readonly(
          "d_validpuppixy",
          [](sutra_wfs_sh &ssh) { return ssh.d_validpuppixy; },
          "Y position of the bottom left corner of each ssp")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("loadarrays", wy::colCast(&sutra_wfs_sh::loadarrays), R"pbdoc(
      Load SH WFS arrays

      Parameters
      ------------
      phasemap: TODO: docstring
      hrmap:
      binmap:
      offsets:
      fluxPerSub: (np.array[ndim=2,dtype=np.float32]): Normalized flux per ssp
      validsubsx: (np.array[ndim=1, dtype=np.int64]): X position of each valid ssp
      validsubsy: (np.array[ndim=1, dtype=np.int64]): Y position of each valid ssp
      istart:
      jstart:
      kernel:
    )pbdoc",
           py::arg("phasemap"), py::arg("hrmap"), py::arg("binmap"),
           py::arg("offsets"), py::arg("fluxPerSub"), py::arg("validsubsx"),
           py::arg("validsubsy"), py::arg("istart"), py::arg("jstart"),
           py::arg("kernel"))

      .def("comp_nphot", &sutra_wfs_sh::comp_nphot, R"pbdoc(
      Compute the currect number of photons for a given system

      Parameters
      ------------
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
      .def("set_bincube", wy::colCast(&sutra_wfs_sh::set_bincube), R"pbdoc(
        Set the bincube of the SH WFS

        Parameters
        ------------
        bincube: (np.array[ndim=3, dtype=np.float32]) : cube of subap. images
        nElem: (int): Number of elements in bincube
      )pbdoc",
           py::arg("bincube"), py::arg("nElem"));
};
