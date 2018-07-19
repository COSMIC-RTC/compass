#include <wyrm>

#include <sutra_wfs_sh.h>

namespace py = pybind11;

void declare_shesha_wfs_sh(py::module &mod) {
  py::class_<sutra_wfs_sh, sutra_wfs>(mod, "SHWFS")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("d_binmap",
                             [](sutra_wfs_sh &ssh) { return ssh.d_binmap; },
                             "TODO: docstring")

      .def_property_readonly("d_istart",
                             [](sutra_wfs_sh &ssh) { return ssh.d_istart; },
                             "X position of the bottom left corner of each ssp")

      .def_property_readonly("d_jstart",
                             [](sutra_wfs_sh &ssh) { return ssh.d_jstart; },
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

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      ;
};
