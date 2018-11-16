#include <wyrm>

#include <sutra_lgs.h>

namespace py = pybind11;

void declare_lgs(py::module &mod) {
  py::class_<sutra_lgs>(mod, "LGS")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("device", [](sutra_lgs &sl) { return sl.device; },
                             "GPU device index")

      .def_property_readonly("nvalid", [](sutra_lgs &sl) { return sl.nvalid; },
                             "TODO: docstring")

      .def_property_readonly("npix", [](sutra_lgs &sl) { return sl.npix; },
                             "TODO: docstring")

      .def_property_readonly("nmaxhr", [](sutra_lgs &sl) { return sl.nmaxhr; },
                             "Size of HR support")

      .def_property_readonly("hg", [](sutra_lgs &sl) { return sl.hg; },
                             "TODO: docstring")

      .def_property_readonly("h0", [](sutra_lgs &sl) { return sl.h0; },
                             "TODO: docstring")

      .def_property_readonly("deltah", [](sutra_lgs &sl) { return sl.deltah; },
                             "TODO: docstring")

      .def_property_readonly("pixsize",
                             [](sutra_lgs &sl) { return sl.pixsize; },
                             "Pixel size on sky[arcsec]")

      .def_property_readonly("nprof", [](sutra_lgs &sl) { return sl.nprof; },
                             "TODO: docstring")

      .def_property_readonly("d_doffaxis",
                             [](sutra_lgs &sl) { return sl.d_doffaxis; },
                             "TODO: docstring")

      .def_property_readonly("d_azimuth",
                             [](sutra_lgs &sl) { return sl.d_azimuth; },
                             "TODO: docstring")

      .def_property_readonly("d_prof1d",
                             [](sutra_lgs &sl) { return sl.d_prof1d; },
                             "TODO: docstring")

      .def_property_readonly("d_profcum",
                             [](sutra_lgs &sl) { return sl.d_profcum; },
                             "TODO: docstring")

      .def_property_readonly("d_prof2d",
                             [](sutra_lgs &sl) { return sl.d_prof2d; },
                             "TODO: docstring")

      .def_property_readonly("d_beam", [](sutra_lgs &sl) { return sl.d_beam; },
                             "TODO: docstring")

      .def_property_readonly("d_ftbeam",
                             [](sutra_lgs &sl) { return sl.d_ftbeam; },
                             "TODO: docstring")

      .def_property_readonly("d_lgskern",
                             [](sutra_lgs &sl) { return sl.d_lgskern; },
                             "TODO: docstring")

      .def_property_readonly("d_ftlgskern",
                             [](sutra_lgs &sl) { return sl.d_ftlgskern; },
                             "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("lgs_init", wy::colCast(&sutra_lgs::lgs_init), R"pbdoc(
        Initialize LGS object

        Parameters
        ------------
        nprof: (int): TODO: docstring
        hg: (float):
        h0: (float):
        deltah: (float):
        pixsize: (float):
        doffaxis:(np.array[ndim= , dtype=np.float32]):
        prof1d:(np.array[ndim= , dtype=np.float32]):
        profcum:(np.array[ndim= , dtype=np.float32]):
        beam:(np.array[ndim= , dtype=np.float32]):
        ftbeam:(np.array[ndim= , dtype=np.complex64]):
        azimuth:(np.array[ndim= , dtype=np.float32]):
    )pbdoc",
           py::arg("nprof"), py::arg("hg"), py::arg("h0"), py::arg("deltah"),
           py::arg("pixsize"), py::arg("doffaxis"), py::arg("prof1d"),
           py::arg("profcum"), py::arg("beam"), py::arg("ftbeam"),
           py::arg("azimuth"))
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      ;
};
