#include <wyrm>

#include <sutra_rtc_brahma.h>

namespace py = pybind11;

std::unique_ptr<sutra_rtc_brahma> rtc_brahma_init(carma_context *context,
                                                  sutra_sensors *wfs,
                                                  sutra_target *target,
                                                  ACE_TCHAR *name) {
  return std::unique_ptr<sutra_rtc_brahma>(
      new sutra_rtc_brahma(context, wfs, target, name));
}

void declare_rtc_brahma(py::module &mod) {
  py::class_<sutra_rtc_brahma, sutra_rtc>(mod, "Rtc_brahma")
      .def(py::init(wy::colCast(rtc_brahma_init)), R"pbdoc(
        Create and initialise a brahma rtc object

        Parameters
        ------------
        context: (carma_context) : current carma context
        wfs:
        target:
        name:
        )pbdoc",
           py::arg("context"), py::arg("wfs"), py::arg("target"),
           py::arg("name"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      // .def_property_readonly("framecounter",
      //                        [](sutra_rtc_brahma &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](sutra_rtc_brahma &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &sutra_rtc_brahma::publish)

      ;
};
