#include <wyrm>

#include <sutra_rtc_brahma.h>

namespace py = pybind11;
using rtc = sutra_rtc<float>;
using rtcH = sutra_rtc<half>;
using rtc_brahma = sutra_rtc_brahma<float>;
using rtc_brahmaH = sutra_rtc_brahma<half>;

std::unique_ptr<rtc_brahma> rtc_brahma_init(carma_context *context,
                                            sutra_sensors *wfs,
                                            sutra_target *target,
                                            ACE_TCHAR *name) {
  return std::unique_ptr<rtc_brahma>(
      new rtc_brahma(context, wfs, target, name));
}

std::unique_ptr<rtc_brahmaH> rtc_brahmaH_init(carma_context *context,
                                              sutra_sensors *wfs,
                                              sutra_target *target,
                                              ACE_TCHAR *name) {
  return std::unique_ptr<rtc_brahmaH>(
      new rtc_brahmaH(context, wfs, target, name));
}

void declare_rtc_brahma(py::module &mod) {
  py::class_<rtc_brahma, rtc>(mod, "Rtc_brahma")
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
      //                        [](rtc_brahma &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](rtc_brahma &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](rtc_brahma &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](rtc_brahma &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](rtc_brahma &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](rtc_brahma &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](rtc_brahma &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](rtc_brahma &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](rtc_brahma &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &rtc_brahma::publish)

      ;

  py::class_<rtc_brahmaH, rtcH>(mod, "Rtc_brahmaH")
      .def(py::init(wy::colCast(rtc_brahmaH_init)), R"pbdoc(
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
      //                        [](rtc_brahmaH &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](rtc_brahmaH &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](rtc_brahmaH &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](rtc_brahmaH &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](rtc_brahmaH &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](rtc_brahmaH &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](rtc_brahmaH &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](rtc_brahmaH &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](rtc_brahmaH &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](rtc_brahmaH &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](rtc_brahmaH &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &rtc_brahmaH::publish)

      ;
};
