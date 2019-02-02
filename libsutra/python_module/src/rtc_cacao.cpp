#include <wyrm>

#include <sutra_rtc_cacao.h>

namespace py = pybind11;

template <typename T>
std::unique_ptr<sutra_rtc_cacao<T>> rtc_cacao_init(
    carma_context *context, std::string iCalFrame_name,
    std::string iLoopFrame_name) {
  return std::unique_ptr<sutra_rtc_cacao<T>>(
      new sutra_rtc_cacao<T>(context, iCalFrame_name, iLoopFrame_name));
}

template <typename T>
void rtc_cacao_impl(py::module &mod, const char *name) {
  using rtc = sutra_rtc;
  using rtc_cacao = sutra_rtc_cacao<T>;

  py::class_<rtc_cacao, rtc>(mod, name)
      .def(py::init(wy::colCast(rtc_cacao_init<T>)), R"pbdoc(
            Create and initialise a cacao rtc object

            Parameters
            ------------
            context: (carma_context) : current carma context
            iCalFrame_name:
            iLoopFrame_name:
            )pbdoc",
           py::arg("iCalFrame_name"), py::arg("iLoopFrame_name"),
           py::arg("name"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      // .def_property_readonly("framecounter",
      //                        [](rtc_cacao &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("wfs_size",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs_size;
      //                        },
      //                        "WFS size")

      // .def_property_readonly("wfs_phase_size",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs_phase_size;
      //                        },
      //                        "WFS phase support size")

      // .def_property_readonly("wfs",
      //                        [](rtc_cacao &st) {
      //                          return st.wfs;
      //                        },
      //                        "WFS object")

      // .def_property_readonly("target_size",
      //                        [](rtc_cacao &st) {
      //                          return st.target_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target_phase_size",
      //                        [](rtc_cacao &st) {
      //                          return st.target_phase_size;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("target",
      //                        [](rtc_cacao &st) {
      //                          return st.target;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nslp",
      //                        [](rtc_cacao &st) {
      //                          return st.nslp;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("ncmd",
      //                        [](rtc_cacao &st) {
      //                          return st.ncmd;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("nvalid",
      //                        [](rtc_cacao &st) {
      //                          return st.nvalid;
      //                        },
      //                        "TODO: docstring")

      // .def_property_readonly("is_initialised",
      //                        [](rtc_cacao &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "TODO: docstring")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &rtc_cacao::publish)

      ;
}

void declare_rtc_cacao(py::module &mod) {
#ifdef CAN_DO_HALF
  // rtc_cacao_impl<half>(mod, "Rtc_cacaoH");
#endif
  rtc_cacao_impl<float>(mod, "Rtc_cacao");
};
