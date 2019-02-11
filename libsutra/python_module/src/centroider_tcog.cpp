#include <wyrm>

#include <sutra_centroider_tcog.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_tcog_impl(py::module &mod, const char *name) {
  using centroider_tcog = sutra_centroider_tcog<Tin, Tcomp>;

  py::class_<centroider_tcog, sutra_centroider<Tin, Tcomp>>(mod, name)
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly("threshold",
                             [](centroider_tcog &sc) { return sc.threshold; },
                             "Threshold value")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_threshold", &centroider_tcog::set_threshold,
           R"pbdoc(
            Set the threshold value of a TCOG centroider

            Parameters
            ------------
            thresh : (float) : threshold value to set

        )pbdoc",
           py::arg("thresh"))

      ;
};
void declare_centroider_tcog(py::module &mod) {
  centroider_tcog_impl<float, float>(mod, "CentroiderTCOG_FF");
  centroider_tcog_impl<uint16_t, float>(mod, "CentroiderTCOG_UF");
}
