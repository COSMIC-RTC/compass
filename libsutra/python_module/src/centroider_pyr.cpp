#include <wyrm>

#include <sutra_centroider_pyr.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_pyr_impl(py::module &mod, const char *name) {
  using centroider_pyr = sutra_centroider_pyr<Tin, Tcomp>;

  py::class_<centroider_pyr, sutra_centroider<Tin, Tcomp>>(mod, name)

      .def_property_readonly(
          "pyr_method", [](centroider_pyr &sc) { return sc.get_method_str(); },
          "Method used for pyramid slopes compuation")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_pyr_method",
           [](centroider_pyr &sc, uint8_t method) {
             return sc.set_method(method);
           },
           R"pbdoc(
            Set the pyramid method for slopes computation

            Parameters
            ------------
            method : (int) : new centroiding method (0: nosinus global
                                                    1: sinus global
                                                    2: nosinus local
                                                    3: sinus local)
                            favor use of shesha_constant.PyrCentroiderMethod

        )pbdoc",
           py::arg("method"))

      .def("set_pyr_thresh",
           [](centroider_pyr &sc, float thresh) {
             return sc.set_valid_thresh(thresh);
           },
           R"pbdoc(
            Set the pyramid threshold value

            Parameters
            ------------
            thresh : (float) : threshold value
        )pbdoc",
           py::arg("thresh"));
};
void declare_centroider_pyr(py::module &mod) {
  centroider_pyr_impl<float, float>(mod, "CentroiderPYR_FF");
  centroider_pyr_impl<uint16_t, float>(mod, "CentroiderPYR_UF");
}
