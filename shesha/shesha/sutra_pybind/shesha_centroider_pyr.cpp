#include <wyrm>

#include <sutra_centroider_pyr.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_shesha_centroider_pyr(py::module &mod) {
  py::class_<sutra_centroider_pyr, sutra_centroider>(mod, "CentroiderPYR")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_pyr_method",
           [](sutra_centroider_pyr &sc, uint8_t method) {
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
           [](sutra_centroider_pyr &sc, float thresh) {
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
