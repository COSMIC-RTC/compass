#include <wyrm>

#include <sutra_centroider_bpcog.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

void declare_centroider_bpcog(py::module &mod) {
  py::class_<sutra_centroider_bpcog, sutra_centroider<float>>(mod,
                                                              "CentroiderBPCOG")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly("nmax",
                             [](sutra_centroider_bpcog &sc) { return sc.nmax; },
                             "Number of brightest pixels")

      .def_property_readonly(
          "d_bpix", [](sutra_centroider_bpcog &sc) { return sc.d_bpix; },
          "Brightest pixels")

      .def_property_readonly(
          "d_bpind", [](sutra_centroider_bpcog &sc) { return sc.d_bpind; },
          "Brightest pixels indices")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_nmax", &sutra_centroider_bpcog::set_nmax,
           R"pbdoc(
            Set the number of brightest pixels considered for COG computation

            Parameters
            ------------
            nmax : (float) : nmax value to set

        )pbdoc",
           py::arg("nmax"))

      ;
};
