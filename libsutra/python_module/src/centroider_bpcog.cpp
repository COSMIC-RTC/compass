#include <wyrm>

#include <sutra_centroider_bpcog.h>

namespace py = pybind11;

template <typename Tin, typename Tcomp>
void centroider_bpcog_impl(py::module &mod, const char *name) {
  using centroider_bpcog = sutra_centroider_bpcog<Tin, Tcomp>;

  py::class_<centroider_bpcog, sutra_centroider<Tin, Tcomp>>(mod, name)
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly("nmax",
                             [](centroider_bpcog &sc) { return sc.nmax; },
                             "Number of brightest pixels")

      .def_property_readonly("d_bpix",
                             [](centroider_bpcog &sc) { return sc.d_bpix; },
                             "Brightest pixels")

      .def_property_readonly("d_bpind",
                             [](centroider_bpcog &sc) { return sc.d_bpind; },
                             "Brightest pixels indices")
      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_nmax", &centroider_bpcog::set_nmax,
           R"pbdoc(
            Set the number of brightest pixels considered for COG computation

            Parameters
            ------------
            nmax : (float) : nmax value to set

        )pbdoc",
           py::arg("nmax"))

      ;
};

void declare_centroider_bpcog(py::module &mod) {
  centroider_bpcog_impl<float, float>(mod, "CentroiderBPCOG_FF");
  centroider_bpcog_impl<uint16_t, float>(mod, "CentroiderBPCOG_UF");
#ifdef CAN_DO_HALF
  centroider_bpcog_impl<float, half>(mod, "CentroiderBPCOG_FH");
  centroider_bpcog_impl<uint16_t, half>(mod, "CentroiderBPCOG_UH");

#endif
}
