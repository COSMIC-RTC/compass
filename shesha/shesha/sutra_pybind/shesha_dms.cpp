#include <wyrm>

#include <sutra_dm.h>

namespace py = pybind11;

std::unique_ptr<sutra_dms> dms_init(int ndms) {
  return std::unique_ptr<sutra_dms>(new sutra_dms(ndms));
};

void declare_shesha_dms(py::module &mod) {
  py::class_<sutra_dms>(mod, "Dms")
      .def(py::init(wy::castParameter(dms_init)), R"pbdoc(
            Create a void DMS object
            Parameters
            ------------
            ndms: (int): number of DM that will be included in DMS
        )pbdoc",
           py::arg("ndms"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly("d_dms",
                             [](sutra_dms &sdms) { return sdms.d_dms; },
                             "Vector of sutra_dm")

      .def_property_readonly("d_dms",
                             [](sutra_dms &sdms) { return sdms.d_dms; },
                             "Vector of sutra_dm")
      .def_property_readonly("ndm",
                             [](sutra_dms &sdms) { return sdms.d_dms.size(); },
                             "Number of sutra_dm in sutra_dms")
      .def_property_readonly("nact_total",
                             [](sutra_dms &sdms) { return sdms.nact_total(); },
                             "Total number of actuators in sutra_dms")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("add_dm",
           [](sutra_dms &sdms, carma_context &context, const char *type,
              float alt, long dim, long ninflu, long influsize, long ninflupos,
              long n_npoints, float push4imat, long nord, int device) {
             return sdms.add_dm(&context, type, alt, dim, ninflu, influsize,
                                ninflupos, n_npoints, push4imat, nord, device);
           },
           R"pbdoc(
        Add a sutra_dm in the sutra_dms vector

        Parameters
        ------------
        context: (carma_context) : current carma context
        type: (str): DM type ("pzt", "kl", or "tt")
        alt: (float): Conjugaison altitude in meters
        dim: (long): Support dimension
        ninflu: (long): Number of actuators
        influsize: (long): Influenction function support size
        ninflupos: (long): Size of _influpos array
        n_npoints: (long): Size of _ninflu array
        push4imat: (float): Voltage to apply for imat computation
        nord: (long): Number of radial order for kl dm (0 if not kl)
        device: (int): Device index
        )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("alt"), py::arg("dim"),
           py::arg("ninflu"), py::arg("influsize"), py::arg("ninflupos"),
           py::arg("n_npoints"), py::arg("push4imat"), py::arg("nord"),
           py::arg("device"))

      .def("remove_dm", wy::castParameter(&sutra_dms::remove_dm),
           R"pbdoc(
        Remove and delete the selected DM from sutra_dms
        Parameters
        ------------
        idx: (int): index of DM
        )pbdoc",
           py::arg("idx"))

      ;
};

void declare_shesha_dm(py::module &mod) { py::class_<sutra_dms>(mod, "Dm"); };
