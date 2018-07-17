#include <wyrm>

#include <sutra_target.h>

namespace py = pybind11;

std::unique_ptr<sutra_target> target_init(carma_context &context,
                                          sutra_telescope *d_tel, int ntargets,
                                          float *xpos, float *ypos,
                                          float *lambda, float *mag,
                                          float zerop, long *sizes, int Npts,
                                          int device) {
  return std::unique_ptr<sutra_target>(
      new sutra_target(&context, d_tel, ntargets, xpos, ypos, lambda, mag,
                       zerop, sizes, Npts, device));
}

void declare_shesha_target(py::module &mod) {
  py::class_<sutra_target>(mod, "Target")
      .def(py::init(wy::colCast(target_init)), R"pbdoc(
        Create and initialise an target object
        Parameters
        ------------
        context: (carma_context) : current carma context
        d_tel: (sutra_telescope) : sutra_telescope object
        ntargets: (int): number of targets
        xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : X positions of each target in arcsec
        ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : Y positions of each target in arcsec
        lambda: (np.ndarray[ndim=1,dtype=np.float32_t]) : Wavelength of each target in µm
        mag: (np.ndarray[ndim=1,dtype=np.float32_t]) : magnitude of each target
        zerop: (float) : Flux at magnitude 0 in photons/m²/s
        sizes: (np.ndarray[ndim=1,dtype=np.int64_t]) : Support size of each target
        Npts : (int): number of points in the pupil
        device: (int): GPU device index
        )pbdoc",
           py::arg("context"), py::arg("d_tel"), py::arg("ntargets"),
           py::arg("xpos"), py::arg("ypos"), py::arg("lambda"), py::arg("mag"),
           py::arg("zerop"), py::arg("sizes"), py::arg("Npts"),
           py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("ntargets",
                             [](sutra_target &st) { return st.ntargets; },
                             "Number of targets")

      .def_property_readonly("d_targets",
                             [](sutra_target &st) -> vector<sutra_source *> & {
                               return st.d_targets;
                             },
                             "Vector of targets")
      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      ;
};
