#include <wyrm>

#include <sutra_target_brahma.h>

namespace py = pybind11;

std::unique_ptr<sutra_target_brahma> target_brahma_init(
    carma_context *context, ACE_TCHAR *name, sutra_telescope *d_tel,
    int subsample, int ntargets, float *xpos, float *ypos, float *lambda,
    float *mag, float zerop, long *sizes, int Npts, int device) {
  return std::unique_ptr<sutra_target_brahma>(
      new sutra_target_brahma(context, name, d_tel, subsample, ntargets, xpos,
                              ypos, lambda, mag, zerop, sizes, Npts, device));
}

void declare_target_brahma(py::module &mod) {
  py::class_<sutra_target_brahma, sutra_target>(mod, "Target_brahma")
      .def(py::init(wy::colCast(target_brahma_init)), R"pbdoc(
        Create and initialise a brahma target object

        Parameters
        ------------
        context: (carma_context) : current carma context
        name:
        subsample:
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
           py::arg("context"), py::arg("name"), py::arg("subsample"),
           py::arg("d_tel"), py::arg("ntargets"), py::arg("xpos"),
           py::arg("ypos"), py::arg("lambda"), py::arg("mag"), py::arg("zerop"),
           py::arg("sizes"), py::arg("Npts"), py::arg("device"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      // .def_property_readonly("framecounter",
      //                        [](sutra_target_brahma &st) { return
      //                        st.framecounter; }, "Frame counter")

      // .def_property_readonly("samplecounter",
      //                        [](sutra_target_brahma &st) {
      //                          return st.samplecounter;
      //                        },
      //                        "Sample counter")

      // .def_property_readonly("subsample",
      //                        [](sutra_target_brahma &st) {
      //                          return st.subsample;
      //                        },
      //                        "Subsample")

      // .def_property_readonly("is_initialised",
      //                        [](sutra_target_brahma &st) {
      //                          return st.is_initialised;
      //                        },
      //                        "is_initialised flag")
      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("publish", &sutra_target_brahma::publish)

      ;
};
