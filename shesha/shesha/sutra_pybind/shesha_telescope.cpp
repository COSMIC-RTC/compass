#include <wyrm>

#include <sutra_telescope.h>

namespace py = pybind11;

std::unique_ptr<sutra_telescope> telescope_init(
    carma_context &context, long n_pup, long npos, float *pupil,
    float *phase_ab_M1, long n_pup_m, float *pupil_m, float *phase_ab_M1_m) {
  return std::unique_ptr<sutra_telescope>(
      new sutra_telescope(&context, n_pup, npos, pupil, phase_ab_M1, n_pup_m,
                          pupil_m, phase_ab_M1_m));
}

void declare_shesha_telescope(py::module &mod) {
  py::class_<sutra_telescope>(mod, "Telescope")
      .def(py::init(wy::colCast(telescope_init)), R"pbdoc(
        Create and initialise a Telescope object
        Parameters
        ------------
        context: (carma_context) : current carma context
        n_pup: (long) : spupil size
        npos : (long): number of points in the pupil
        pupil: (np.ndarray[ndim=2, dtype=np.float32_t]) : spupil
        phase_ab_M1: (np.ndarray[ndim=2, dtype=np.float32_t]) : M1 aberrations on spupil support
        n_pup_m: (long) : mpupil size
        pupil_m: (np.ndarray[ndim=2, dtype=np.float32_t]) : mpupil
        phase_ab_M1_m: (np.ndarray[ndim=2, dtype=np.float32_t]) : M1 aberrations on mpupil support
        )pbdoc",
           py::arg("context"), py::arg("n_pup"), py::arg("npos"),
           py::arg("pupil"), py::arg("phase_ab_M1"), py::arg("n_pup_m"),
           py::arg("pupil_m"), py::arg("phase_ab_M1_m"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("device",
                             [](sutra_telescope &sp) { return sp.device; },
                             "Device number")

      .def_property_readonly("pup_size",
                             [](sutra_telescope &sp) { return sp.pup_size; },
                             "Small Pupil size")

      .def_property_readonly("pup_size_m",
                             [](sutra_telescope &sp) { return sp.pup_size_m; },
                             "Medium Pupil size")

      .def_property_readonly(
          "num_eleme_pup", [](sutra_telescope &sp) { return sp.num_eleme_pup; },
          "number of points in the pupil")

      .def_property_readonly("d_pupil",
                             [](sutra_telescope &sp) { return sp.d_pupil; },
                             "Small pupil of the Telescope")

      .def_property_readonly("d_pupil_m",
                             [](sutra_telescope &sp) { return sp.d_pupil_m; },
                             "Medium pupil of the Telescope")

      .def_property_readonly(
          "d_phase_ab_M1", [](sutra_telescope &sp) { return sp.d_phase_ab_M1; },
          "M1 aberrations on the small pupil")

      .def_property_readonly(
          "d_phase_ab_M1_m",
          [](sutra_telescope &sp) { return sp.d_phase_ab_M1_m; },
          "M1 aberrations on the medium pupil")

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_pupil",
           [](sutra_telescope &sp,
              py::array_t<float, py::array::f_style | py::array::forcecast>
                  data) {
             if (sp.d_pupil->getNbElem() == data.size())
               sp.d_pupil->host2device(data.mutable_data());
             else
               DEBUG_TRACE("Wrong dimensions");
           },
           R"pbdoc(
        Set the small pupil
        Parameters
        ------------
        pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  small pupil
      )pbdoc",
           py::arg("pup"))

      .def("set_pupil_m",
           [](sutra_telescope &sp,
              py::array_t<float, py::array::f_style | py::array::forcecast>
                  data) {
             if (sp.d_pupil_m->getNbElem() == data.size())
               sp.d_pupil_m->host2device(data.mutable_data());
             else
               DEBUG_TRACE("Wrong dimensions");
           },
           R"pbdoc(
        Set the medium pupil
        Parameters
        ------------
        pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  medium pupil
      )pbdoc",
           py::arg("pup"))

      .def("set_phase_ab_M1",
           [](sutra_telescope &sp,
              py::array_t<float, py::array::f_style | py::array::forcecast>
                  data) {
             if (sp.d_phase_ab_M1->getNbElem() == data.size())
               sp.d_phase_ab_M1->host2device(data.mutable_data());
             else
               DEBUG_TRACE("Wrong dimensions");
           },
           R"pbdoc(
        Set the M1 phase aberration in the small pupil
        Parameters
        ------------
        phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the small pupil
      )pbdoc",
           py::arg("phase_ab"))

      .def("set_phase_ab_M1_m",
           [](sutra_telescope &sp,
              py::array_t<float, py::array::f_style | py::array::forcecast>
                  data) {
             if (sp.d_phase_ab_M1_m->getNbElem() == data.size())
               sp.d_phase_ab_M1_m->host2device(data.mutable_data());
             else
               DEBUG_TRACE("Wrong dimensions");
           },
           R"pbdoc(
        Set the M1 phase aberration in the medium pupil
        Parameters
        ------------
        phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the medium pupil
      )pbdoc",
           py::arg("phase_ab"))

      ;
};
