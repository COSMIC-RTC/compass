#include <wyrm>

#include <sutra_turbu.h>

namespace py = pybind11;

std::unique_ptr<sutra_atmos> atmos_init(carma_context &context, int nscreens,
                                        float global_r0, float *r0_per_layers,
                                        long *dim_screens, long *stencil_size,
                                        float *altitude, float *windspeed,
                                        float *winddir, float *deltax,
                                        float *deltay, int device) {
  return std::unique_ptr<sutra_atmos>(new sutra_atmos(
      &context, nscreens, global_r0, r0_per_layers, dim_screens, stencil_size,
      altitude, windspeed, winddir, deltax, deltay, device));
}

void declare_shesha_atmos(py::module &mod) {
  py::class_<sutra_atmos>(mod, "Atmos")
      .def(py::init(wy::castParameter(atmos_init)), R"pbdoc(
        Create and initialise an atmos object
        Parameters
        ------------
        context: (carma_context) : current carma context
        nscreens: (float) : number of turbulent layers
        r0_per_layer: (float) : global r0
        size1: (np.ndarray[ndim=1, dtype=np.int64_t]) : First size of screens
        size2: (np.ndarray[ndim=1, dtype=np.int64_t]) : Second size of screens
        altitude: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes [m]
        windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind speed [m/s]
        winddir: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind direction [deg]
        deltax: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltax pixels in the x-direction at each iteration
        deltay: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltay pixels in the y-direction at each iteration
        device: (int): GPU device index
        )pbdoc",
           py::arg("context"), py::arg("nscreens"), py::arg("global_r0"),
           py::arg("r0_per_layer"), py::arg("size1"), py::arg("size2"),
           py::arg("altitude"), py::arg("windspeed"), py::arg("winddir"),
           py::arg("deltax"), py::arg("deltay"), py::arg("device"))

      .def("move_atmos", &sutra_atmos::move_atmos, R"pbdoc(
        Move the turbulence in the atmos screen following loaded
        parameters such as windspeed and wind direction
        )pbdoc")

      .def("set_global_r0", wy::castParameter(&sutra_atmos::set_global_r0),
           R"pbdoc(
        Change the current global r0 of all layers

        Parameters
        ------------
        r0: (float): r0 @ 0.5 microns
        )pbdoc",
           py::arg("r0"))

      .def("init_screen", wy::castParameter(&sutra_atmos::init_screen), R"pbdoc(
        Initialize an newly allocated screen

        Parameters
        ------------
        altitude: (float): altitude of the screen
        A: (np.ndarray[ndim=2, dtype=np.float32]): A matrix (cf. Assemat)
        B: (np.ndarray[ndim=2, dtype=np.float32]): B matrix (cf. Assemat)
        istencilx: (np.ndarray[ndim=2, dtype=int32]): X stencil index
        istencily: (np.ndarray[ndim=2, dtype=int32]): Y stencil index
        seed: (int): seed for RNG
        )pbdoc",
           py::arg("altitude"), py::arg("A"), py::arg("B"),
           py::arg("istencilx"), py::arg("istencily"), py::arg("seed"))

      .def("add_screen", wy::castParameter(&sutra_atmos::add_screen), R"pbdoc(
        Add a screen to the atmos object.

        Parameters
        ------------
        altitude: (float) : altitude of the screen in meters
        size: (long) : dimension of the screen (size x size)
        stencil_size: (long): dimension of the stencil
        r0: (float) : frac of r0**(5/3)
        windspeed: (float) : windspeed of the screen [m/s]
        winddir: (float) : wind direction (deg)
        deltax: (float) : extrude deltax pixels in the x-direction at each iteration
        deltay: (float) : extrude deltay pixels in the y-direction at each iteration
        device: (int) : device number
        )pbdoc",
           py::arg("altitude"), py::arg("size"), py::arg("stencil_size"),
           py::arg("r0"), py::arg("windspeed"), py::arg("winddir"),
           py::arg("deltax"), py::arg("deltay"), py::arg("device"))

      .def("refresh_screen", wy::castParameter(&sutra_atmos::refresh_screen),
           R"pbdoc(
        Refresh the selected screen by extrusion
        Parameters
        ------------
        alt: (float): altitude of the screen to refresh
        )pbdoc",
           py::arg("alt"))

      .def("del_screen", wy::castParameter(&sutra_atmos::del_screen), R"pbdoc(
        Delete the selected screen
        Parameters
        ------------
        alt: (float): altitude of the screen to delete
        )pbdoc",
           py::arg("alt"))

      .def("set_seed", wy::castParameter(&sutra_atmos::set_seed), R"pbdoc(
        Set the seed of the selected screen RNG

        Parameters
        ------------
        alt: (float) :altitude of the screen to modify
        seed: (int) :new seed
        )pbdoc",
           py::arg("alt"), py::arg("seed"))

      .def("get_screen",
           [](sutra_atmos &sa, float alt) {
             sutra_tscreen *tscreen = sa.d_screens[alt];
             // py::array<float> screen(tscreen->d_screen->getNbElem());
             return tscreen->d_tscreen->d_screen;
           },
           py::return_value_policy::reference_internal);
};
