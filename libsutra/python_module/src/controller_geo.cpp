#include <wyrm>

#include <sutra_controller_geo.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;
using controller_geo = sutra_controller_geo<float>;

void declare_controller_geo(py::module &mod) {
  py::class_<controller_geo, sutra_controller<float>>(mod, "ControllerGEO")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("gain", [](controller_geo &sc) { return sc.gain; },
                             "Controller gain")

      .def_property_readonly("Nphi", [](controller_geo &sc) { return sc.Nphi; },
                             "Number of points in the pupil")

      .def_property_readonly("Ntt", [](controller_geo &sc) { return sc.Ntt; },
                             "Number of tip-tilt mirror")

      .def_property_readonly("d_gain",
                             [](controller_geo &sc) { return sc.d_gain; },
                             "vector of modal gains")

      .def_property_readonly("d_phi",
                             [](controller_geo &sc) { return sc.d_phi; },
                             "Phase in the pupil without piston (double)")

      .def_property_readonly("d_phif",
                             [](controller_geo &sc) { return sc.d_phif; },
                             "Phase in the pupil without piston (float)")

      .def_property_readonly("d_indx_pup",
                             [](controller_geo &sc) { return sc.d_indx_pup; },
                             "Indices of the valid pixels in spupil")

      .def_property_readonly("d_indx_mpup",
                             [](controller_geo &sc) { return sc.d_indx_mpup; },
                             "Indices of the valid pixels in mpupil")

      .def_property_readonly(
          "d_IFsparse", [](controller_geo &sc) { return sc.d_IFsparse; },
          "Influence functions in the pupil (sparse representation")

      .def_property_readonly("d_geocov",
                             [](controller_geo &sc) { return sc.d_geocov; },
                             "Geometric covariance matrix")

      .def_property_readonly("d_compdouble",
                             [](controller_geo &sc) { return sc.d_compdouble; },
                             "Buffer for computation (double precision)")

      .def_property_readonly("d_compfloat",
                             [](controller_geo &sc) { return sc.d_compfloat; },
                             "Buffer for computation (simple precision)")

      .def_property_readonly("d_TT", [](controller_geo &sc) { return sc.d_TT; },
                             "Tip-tilt influence functions")

      .def_property_readonly("d_geocovTT",
                             [](controller_geo &sc) { return sc.d_geocovTT; },
                             "Geometric covariance matrix for TT mirror")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("load_Btt", wy::colCast(&controller_geo::load_Btt),
           R"pbdoc(
               Load the Btt modal basis in the geo controller for ROKET

               Parameters
               ------------
               Btt_pzt: (np.array[ndim=2,dtype=np.float32]) : PZT DM component
               Btt_tt: (np.array[ndim=2,dtype=np.float32]) : TT mirror component
               )pbdoc",
           py::arg("Btt_pzt"), py::arg("Btt_tt"))

      .def("init_proj_sparse", wy::colCast(&controller_geo::init_proj_sparse),
           R"pbdoc(
        Initializes projection matrices

        Parameters
        ------------
        dms: (sutra_dms): sutra_dms object
        indx_dm: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the pupil in the DM support
        unitpervolt: (np.array[ndim=1,dtype=np.float32]): Unit per volt of each DM
        indx_pup: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the small pupil
        indx_mpup: (np.array[ndim=1,dtype=np.int64]): Indices of valid pixels of the medium pupil
        roket: (bool): ROKET flag
    )pbdoc",
           py::arg("dms"), py::arg("indx_dm"), py::arg("unitpervolt"),
           py::arg("indx_pup"), py::arg("indx_mpup"), py::arg("roket"))

      .def("comp_dphi", wy::colCast(&controller_geo::comp_dphi),
           R"pbdoc(
      Get the pupil phase and remove piston before projection

      Parameters
      ------------
      source: (sutra_source): Phase source
      wfs_direction: (bool): Must be True if the source is a WFS GS
    )pbdoc",
           py::arg("source"), py::arg("wfs_direction"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def("set_gain", wy::colCast(&controller_geo::set_gain), R"pbdoc(
      Set the controller gain

      Parameters
      ------------
      gain: (float): gain to set
    )pbdoc",
           py::arg("gain"))

      .def("load_mgain", wy::colCast(&controller_geo::load_mgain),
           R"pbdoc(
      Set the controller modal gains

      Parameters
      ------------
      mgain: (np.array[ndim1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("mgain"))

      ;
};
