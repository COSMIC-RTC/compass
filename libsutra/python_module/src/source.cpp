#include <wyrm>

#include <sutra_source.h>
#include "declare_name.hpp"

namespace py = pybind11;

void declare_source(py::module &mod) {
  py::class_<sutra_source>(mod, "Source")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](sutra_source &ss) { return ss.device; },
          "GPU device index")

      .def_property_readonly(
          "posx", [](sutra_source &ss) { return ss.posx; },
          "X position of the source")

      .def_property_readonly(
          "posy", [](sutra_source &ss) { return ss.posy; },
          "Y position of the source")

      .def_property_readonly(
          "npts", [](sutra_source &ss) { return ss.npts; },
          "Number of points in the pupil")

      .def_property_readonly(
          "mag", [](sutra_source &ss) { return ss.mag; },
          "Magnitude of the source")

      .def_property_readonly(
          "lambda", [](sutra_source &ss) { return ss.lambda; },
          "Wavelength of the source")

      .def_property_readonly(
          "zp", [](sutra_source &ss) { return ss.zp; }, "Flux at magnitude 0")

      .def_property_readonly(
          "scale", [](sutra_source &ss) { return ss.scale; },
          "Phase scale factor (2*pi/lambda)")

      .def_property_readonly(
          "lgs", [](sutra_source &ss) { return ss.lgs; }, "Boolean for LGS")

      .def_property_readonly(
          "G", [](sutra_source &ss) { return ss.G; },
          "Magnifying factor for WFS misalignment")

      .def_property_readonly(
          "thetaML", [](sutra_source &ss) { return ss.thetaML; },
          "Pupil rotation angle for WFS misalignment")

      .def_property_readonly(
          "dx", [](sutra_source &ss) { return ss.dx; },
          "X axis WFS misalignment [pixels]")

      .def_property_readonly(
          "dy", [](sutra_source &ss) { return ss.dy; },
          "Y axis WFS misalignment [pixels]")

      .def_property_readonly(
          "type", [](sutra_source &ss) { return ss.type; },
          "Type of source (Target or WFS GS)")

      .def_property_readonly(
          "block_size", [](sutra_source &ss) { return ss.block_size; },
          "Optimum block size of device")

      .def_property_readonly(
          "strehl_se", [](sutra_source &ss) { return ss.strehl_se; },
          "Short exposure Strehl ratio")

      .def_property_readonly(
          "strehl_le", [](sutra_source &ss) { return ss.strehl_le; },
          "Long exposure Strehl ratio")

      .def_property_readonly(
          "ref_strehl", [](sutra_source &ss) { return ss.ref_strehl; },
          "reference for Strehl computation (Airy)")

      .def_property_readonly(
          "strehl_counter", [](sutra_source &ss) { return ss.strehl_counter; },
          "Counter for LE Strehl computation")

      .def_property_readonly(
          "phase_var", [](sutra_source &ss) { return ss.phase_var; },
          "Short exposure variance in the pupil [µm²]")

      .def_property_readonly(
          "phase_var_avg", [](sutra_source &ss) { return ss.phase_var_avg; },
          "Long exposure variance in the pupil [µm²]")

      .def_property_readonly(
          "phase_var_count",
          [](sutra_source &ss) { return ss.phase_var_count; },
          "Counter fo long exposure variance computation")

      .def_property_readonly(
          "d_phase", [](sutra_source &ss) { return ss.d_phase->d_screen; },
          "Phase screen of the source")

      .def_property_readonly(
          "phase_telemetry",
          [](sutra_source &ss) { return ss.phase_telemetry; },
          "TODO: docstring")

      .def_property_readonly(
          "d_lgs", [](sutra_source &ss) { return ss.d_lgs; },
          "LGS structure of WFS")

      .def_property_readonly(
          "d_pupil", [](sutra_source &ss) { return ss.d_pupil; }, "Pupil mask")

      .def_property_readonly(
          "d_image_se", [](sutra_source &ss) { return ss.d_image_se; },
          "Short exposure image of the source")

      .def_property_readonly(
          "d_image_le", [](sutra_source &ss) { return ss.d_image_le; },
          "Long exposure image of the source")

      .def_property_readonly(
          "d_amplipup", [](sutra_source &ss) { return ss.d_amplipup; },
          "Complex amplitude in the pupil plane")

      .def_property_readonly(
          "d_phasepts", [](sutra_source &ss) { return ss.d_phasepts; },
          "Phase on the valid pixels of the pupil plane")

      .def_property_readonly(
          "d_wherephase", [](sutra_source &ss) { return ss.d_wherephase; },
          "Indices of the valid pixels of the pupil")

      .def_property_readonly(
          "d_ncpa_phase", [](sutra_source &ss) { return ss.d_ncpa_phase; },
          "NCPA phase")

      .def_property_readonly(
          "xoff",
          [](sutra_source &ss) -> map<type_screen, float> & { return ss.xoff; },
          "X offset for raytracing")

      .def_property_readonly(
          "yoff",
          [](sutra_source &ss) -> map<type_screen, float> & { return ss.yoff; },
          "Y offset for raytracing")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("add_layer", wy::colCast(&sutra_source::add_layer), R"pbdoc(
        Add a phase screen "dm" or "atmos" as layers to consider for raytracing

        Parameters
        ------------
        context: (carma_context) : carma or carmaWrap context
        type: (str) : "atmos" or "dm"
        xoff: (float) : x-offset for raytracing
        yoff: (float) : y-offset for raytracing
        )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("xoff"),
           py::arg("yoff"))

      .def("remove_layer", wy::colCast(&sutra_source::remove_layer), R"pbdoc(
        Remove a phase screen for raytracing

        Parameters
        ------------
        type: (str) : "atmos" or "dm"
        idx: (int) : index of the DM or turbulent layer to remove
        )pbdoc",
           py::arg("type"), py::arg("idx"))

      .def("comp_image", &sutra_source::comp_image, R"pbdoc(
        Compute short and long exposure images

        Parameters
        ------------
        puponly: (int) : Airy computation
        comp_le: (bool) : Flag for computing LE image
        )pbdoc",
           py::arg("puponly") = 0, py::arg("comp_le") = true)

      .def("init_strehlmeter", &sutra_source::init_strehlmeter,
           "Initialize Strehl ratio computation")

      .def("reset_strehlmeter", &sutra_source::reset_strehlmeter,
           "Reset Strehl ratio")

      .def("reset_phase", &sutra_source::reset_phase, "Reset the phase screen")

      .def("comp_strehl", &sutra_source::comp_strehl, "Compute Strehl ratio",
           py::arg("do_fit") = true)

      .def("raytrace", (int (sutra_source::*)(bool)) & sutra_source::raytrace,
           R"pbdoc(
        Raytrace through ncpa layers

        Parameters
        ------------
        rst: (bool): reset screen phase before raytracing
    )pbdoc",
           py::arg("rst") = false)

      .def("raytrace",
           (int (sutra_source::*)(sutra_telescope * tel, bool)) &
               sutra_source::raytrace,
           R"pbdoc(
        Raytrace through telescope aberrations

        Parameters
        ------------
        tel: (sutra_telescope): sutra_telescope object
        rst: (bool): reset screen phase before raytracing
    )pbdoc",
           py::arg("tel"), py::arg("rst") = false)

      .def("raytrace",
           (int (sutra_source::*)(sutra_atmos * atmos, bool)) &
               sutra_source::raytrace,
           R"pbdoc(
        Raytrace through turbulent layers. Calling this function will automatically reset the screen phase before raytracing.

        Parameters
        ------------
        atmos: (sutra_atmos): sutra_atmos object
        async: (bool): asynchronous mode
    )pbdoc",
           py::arg("atmos"), py::arg("async") = false)

      .def("raytrace",

           (int (sutra_source::*)(sutra_dms * dms, bool, bool, bool)) &
               sutra_source::raytrace,
           R"pbdoc(
        Raytrace through DMs

        Parameters
        ------------
        dms: (sutra_dms): sutra_dms object
        rst: (bool): reset phase screen before raytracing
        do_phase_var: (bool): compute the residual phase variance
        async: (bool): asynchronous mode
    )pbdoc",
           py::arg("dms"), py::arg("rst") = false,
           py::arg("do_phase_var") = true, py::arg("async") = false)

      .def("raytrace",

           (int (sutra_source::*)(sutra_telescope * tel, sutra_atmos * atm,
                                  sutra_dms * dms, bool, bool)) &
               sutra_source::raytrace,
           R"pbdoc(
        Raytrace through all layers (turbu, dms, telescope, ncpa)

        Parameters
        ------------
        tel: (sutra_tel): sutra_telescope object
        atm: (sutra_atmos): sutra_atmos object
        dms: (sutra_dms): sutra_dms object
        do_phase_var: (bool): compute the residual phase variance
        async: (bool): asynchronous mode
    )pbdoc",
           py::arg("dms"), py::arg("atm"), py::arg("tel"),
           py::arg("do_phase_var") = true, py::arg("async") = false)

      .def("__str__",
           [](sutra_source &ss) {
             std::cout << "Type | LGS | position(\") |  Mag | Lambda (mic.)"
                       << std::endl;
             std::cout << ss.type << " | " << ss.lgs << " | "
                       << "(" << ss.posx << "," << ss.posy << ") | " << ss.mag
                       << " | " << ss.lambda << std::endl;

             return "";
           })

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def(
          "set_ncpa",
          [](sutra_source &ss,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (data.size() == ss.d_phase->d_screen->getNbElem()) {
              if (ss.d_ncpa_phase == nullptr)
                ss.d_ncpa_phase = new carma_obj<float>(ss.d_phase->d_screen);
              ss.d_ncpa_phase->host2device(data.mutable_data());
            } else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
                      Set the NCPA phase

                      Parameters
                      ------------
                      data: (np.array(ndim=2,dtype=np.float32)): NCPA phase to set
                  )pbdoc",
          py::arg("data"))

      .def(
          "set_phase",
          [](sutra_source &ss,
             py::array_t<float, py::array::f_style | py::array::forcecast>
                 data) {
            if (data.size() == ss.d_phase->d_screen->getNbElem()) {
              ss.d_phase->d_screen->host2device(data.mutable_data());
            } else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
                      Set the target screen phase

                      Parameters
                      ------------
                      data: (np.array(ndim=2,dtype=np.float32)): target phase to set
                  )pbdoc",
          py::arg("data"))

      ;
};
