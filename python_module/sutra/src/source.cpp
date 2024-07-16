// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      source.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraSource
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_source.hpp>

namespace py = pybind11;

void declare_source(py::module &mod) {
  py::class_<SutraSource>(mod, "Source")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

      .def_property_readonly(
          "device", [](SutraSource &ss) { return ss.device; },
          "GPU device index")

      .def_property_readonly(
          "posx", [](SutraSource &ss) { return ss.posx; },
          "X position of the source")

      .def_property_readonly(
          "posy", [](SutraSource &ss) { return ss.posy; },
          "Y position of the source")

      .def_property_readonly(
          "npts", [](SutraSource &ss) { return ss.npts; },
          "Number of points in the pupil")

      .def_property_readonly(
          "mag", [](SutraSource &ss) { return ss.mag; },
          "Magnitude of the source")

      .def_property_readonly(
          "lambda_um", [](SutraSource &ss) { return ss.lambda; },
          "Wavelength of the source in µm")

      .def_property_readonly(
          "zp", [](SutraSource &ss) { return ss.zp; }, "Flux at magnitude 0")

      .def_property_readonly(
          "scale", [](SutraSource &ss) { return ss.scale; },
          "Phase scale factor (2*pi/lambda)")

      .def_property_readonly(
          "lgs", [](SutraSource &ss) { return ss.lgs; }, "Boolean for LGS")

      .def_property_readonly(
          "G", [](SutraSource &ss) { return ss.G; },
          "Magnifying factor for WFS misalignment")

      .def_property_readonly(
          "thetaML", [](SutraSource &ss) { return ss.thetaML; },
          "Pupil rotation angle for WFS misalignment")

      .def_property_readonly(
          "dx", [](SutraSource &ss) { return ss.dx; },
          "X axis WFS misalignment [pixels]")

      .def_property_readonly(
          "dy", [](SutraSource &ss) { return ss.dy; },
          "Y axis WFS misalignment [pixels]")

      .def_property_readonly(
          "type", [](SutraSource &ss) { return ss.type; },
          "Type of source (Target or WFS GS)")

      .def_property_readonly(
          "block_size", [](SutraSource &ss) { return ss.block_size; },
          "Optimum block size of device")

      .def_property_readonly(
          "strehl_se", [](SutraSource &ss) { return ss.strehl_se; },
          "Short exposure Strehl ratio")

      .def_property_readonly(
          "strehl_le", [](SutraSource &ss) { return ss.strehl_le; },
          "Long exposure Strehl ratio")

      .def_property_readonly(
          "ref_strehl", [](SutraSource &ss) { return ss.ref_strehl; },
          "reference for Strehl computation (Airy)")

      .def_property_readonly(
          "strehl_counter", [](SutraSource &ss) { return ss.strehl_counter; },
          "Counter for LE Strehl computation")

      .def_property_readonly(
          "phase_var", [](SutraSource &ss) { return ss.phase_var; },
          "Short exposure variance in the pupil [µm²]")

      .def_property_readonly(
          "phase_var_avg", [](SutraSource &ss) { return ss.phase_var_avg; },
          "Long exposure variance in the pupil [µm²]")

      .def_property_readonly(
          "phase_var_count", [](SutraSource &ss) { return ss.phase_var_count; },
          "Counter fo int64_t exposure variance computation")

      .def_property_readonly(
          "d_phase", [](SutraSource &ss) { return ss.d_phase->d_screen; },
          "Phase screen of the source")

      .def_property_readonly(
          "phase_telemetry", [](SutraSource &ss) { return ss.phase_telemetry; },
          "TODO: docstring")

      .def_property_readonly(
          "d_lgs", [](SutraSource &ss) { return ss.d_lgs; },
          "LGS structure of WFS")

      .def_property_readonly(
          "d_pupil", [](SutraSource &ss) { return ss.d_pupil; }, "Pupil mask")

      .def_property_readonly(
          "d_image_se", [](SutraSource &ss) { return ss.d_image_se; },
          "Short exposure image of the source")

      .def_property_readonly(
          "d_image_le", [](SutraSource &ss) { return ss.d_image_le; },
          "Long exposure image of the source")

      .def_property_readonly(
          "d_amplipup", [](SutraSource &ss) { return ss.d_amplipup; },
          "Complex amplitude in the pupil plane")

      .def_property_readonly(
          "d_phasepts", [](SutraSource &ss) { return ss.d_phasepts; },
          "Phase on the valid pixels of the pupil plane")

      .def_property_readonly(
          "d_wherephase", [](SutraSource &ss) { return ss.d_wherephase; },
          "Indices of the valid pixels of the pupil")

      .def_property_readonly(
          "d_ncpa_phase", [](SutraSource &ss) { return ss.d_ncpa_phase; },
          "NCPA phase")

      .def_property_readonly(
          "xoff",
          [](SutraSource &ss) -> map<type_screen, float> & { return ss.xoff; },
          "X offset for raytracing")

      .def_property_readonly(
          "yoff",
          [](SutraSource &ss) -> map<type_screen, float> & { return ss.yoff; },
          "Y offset for raytracing")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("add_layer", &SutraSource::add_layer, R"pbdoc(
    Add a phase screen "dm" or "atmos" as layers to consider for raytracing

    Args:
        context: (CarmaContext) : carma or carma context

        type: (str) : "atmos" or "dm"

        xoff: (float) : x-offset for raytracing

        yoff: (float) : y-offset for raytracing
        )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("xoff"),
           py::arg("yoff"))

      .def("remove_layer", &SutraSource::remove_layer, R"pbdoc(
    Remove a phase screen for raytracing

    Args:
        type: (str) : "atmos" or "dm"

        idx: (int32_t) : index of the DM or turbulent layer to remove
        )pbdoc",
           py::arg("type"), py::arg("idx"))

      .def("comp_image", &SutraSource::comp_image, R"pbdoc(
    Compute short and int64_t exposure images

    Args:
        puponly: (int32_t) : Airy computation

        comp_le: (bool) : Flag for computing LE image
        )pbdoc",
           py::arg("puponly") = 0, py::arg("comp_le") = true)

      .def("init_strehlmeter", &SutraSource::init_strehlmeter,
           "Initialize Strehl ratio computation")

      .def("reset_strehlmeter", &SutraSource::reset_strehlmeter,
           "Reset Strehl ratio")

      .def("reset_phase", &SutraSource::reset_phase, "Reset the phase screen")

      .def("comp_strehl", &SutraSource::comp_strehl, "Compute Strehl ratio",
           py::arg("do_fit") = true)

      .def("raytrace", (int32_t (SutraSource::*)(bool)) & SutraSource::raytrace,
           R"pbdoc(
    Raytrace through ncpa layers

    Args:
        rst: (bool): reset screen phase before raytracing
    )pbdoc",
           py::arg("rst") = false)

      .def("raytrace",
           (int32_t (SutraSource::*)(SutraTelescope * tel, bool)) &
               SutraSource::raytrace,
           R"pbdoc(
    Raytrace through telescope aberrations

    Args:
        tel: (SutraTelescope): SutraTelescope object

        rst: (bool): reset screen phase before raytracing
    )pbdoc",
           py::arg("tel"), py::arg("rst") = false)

      .def("raytrace",
           (int32_t (SutraSource::*)(SutraAtmos * atmos, bool)) &
               SutraSource::raytrace,
           R"pbdoc(
    Raytrace through turbulent layers. Calling this function will automatically reset the screen phase before raytracing.

    Args:
        atmos: (SutraAtmos): SutraAtmos object

        do_async: (bool): asynchronous mode
    )pbdoc",
           py::arg("atmos"), py::arg("do_async") = false)

      .def("raytrace",

           (int32_t (SutraSource::*)(SutraDms * dms, bool, bool, bool)) &
               SutraSource::raytrace,
           R"pbdoc(
    Raytrace through DMs

    Args:
        dms: (SutraDms): SutraDms object

        rst: (bool): reset phase screen before raytracing

        do_phase_var: (bool): compute the residual phase variance

        do_async: (bool): asynchronous mode
    )pbdoc",
           py::arg("dms"), py::arg("rst") = false,
           py::arg("do_phase_var") = true, py::arg("do_async") = false)

      .def("raytrace",

           (int32_t (SutraSource::*)(SutraTelescope * tel, SutraAtmos * atm,
                                 SutraDms * dms, bool, bool)) &
               SutraSource::raytrace,
           R"pbdoc(
    Raytrace through all layers (turbu, dms, telescope, ncpa)

    Args:
        tel: (sutra_tel): SutraTelescope object

        atm: (SutraAtmos): SutraAtmos object

        dms: (SutraDms): SutraDms object

        do_phase_var: (bool): compute the residual phase variance

        do_async: (bool): asynchronous mode
    )pbdoc",
           py::arg("dms"), py::arg("atm"), py::arg("tel"),
           py::arg("do_phase_var") = true, py::arg("do_async") = false)

      .def("__str__",
           [](SutraSource &ss) {
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

      .def(
          "set_ncpa",
          [](SutraSource &ss,
             ArrayFStyle<float>
                 data) {
            if (data.size() == ss.d_phase->d_screen->get_nb_elements()) {
              if (ss.d_ncpa_phase == nullptr)
                ss.d_ncpa_phase = new CarmaObj<float>(ss.d_phase->d_screen);
              ss.d_ncpa_phase->host2device(data.mutable_data());
            } else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the NCPA phase

    Args:
        data: (np.array(ndim=2,dtype=np.float32)): NCPA phase to set
                  )pbdoc",
          py::arg("data"))

      .def(
          "set_phase",
          [](SutraSource &ss,
             ArrayFStyle<float>
                 data) {
            if (data.size() == ss.d_phase->d_screen->get_nb_elements()) {
              ss.d_phase->d_screen->host2device(data.mutable_data());
            } else
              DEBUG_TRACE("Wrong dimensions");
          },
          R"pbdoc(
    Set the target screen phase

    Args:
        data: (np.array(ndim=2,dtype=np.float32)): target phase to set
    )pbdoc",
          py::arg("data"))

      ;
};
