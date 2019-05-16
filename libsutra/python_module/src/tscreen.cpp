#include <wyrm>

#include <sutra_tscreen.h>
#include "declare_name.hpp"

namespace py = pybind11;

void declare_tscreen(py::module &mod) {
  py::class_<sutra_tscreen>(mod, "Tscreen")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](sutra_tscreen &st) { return st.device; }, "Device index")

      .def_property_readonly(
          "d_A", [](sutra_tscreen &st) { return st.d_A; },
          "A matrix for extrusion")

      .def_property_readonly(
          "d_B", [](sutra_tscreen &st) { return st.d_B; },
          "B matrix for extrusion")

      .def_property_readonly(
          "d_istencilx", [](sutra_tscreen &st) { return st.d_istencilx; },
          "stencil for row extrusion")

      .def_property_readonly(
          "d_istencily", [](sutra_tscreen &st) { return st.d_istencily; },
          "stencil for column extrusion")

      .def_property_readonly(
          "d_z", [](sutra_tscreen &st) { return st.d_z; },
          "tmp array for extrusion process")

      .def_property_readonly(
          "d_noise", [](sutra_tscreen &st) { return st.d_noise; },
          "random numbers for extrusion")

      .def_property_readonly(
          "d_ytmp", [](sutra_tscreen &st) { return st.d_ytmp; },
          "contains the extrude update")

      .def_property_readonly(
          "screen_size", [](sutra_tscreen &st) { return st.screen_size; },
          "size of phase screen")

      .def_property_readonly(
          "r0", [](sutra_tscreen &st) { return st.r0; }, "layer r0 in pixels")

      .def_property_readonly(
          "amplitude", [](sutra_tscreen &st) { return st.amplitude; },
          "amplitude for extrusion (r0**(-5/6)")

      .def_property_readonly(
          "altitude", [](sutra_tscreen &st) { return st.altitude; },
          "altitude of the phase screen")

      .def_property_readonly(
          "windspeed", [](sutra_tscreen &st) { return st.windspeed; },
          "wind speed of phase screen")

      .def_property_readonly(
          "winddir", [](sutra_tscreen &st) { return st.winddir; },
          "wind direction of phase screen")

      .def_property_readonly(
          "deltax", [](sutra_tscreen &st) { return st.deltax; },
          "number of columns to extrude per iteration")

      .def_property_readonly(
          "deltay", [](sutra_tscreen &st) { return st.deltay; },
          "number of rows to extrude per iteration")

      .def_property_readonly(
          "accumx", [](sutra_tscreen &st) { return st.accumx; },
          "accumulate columns to extrude")

      .def_property_readonly(
          "accumy", [](sutra_tscreen &st) { return st.accumy; },
          "accumulate rows to extrude")

      .def_property_readonly(
          "d_screen", [](sutra_tscreen &st) { return st.d_tscreen->d_screen; },
          "Turbulent phase screen");
};
