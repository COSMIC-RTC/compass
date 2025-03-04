// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      controller.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraController
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "sutraUtils.hpp"

#include <sutra_controller.hpp>

namespace py = pybind11;

template <typename Tcomp, typename Tout>
int32_t 
add_perturb_voltage(SutraController<Tcomp, Tout> &sc, string name, ArrayFStyle<Tcomp> &perturb, int32_t N)
{
    return sc.add_perturb_voltage(name, perturb.mutable_data(), N);
}

template <typename Tcomp, typename Tout>
int32_t 
set_perturb_voltage(SutraController<Tcomp, Tout> &sc, string name, ArrayFStyle<Tcomp> &perturb, int32_t N)
{
    return sc.set_perturb_voltage(name, perturb.mutable_data(), N);
}

template <typename Tcomp, typename Tout>
int32_t 
set_com(SutraController<Tcomp, Tout> &sc, ArrayFStyle<Tcomp> &perturb, int32_t N)
{
    return sc.set_com(perturb.mutable_data(), N);
}


template <typename Tcomp, typename Tout>
void controller_impl(py::module &mod, const char *name)
{
    using controller = SutraController<Tcomp, Tout>;

    py::class_<controller>(mod, name)

        //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
        //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
        //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
        //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
        //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
        //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝

        .def_property_readonly(
            "context", [](controller &sc)
            { return sc.current_context; },
            "GPU context")

        .def_property_readonly(
            "device", [](controller &sc)
            { return sc.device; },
            "GPU device index")

        .def_property_readonly(
            "type", [](controller &sc)
            { return sc.get_type(); },
            "Controller type")

        .def_property_readonly(
            "nactu", [](controller &sc)
            { return sc.nactu(); },
            "Number of actuators to control")

        .def_property_readonly(
            "nslope", [](controller &sc)
            { return sc.nslope(); },
            "Number of slopes")

        .def_property_readonly(
            "gain", [](controller &sc)
            { return sc.gain; },
            "Controller gain")
        .def_property_readonly(
            "open_loop", [](controller &sc)
            { return sc.open_loop; },
            "Open loop flag")

        .def_property_readonly(
            "delay", [](controller &sc)
            { return sc.delay; },
            "Loop delay")

        .def_property_readonly(
            "d_dmseen", [](controller &sc)
            { return sc.d_dmseen; },
            "Vector of SutraDm commanded")

        .def_property_readonly(
            "centro_idx", [](controller &sc)
            { return sc.centro_idx; },
            "Indices of the handled centroiders")

        .def_property_readonly(
            "d_centroids", [](controller &sc)
            { return sc.d_centroids; },
            "Slopes vector")

        .def_property_readonly(
            "d_com", [](controller &sc)
            { return sc.d_com; },
            "Current command vector")

        .def_property_readonly(
            "d_com1", [](controller &sc)
            { return sc.d_com1; },
            "Command vector at iteration k-1")

        .def_property_readonly(
            "d_circularComs0",
            [](controller &sc)
            { return sc.d_circular_coms[0]; },
            "Oldest command vector in the circular buffer")

        .def_property_readonly(
            "d_circularComs1",
            [](controller &sc)
            { return sc.d_circular_coms[1]; },
            "Second oldest Command vector in the circular buffer")

        .def_property_readonly(
            "comRange",
            [](controller &sc)
            {
                return std::make_tuple(sc.volt_min, sc.volt_max);
            },
            "Tuple (volt_min, volt_max) used for command clipping")
        .def_property_readonly(
            "val_max", [](controller &sc)
            { return sc.val_max; },
            "Maximum value for d_voltage (ADU). Only used if "
            "output is expected in uint16")

        .def_property_readonly(
            "d_perturb_map",
            [](controller &sc)
                -> map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>> &
            {
                return sc.d_perturb_map;
            },
            "Perturbation voltage buffers")

        .def_property_readonly(
            "d_voltage", [](controller &sc)
            { return sc.d_voltage; },
            "Total voltage to apply on the DMs")

        .def_property_readonly(
            "d_com_clipped", [](controller &sc)
            { return sc.d_com_clipped; },
            "Delayed commands")

        //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
        //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
        //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
        //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
        //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
        //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

        .def("add_perturb", &controller::add_perturb,
             "Add the perturbation voltage to the command")

        .def("command_delay", &controller::command_delay, "Delay the command")

        .def("reset_coms", &controller::reset_coms,
             "Reset the commands circular buffer")

        .def("enable_perturb_voltage",
             &controller::enable_perturb_voltage,
             R"pbdoc(
    Enable a perturbation voltage buffer

    Args:
        name: (str): name of the buffer to enable
     )pbdoc",
             py::arg("name"))

        .def("disable_perturb_voltage",
             &controller::disable_perturb_voltage,
             R"pbdoc(
    Disable a perturbation voltage buffer

    Args:
        name: (str): name of the buffer to enable
     )pbdoc",
             py::arg("name"))

        .def("comp_voltage", &controller::comp_voltage,
             "Computes the final voltage to send to the DM")

        .def("clip_commands", &controller::clip_commands,
             "Clip the commands between volt_min and volt_max (values set in the "
             "controller)")

        //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
        //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
        //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
        //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
        //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
        //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝

        .def("add_perturb_voltage", &add_perturb_voltage<Tcomp, Tout>,
             R"pbdoc(
    Add a new perturbation voltage buffer

    Args:
        name: (str): name of the new buffer

        perturb: (np.array[ndim=2,dtype=np.float32]): perturbation voltage to set

        N: (int): Number of perturb voltage vectors in the buffer
     )pbdoc",
             py::arg("name"), py::arg("perturb"), py::arg("N"))

        .def("set_perturb_voltage", &set_perturb_voltage<Tcomp, Tout>,
             R"pbdoc(
    Set an existing perturbation voltage buffer

    Args:
        name: (str): name of the buffer

        perturb: (np.array[ndim=2,dtype=np.float32]): perturbation voltage to set

        N: (int): Number of perturb voltage vectors in the buffer
     )pbdoc",
             py::arg("name"), py::arg("perturb"), py::arg("N"))

        .def("remove_perturb_voltage",
             &controller::remove_perturb_voltage,
             R"pbdoc(
    Remove a perturbation voltage buffer

    Args:
        name: (str): name of the buffer to remove
     )pbdoc",
             py::arg("name"))

        .def("reset_perturb_voltage", &controller::reset_perturb_voltage,
             R"pbdoc(
    Remove all perturbation voltage buffers
     )pbdoc")

        .def("set_open_loop", &controller::set_open_loop,
             R"pbdoc(
    Open (1) or close (0) the loop

    Args:
        status: (int): open loop status

        rst: (bool): reset integrator if True
     )pbdoc",
             py::arg("status"), py::arg("rst") = true)

        .def("set_delay", &controller::set_delay,
             R"pbdoc(
    Set the delay

    Args:
        delay: (float): loop delay in frames
     )pbdoc",
             py::arg("delay"))

        .def("set_gain", &controller::set_gain,
             R"pbdoc(
    Set the gain

    Args:
        gain: (float): loop gain
     )pbdoc",
             py::arg("gain"))

        .def(
            "set_comRange",
            [](controller &sc, float volt_min, float volt_max)
            {
                sc.set_volt_max(volt_max);
                sc.set_volt_min(volt_min);
            },
            R"pbdoc(
    Set the volt_min and volt_max value for command clipping

    Args:
        volt_min: (float): volt_min value for clipping

        volt_max: (float): volt_max value for clipping
     )pbdoc",
            py::arg("volt_min"), py::arg("volt_max"))

        .def("set_volt_max", &controller::set_volt_max,
             R"pbdoc(
    Set the volt_max value for command clipping

    Args:
        volt_max: (float): volt_max value for clipping
     )pbdoc",
             py::arg("volt_max"))

        .def("set_val_max", &controller::set_val_max,
             R"pbdoc(
    Set the val_max value for command conversion

    Args:
        val_max: (float): val_max value for conversion
     )pbdoc",
             py::arg("val_max"))

        .def("set_com", &set_com<Tcomp, Tout>, R"pbdoc(
    Set the command vector of the controller

    Args:
        com: (np.array[ndim=3, dtype=np.float32]) : command vector

        nElem: (int): Number of elements in com
        )pbdoc",
             py::arg("com"), py::arg("nElem"));
}

void declare_controller(py::module &mod)
{
    controller_impl<float, float>(mod, "Controller_FF");
    controller_impl<float, uint16_t>(mod, "Controller_FU");
}
