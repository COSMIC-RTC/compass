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

//! \file      context.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaContext
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.hpp>

namespace py = pybind11;

void declare_carma_context(py::module &mod) {
  py::class_<CarmaDevice>(mod, "device")
      .def_property_readonly("id", &CarmaDevice::get_properties)
      .def_property_readonly("compute_perf", &CarmaDevice::get_compute_perf)
      .def_property_readonly("cores_per_sm", &CarmaDevice::get_cores_per_sm)
      .def_property_readonly("name", &CarmaDevice::get_name)
      .def_property_readonly("total_mem", &CarmaDevice::get_total_mem)
      .def_property_readonly("total_mem", &CarmaDevice::get_free_mem);

  py::class_<CarmaContext>(mod, "context")
      .def_property_readonly("ndevices", &CarmaContext::get_ndevice)

      .def_static("get_instance", &CarmaContext::instance,
                  py::return_value_policy::reference)
      .def_static("get_instance_1gpu", &CarmaContext::instance_1gpu,
                  py::return_value_policy::reference)
      .def_static("get_instance_ngpu",
                  [](int32_t nb_devices, const py::array_t<int32_t> &devices_id) -> CarmaContext & {
                      return CarmaContext::instance_ngpu(nb_devices, devices_id.data());
                  },
                  py::return_value_policy::reference)

      .def_property_readonly("ndevice", &CarmaContext::get_ndevice)
      .def_property_readonly("active_device", &CarmaContext::get_active_device)
      .def_property_readonly("activeRealDevice",
                             &CarmaContext::get_active_real_device)
      .def_property_readonly("cudaRuntimeGetVersion",
                             &CarmaContext::get_cuda_runtime_get_version)
      .def_property_readonly("driverVersion",
                             &CarmaContext::get_cuda_driver_get_version)
      .def("get_device", &CarmaContext::get_device,
           py::return_value_policy::reference)
      .def("set_active_device",
           [](CarmaContext &cc, int32_t new_device) {
             return cc._set_active_device(new_device, 1, __FILE__, __LINE__);
           })
      .def("set_active_device_force",
           [](CarmaContext &cc, int32_t new_device) {
             return cc.set_active_device_force(new_device, 0);
           })
      // .def("set_active_deviceForCpy", &CarmaContext::set_active_deviceForCpy);
      .def(
          "activate_tensor_cores",
          [](CarmaContext &cc, bool flag) {
            int32_t ndevices = cc.get_ndevice();
            for (int32_t i = 0; i < ndevices; i++) {
              cc.get_device(i)->set_cublas_math_mode(flag);
            }
          },
          "Set the cublas math mode using tensor cores or not",
          py::arg("flag"));

  mod.def("deviceSync", &__carma_safe_device_synchronize);
}
