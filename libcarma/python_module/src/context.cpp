// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      context.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaContext
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>
namespace py = pybind11;

void declare_carmaWrap_context(py::module &mod) {
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
                  wy::colCast(CarmaContext::instance_ngpu),
                  py::return_value_policy::reference)

      // .def(py::init([](py::buffer data, ::CORBA::Boolean copy) {
      //   py::buffer_info info = data.request();
      //   if(info.format != py::format_descriptor<::CORBA::ULong>::format())
      //     throw invalid_argument("Buffer given has an incorrect format,
      //     expected: uint32");
      //   ::CORBA::ULong *ptr;
      //   if (copy == 0) {
      //     ptr = reinterpret_cast<::CORBA::ULong *>(info.ptr);
      //   } else {
      //     ssize_t size = info.itemsize * info.size;
      //     ptr = new ::CORBA::ULong[info.size];
      //     memcpy(ptr, info.ptr, size);
      //   }
      //   return unique_ptr<B::Dims>(new B::Dims(info.size, info.size, ptr,
      //   copy));
      // }))

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
           [](CarmaContext &cc, int new_device) {
             return cc._set_active_device(new_device, 1, __FILE__, __LINE__);
           })
      .def("set_active_device_force",
           [](CarmaContext &cc, int new_device) {
             return cc.set_active_device_force(new_device, 0);
           })
      // .def("set_active_deviceForCpy", &CarmaContext::set_active_deviceForCpy);
      .def(
          "activate_tensor_cores",
          [](CarmaContext &cc, bool flag) {
            int ndevices = cc.get_ndevice();
            for (int i = 0; i < ndevices; i++) {
              cc.get_device(i)->set_cublas_math_mode(flag);
            }
          },
          "Set the cublas math mode using tensor cores or not",
          py::arg("flag"));

  mod.def("deviceSync", &__carma_safe_device_synchronize);
}
