#include <wyrm>

#include <carma.h>

namespace py = pybind11;

void declare_naga_context(py::module &mod) {
  py::class_<carma_device>(mod, "naga_device")
      .def_property_readonly("id", &carma_device::get_properties)
      .def_property_readonly("compute_perf", &carma_device::get_compute_perf)
      .def_property_readonly("cores_per_sm", &carma_device::get_cores_per_sm)
      .def_property_readonly("p2p_activate", &carma_device::isP2P_active)
      .def_property_readonly("name", &carma_device::getName)
      .def_property_readonly("totalMem", &carma_device::getTotalMem)
      .def_property_readonly("freeMem", &carma_device::getFreeMem);

  py::class_<carma_context>(mod, "naga_context")
      .def_property_readonly("ndevices", &carma_context::get_ndevice)

      .def_static("get_instance", &carma_context::instance,
                  py::return_value_policy::reference)
      .def_static("get_instance_1gpu", &carma_context::instance_1gpu,
                  py::return_value_policy::reference)
      .def_static("get_instance_ngpu",
                  wy::colCast(carma_context::instance_ngpu),
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

      .def_property_readonly("ndevice", &carma_context::get_ndevice)
      .def_property_readonly("activeDevice", &carma_context::get_activeDevice)
      .def_property_readonly("activeRealDevice",
                             &carma_context::get_activeRealDevice)
      .def_property_readonly("cudaRuntimeGetVersion",
                             &carma_context::get_cudaRuntimeGetVersion)
      .def_property_readonly("driverVersion",
                             &carma_context::get_cudaDriverGetVersion)
      .def("get_device", &carma_context::get_device)
      .def("set_activeDevice",
           [](carma_context &cc, int newDevice) {
             return cc._set_activeDevice(newDevice, 1, __FILE__, __LINE__);
           })
      .def("set_activeDeviceForce",
           [](carma_context &cc, int newDevice) {
             return cc._set_activeDeviceForce(newDevice, 1, __FILE__, __LINE__);
           })
      // .def("set_activeDeviceForCpy", &carma_context::set_activeDeviceForCpy);
      ;

  mod.def("deviceSync", &__carmaSafeDeviceSynchronize);
}
