#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <carma.h>
#include "host_obj.hpp"
#include "type_list.hpp"

// #ifdef CAN_DO_HALF
// using TypeListHostObj =
//     GenericTypeList<int, float, double, half, cuFloatComplex,
//     cuDoubleComplex>;
// #else
using TypeListHostObj =
    GenericTypeList<int, float, double, cuFloatComplex, cuDoubleComplex>;
// #endif

void declare_carmaWrap_host_obj(py::module &mod) {
  py::enum_<MemAlloc>(mod, "MemAlloc")
      .value("MA_MALLOC", MemAlloc::MA_MALLOC)
      .value("MA_PAGELOCK", MemAlloc::MA_PAGELOCK)
      .value("MA_ZEROCPY", MemAlloc::MA_ZEROCPY)
      .value("MA_PORTABLE", MemAlloc::MA_PORTABLE)
      .value("MA_WRICOMB", MemAlloc::MA_WRICOMB)
      .value("MA_GENEPIN", MemAlloc::MA_GENEPIN)
      .export_values();

  apply<CarmaHostObjInterfacer, TypeListHostObj>(mod);
}
