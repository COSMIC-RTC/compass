#include <carma.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void declare_carmaWrap_obj(py::module &);
void declare_carmaWrap_sparse_obj(py::module &);
void declare_carmaWrap_context(py::module &);
void declare_carmaWrap_timer(py::module &);

#ifdef CAN_DO_HALF
void declare_half_setter_getter(py::module &mod);
#endif

// Expose classes and methods to Python
PYBIND11_MODULE(carmaWrap, mod) {
  mod.doc() = "";
  PYBIND11_NUMPY_DTYPE(float2, x, y);
  PYBIND11_NUMPY_DTYPE(double2, x, y);
  declare_carmaWrap_context(mod);
  declare_carmaWrap_obj(mod);
  declare_carmaWrap_sparse_obj(mod);
#ifdef CAN_DO_HALF
  declare_half_setter_getter(mod);
#endif
  declare_carmaWrap_timer(mod);

  // declare_carmaWrap_obj<cuDoubleComplex>(mod, "double_complex");

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
