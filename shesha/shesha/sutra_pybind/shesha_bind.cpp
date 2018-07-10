#include <pybind11/pybind11.h>

namespace py = pybind11;

void declare_shesha_atmos(py::module &);

// Expose classes and methods to Python
PYBIND11_MODULE(shesha_bind, mod) {
  mod.doc() = "Binding module for libsutra into shesha";

  declare_shesha_atmos(mod);
}
