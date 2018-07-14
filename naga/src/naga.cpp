#include <carma.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename T>
void declare_naga_obj(py::module &, std::string);

void declare_naga_context(py::module &);

// Expose classes and methods to Python
PYBIND11_MODULE(naga, mod) {
  mod.doc() = "";

  declare_naga_context(mod);
  // declare_naga_obj<int>(mod, "int");
  // declare_naga_obj<unsigned int>(mod, "uint");
  declare_naga_obj<float>(mod, "float");
  declare_naga_obj<double>(mod, "double");
  // declare_naga_obj<float2>(mod, "float2");
  // declare_naga_obj<double2>(mod, "double2");
  // declare_naga_obj<cuFloatComplex>(mod, "float_complex");
  // declare_naga_obj<cuDoubleComplex>(mod, "double_complex");
}
