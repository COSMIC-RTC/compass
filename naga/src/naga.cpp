#include <carma.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename T>
void declare_naga_obj(py::module &, std::string);
template <typename T>
void declare_naga_sparse_obj(py::module &, std::string);

void declare_naga_context(py::module &);
void declare_half_setter_getter(py::module &mod);

// Expose classes and methods to Python
PYBIND11_MODULE(naga, mod) {
  mod.doc() = "";
  PYBIND11_NUMPY_DTYPE(float2, x, y);
  declare_naga_context(mod);
  declare_naga_obj<int>(mod, "int");
  // declare_naga_obj<unsigned int>(mod, "uint");
  declare_naga_obj<float>(mod, "float");
  declare_naga_obj<double>(mod, "double");
  declare_naga_sparse_obj<float>(mod, "float");
  declare_naga_sparse_obj<double>(mod, "double");
  // declare_naga_obj<float2>(mod, "float2");
  // declare_naga_obj<double2>(mod, "double2");
  declare_naga_obj<cuFloatComplex>(mod, "float_complex");
  declare_naga_obj<half>(mod, "half");
  declare_half_setter_getter(mod);

  // declare_naga_obj<cuDoubleComplex>(mod, "double_complex");

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
