#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

namespace py = pybind11;

void declare_carmaWrap_host_obj(py::module &);
void declare_carmaWrap_obj(py::module &);
void declare_carmaWrap_sparse_obj(py::module &);
void declare_carmaWrap_context(py::module &);
void declare_carmaWrap_timer(py::module &);

// Expose classes and methods to Python
PYBIND11_MODULE(carmaWrap, mod) {
  mod.doc() = "";
  declare_carmaWrap_context(mod);
  declare_carmaWrap_obj(mod);
  declare_carmaWrap_host_obj(mod);
  declare_carmaWrap_sparse_obj(mod);
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
