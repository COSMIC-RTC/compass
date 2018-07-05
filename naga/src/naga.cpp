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
 
} 
