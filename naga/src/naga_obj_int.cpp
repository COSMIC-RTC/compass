#include "naga_obj.hpp"

namespace py = pybind11;

template void declare_naga_obj<int>(py::module &mod, std::string suffix);
