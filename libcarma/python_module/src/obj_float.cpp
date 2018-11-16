#include "obj.hpp"

namespace py = pybind11;

template void declare_carmaWrap_obj<float>(py::module &mod, std::string suffix);
