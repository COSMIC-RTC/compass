#include "obj.hpp"

namespace py = pybind11;

template void declare_carmaWrap_obj<cuFloatComplex>(py::module &mod,
                                               std::string suffix);
