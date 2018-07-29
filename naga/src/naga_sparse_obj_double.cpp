#include "naga_sparse_obj.hpp"

namespace py = pybind11;

template void declare_naga_sparse_obj<double>(py::module &mod,
                                              std::string suffix);
