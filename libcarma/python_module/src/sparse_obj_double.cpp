#include "sparse_obj.hpp"

namespace py = pybind11;

template void declare_carmaWrap_sparse_obj<double>(py::module &mod,
                                              std::string suffix);
