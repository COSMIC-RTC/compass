#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <carma.h>
#include <carma_sparse_obj.h>

#include "sparse_obj.hpp"
#include "type_list.hpp"

using TypeListSparseObj = GenericTypeList<float, double>;

void declare_carmaWrap_sparse_obj(py::module &mod) {
  apply<CarmaSparseObjInterfacer, TypeListSparseObj>(mod);
}
