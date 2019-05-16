#include <carma.h>
#include <carma_sparse_obj.h>

#include "sparse_obj.hpp"

using TypeListSparseObj = GenericTypeList<float, double>;

void declare_carmaWrap_sparse_obj(py::module &mod) {
  apply<CarmaSparseObjInterfacer, TypeListSparseObj>(mod);
}
