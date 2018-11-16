#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <carma.h>
#include <carma_sparse_obj.h>

namespace py = pybind11;

template <typename T>
void declare_carmaWrap_sparse_obj(py::module &mod, std::string suffix)
{

  using Class = carma_sparse_obj<T>;

  py::class_<Class>(mod, ("sparse_obj_" + suffix).c_str())
    .def("get_csr",[](Class &frame){
        py::object CSR = py::module::import("scipy.sparse").attr("csr_matrix");
        int dim1 = frame.getDims(1);
        int dim2 = frame.getDims(2);
        int nnz = frame.nz_elem;

        std::vector<int> rowind = std::vector<int>(dim1 + 1);
        std::vector<int> colind = std::vector<int>(nnz);
        std::vector<T> data = std::vector<T>(nnz);

        frame.sparse_to_host(rowind.data(), colind.data(), data.data());
        py::tuple csrMat = py::make_tuple(py::array(data.size(),data.data()), py::array(colind.size(),colind.data()),py::array(rowind.size(), rowind.data()));
        py::tuple shape = py::make_tuple(dim1,dim2);


        return CSR(csrMat,py::arg("shape")=shape);
        })
;
};
