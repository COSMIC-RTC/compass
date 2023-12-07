// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sparse_obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaSparseObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _WRAP_SPARSE_OBJ_H_
#define _WRAP_SPARSE_OBJ_H_

#include "declare_name.hpp"

#include <carma.h>
#include <carma_sparse_obj.h>

#include <type_list.hpp>

namespace py = pybind11;

struct CarmaSparseObjInterfacer {
  template <typename T> static void call(py::module &mod) {
    auto name = appendName<T>("sparse_obj_");
    using Class = CarmaSparseObj<T>;

    py::class_<Class>(mod, name.data(), py::buffer_protocol())
    .def("get_csr",[](Class &frame){
        py::object CSR = py::module::import("scipy.sparse").attr("csr_matrix");
        int32_t dim1 = frame.get_dims(1);
        int32_t dim2 = frame.get_dims(2);
        int32_t nnz = frame.nz_elem;

        std::vector<int32_t> rowind = std::vector<int32_t>(dim1 + 1);
        std::vector<int32_t> colind = std::vector<int32_t>(nnz);
        std::vector<T> data = std::vector<T>(nnz);

        frame.sparse_to_host(rowind.data(), colind.data(), data.data());
        py::tuple csrMat = py::make_tuple(py::array(data.size(),data.data()), py::array(colind.size(),colind.data()),py::array(rowind.size(), rowind.data()));
        py::tuple shape = py::make_tuple(dim1,dim2);


        return CSR(csrMat,py::arg("shape")=shape);
        });
}
};
#endif