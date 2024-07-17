// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sparse_obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaSparseObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _WRAP_SPARSE_OBJ_H_
#define _WRAP_SPARSE_OBJ_H_

#include "declare_name.hpp"

#include <carma.hpp>
#include <carma_sparse_obj.hpp>

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
