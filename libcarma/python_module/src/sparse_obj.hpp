// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sparse_obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaSparseObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#ifndef _WRAP_SPARSE_OBJ_H_
#define _WRAP_SPARSE_OBJ_H_

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

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
        int dim1 = frame.get_dims(1);
        int dim2 = frame.get_dims(2);
        int nnz = frame.nz_elem;

        std::vector<int> rowind = std::vector<int>(dim1 + 1);
        std::vector<int> colind = std::vector<int>(nnz);
        std::vector<T> data = std::vector<T>(nnz);

        frame.sparse_to_host(rowind.data(), colind.data(), data.data());
        py::tuple csrMat = py::make_tuple(py::array(data.size(),data.data()), py::array(colind.size(),colind.data()),py::array(rowind.size(), rowind.data()));
        py::tuple shape = py::make_tuple(dim1,dim2);


        return CSR(csrMat,py::arg("shape")=shape);
        });
}
};
#endif
