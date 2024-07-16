// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      host_obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaHostObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _WRAP_HOST_OBJ_H_
#define _WRAP_HOST_OBJ_H_

#include "declare_name.hpp"

#include <carma.hpp>

#include <type_list.hpp>

namespace py = pybind11;


struct CarmaHostObjInterfacer {
  template <typename T> static void call(py::module &mod) {
    auto name = appendName<T>("host_obj_");
    using Class = CarmaHostObj<T>;

    py::class_<Class>(mod, name.data(), py::buffer_protocol())
        .def(py::init([](const py::array_t<T, py::array::f_style |
                                                  py::array::forcecast> &data) {
               int32_t ndim = data.ndim() + 1;
               std::vector<int64_t> data_dims(ndim);
               data_dims[0] = data.ndim();
               copy(data.shape(), data.shape() + data.ndim(),
                    begin(data_dims) + 1);
               return std::unique_ptr<Class>(
                   new Class(data_dims.data(), (const T *)data.data()));
             }),
             "TODO", // TODO do the documentation...
             py::arg("h_data").none(false))

        .def(py::init([](const Class &data) {
               return std::unique_ptr<Class>(new Class(&data));
             }),
             "TODO", // TODO do the documentation...
             py::arg("d_data").none(false))

        .def_buffer([](Class &frame) -> py::buffer_info {

          const int64_t *dims = frame.get_dims();
          std::vector<ssize_t> shape(dims[0]);
          std::vector<ssize_t> strides(dims[0]);
          ssize_t stride = sizeof(T);

          // C-style
          // for (ssize_t dim(dims[0] - 1); dim >= 0; --dim)
          // {
          //   shape[dim] = dims[dim + 1];
          //   strides[dim] = stride;
          //   stride *= shape[dim];
          // }

          // F-style
          for (ssize_t dim(0); dim < dims[0]; ++dim) {
            shape[dim] = dims[dim + 1];
            strides[dim] = stride;
            stride *= shape[dim];
          }

          return py::buffer_info(frame.get_data(), sizeof(T),
                                 py::format_descriptor<T>::format(), dims[0],
                                 shape, strides);
        })

        // .def("__repr__", &std::string)

        // const int64_t *get_dims()
        .def_property_readonly("shape",
                               [](Class &frame) -> py::array_t<int64_t> {
                                 int64_t nb_dim = frame.get_dims(0);
                                 const int64_t *c_dim = frame.get_dims() + 1;
                                 return py::array_t<int64_t>(nb_dim, c_dim);
                               },
                               "TODO") // TODO do the documentation...

        // int32_t get_nb_elements()
        .def_property_readonly("nbElem", &Class::get_nb_elements,
                               "TODO") // TODO do the documentation...
    ;

  }
};
#endif
