// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      obj_half.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for half float CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _WRAP_OBJ_HALF_H_
#define _WRAP_OBJ_HALF_H_

#include "declare_name.hpp"

#include <carma.hpp>

#include <type_list.hpp>

namespace py = pybind11;

struct CarmaObjHalfInterfacer {
  template <typename T>
  static void call(py::module &mod) {
    auto name = appendName<T>("obj_");
    using Class = CarmaObj<T>;
    using ClassHost = CarmaHostObj<T>;

    py::class_<Class> carmaWrapObj(mod, name.data(), py::buffer_protocol());
    carmaWrapObj.def(
        py::init([](CarmaContext &c,
                    const py::array_t<T, py::array::f_style |
                                             py::array::forcecast> &data) {
          int32_t ndim = data.ndim() + 1;
          std::vector<int64_t> data_dims(ndim);
          data_dims[0] = data.ndim();
          copy(data.shape(), data.shape() + data.ndim(), begin(data_dims) + 1);
          return std::unique_ptr<Class>(
              new Class(&c, data_dims.data(), (const T *)data.data()));
        }),
        "TODO",  // TODO do the documentation...
        py::arg("context").none(false), py::arg("h_data").none(false));

    carmaWrapObj.def(py::init([](CarmaContext &c, const Class &data) {
                       return std::unique_ptr<Class>(new Class(&c, &data));
                     }),
                     "TODO",  // TODO do the documentation...
                     py::arg("context").none(false),
                     py::arg("d_data").none(false));

    carmaWrapObj.def_buffer([](Class &frame) -> py::buffer_info {
      frame.sync_h_data();

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

      return py::buffer_info(frame.get_h_data(), sizeof(T),
                             py::format_descriptor<T>::format(), dims[0], shape,
                             strides);
    });

    carmaWrapObj.def("__repr__", &Class::to_string);

    // int32_t get_nb_streams()
    carmaWrapObj.def_property_readonly("nb_streams", &Class::get_nb_streams,
                                       "TODO");  // TODO do the documentation...
    // int32_t add_stream()
    carmaWrapObj.def("add_stream", (int32_t (Class::*)()) & Class::add_stream,
                     "TODO");  // TODO do the documentation...
    // int32_t add_stream(int32_t nb)
    carmaWrapObj.def("add_stream", (int32_t (Class::*)(int32_t)) & Class::add_stream,
                     "TODO",
                     py::arg("np"));  // TODO do the documentation...
    // int32_t del_stream()
    carmaWrapObj.def("del_stream", (int32_t (Class::*)()) & Class::del_stream,
                     "TODO");  // TODO do the documentation...
    // int32_t del_stream(int32_t nb)
    carmaWrapObj.def("del_stream", (int32_t (Class::*)(int32_t)) & Class::del_stream,
                     "TODO",
                     py::arg("np"));  // TODO do the documentation...
    // int32_t wait_stream(int32_t stream)
    carmaWrapObj.def("wait_stream", &Class::wait_stream, "TODO",
                     py::arg("steam"));  // TODO do the documentation...
    carmaWrapObj.def(
        "swap_ptr", [](Class &obj, Class &obj2) { obj.swap_ptr(obj2.get_data()); },
        "TODO",
        py::arg("ptr"));  // TODO do the documentation...
    // int32_t wait_all_streams()
    carmaWrapObj.def("wait_all_streams", &Class::wait_all_streams,
                     "TODO");  // TODO do the documentation...

    // const int64_t *get_dims()
    carmaWrapObj.def_property_readonly(
        "shape",
        [](Class &frame) -> py::array_t<int64_t> {
          int64_t nb_dim = frame.get_dims(0);
          const int64_t *c_dim = frame.get_dims() + 1;
          return py::array_t<int64_t>(nb_dim, c_dim);
        },
        "TODO");  // TODO do the documentation...

    // int32_t get_nb_elements()
    carmaWrapObj.def_property_readonly("nbElem", &Class::get_nb_elements,
                                       "TODO");  // TODO do the documentation...
    // CarmaContext* get_context()
    carmaWrapObj.def_property_readonly("context", &Class::get_context,
                                       "TODO");  // TODO do the documentation...
    // int32_t get_device()
    carmaWrapObj.def_property_readonly("device", &Class::get_device,
                                       "TODO");  // TODO do the documentation...
    // int32_t get_o_data()
    // carmaWrapObj.def_property_readonly("o_data", &Class::get_o_data_value,
    //                                    "TODO");  // TODO do the
    //                                    documentation...

    // int32_t host2device(T_data *data);
    carmaWrapObj.def(
        "host2device",
        [](Class &c,
           py::array_t<T, py::array::f_style | py::array::forcecast> &data) {
          c.host2device(data.data());
        },
        "TODO",
        py::arg("data").none(false));  // TODO do the documentation...
    // int32_t device2host(T_data *data);
    carmaWrapObj.def(
        "device2host",
        [](Class &c,
           py::array_t<T, py::array::f_style | py::array::forcecast> &data) {
          c.device2host(data.mutable_data());
        },
        "TODO",
        py::arg("data").none(false));  // TODO do the documentation...

    // int32_t copy_into(T_data *data, int32_t nb_elem);
    carmaWrapObj.def(
        "copy_into",
        [](Class &src, Class &dest, int64_t nb_elem) {
          if (nb_elem < 0) {
            nb_elem = src.get_nb_elements();
          }
          src.copy_into(dest, nb_elem);
        },
        "TODO", py::arg("dest"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
    // int32_t copy_from(T_data *data, int32_t nb_elem);
    carmaWrapObj.def(
        "copy_from",
        [](Class &dest, Class &src, int64_t nb_elem) {
          if (nb_elem < 0) {
            nb_elem = dest.get_nb_elements();
          }
          dest.copy_from(src, nb_elem);
        },
        "TODO", py::arg("data"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
#ifdef USE_OCTOPUS
    carmaWrapObj.def("copy_into",
                     (int32_t (Class::*)(ipc::Cacao<T> *)) & Class::copy_into);
    carmaWrapObj.def("copy_from",
                     (int32_t (Class::*)(ipc::Cacao<T> *)) & Class::copy_from);
#endif
    // inline int32_t reset()
    carmaWrapObj.def("reset", (int32_t (Class::*)(void)) &Class::reset,
                     "TODO");  // TODO do the documentation...
  }
};
#endif
