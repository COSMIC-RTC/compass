// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      obj_half.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for half float carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _WRAP_OBJ_HALF_H_
#define _WRAP_OBJ_HALF_H_

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

#include <type_list.hpp>

namespace py = pybind11;

struct CarmaObjHalfInterfacer {
  template <typename T>
  static void call(py::module &mod) {
    auto name = appendName<T>("obj_");
    using Class = carma_obj<T>;
    using ClassHost = carma_host_obj<T>;

    py::class_<Class> carmaWrapObj(mod, name.data(), py::buffer_protocol());
    carmaWrapObj.def(
        py::init([](carma_context &c,
                    const py::array_t<T, py::array::f_style |
                                             py::array::forcecast> &data) {
          int ndim = data.ndim() + 1;
          std::vector<long> data_dims(ndim);
          data_dims[0] = data.ndim();
          copy(data.shape(), data.shape() + data.ndim(), begin(data_dims) + 1);
          return std::unique_ptr<Class>(
              new Class(&c, data_dims.data(), (const T *)data.data()));
        }),
        "TODO",  // TODO do the documentation...
        py::arg("context").none(false), py::arg("h_data").none(false));

    carmaWrapObj.def(py::init([](carma_context &c, const Class &data) {
                       return std::unique_ptr<Class>(new Class(&c, &data));
                     }),
                     "TODO",  // TODO do the documentation...
                     py::arg("context").none(false),
                     py::arg("d_data").none(false));

    carmaWrapObj.def_buffer([](Class &frame) -> py::buffer_info {
      frame.sync_h_data();

      const long *dims = frame.getDims();
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

    carmaWrapObj.def("__repr__", &Class::toString);

    // int get_nbStreams()
    carmaWrapObj.def_property_readonly("nbStreams", &Class::get_nbStreams,
                                       "TODO");  // TODO do the documentation...
    // int add_stream()
    carmaWrapObj.def("add_stream", (int (Class::*)()) & Class::add_stream,
                     "TODO");  // TODO do the documentation...
    // int add_stream(int nb)
    carmaWrapObj.def("add_stream", (int (Class::*)(int)) & Class::add_stream,
                     "TODO",
                     py::arg("np"));  // TODO do the documentation...
    // int del_stream()
    carmaWrapObj.def("del_stream", (int (Class::*)()) & Class::del_stream,
                     "TODO");  // TODO do the documentation...
    // int del_stream(int nb)
    carmaWrapObj.def("del_stream", (int (Class::*)(int)) & Class::del_stream,
                     "TODO",
                     py::arg("np"));  // TODO do the documentation...
    // int wait_stream(int stream)
    carmaWrapObj.def("wait_stream", &Class::wait_stream, "TODO",
                     py::arg("steam"));  // TODO do the documentation...
    carmaWrapObj.def(
        "swapPtr", [](Class &obj, Class &obj2) { obj.swapPtr(obj2.getData()); },
        "TODO",
        py::arg("ptr"));  // TODO do the documentation...
    // int wait_all_streams()
    carmaWrapObj.def("wait_all_streams", &Class::wait_all_streams,
                     "TODO");  // TODO do the documentation...

    // const long *getDims()
    carmaWrapObj.def_property_readonly(
        "shape",
        [](Class &frame) -> py::array_t<long> {
          long nb_dim = frame.getDims(0);
          const long *c_dim = frame.getDims() + 1;
          return py::array_t<long>(nb_dim, c_dim);
        },
        "TODO");  // TODO do the documentation...

    // int getNbElem()
    carmaWrapObj.def_property_readonly("nbElem", &Class::getNbElem,
                                       "TODO");  // TODO do the documentation...
    // carma_context* getContext()
    carmaWrapObj.def_property_readonly("context", &Class::getContext,
                                       "TODO");  // TODO do the documentation...
    // int getDevice()
    carmaWrapObj.def_property_readonly("device", &Class::getDevice,
                                       "TODO");  // TODO do the documentation...
    // int getOData()
    // carmaWrapObj.def_property_readonly("o_data", &Class::getODataValue,
    //                                    "TODO");  // TODO do the
    //                                    documentation...

    // int host2device(T_data *data);
    carmaWrapObj.def(
        "host2device",
        [](Class &c,
           py::array_t<T, py::array::f_style | py::array::forcecast> &data) {
          c.host2device(data.data());
        },
        "TODO",
        py::arg("data").none(false));  // TODO do the documentation...
    // int device2host(T_data *data);
    carmaWrapObj.def(
        "device2host",
        [](Class &c,
           py::array_t<T, py::array::f_style | py::array::forcecast> &data) {
          c.device2host(data.mutable_data());
        },
        "TODO",
        py::arg("data").none(false));  // TODO do the documentation...

    // int copyInto(T_data *data, int nb_elem);
    carmaWrapObj.def(
        "copyInto",
        [](Class &src, Class &dest, long nb_elem) {
          if (nb_elem < 0) {
            nb_elem = src.getNbElem();
          }
          src.copyInto(dest, nb_elem);
        },
        "TODO", py::arg("dest"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
    // int copyFrom(T_data *data, int nb_elem);
    carmaWrapObj.def(
        "copyFrom",
        [](Class &dest, Class &src, long nb_elem) {
          if (nb_elem < 0) {
            nb_elem = dest.getNbElem();
          }
          dest.copyFrom(src, nb_elem);
        },
        "TODO", py::arg("data"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
#ifdef USE_OCTOPUS
    carmaWrapObj.def("copyInto",
                     (int (Class::*)(ipc::Cacao<T> *)) & Class::copyInto);
    carmaWrapObj.def("copyFrom",
                     (int (Class::*)(ipc::Cacao<T> *)) & Class::copyFrom);
#endif
    // inline int reset()
    carmaWrapObj.def("reset", &Class::reset,
                     "TODO");  // TODO do the documentation...
  }
};
#endif
