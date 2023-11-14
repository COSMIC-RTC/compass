// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      host_obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaHostObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _WRAP_HOST_OBJ_H_
#define _WRAP_HOST_OBJ_H_

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

#include <type_list.hpp>

namespace py = pybind11;



struct CarmaHostObjInterfacer {
  template <typename T> static void call(py::module &mod) {
    auto name = appendName<T>("host_obj_");
    using Class = CarmaHostObj<T>;

    py::class_<Class>(mod, name.data(), py::buffer_protocol())
        .def(py::init([](const py::array_t<T, py::array::f_style |
                                                  py::array::forcecast> &data) {
               int ndim = data.ndim() + 1;
               std::vector<long> data_dims(ndim);
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

          const long *dims = frame.get_dims();
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

        // const long *get_dims()
        .def_property_readonly("shape",
                               [](Class &frame) -> py::array_t<long> {
                                 long nb_dim = frame.get_dims(0);
                                 const long *c_dim = frame.get_dims() + 1;
                                 return py::array_t<long>(nb_dim, c_dim);
                               },
                               "TODO") // TODO do the documentation...

        // int get_nb_elements()
        .def_property_readonly("nbElem", &Class::get_nb_elements,
                               "TODO") // TODO do the documentation...
    ;

#ifdef USE_MAGMA
    // MAGMA functions

    // template<class T>
    // int carma_svd(CarmaObj<T> *imat, CarmaObj<T> *eigenvals,
    //               CarmaObj<T> *mod2act, CarmaObj<T> *mes2mod);
    mod.def(appendName<T>("magma_svd_cpu_").data(), py::overload_cast<Class *,
                  Class *, Class *, Class *>(&carma_magma_svd_cpu<T>));

    // TODO after CarmaHostObj
    // template<class T>
    // int carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
    // *eigenvals);
    // mod.def(appendName<T>("syevd_").data(), &carma_magma_syevd<T>);
    mod.def(appendName<T>("magma_syevd_cpu_").data(), py::overload_cast<char, Class *, Class *>(&carma_magma_syevd_cpu<T>));

    // template<class T, int method>
    // int carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
    // *eigenvals); template<class T> int carma_magma_syevd_m(long ngpu, char jobz,
    // long N, T *mat, T *eigenvals); template<class T> int carma_magma_syevd_m(long
    // ngpu, char jobz, CarmaHostObj<T> *mat,
    //                   CarmaHostObj<T> *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
    //                   CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
    // template<class T>
    // int carma_magma_getri(CarmaObj<T> *d_iA);
    mod.def(appendName<T>("magma_getri_cpu_").data(), py::overload_cast<Class *>(&carma_magma_getri_cpu<T>));

    // template<class T>
    // int carma_magma_potr_inv(CarmaObj<T> *d_iA);
    mod.def(appendName<T>("magma_potri_cpu_").data(), py::overload_cast<Class *>(&carma_magma_potr_inv_cpu<T>));

    // TODO after CarmaHostObj
    // template<class T>
    // int carma_magma_potr_inv_m(long num_gpus, CarmaHostObj<T> *h_A, CarmaObj<T>
    // *d_iA);

    // MAGMA functions (direct access)
    // template<class T>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T, int method>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
    // template<class T>
    // int carma_magma_potr_inv_m(long num_gpus, long N, T *h_A, T *d_iA);
#endif

  }
};
#endif
