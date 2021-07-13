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

//! \file      obj_complex.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for complex CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _WRAP_OBJ_COMPLEX_H_
#define _WRAP_OBJ_COMPLEX_H_

#include "declare_name.hpp"

#include <carma.h>

#include <wyrm>

#include <type_list.hpp>

namespace py = pybind11;

struct CarmaObjComplexInterfacer {
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
          int ndim = data.ndim() + 1;
          std::vector<long> data_dims(ndim);
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

      return py::buffer_info(frame.get_h_data(), sizeof(T),
                             py::format_descriptor<T>::format(), dims[0], shape,
                             strides);
    });

    carmaWrapObj.def("__repr__", &Class::to_string);

    // int get_nb_streams()
    carmaWrapObj.def_property_readonly("nb_streams", &Class::get_nb_streams,
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
        "swap_ptr", [](Class &obj, Class &obj2) { obj.swap_ptr(obj2.get_data()); },
        "TODO",
        py::arg("ptr"));  // TODO do the documentation...
    // int wait_all_streams()
    carmaWrapObj.def("wait_all_streams", &Class::wait_all_streams,
                     "TODO");  // TODO do the documentation...

    // const long *get_dims()
    carmaWrapObj.def_property_readonly(
        "shape",
        [](Class &frame) -> py::array_t<long> {
          long nb_dim = frame.get_dims(0);
          const long *c_dim = frame.get_dims() + 1;
          return py::array_t<long>(nb_dim, c_dim);
        },
        "TODO");  // TODO do the documentation...

    // int get_nb_elements()
    carmaWrapObj.def_property_readonly("nbElem", &Class::get_nb_elements,
                                       "TODO");  // TODO do the documentation...
    // CarmaContext* get_context()
    carmaWrapObj.def_property_readonly("context", &Class::get_context,
                                       "TODO");  // TODO do the documentation...
    // int get_device()
    carmaWrapObj.def_property_readonly("device", &Class::get_device,
                                       "TODO");  // TODO do the documentation...
    // int get_o_data()
    // carmaWrapObj.def_property_readonly("o_data", &Class::get_o_data_value,
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

    // int copy_into(T_data *data, int nb_elem);
    carmaWrapObj.def(
        "copy_into",
        [](Class &src, Class &dest, long nb_elem) {
          if (nb_elem < 0) {
            nb_elem = src.get_nb_elements();
          }
          src.copy_into(dest, nb_elem);
        },
        "TODO", py::arg("dest"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
    // int copy_from(T_data *data, int nb_elem);
    carmaWrapObj.def(
        "copy_from",
        [](Class &dest, Class &src, long nb_elem) {
          if (nb_elem < 0) {
            nb_elem = dest.get_nb_elements();
          }
          dest.copy_from(src, nb_elem);
        },
        "TODO", py::arg("data"),
        py::arg("nb_elem") = -1);  // TODO do the documentation...
#ifdef USE_OCTOPUS
    carmaWrapObj.def("copy_into",
                     (int (Class::*)(ipc::Cacao<T> *)) & Class::copy_into);
    carmaWrapObj.def("copy_from",
                     (int (Class::*)(ipc::Cacao<T> *)) & Class::copy_from);
#endif
    // inline int reset()
    carmaWrapObj.def("reset", (int (Class::*)(void)) &Class::reset,
                     "TODO");  // TODO do the documentation...

    /**< sum */
    // T_data sum();
    // carmaWrapObj.def("sum", &Class::sum,
    //                  "TODO");  // TODO do the documentation...
    // void init_reduceCub();
    carmaWrapObj.def("init_reduceCub", &Class::init_reduceCub,
                     "TODO");  // TODO do the documentation...
    // void reduceCub();
    carmaWrapObj.def("reduceCub", (void (Class::*)(void)) &Class::reduceCub,
                     "TODO");  // TODO do the documentation...
    // void clip(T_data min, T_data max);
    // carmaWrapObj.def(
    //     "clip", &Class::clip, "TODO", py::arg("data_min").none(false),
    //     py::arg("data_max").none(false));  // TODO do the documentation...

    // /**< transpose */
    // int transpose(CarmaObj<T_data> *source);
    carmaWrapObj.def(
        "transpose", &Class::transpose, "TODO",
        py::arg("source").none(false));  // TODO do the documentation...
    // //CarmaObj<T_data>& operator= (const CarmaObj<T_data>& obj);

    // /**< Cublas V2 */
    // int imax(int incx);
    carmaWrapObj.def("aimax", &Class::aimax, "TODO",
                     py::arg("incx") = 1);  // TODO do the documentation...
    // int imin(int incx);
    carmaWrapObj.def("aimin", &Class::aimin, "TODO",
                     py::arg("incx") = 1);  // TODO do the documentation...
    // T_data asum(int incx);
    // carmaWrapObj.def("asum", &Class::asum, "TODO",
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // T_data nrm2(int incx);
    // carmaWrapObj.def("nrm2", &Class::nrm2, "TODO",
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // T_data dot(CarmaObj<T_data> *source, int incx, int incy);
    // carmaWrapObj.def("dot", &Class::dot, "TODO", py::arg("source").none(false),
    //                  py::arg("incx") = 1,
    //                  py::arg("incy") = 1);  // TODO do the documentation...
    // void scale(T_data alpha, int incx);
    // carmaWrapObj.def("scale", &Class::scale, "TODO",
    //                  py::arg("scale").none(false),
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // void swap(CarmaObj<T_data> *source, int incx, int incy);
    carmaWrapObj.def("swap", &Class::swap, "TODO",
                     py::arg("source").none(false), py::arg("incx") = 1,
                     py::arg("incy") = 1);  // TODO do the documentation...
    // void copy(CarmaObj<T_data> *source, int incx, int incy);
    carmaWrapObj.def("copy", &Class::copy,
                     "TODO");  // TODO do the documentation...

    // /**< Curand */
    carmaWrapObj.def("is_rng_init", &Class::is_rng_init);
    // int init_prng();
    carmaWrapObj.def("init_prng", (int (Class::*)()) & Class::init_prng);
    // int init_prng(long seed);
    carmaWrapObj.def("init_prng", (int (Class::*)(long)) & Class::init_prng);
    // int destroy_prng();
    carmaWrapObj.def("destroy_prng", &Class::destroy_prng);
    // int prng(T_data *output, char gtype, float alpha, float beta);
    // carmaWrapObj.def("prng",
    //                  (int (Class::*)(T *, char, float, float)) & Class::prng);
    // int prng(T_data *output, char gtype, float alpha);
    // carmaWrapObj.def("prng", (int (Class::*)(T *, char, float)) & Class::prng);
    // int prng(char gtype, float alpha, float beta);
    carmaWrapObj.def("prng",
                     (int (Class::*)(char, float, float)) & Class::prng);
    // int prng(char gtype, float alpha);
    carmaWrapObj.def("prng", (int (Class::*)(char, float)) & Class::prng);
    // int prng(char gtype);
    carmaWrapObj.def("prng", (int (Class::*)(char)) & Class::prng);

    carmaWrapObj.def(
        "random",
        [](Class &data, int seed, char gtype) {
          data.init_prng(seed);
          data.prng(gtype);
        },
        py::arg("seed") = 1234, py::arg("j") = 'U');

    carmaWrapObj.def(
        "random_host",
        [](Class &data, int seed, char gtype) {
          data.init_prng_host(seed);
          data.prng_host(gtype);
        },
        py::arg("seed") = 1234, py::arg("j") = 'U');

    // // int prng_montagn( float init_montagn );
    // carmaWrapObj.def("prng_montagn", &Class::prng_montagn);

    // // int init_prng_host(int seed);
    // carmaWrapObj.def("init_prng_host",
    //                  (int (Class::*)(int)) & Class::init_prng_host);
    // // int prng_host(char gtype);
    // carmaWrapObj.def("prng_host", (int (Class::*)(char)) & Class::prng_host);
    // // int prng_host(char gtype, T_data stddev);
    // carmaWrapObj.def("prng_host", (int (Class::*)(char, T)) & Class::prng_host);
    // // int prng_host(char gtype, T_data stddev, T_data alpha);
    // carmaWrapObj.def("prng_host",
    //                  (int (Class::*)(char, T, T)) & Class::prng_host);
    // // int destroy_prng_host();
    // carmaWrapObj.def("destroy_prng_host", &Class::destroy_prng_host);

    carmaWrapObj.def(
        "fft",
        [](Class &data, Class &dest, int direction) {
          throw std::runtime_error("not implemented");
          //  const long *dims = data.get_dims();
          //  cufftHandle *handle = data.get_plan();
          //  if(dest == nullptr) {
          //    dest = Class(data.get_context(), dims);
          //  }
          //  carma_initfft(dims, handle, carma_select_plan<T,T>());
          //  CarmaFFT(data.get_data(), dest.get_data(), direction, handle);
        },
        py::arg("dest") = nullptr, py::arg("direction") = 1);
    // CU functions clip
    // template<class T_data>
    // void clip_array(T_data *d_data, T_data min, T_data max, int N,
    // CarmaDevice *device);

    // CU functions sum
    // template<class T_data>
    // void reduce(int size, int threads, int blocks, T_data *d_idata,
    //             T_data *d_odata);
    // template<class T_data>
    // T_data reduce(T_data * data, int N);

    // CU functions transpose
    // template<class T_data>
    // int transposeCU(T_data *d_idata, T_data *d_odata, long N1, long N2);

    // CU functions generic
    // template<class T_data>
    // int launch_generic1d(T_data *d_idata, T_data *d_odata, int N,
    //                     CarmaDevice *device);
    // template<class T_data>
    // int launch_generic2d(T_data *d_odata, T_data *d_idata, int N1, int N2);

    // CU functions curand
    // int carma_prng_init(int *seed, const int nb_threads, const int nb_blocks,
    //                     curandState *state);
    // template<class T>
    // int carma_prng_cu(T *results, const int nb_threads, const int nb_blocks,
    //                   curandState *state, char gtype, int n, float alpha,
    //                   float beta);
    // template<class T>
    // int carma_curand_montagn(curandState *state, T *d_odata, int N,
    // CarmaDevice *device);

    // CU functions fft
    // template<class T_in, class T_out>
    // cufftType carma_select_plan();
    // template<class T_in, class T_out>
    // void carma_initfft(const long *dims_data, cufftHandle *plan, cufftType
    // type_plan); template<class T_in, class T_out> int CarmaFFT(T_in *input,
    // T_out *output, int dir, cufftHandle plan);

    // CU functions generic
    // template<class T_data>
    // int fillindex(T_data *d_odata, T_data *d_idata, int *indx, int N,
    //               CarmaDevice *device);
    // template<class T_data>
    // int fillvalues(T_data *d_odata, T_data val, int N,
    //               CarmaDevice *device);
    // template<class T>
    // int getarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
    //               CarmaDevice *device);
    // template<class T>
    // int fillarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int
    // N,
    //                 CarmaDevice *device);
    // template<class T>
    // int fillarray2d2(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int
    // N,
    //                 CarmaDevice *device);
    // template<class T>
    // int fill_sym_matrix(char src_uplo, T *d_data, int Ncol, int N,
    //                     CarmaDevice *device);
    // template<class T>
    // int carma_plus(T *d_odata, T elpha, int N, CarmaDevice *device);
    // template<class T>
    // int carma_plusai(T *d_odata, T *i_data, int i, int sgn, int N,
    //                 CarmaDevice *device);

    // CU functions fftconv
    // int fftconv_unpad(float *d_odata, float *d_idata, int fftW, int dataH,
    //                   int dataW, int N, int n, int nim);
    // int carma_initfftconv(CarmaObjS *data_in, CarmaObjS *kernel_in, CarmaObjS
    // *padded_data, CarmaObjC *padded_spectrum, int kernelY, int kernelX);

    // CPP functions fftconv
    // int carma_fftconv(CarmaObjS *data_out, CarmaObjS *padded_data,
    //                   CarmaObjC *padded_spectrum, int kernelY, int kernelX);

    mod.def(
        appendName<T>("syevd_").data(),
        [](Class &d_mat_a, Class &eigenvals, Class *d_U, bool computeU) {
          if (d_U == nullptr) {
            if (computeU) {
              carma_syevd(SOLVER_EIG_MODE_VECTOR, &d_mat_a, &eigenvals);
            } else {
              carma_syevd(SOLVER_EIG_MODE_NOVECTOR, &d_mat_a, &eigenvals);
            }
          } else {
            d_U->copy_from(d_mat_a, d_mat_a.get_nb_elements());
            if (computeU) {
              carma_syevd(SOLVER_EIG_MODE_VECTOR, d_U, &eigenvals);
            } else {
              carma_syevd(SOLVER_EIG_MODE_NOVECTOR, d_U, &eigenvals);
            }
          }
        },
        py::arg("d_mat_a"), py::arg("eigenvals"), py::arg("d_U") = nullptr,
        py::arg("computeU") = true);

#ifdef USE_MAGMA
    // MAGMA functions

    // template<class T>
    // int carma_svd(CarmaObj<T> *imat, CarmaObj<T> *eigenvals,
    //               CarmaObj<T> *mod2act, CarmaObj<T> *mes2mod);
    // mod.def(appendName<T>("svd_").data(), &carma_svd<T>);

    // TODO after CarmaHostObj
    // template<class T>
    // int carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
    // *eigenvals);
    // mod.def( appendName<T>("syevd_").data(), py::overload_cast<char, Class
    // *, ClassHost *>(&carma_magma_syevd<T>)); mod.def(
    // appendName<T>("syevd_").data(), py::overload_cast<char, long, T *, T
    // *>(&carma_magma_syevd<T>));
    mod.def(
        appendName<T>("magma_syevd_").data(),
        [](Class &d_mat_a, ClassHost &eigenvals, Class *d_U, bool computeU) {
          if (d_U == nullptr) {
            if (computeU) {
              carma_magma_syevd(SOLVER_EIG_MODE_VECTOR, &d_mat_a, &eigenvals);
            } else {
              carma_magma_syevd(SOLVER_EIG_MODE_NOVECTOR, &d_mat_a, &eigenvals);
            }
          } else {
            d_U->copy_from(d_mat_a, d_mat_a.get_nb_elements());
            if (computeU) {
              carma_magma_syevd(SOLVER_EIG_MODE_VECTOR, d_U, &eigenvals);
            } else {
              carma_magma_syevd(SOLVER_EIG_MODE_NOVECTOR, d_U, &eigenvals);
            }
          }
        },
        py::arg("d_mat_a"), py::arg("eigenvals"), py::arg("d_U") = nullptr,
        py::arg("computeU") = true);

    // template<class T, int method>
    // int carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
    // *eigenvals); template<class T> int carma_magma_syevd_m(long ngpu, char
    // jobz, long N, T *mat, T *eigenvals); template<class T> int
    // carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
    //                   CarmaHostObj<T> *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
    //                   CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
    // template<class T>
    // int carma_magma_getri(CarmaObj<T> *d_iA);
    mod.def(appendName<T>("magma_getri_").data(), &carma_magma_getri<T>);

    // template<class T>
    // int carma_magma_potr_inv(CarmaObj<T> *d_iA);
    mod.def(appendName<T>("magma_potri_").data(), &carma_magma_potr_inv<T>);

    // TODO after CarmaHostObj
    // template<class T>
    // int carma_magma_potr_inv_m(long num_gpus, CarmaHostObj<T> *h_A,
    // CarmaObj<T> *d_iA);

    // MAGMA functions (direct access)
    // template<class T>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T, int method>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, long N, T *mat, T
    // *eigenvals); template<class T> int carma_magma_potr_inv_m(long num_gpus,
    // long N, T *h_A, T *d_iA);
#endif

  }
};
#endif
