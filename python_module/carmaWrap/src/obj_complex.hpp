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

//! \file      obj_complex.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for complex CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _WRAP_OBJ_COMPLEX_H_
#define _WRAP_OBJ_COMPLEX_H_

#include "declare_name.hpp"

#include <carma.hpp>

#include <type_list.hpp>

namespace py = pybind11;

template <typename T>
struct ToCuComplex
{
  using Type = T;
  using CarmaObjType = CarmaObj<T>;
  using CarmaHostObjType = CarmaHostObj<T>;
};

template <>
struct ToCuComplex<float>
{
  using Type = cuFloatComplex;
  using CarmaObjType = CarmaObj<cuFloatComplex>;
  using CarmaHostObjType = CarmaHostObj<cuFloatComplex>;
};
template <>
struct ToCuComplex<double>
{
  using Type = cuDoubleComplex;
  using CarmaObjType = CarmaObj<cuDoubleComplex>;
  using CarmaHostObjType = CarmaHostObj<cuDoubleComplex>;
};

struct CarmaObjComplexInterfacer
{
  template <typename T>
  static void call(py::module &mod)
  {
    auto name = appendName<T>("obj_") + "_complex";
    using Type = typename ToCuComplex<T>::Type;
    using Class = typename ToCuComplex<T>::CarmaObjType;
    using ClassHost = typename ToCuComplex<T>::CarmaHostObjType;

    py::class_<Class> carmaWrapObj(mod, name.data(), py::buffer_protocol());
    carmaWrapObj.def(
        py::init([](CarmaContext &c,
                    const py::array_t<std::complex<T>, py::array::f_style |
                                                           py::array::forcecast> &data)
                 {
          int32_t ndim = data.ndim() + 1;
          std::vector<int64_t> data_dims(ndim);
          data_dims[0] = data.ndim();
          copy(data.shape(), data.shape() + data.ndim(), begin(data_dims) + 1);
          return std::unique_ptr<Class>(
              new Class(&c, data_dims.data(), (const Type *)data.data())); }),
        "TODO", // TODO do the documentation...
        py::arg("context").none(false), py::arg("h_data").none(false));

    carmaWrapObj.def(py::init([](CarmaContext &c, const Class &data)
                              { return std::unique_ptr<Class>(new Class(&c, &data)); }),
                     "TODO", // TODO do the documentation...
                     py::arg("context").none(false),
                     py::arg("d_data").none(false));

    carmaWrapObj.def_buffer([](Class &frame) -> py::buffer_info
                            {
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

      return py::buffer_info(frame.get_h_data(), sizeof(std::complex<T>),
                             py::format_descriptor<std::complex<T>>::format(), dims[0], shape,
                             strides); });

    carmaWrapObj.def("__repr__", &Class::to_string);

    // int32_t get_nb_streams()
    carmaWrapObj.def_property_readonly("nb_streams", &Class::get_nb_streams,
                                       "TODO"); // TODO do the documentation...
    // int32_t add_stream()
    carmaWrapObj.def("add_stream", (int32_t(Class::*)()) & Class::add_stream,
                     "TODO"); // TODO do the documentation...
    // int32_t add_stream(int32_t nb)
    carmaWrapObj.def("add_stream", (int32_t(Class::*)(int32_t)) & Class::add_stream,
                     "TODO",
                     py::arg("np")); // TODO do the documentation...
    // int32_t del_stream()
    carmaWrapObj.def("del_stream", (int32_t(Class::*)()) & Class::del_stream,
                     "TODO"); // TODO do the documentation...
    // int32_t del_stream(int32_t nb)
    carmaWrapObj.def("del_stream", (int32_t(Class::*)(int32_t)) & Class::del_stream,
                     "TODO",
                     py::arg("np")); // TODO do the documentation...
    // int32_t wait_stream(int32_t stream)
    carmaWrapObj.def("wait_stream", &Class::wait_stream, "TODO",
                     py::arg("steam")); // TODO do the documentation...
    carmaWrapObj.def(
        "swap_ptr", [](Class &obj, Class &obj2)
        { obj.swap_ptr(obj2.get_data()); },
        "TODO",
        py::arg("ptr")); // TODO do the documentation...
    // int32_t wait_all_streams()
    carmaWrapObj.def("wait_all_streams", &Class::wait_all_streams,
                     "TODO"); // TODO do the documentation...

    // const int64_t *get_dims()
    carmaWrapObj.def_property_readonly(
        "shape",
        [](Class &frame) -> py::array_t<int64_t>
        {
          int64_t nb_dim = frame.get_dims(0);
          const int64_t *c_dim = frame.get_dims() + 1;
          return py::array_t<int64_t>(nb_dim, c_dim);
        },
        "TODO"); // TODO do the documentation...

    // int32_t get_nb_elements()
    carmaWrapObj.def_property_readonly("nbElem", &Class::get_nb_elements,
                                       "TODO"); // TODO do the documentation...
    // CarmaContext* get_context()
    carmaWrapObj.def_property_readonly("context", &Class::get_context,
                                       "TODO"); // TODO do the documentation...
    // int32_t get_device()
    carmaWrapObj.def_property_readonly("device", &Class::get_device,
                                       "TODO"); // TODO do the documentation...
    // int32_t get_o_data()
    // carmaWrapObj.def_property_readonly("o_data", &Class::get_o_data_value,
    //                                    "TODO");  // TODO do the
    //                                    documentation...

    // int32_t host2device(T_data *data);
    carmaWrapObj.def(
        "host2device",
        [](Class &c,
           py::array_t<std::complex<T>, py::array::f_style | py::array::forcecast> &data)
        {
          c.host2device((Type*)data.data());
        },
        "TODO",
        py::arg("data").none(false)); // TODO do the documentation...
    // int32_t device2host(T_data *data);
    carmaWrapObj.def(
        "device2host",
        [](Class &c,
           py::array_t<std::complex<T>, py::array::f_style | py::array::forcecast> &data)
        {
          c.device2host((Type*)data.mutable_data());
        },
        "TODO",
        py::arg("data").none(false)); // TODO do the documentation...

    // int32_t copy_into(T_data *data, int32_t nb_elem);
    carmaWrapObj.def(
        "copy_into",
        [](Class &src, Class &dest, int64_t nb_elem)
        {
          if (nb_elem < 0)
          {
            nb_elem = src.get_nb_elements();
          }
          src.copy_into(dest, nb_elem);
        },
        "TODO", py::arg("dest"),
        py::arg("nb_elem") = -1); // TODO do the documentation...
    // int32_t copy_from(T_data *data, int32_t nb_elem);
    carmaWrapObj.def(
        "copy_from",
        [](Class &dest, Class &src, int64_t nb_elem)
        {
          if (nb_elem < 0)
          {
            nb_elem = dest.get_nb_elements();
          }
          dest.copy_from(src, nb_elem);
        },
        "TODO", py::arg("data"),
        py::arg("nb_elem") = -1); // TODO do the documentation...
#ifdef USE_OCTOPUS
    carmaWrapObj.def("copy_into",
                     (int32_t(Class::*)(ipc::Cacao<Type> *)) & Class::copy_into);
    carmaWrapObj.def("copy_from",
                     (int32_t(Class::*)(ipc::Cacao<Type> *)) & Class::copy_from);
#endif
    // inline int32_t reset()
    carmaWrapObj.def("reset", (int32_t(Class::*)(void)) & Class::reset,
                     "TODO"); // TODO do the documentation...

    /**< sum */
    // T_data sum();
    // carmaWrapObj.def("sum", &Class::sum,
    //                  "TODO");  // TODO do the documentation...
    // void init_reduceCub();
    carmaWrapObj.def("init_reduceCub", &Class::init_reduceCub,
                     "TODO"); // TODO do the documentation...
    // void reduceCub();
    carmaWrapObj.def("reduceCub", (void(Class::*)(void)) & Class::reduceCub,
                     "TODO"); // TODO do the documentation...
    // void clip(T_data min, T_data max);
    // carmaWrapObj.def(
    //     "clip", &Class::clip, "TODO", py::arg("data_min").none(false),
    //     py::arg("data_max").none(false));  // TODO do the documentation...

    // /**< transpose */
    // int32_t transpose(CarmaObj<T_data> *source);
    carmaWrapObj.def(
        "transpose", &Class::transpose, "TODO",
        py::arg("source").none(false)); // TODO do the documentation...
    // //CarmaObj<T_data>& operator= (const CarmaObj<T_data>& obj);

    // /**< Cublas V2 */
    // int32_t imax(int32_t incx);
    carmaWrapObj.def("aimax", &Class::aimax, "TODO",
                     py::arg("incx") = 1); // TODO do the documentation...
    // int32_t imin(int32_t incx);
    carmaWrapObj.def("aimin", &Class::aimin, "TODO",
                     py::arg("incx") = 1); // TODO do the documentation...
    // T_data asum(int32_t incx);
    // carmaWrapObj.def("asum", &Class::asum, "TODO",
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // T_data nrm2(int32_t incx);
    // carmaWrapObj.def("nrm2", &Class::nrm2, "TODO",
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // T_data dot(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
    // carmaWrapObj.def("dot", &Class::dot, "TODO", py::arg("source").none(false),
    //                  py::arg("incx") = 1,
    //                  py::arg("incy") = 1);  // TODO do the documentation...
    // void scale(T_data alpha, int32_t incx);
    // carmaWrapObj.def("scale", &Class::scale, "TODO",
    //                  py::arg("scale").none(false),
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // void swap(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
    carmaWrapObj.def("swap", &Class::swap, "TODO",
                     py::arg("source").none(false), py::arg("incx") = 1,
                     py::arg("incy") = 1); // TODO do the documentation...
    // void copy(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
    carmaWrapObj.def("copy", &Class::copy,
                     "TODO"); // TODO do the documentation...

    // /**< Curand */
    carmaWrapObj.def("is_rng_init", &Class::is_rng_init);
    // int32_t init_prng();
    carmaWrapObj.def("init_prng", (int32_t(Class::*)()) & Class::init_prng);
    // int32_t init_prng(int64_t seed);
    carmaWrapObj.def("init_prng", (int32_t(Class::*)(int64_t)) & Class::init_prng);
    // int32_t destroy_prng();
    carmaWrapObj.def("destroy_prng", &Class::destroy_prng);
    // int32_t prng(T_data *output, char gtype, float alpha, float beta);
    // carmaWrapObj.def("prng",
    //                  (int32_t (Class::*)(T *, char, float, float)) & Class::prng);
    // int32_t prng(T_data *output, char gtype, float alpha);
    // carmaWrapObj.def("prng", (int32_t (Class::*)(T *, char, float)) & Class::prng);
    // int32_t prng(char gtype, float alpha, float beta);
    carmaWrapObj.def("prng",
                     (int32_t(Class::*)(char, float, float)) & Class::prng);
    // int32_t prng(char gtype, float alpha);
    carmaWrapObj.def("prng", (int32_t(Class::*)(char, float)) & Class::prng);
    // int32_t prng(char gtype);
    carmaWrapObj.def("prng", (int32_t(Class::*)(char)) & Class::prng);

    carmaWrapObj.def(
        "random",
        [](Class &data, int32_t seed, char gtype)
        {
          data.init_prng(seed);
          data.prng(gtype);
        },
        py::arg("seed") = 1234, py::arg("j") = 'U');

    carmaWrapObj.def(
        "random_host",
        [](Class &data, int32_t seed, char gtype)
        {
          data.init_prng_host(seed);
          data.prng_host(gtype);
        },
        py::arg("seed") = 1234, py::arg("j") = 'U');

    // // int32_t prng_montagn( float init_montagn );
    // carmaWrapObj.def("prng_montagn", &Class::prng_montagn);

    // // int32_t init_prng_host(int32_t seed);
    // carmaWrapObj.def("init_prng_host",
    //                  (int32_t (Class::*)(int32_t)) & Class::init_prng_host);
    // // int32_t prng_host(char gtype);
    // carmaWrapObj.def("prng_host", (int32_t (Class::*)(char)) & Class::prng_host);
    // // int32_t prng_host(char gtype, T_data stddev);
    // carmaWrapObj.def("prng_host", (int32_t (Class::*)(char, T)) & Class::prng_host);
    // // int32_t prng_host(char gtype, T_data stddev, T_data alpha);
    // carmaWrapObj.def("prng_host",
    //                  (int32_t (Class::*)(char, T, T)) & Class::prng_host);
    // // int32_t destroy_prng_host();
    // carmaWrapObj.def("destroy_prng_host", &Class::destroy_prng_host);

    carmaWrapObj.def(
        "fft",
        [](Class &data, Class &dest, int32_t direction)
        {
          throw std::runtime_error("not implemented");
          //  const int64_t *dims = data.get_dims();
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
    // void clip_array(T_data *d_data, T_data min, T_data max, int32_t N,
    // CarmaDevice *device);

    // CU functions sum
    // template<class T_data>
    // void reduce(int32_t size, int32_t threads, int32_t blocks, T_data *d_idata,
    //             T_data *d_odata);
    // template<class T_data>
    // T_data reduce(T_data * data, int32_t N);

    // CU functions transpose
    // template<class T_data>
    // int32_t transposeCU(T_data *d_idata, T_data *d_odata, int64_t N1, int64_t N2);

    // CU functions generic
    // template<class T_data>
    // int32_t launch_generic1d(T_data *d_idata, T_data *d_odata, int32_t N,
    //                     CarmaDevice *device);
    // template<class T_data>
    // int32_t launch_generic2d(T_data *d_odata, T_data *d_idata, int32_t N1, int32_t N2);

    // CU functions curand
    // int32_t carma_prng_init(int32_t *seed, const int32_t nb_threads, const int32_t nb_blocks,
    //                     curandState *state);
    // template<class T>
    // int32_t carma_prng_cu(T *results, const int32_t nb_threads, const int32_t nb_blocks,
    //                   curandState *state, char gtype, int32_t n, float alpha,
    //                   float beta);
    // template<class T>
    // int32_t carma_curand_montagn(curandState *state, T *d_odata, int32_t N,
    // CarmaDevice *device);

    // CU functions fft
    // template<class T_in, class T_out>
    // cufftType carma_select_plan();
    // template<class T_in, class T_out>
    // void carma_initfft(const int64_t *dims_data, cufftHandle *plan, cufftType
    // type_plan); template<class T_in, class T_out> int32_t CarmaFFT(T_in *input,
    // T_out *output, int32_t dir, cufftHandle plan);

    // CU functions generic
    // template<class T_data>
    // int32_t fillindex(T_data *d_odata, T_data *d_idata, int32_t *indx, int32_t N,
    //               CarmaDevice *device);
    // template<class T_data>
    // int32_t fillvalues(T_data *d_odata, T_data val, int32_t N,
    //               CarmaDevice *device);
    // template<class T>
    // int32_t getarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
    //               CarmaDevice *device);
    // template<class T>
    // int32_t fillarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t
    // N,
    //                 CarmaDevice *device);
    // template<class T>
    // int32_t fillarray2d2(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t
    // N,
    //                 CarmaDevice *device);
    // template<class T>
    // int32_t fill_sym_matrix(char src_uplo, T *d_data, int32_t Ncol, int32_t N,
    //                     CarmaDevice *device);
    // template<class T>
    // int32_t carma_plus(T *d_odata, T elpha, int32_t N, CarmaDevice *device);
    // template<class T>
    // int32_t carma_plusai(T *d_odata, T *i_data, int32_t i, int32_t sgn, int32_t N,
    //                 CarmaDevice *device);

    // CU functions fftconv
    // int32_t fftconv_unpad(float *d_odata, float *d_idata, int32_t fftW, int32_t dataH,
    //                   int32_t dataW, int32_t N, int32_t n, int32_t nim);
    // int32_t carma_initfftconv(CarmaObjS *data_in, CarmaObjS *kernel_in, CarmaObjS
    // *padded_data, CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX);

    // CPP functions fftconv
    // int32_t carma_fftconv(CarmaObjS *data_out, CarmaObjS *padded_data,
    //                   CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX);

    mod.def(
        appendName<T>("syevd_").data(),
        [](Class &d_mat_a, Class &eigenvals, Class *d_U, bool computeU)
        {
          if (d_U == nullptr)
          {
            if (computeU)
            {
              carma_syevd(SOLVER_EIG_MODE_VECTOR, &d_mat_a, &eigenvals);
            }
            else
            {
              carma_syevd(SOLVER_EIG_MODE_NOVECTOR, &d_mat_a, &eigenvals);
            }
          }
          else
          {
            d_U->copy_from(d_mat_a, d_mat_a.get_nb_elements());
            if (computeU)
            {
              carma_syevd(SOLVER_EIG_MODE_VECTOR, d_U, &eigenvals);
            }
            else
            {
              carma_syevd(SOLVER_EIG_MODE_NOVECTOR, d_U, &eigenvals);
            }
          }
        },
        py::arg("d_mat_a"), py::arg("eigenvals"), py::arg("d_U") = nullptr,
        py::arg("computeU") = true);
  }
};
#endif
