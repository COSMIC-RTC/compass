// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COSMIC Team <https://github.com/COSMIC-RTC/compass>
//  All rights reserved.
// -----------------------------------------------------------------------------

//! \file      obj.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _WRAP_OBJ_H_
#define _WRAP_OBJ_H_

#include "declare_name.hpp"

#include <carma.hpp>

#include <type_list.hpp>

namespace py = pybind11;

struct CarmaObjInterfacer {
  template <typename T>
  static void call(py::module &mod) {
    auto name = appendName<T>("obj_");
    using Class = CarmaObj<T>;
    using ClassHost = CarmaHostObj<T>;

    py::class_<Class>(mod, name.data(), py::buffer_protocol())
        .def(py::init([](CarmaContext &c,
                         const py::array_t<T, py::array::f_style |
                                                  py::array::forcecast> &data) {
               int32_t ndim = data.ndim() + 1;
               std::vector<int64_t> data_dims(ndim);
               data_dims[0] = data.ndim();
               copy(data.shape(), data.shape() + data.ndim(),
                    begin(data_dims) + 1);
               return std::unique_ptr<Class>(
                   new Class(&c, data_dims.data(), (const T *)data.data()));
             }),
             "TODO",  // TODO do the documentation...
             py::arg("context").none(false), py::arg("h_data").none(false))

        .def(py::init([](CarmaContext &c, const Class &data) {
               return std::unique_ptr<Class>(new Class(&c, &data));
             }),
             "TODO",  // TODO do the documentation...
             py::arg("context").none(false), py::arg("d_data").none(false))

        .def_buffer([](Class &frame) -> py::buffer_info {
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
                                 py::format_descriptor<T>::format(), dims[0],
                                 shape, strides);
        })

        .def("to_cupy", [](Class &cls) {
          py::object MemoryPointer = py::module::import("cupy.cuda.memory").attr("MemoryPointer");
          py::object UnownedMemory = py::module::import("cupy.cuda.memory").attr("UnownedMemory");
          py::object ndarray = py::module::import("cupy").attr("ndarray");
          const int64_t *dims = cls.get_dims();
          std::vector<ssize_t> shape(dims[0]);
          std::vector<ssize_t> strides(dims[0]);
          ssize_t stride = sizeof(T);

          for (ssize_t dim(0); dim < dims[0]; ++dim) {
            shape[dim] = dims[dim + 1];
            strides[dim] = stride;
            stride *= shape[dim];
          }
          return ndarray(shape,
                        py::arg("dtype")=py::format_descriptor<T>::format(),
                        py::arg("memptr")=MemoryPointer(UnownedMemory(reinterpret_cast<intptr_t>(cls.get_data()), stride, py::cast<py::none>(Py_None)), 0),
                        py::arg("strides")=strides);
        }, "TODO")

        .def("__repr__", &Class::to_string)

        // int32_t get_nb_streams()
        .def_property_readonly("nb_streams", &Class::get_nb_streams,
                               "TODO")  // TODO do the documentation...
        // int32_t add_stream()
        .def("add_stream", (int32_t (Class::*)()) & Class::add_stream,
             "TODO")  // TODO do the documentation...
        // int32_t add_stream(int32_t nb)
        .def("add_stream", (int32_t (Class::*)(int32_t)) & Class::add_stream, "TODO",
             py::arg("np"))  // TODO do the documentation...
        // int32_t del_stream()
        .def("del_stream", (int32_t (Class::*)()) & Class::del_stream,
             "TODO")  // TODO do the documentation...
        // int32_t del_stream(int32_t nb)
        .def("del_stream", (int32_t (Class::*)(int32_t)) & Class::del_stream, "TODO",
             py::arg("np"))  // TODO do the documentation...
        // int32_t wait_stream(int32_t stream)
        .def("wait_stream", &Class::wait_stream, "TODO",
             py::arg("steam"))  // TODO do the documentation...
        .def("swap_ptr", [](Class &obj, Class &obj2){
          obj.swap_ptr(obj2.get_data());
        }, "TODO",
             py::arg("ptr"))  // TODO do the documentation...
    // int32_t wait_all_streams()
        .def("wait_all_streams", &Class::wait_all_streams,
             "TODO")  // TODO do the documentation...

        // const int64_t *get_dims()
        .def_property_readonly("shape",
                               [](Class &frame) -> py::array_t<int64_t> {
                                 int64_t nb_dim = frame.get_dims(0);
                                 const int64_t *c_dim = frame.get_dims() + 1;
                                 return py::array_t<int64_t>(nb_dim, c_dim);
                               },
                               "TODO")  // TODO do the documentation...

        // int32_t get_nb_elements()
        .def_property_readonly("nbElem", &Class::get_nb_elements,
                               "TODO")  // TODO do the documentation...
        // CarmaContext* get_context()
        .def_property_readonly("context", &Class::get_context,
                               "TODO")  // TODO do the documentation...
        // int32_t get_device()
        .def_property_readonly("device", &Class::get_device,
                               "TODO")  // TODO do the documentation...
        // int32_t get_o_data()
        .def_property_readonly("o_data", &Class::get_o_data_value,
                               "TODO")  // TODO do the documentation...
        .def_property_readonly("d_ptr", [](Class &cls) {return reinterpret_cast<intptr_t>(cls.get_data());}, "TODO")
        // int32_t host2device(T_data *data);
        .def("host2device",
             [](Class &c,
                py::array_t<T, py::array::f_style | py::array::forcecast>
                    &data) { c.host2device((const T *)data.data()); },
             "TODO",
             py::arg("data").none(false))  // TODO do the documentation...
        // int32_t device2host(T_data *data);
        .def("device2host",
             [](Class &c,
                py::array_t<T, py::array::f_style | py::array::forcecast>
                    &data) { c.device2host((T *)data.mutable_data()); },
             "TODO",
             py::arg("data").none(false))  // TODO do the documentation...

        // int32_t copy_into(T_data *data, int32_t nb_elem);
        .def("copy_into",
             [](Class &src, Class &dest, int64_t nb_elem) {
               if (nb_elem < 0) {
                 nb_elem = src.get_nb_elements();
               }
               src.copy_into(dest, nb_elem);
             },
             "TODO", py::arg("dest"),
             py::arg("nb_elem") = -1)  // TODO do the documentation...
        // int32_t copy_from(T_data *data, int32_t nb_elem);
        .def("copy_from",
             [](Class &dest, Class &src, int64_t nb_elem) {
               if (nb_elem < 0) {
                 nb_elem = dest.get_nb_elements();
               }
               dest.copy_from(src, nb_elem);
             },
             "TODO", py::arg("data"),
             py::arg("nb_elem") = -1)  // TODO do the documentation...
        // inline int32_t reset()
        .def("reset", (int32_t (Class::*)(void)) &Class::reset, "TODO")  // TODO do the documentation...

        /**< sum */
        // T_data sum();
        .def("sum", &Class::sum, "TODO")  // TODO do the documentation...
        // void init_reduceCub();
        .def("init_reduceCub", &Class::init_reduceCub,
             "TODO")  // TODO do the documentation...
        // void reduceCub();
        .def("reduceCub", (void (Class::*)(void)) &Class::reduceCub,
             "TODO")  // TODO do the documentation...
        // void clip(T_data min, T_data max);
        .def("clip", (void (Class::*)(T, T)) &Class::clip, "TODO", py::arg("data_min").none(false),
             py::arg("data_max").none(false))  // TODO do the documentation...

        // /**< transpose */
        // int32_t transpose(CarmaObj<T_data> *source);
        .def("transpose", &Class::transpose, "TODO",
             py::arg("source").none(false))  // TODO do the documentation...
        // //CarmaObj<T_data>& operator= (const CarmaObj<T_data>& obj);

        // /**< Cublas V2 */
        // int32_t imax(int32_t incx);
        .def("aimax", &Class::aimax, "TODO",
             py::arg("incx") = 1)  // TODO do the documentation...
        // int32_t imin(int32_t incx);
        .def("aimin", &Class::aimin, "TODO",
             py::arg("incx") = 1)  // TODO do the documentation...
        // T_data asum(int32_t incx);
        .def("asum", &Class::asum, "TODO",
             py::arg("incx") = 1)  // TODO do the documentation...
        // T_data nrm2(int32_t incx);
        .def("nrm2", &Class::nrm2, "TODO",
             py::arg("incx") = 1)  // TODO do the documentation...
        // T_data dot(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
        .def("dot", &Class::dot, "TODO", py::arg("source").none(false),
             py::arg("incx") = 1,
             py::arg("incy") = 1)  // TODO do the documentation...
        // void scale(T_data alpha, int32_t incx);
        .def("scale", &Class::scale, "TODO", py::arg("scale").none(false),
             py::arg("incx") = 1)  // TODO do the documentation...
        // void swap(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
        .def("swap", &Class::swap, "TODO", py::arg("source").none(false),
             py::arg("incx") = 1,
             py::arg("incy") = 1)  // TODO do the documentation...
        // void copy(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
        .def("copy", &Class::copy, "TODO")  // TODO do the documentation...
        // void axpy(T_data alpha, CarmaObj<T_data> *source, int32_t incx, int32_t
        // incy);
        .def("axpy", &Class::axpy, "TODO", py::arg("alpha"),
             py::arg("source").none(false), py::arg("incx") = 1,
             py::arg("incy") = 1, py::arg("offset") = 0) // TODO do the documentation...
        // void rot(CarmaObj<T_data> *source, int32_t incx, int32_t incy, T_data sc,
        //          T_data ss);
        .def("rot", &Class::rot, "TODO")  // TODO do the documentation...

        // void gemv(char trans, T_data alpha, CarmaObj<T_data> *matA, int32_t lda,
        //           CarmaObj<T_data> *vectx, int32_t incx, T_data beta, int32_t incy);
        .def("gemv",
             [](Class &mat, Class &vectx, T alpha, char op, Class *vecty,
                T beta) {
               if (vecty == nullptr) {
                 int64_t dims[] = {1, 0};
                 if (op == 'N' || op == 'n') {
                   dims[1] = mat.get_dims(1);
                 } else {
                   dims[1] = mat.get_dims(2);
                 }
                 vecty = new Class(mat.get_context(), dims);
                 vecty->reset();
               }
               vecty->gemv(op, alpha, &mat, mat.get_dims(1), &vectx, 1, beta, 1);
               return vecty;
             },
             "this method performs one of the matrix‐vector operations vecty = "
             "alpha * op(mat) * vectx + beta * vecty",
             py::arg("vectx"), py::arg("alpha") = 1, py::arg("op") = 'N',
             py::arg("vecty") = nullptr, py::arg("beta") = 0)  // &Class::gemv)
    // void ger(T_data alpha, CarmaObj<T_data> *vectx, int32_t incx,
    //          CarmaObj<T_data> *vecty, int32_t incy, int32_t lda);

        // void ger(T_data alpha, CarmaObj<T_data> *vectx, int32_t incx,
        //          CarmaObj<T_data> *vecty, int32_t incy, int32_t lda);
        .def("ger",
             [](Class &vectx, Class &vecty, Class *mat, T alpha) {
               std::unique_ptr<Class> ptr_res;
               if (mat == nullptr) {
                 int64_t dims[] = {2, vectx.get_nb_elements(), vecty.get_nb_elements()};
                 mat = new Class(vectx.get_context(), dims);
                 mat->reset();
               }
               mat->ger(alpha, &vectx, 1, &vecty, 1, vectx.get_nb_elements());
               return mat;
             },
             "this method performs the symmetric rank 1 operation A = alpha * "
             "x * y T + A",
             py::arg("vecty"), py::arg("mat") = nullptr,
             py::arg("alpha") = 1)  // &Class::ger)
        // void symv(char uplo, T_data alpha, CarmaObj<T_data> *matA,
        //           int32_t lda, CarmaObj<T_data> *vectx, int32_t incx, T_data beta,
        //           int32_t incy);
        .def("symv",
             [](Class &mat, Class &vectx, T alpha, char uplo, Class *vecty,
                T beta) {
               int32_t lda = mat.get_dims(2);
               if (vecty == nullptr) {
                 int64_t dims[] = {1, lda};
                 vecty = new Class(mat.get_context(), dims);
                 vecty->reset();
               }
               vecty->symv(uplo, alpha, &mat, lda, &vectx, 1, beta, 1);
               return vecty;
             },
             "this method performs one of the matrix‐vector operations vecty = "
             "alpha * mat * vectx + beta * vecty",
             py::arg("vectx"), py::arg("alpha") = 1, py::arg("uplo") = 'l',
             py::arg("vecty") = nullptr, py::arg("beta") = 0)  // &Class::gemv)
        // void gemm(char transa, char transb, T_data alpha, CarmaObj<T_data>
        // *matA,
        //           int32_t lda, CarmaObj<T_data> *matB, int32_t ldb, T_data beta, int32_t
        //           ldc);
        .def("gemm",
             [](Class &matA, Class &matB, char op_a, char op_b, T alpha,
                Class *matC, T beta) /*-> std::unique_ptr<Class>*/ {
               int32_t lda, ldb, ldc, m, n, k;
               if (op_a == 'N' || op_a == 'n') {
                 m = matA.get_dims(1);
                 k = matA.get_dims(2);
               } else {
                 m = matA.get_dims(2);
                 k = matA.get_dims(1);
               }

               if (op_b == 'N' || op_b == 'n') {
                 k = matB.get_dims(1);
                 n = matB.get_dims(2);
               } else {
                 k = matB.get_dims(2);
                 n = matB.get_dims(1);
               }

               if (matC == nullptr) {
                 int64_t dims[] = {2, m, n};
                 matC = new Class(matA.get_context(), dims);
                 matC->reset();
               }

               carma_gemm<T>(matA.get_context()->get_cublas_handle(), op_a, op_b,
                             m, n, k, alpha, matA, matA.get_dims(1), matB,
                             matB.get_dims(1), beta, *matC, matC->get_dims(1));
               return matC;
             },
             "this method performs one of the matrix‐marix operations matC = "
             "alpha * op_a(matA) * op_b(matB) + beta * matC",
             py::arg("matB"), py::arg("op_a") = 'N', py::arg("op_b") = "N",
             py::arg("alpha") = 1, py::arg("matC") = nullptr,
             py::arg("beta") = 0)

        // void symm(char side, char uplo, T_data alpha,
        //           CarmaObj<T_data> *matA, int32_t lda, CarmaObj<T_data> *matB,
        //           int32_t ldb, T_data beta, int32_t ldc);
        .def("symm",
             [](Class &matA, Class &matB, T alpha, Class *matC, T beta,
                char side, char uplo) {
               if (matC == nullptr) {
                 int64_t dims[] = {2, matB.get_dims(1), matB.get_dims(2)};
                 matC = new Class(matA.get_context(), dims);
               }
               carma_symm<T>(matA.get_context()->get_cublas_handle(), side, uplo,
                             matB.get_dims(1), matB.get_dims(2), alpha, matA,
                             matA.get_dims(1), matB, matB.get_dims(1), beta,
                             *matC, matC->get_dims(1));
               return matC;
               // matA.symm(uplo, alpha, &matB, lda, &vectx, 1, beta, 1);
             },
             "this method performs one of the matrix‐marix operations matC = "
             "alpha * matA * matB + beta * C",
             py::arg("matB"), py::arg("alpha") = 1, py::arg("matC") = nullptr,
             py::arg("beta") = 0, py::arg("side") = "l", py::arg("uplo") = "u")

        /*
        template <class T_data>
        cublasStatus_t carma_symm(cublasHandle_t cublas_handle, char side,
                                  char uplo, int32_t m, int32_t n, T_data alpha,
                                  T_data *matA, int32_t lda, T_data *matB, int32_t ldb,
                                  T_data beta, T_data *matC, int32_t ldc);
        */

        // void syrk(char uplo, char transa, T_data alpha,
        //           CarmaObj<T_data> *matA, int32_t lda, T_data beta, int32_t ldc);
        .def("syrk",
             [](Class &matA, char fill, char op, T alpha, Class *matC, T beta) {
               int32_t n, k;
               if (op == 'N' || op == 'n') {
                 n = matA.get_dims(1);
                 k = matA.get_dims(2);
               } else {
                 n = matA.get_dims(2);
                 k = matA.get_dims(1);
               }
               if (matC == nullptr) {
                 int64_t dims[] = {2, n, n};
                 matC = new Class(matA.get_context(), dims);
                 matC->reset();
               }
               carma_syrk<T>(matA.get_context()->get_cublas_handle(), fill, op, n,
                             k, alpha, matA, matA.get_dims(1), beta, *matC,
                             matC->get_dims(1));
               return matC;
             },
             "this method performs the symmetric rank- k update",
             py::arg("fill") = "U", py::arg("op") = 'N', py::arg("alpha") = 1,
             py::arg("matC") = nullptr, py::arg("beta") = 0)
        // void syrkx(char uplo, char transa, T_data alpha,
        //            CarmaObj<T_data> *matA, int32_t lda, CarmaObj<T_data> *matB,
        //            int32_t ldb, T_data beta, int32_t ldc);
        .def("syrkx",
             [](Class &matA, Class &matB, char fill, char op, T alpha,
                Class *matC, T beta) {
               int32_t n, k;
               if (op == 'N' || op == 'n') {
                 n = matA.get_dims(1);
                 k = matA.get_dims(2);
               } else {
                 n = matA.get_dims(2);
                 k = matA.get_dims(1);
               }
               if (matC == nullptr) {
                 int64_t dims[] = {2, n, n};
                 matC = new Class(matA.get_context(), dims);
                 matC->reset();
               }
               carma_syrkx<T>(matA.get_context()->get_cublas_handle(), fill, op,
                              n, k, alpha, matA, matA.get_dims(1), matB,
                              matB.get_dims(1), beta, *matC, matC->get_dims(1));
               return matC;
             },
             "this method performs the symmetric rank- k update",
             py::arg("matB"), py::arg("fill") = "U", py::arg("op") = 'N',
             py::arg("alpha") = 1, py::arg("matC") = nullptr,
             py::arg("beta") = 0)
        // void geam(char transa, char transb, T_data alpha, CarmaObj<T_data>
        // *matA,
        //           int32_t lda, T_data beta, CarmaObj<T_data> *matB, int32_t ldb, int32_t
        //           ldc);
        .def("geam",
             [](Class &matA, Class &matB, char opA, char opB, T alpha,
                Class *matC, T beta) {
               int32_t m, n;
               if (opA == 'N' || opA == 'n') {
                 m = matA.get_dims(1);
                 n = matA.get_dims(2);
               } else {
                 m = matA.get_dims(2);
                 n = matA.get_dims(1);
               }
               if (matC == nullptr) {
                 int64_t dims[] = {2, m, n};
                 matC = new Class(matA.get_context(), dims);
                 matC->reset();
               }
               carma_geam<T>(matA.get_context()->get_cublas_handle(), opA, opB, m,
                             n, alpha, matA, matA.get_dims(1), beta, matB,
                             matB.get_dims(1), *matC, matC->get_dims(1));
               return matC;
             },
             "this method performs the symmetric rank- k update",
             py::arg("matB"), py::arg("opA") = 'N', py::arg("opB") = 'N',
             py::arg("alpha") = 1, py::arg("matC") = nullptr,
             py::arg("beta") = 0)
        .def("dgmm",
             [](Class &matA, Class &vectX, T alpha, char side, Class *matC,
                int32_t incx) {
               if (matC == nullptr) {
                 int64_t dims[] = {2, matA.get_dims(1), matA.get_dims(2)};
                 matC = new Class(matA.get_context(), dims);
                 matC->reset();
               }
               carma_dgmm<T>(matA.get_context()->get_cublas_handle(), side,
                             matA.get_dims(1), matA.get_dims(2), matA,
                             matA.get_dims(1), vectX, incx, *matC,
                             matC->get_dims(1));
               return matC;
             },
             "this method performs one of the matrix‐marix operations matC = "
             "diag(vectX)*matA if side='l'",
             py::arg("vectX"), py::arg("alpha") = 1, py::arg("side") = "r",
             py::arg("matC") = nullptr, py::arg("incx") = 1)

        // /**< Curand */
        .def("is_rng_init", &Class::is_rng_init)
        // int32_t init_prng();
        .def("init_prng", (int32_t (Class::*)()) & Class::init_prng)
        // int32_t init_prng(int64_t seed);
        .def("init_prng", (int32_t (Class::*)(int64_t)) & Class::init_prng)
        // int32_t destroy_prng();
        .def("destroy_prng", &Class::destroy_prng)
        // int32_t prng(T_data *output, char gtype, float alpha, float beta);
        .def("prng", (int32_t (Class::*)(T *, char, float, float)) & Class::prng)
        // int32_t prng(T_data *output, char gtype, float alpha);
        .def("prng", (int32_t (Class::*)(T *, char, float)) & Class::prng)
        // int32_t prng(char gtype, float alpha, float beta);
        .def("prng", (int32_t (Class::*)(char, float, float)) & Class::prng)
        // int32_t prng(char gtype, float alpha);
        .def("prng", (int32_t (Class::*)(char, float)) & Class::prng)
        // int32_t prng(char gtype);
        .def("prng", (int32_t (Class::*)(char)) & Class::prng)

        .def("random",
             [](Class &data, int32_t seed, char gtype) {
               data.init_prng(seed);
               data.prng(gtype);
             },
             py::arg("seed") = 1234, py::arg("j") = 'U')

        .def("random_host",
             [](Class &data, int32_t seed, char gtype) {
               data.init_prng_host(seed);
               data.prng_host(gtype);
             },
             py::arg("seed") = 1234, py::arg("j") = 'U')

        // int32_t prng_montagn( float init_montagn );
        .def("prng_montagn", &Class::prng_montagn)

        // int32_t init_prng_host(int32_t seed);
        .def("init_prng_host", (int32_t (Class::*)(int32_t)) & Class::init_prng_host)
        // int32_t prng_host(char gtype);
        .def("prng_host", (int32_t (Class::*)(char)) & Class::prng_host)
        // int32_t prng_host(char gtype, T_data stddev);
        .def("prng_host", (int32_t (Class::*)(char, T)) & Class::prng_host)
        // int32_t prng_host(char gtype, T_data stddev, T_data alpha);
        .def("prng_host", (int32_t (Class::*)(char, T, T)) & Class::prng_host)
        // int32_t destroy_prng_host();
        .def("destroy_prng_host", &Class::destroy_prng_host)

        .def("fft",
             [](Class &data, Class &dest, int32_t direction) {
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
    // int32_t fillarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
    //                 CarmaDevice *device);
    // template<class T>
    // int32_t fillarray2d2(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
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

    mod.def(
        appendName<T>("potri_").data(),
        [](Class &d_A, Class *d_res) {
          if (d_res == nullptr) {
            carma_potr_inv(&d_A);
          } else {
            d_res->copy_from(d_A, d_A.get_nb_elements());
            carma_potr_inv(d_res);
          }
        },
        py::arg("d_A"), py::arg("d_res") = nullptr);

  }
};
#endif
