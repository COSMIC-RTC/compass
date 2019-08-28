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

    /**< sum */
    // T_data sum();
    // carmaWrapObj.def("sum", &Class::sum,
    //                  "TODO");  // TODO do the documentation...
    // void init_reduceCub();
    carmaWrapObj.def("init_reduceCub", &Class::init_reduceCub,
                     "TODO");  // TODO do the documentation...
    // void reduceCub();
    carmaWrapObj.def("reduceCub", &Class::reduceCub,
                     "TODO");  // TODO do the documentation...
    // void clip(T_data min, T_data max);
    // carmaWrapObj.def(
    //     "clip", &Class::clip, "TODO", py::arg("data_min").none(false),
    //     py::arg("data_max").none(false));  // TODO do the documentation...

    // /**< transpose */
    // int transpose(carma_obj<T_data> *source);
    carmaWrapObj.def(
        "transpose", &Class::transpose, "TODO",
        py::arg("source").none(false));  // TODO do the documentation...
    // //carma_obj<T_data>& operator= (const carma_obj<T_data>& obj);

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
    // T_data dot(carma_obj<T_data> *source, int incx, int incy);
    // carmaWrapObj.def("dot", &Class::dot, "TODO", py::arg("source").none(false),
    //                  py::arg("incx") = 1,
    //                  py::arg("incy") = 1);  // TODO do the documentation...
    // void scale(T_data alpha, int incx);
    // carmaWrapObj.def("scale", &Class::scale, "TODO",
    //                  py::arg("scale").none(false),
    //                  py::arg("incx") = 1);  // TODO do the documentation...
    // void swap(carma_obj<T_data> *source, int incx, int incy);
    carmaWrapObj.def("swap", &Class::swap, "TODO",
                     py::arg("source").none(false), py::arg("incx") = 1,
                     py::arg("incy") = 1);  // TODO do the documentation...
    // void copy(carma_obj<T_data> *source, int incx, int incy);
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
          //  const long *dims = data.getDims();
          //  cufftHandle *handle = data.getPlan();
          //  if(dest == nullptr) {
          //    dest = Class(data.getContext(), dims);
          //  }
          //  carma_initfft(dims, handle, carma_select_plan<T,T>());
          //  carma_fft(data.getData(), dest.getData(), direction, handle);
        },
        py::arg("dest") = nullptr, py::arg("direction") = 1);
    // CU functions clip
    // template<class T_data>
    // void clip_array(T_data *d_data, T_data min, T_data max, int N,
    // carma_device *device);

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
    //                     carma_device *device);
    // template<class T_data>
    // int launch_generic2d(T_data *d_odata, T_data *d_idata, int N1, int N2);

    // CU functions curand
    // int carma_prng_init(int *seed, const int nThreads, const int nBlocks,
    //                     curandState *state);
    // template<class T>
    // int carma_prng_cu(T *results, const int nThreads, const int nBlocks,
    //                   curandState *state, char gtype, int n, float alpha,
    //                   float beta);
    // template<class T>
    // int carma_curand_montagn(curandState *state, T *d_odata, int N,
    // carma_device *device);

    // CU functions fft
    // template<class T_in, class T_out>
    // cufftType carma_select_plan();
    // template<class T_in, class T_out>
    // void carma_initfft(const long *dims_data, cufftHandle *plan, cufftType
    // tPlan); template<class T_in, class T_out> int carma_fft(T_in *input,
    // T_out *output, int dir, cufftHandle plan);

    // CU functions generic
    // template<class T_data>
    // int fillindex(T_data *d_odata, T_data *d_idata, int *indx, int N,
    //               carma_device *device);
    // template<class T_data>
    // int fillvalues(T_data *d_odata, T_data val, int N,
    //               carma_device *device);
    // template<class T>
    // int getarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
    //               carma_device *device);
    // template<class T>
    // int fillarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int
    // N,
    //                 carma_device *device);
    // template<class T>
    // int fillarray2d2(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int
    // N,
    //                 carma_device *device);
    // template<class T>
    // int fill_sym_matrix(char src_uplo, T *d_data, int Ncol, int N,
    //                     carma_device *device);
    // template<class T>
    // int carma_plus(T *d_odata, T elpha, int N, carma_device *device);
    // template<class T>
    // int carma_plusai(T *d_odata, T *i_data, int i, int sgn, int N,
    //                 carma_device *device);

    // CU functions fftconv
    // int fftconv_unpad(float *d_odata, float *d_idata, int fftW, int dataH,
    //                   int dataW, int N, int n, int nim);
    // int carma_initfftconv(caObjS *data_in, caObjS *kernel_in, caObjS
    // *padded_data, caObjC *padded_spectrum, int kernelY, int kernelX);

    // CPP functions fftconv
    // int carma_fftconv(caObjS *data_out, caObjS *padded_data,
    //                   caObjC *padded_spectrum, int kernelY, int kernelX);

    // MAGMA functions

    // template<class T>
    // int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
    //               carma_obj<T> *mod2act, carma_obj<T> *mes2mod);
    // mod.def(appendName<T>("svd_").data(), &carma_svd<T>);

    // TODO after carma_host_obj
    // template<class T>
    // int carma_magma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T>
    // *eigenvals);
    // mod.def( appendName<T>("syevd_").data(), py::overload_cast<char, Class
    // *, ClassHost *>(&carma_magma_syevd<T>)); mod.def(
    // appendName<T>("syevd_").data(), py::overload_cast<char, long, T *, T
    // *>(&carma_magma_syevd<T>));
    mod.def(
        appendName<T>("magma_syevd_").data(),
        [](Class &d_A, ClassHost &eigenvals, Class *d_U, bool computeU) {
          if (d_U == nullptr) {
            if (computeU) {
              carma_magma_syevd('V', &d_A, &eigenvals);
            } else {
              carma_magma_syevd('N', &d_A, &eigenvals);
            }
          } else {
            d_U->copyFrom(d_A, d_A.getNbElem());
            if (computeU) {
              carma_magma_syevd('V', d_U, &eigenvals);
            } else {
              carma_magma_syevd('N', d_U, &eigenvals);
            }
          }
        },
        py::arg("d_A"), py::arg("eigenvals"), py::arg("d_U") = nullptr,
        py::arg("computeU") = true);

    mod.def(
        appendName<T>("syevd_").data(),
        [](Class &d_A, Class &eigenvals, Class *d_U, bool computeU) {
          if (d_U == nullptr) {
            if (computeU) {
              carma_syevd(CUSOLVER_EIG_MODE_VECTOR, &d_A, &eigenvals);
            } else {
              carma_syevd(CUSOLVER_EIG_MODE_NOVECTOR, &d_A, &eigenvals);
            }
          } else {
            d_U->copyFrom(d_A, d_A.getNbElem());
            if (computeU) {
              carma_syevd(CUSOLVER_EIG_MODE_VECTOR, d_U, &eigenvals);
            } else {
              carma_syevd(CUSOLVER_EIG_MODE_NOVECTOR, d_U, &eigenvals);
            }
          }
        },
        py::arg("d_A"), py::arg("eigenvals"), py::arg("d_U") = nullptr,
        py::arg("computeU") = true);
    // template<class T, int method>
    // int carma_magma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T>
    // *eigenvals); template<class T> int carma_magma_syevd_m(long ngpu, char
    // jobz, long N, T *mat, T *eigenvals); template<class T> int
    // carma_magma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
    //                   carma_host_obj<T> *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
    //                   carma_host_obj<T> *eigenvals, carma_host_obj<T> *U);
    // template<class T>
    // int carma_magma_getri(carma_obj<T> *d_iA);
    mod.def(appendName<T>("magma_getri_").data(), &carma_magma_getri<T>);

    // template<class T>
    // int carma_magma_potri(carma_obj<T> *d_iA);
    mod.def(appendName<T>("magma_potri_").data(), &carma_magma_potri<T>);

    // TODO after carma_host_obj
    // template<class T>
    // int carma_magma_potri_m(long num_gpus, carma_host_obj<T> *h_A,
    // carma_obj<T> *d_iA);

    // MAGMA functions (direct access)
    // template<class T>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T, int method>
    // int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
    // template<class T>
    // int carma_magma_syevd_m(long ngpu, char jobz, long N, T *mat, T
    // *eigenvals); template<class T> int carma_magma_potri_m(long num_gpus,
    // long N, T *h_A, T *d_iA);

    // CULA functions
    // template<class T>
    // int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
    //                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod);

    // int snapTransformSize(unsigned int dataSize);
  }
};
#endif
