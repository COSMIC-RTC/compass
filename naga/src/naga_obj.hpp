#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <carma.h>

namespace py = pybind11;

template <typename T>
void declare_naga_obj(py::module &mod, std::string suffix)
{

  using Class = carma_obj<T>;

  py::class_<Class>(mod, ("naga_obj_" + suffix).c_str(), py::buffer_protocol())
      .def(py::init([](carma_context &c, const py::array_t<T> &data) {
        int ndim = data.ndim() + 1;
        std::vector<long> data_dims(ndim);
        data_dims[0] = data.ndim();
        copy(data.shape(), data.shape()+data.ndim(), begin(data_dims) + 1);
        return std::unique_ptr<Class>(new Class(&c, data_dims.data(), (const T *)data.data()));
      }),
      py::arg("context").none(false),
      py::arg("data").none(false))

      // .def(py::init([](carma_context &c, Class &data) {
      //   return std::unique_ptr<Class>(new Class(&c, &data));
      // }))

      // .def_buffer([](Class &frame) -> py::array_t<T> {

      //   py::array_t<T> array(frame.getNbElem());
      //   auto info = array.request();

      //   const long *dims = frame.getDims();
      //   std::vector<ssize_t> shape(dims[0]);
      //   std::vector<ssize_t> strides(dims[0]);
      //   ssize_t stride = sizeof(T);
      //   for (ssize_t dim(dims[0] - 1); dim >= 0; --dim)
      //   {
      //     // cerr << dim <<  endl;
      //     shape[dim] = dims[dim + 1];
      //     strides[dim] = stride;
      //     stride *= shape[dim];
      //   }

      //   info.ndim = dims[0];
      //   info.shape = shape;
      //   info.strides = strides;

      //   return array;
      // })

      .def_buffer([](Class &frame) -> py::buffer_info {
        frame.sync_h_data();
        const long *dims = frame.getDims();
        std::vector<ssize_t> shape(dims[0]);
        std::vector<ssize_t> strides(dims[0]);
        ssize_t stride = sizeof(T);
        for (ssize_t dim(dims[0] - 1); dim >= 0; --dim)
        {
          // cerr << dim <<  endl;
          shape[dim] = dims[dim + 1];
          strides[dim] = stride;
          stride *= shape[dim];
        }

        return py::buffer_info(frame.get_h_data(),
                               sizeof(T),
                               py::format_descriptor<T>::format(),
                               dims[0],
                               shape,
                               strides);
      })

      .def("__repr__", &Class::toString)

      // int get_nbStreams()
      .def("get_nbStreams", &Class::get_nbStreams)
      // int add_stream()
      .def("add_stream", (int (Class::*)()) & Class::add_stream)
      // int add_stream(int nb)
      .def("add_stream", (int (Class::*)(int)) & Class::add_stream)
      // int del_stream()
      .def("del_stream", (int (Class::*)()) & Class::del_stream)
      // int del_stream(int nb)
      .def("del_stream", (int (Class::*)(int)) & Class::del_stream)
      // int wait_stream(int stream)
      .def("wait_stream", &Class::wait_stream)
      // int wait_all_streams()
      .def("wait_all_streams", &Class::wait_all_streams)

      // const long *getDims()
      // .def("getDims", [](Class &frame) -> py::buffer_info {
      //   long ndim = frame.getDims(0);
      //   std::vector<ssize_t> shape(1);
      //   std::vector<ssize_t> strides(1);
      //   shape[0] = ndim;
      //   strides[0] = sizeof(T);
      //   return py::buffer_info(&frame.getDims(1),
      //                       sizeof(long),
      //                       py::format_descriptor<long>::format(),
      //                       1,
      //                       shape,
      //                       strides);

      // })
      // int getNbElem()
      .def("getNbElem", &Class::getNbElem)
      // carma_context* getContext()
      .def("getContext", &Class::getContext)
      // int getDevice()
      .def("getDevice", &Class::getDevice)

      // int host2device(T_data *data);
      .def("host2device", [](Class &c, py::array_t<T> &data) {
        c.host2device((const T*)data.data());
      })
      // int device2host(T_data *data);
      .def("device2host", [](Class &c, py::array_t<T> &data) {
        c.device2host((T*)data.mutable_data());
      })

      // int copyInto(T_data *data, int nb_elem);
      // int copyFrom(T_data *data, int nb_elem);

      // inline int reset()
      .def("reset", &Class::reset)

      /**< sum */
      // T_data sum();
      .def("sum", &Class::sum)
      // void clip(T_data min, T_data max);
      .def("clip", &Class::clip)

      // /**< transpose */
      // int transpose(carma_obj<T_data> *source);
      .def("transpose", &Class::transpose)
      // //carma_obj<T_data>& operator= (const carma_obj<T_data>& obj);

      // /**< Cublas V2 */
      // int imax(int incx);
      .def("imax", &Class::imax)
      // int imin(int incx);
      .def("imin", &Class::imin)
      // T_data asum(int incx);
      .def("asum", &Class::asum)
      // T_data nrm2(int incx);
      .def("nrm2", &Class::nrm2)
      // T_data dot(carma_obj<T_data> *source, int incx, int incy);
      .def("dot", &Class::dot)
      // void scale(T_data alpha, int incx);
      .def("scale", &Class::scale)
      // void swap(carma_obj<T_data> *source, int incx, int incy);
      .def("swap", &Class::swap)
      // void copy(carma_obj<T_data> *source, int incx, int incy);
      .def("copy", &Class::copy)
      // void axpy(T_data alpha, carma_obj<T_data> *source, int incx, int incy);
      .def("axpy", &Class::axpy)
      // void rot(carma_obj<T_data> *source, int incx, int incy, T_data sc,
      //          T_data ss);
      .def("rot", &Class::rot)

      // void gemv(char trans, T_data alpha, carma_obj<T_data> *matA, int lda,
      //           carma_obj<T_data> *vectx, int incx, T_data beta, int incy);
      .def("gemv", [](Class &mat, Class &vectx, Class &vecty, char op, T alpha, T beta) {
        int lda;
        if (op == 'N' || op == 'n')
        {
          lda = mat.getDims(2);
        }
        else
        {
          lda = mat.getDims(1);
        }
        vecty.gemv(op, alpha, &mat, lda, &vectx, 1, beta, 1);
      },
           "this method performs one of the matrix‐vector operations vecty = alpha * op(mat) * vectx + beta * vecty", py::arg("vectx"), py::arg("vecty"), py::arg("op") = 'N', py::arg("alpha") = 1, py::arg("beta") = 0) // &Class::gemv)
      .def("gemv", [](Class &mat, Class &vectx, char op, T alpha) -> std::unique_ptr<Class> {
        int lda, m, n;
        char op_cublas;
        if (op == 'N' || op == 'n')
        {
          m = mat.getDims(2);
          n = mat.getDims(1);
          lda = m;
          op_cublas = 't';
        }
        else
        {
          m = mat.getDims(1);
          n = mat.getDims(2);
          lda = m;
          op_cublas = 'n';
        }
        long dims[] = {1, n};
        std::unique_ptr<Class> vecty(new Class(mat.getContext(), dims));
        T beta = 0;
        carma_gemv(mat.getContext()->get_cublasHandle(), op_cublas, m, n, alpha, mat.getData(), lda, vectx.getData(), 1, beta, vecty->getData(), 1);
        // vecty->gemv(op, alpha, &mat, lda, &vectx, 1, beta, 1);
        return vecty;
      },
           "this method performs one of the matrix‐vector operations vecty = alpha * op(mat) * vectx", py::arg("vectx"), py::arg("op") = 'N', py::arg("alpha") = 1) // &Class::gemv)
      // void ger(T_data alpha, carma_obj<T_data> *vectx, int incx,
      //          carma_obj<T_data> *vecty, int incy, int lda);
      .def("ger", [](Class &vectx, Class &vecty, Class *mat, T alpha) -> std::unique_ptr<Class> {
        std::unique_ptr<Class> ptr_res;
        if (mat == nullptr)
        {
          long dims[] = {2, vectx.getNbElem(), vecty.getNbElem()};
          ptr_res = std::unique_ptr<Class>(new Class(vectx.getContext(), dims));
        }
        else
        {
          ptr_res = std::unique_ptr<Class>(new Class(mat));
        }
        ptr_res->ger(alpha, &vectx, 1, &vecty, 1, vectx.getNbElem());
        return ptr_res;
      },
           "this method performs the symmetric rank 1 operation A = alpha * x * y T + A", py::arg("vecty"), py::arg("mat") = nullptr, py::arg("alpha") = 1) // &Class::ger)
      // void symv(cublasFillMode_t uplo, T_data alpha, carma_obj<T_data> *matA,
      //           int lda, carma_obj<T_data> *vectx, int incx, T_data beta, int incy);
      .def("symv", [](Class &mat, Class &vectx, Class &vecty, T alpha, T beta) {
        int lda = mat.getDims(2);
        vecty.symv(CUBLAS_FILL_MODE_LOWER, alpha, &mat, lda, &vectx, 1, beta, 1);
      },
           "this method performs one of the matrix‐vector operations vecty = alpha * mat * vectx + beta * vecty", py::arg("vectx"), py::arg("vecty"), py::arg("alpha") = 1, py::arg("beta") = 0) // &Class::gemv)
      .def("symv", [](Class &mat, Class &vectx, T alpha) -> std::unique_ptr<Class> {
        int lda = mat.getDims(2);
        long dims[] = {1, lda};
        std::unique_ptr<Class> vecty(new Class(mat.getContext(), dims));
        T beta = 0;
        vecty->symv(CUBLAS_FILL_MODE_LOWER, alpha, &mat, lda, &vectx, 1, beta, 1);
        return vecty;
      },
           "this method performs one of the matrix‐vector operations vecty = alpha * mat * vectx", py::arg("vectx"), py::arg("alpha") = 1) // &Class::gemv)

      // void gemm(char transa, char transb, T_data alpha, carma_obj<T_data> *matA,
      //           int lda, carma_obj<T_data> *matB, int ldb, T_data beta, int ldc);
      .def("gemm",  [](Class &matA, Class &matB, char op_a, char op_b, T alpha, Class &matC, T beta
      // , int m, int n, int k, int lda, int ldb, int ldc
      ) /*-> std::unique_ptr<Class>*/  {
        char op_cublas_A='t', op_cublas_B='t';

// 1 2 4 4 2 1  seems GOOD!
// 128 256 512 512 256 128

        int lda, ldb, ldc, m, n, k;
        if (op_a == 'N' || op_a == 'n')
        {
          m = matA.getDims(1); // 1
          n = matB.getDims(2); // 2
          k = matA.getDims(2); // 2
          lda = k; // m
          ldb = n; // k
          ldc = m; // m
          op_cublas_A = 't';
          op_cublas_B = 't';
        }
        else
        {
          m = matA.getDims(1);
          n = matA.getDims(2);
          lda = matA.getDims(1);
          op_cublas_A = 'n';
        }

        long dims[] = {2, n, m};
        std::unique_ptr<Class> matTMP(new Class(matA.getContext(),dims));
        carma_geam<T>(matA.getContext()->get_cublasHandle(), 'n', 't',m, n, 0, matTMP->getData(), m, 1, matC.getData(), n, matTMP->getData(), m);
        // T beta = 0;
        carma_gemm<T>(matA.getContext()->get_cublasHandle(), op_cublas_B, op_cublas_A, m, n, k, alpha, matA.getData(), lda, matB.getData(), ldb, beta, matTMP->getData(), ldc);
        carma_geam<T>(matA.getContext()->get_cublasHandle(), 'n', 't',n,m, 0, matC.getData(), n, 1, matTMP->getData(), m, matC.getData(), n);
        // return matC;
      })

// cublasStatus_t carma_gemm<float>(cublasHandle_t cublas_handle, char transa,
//                                  char transb, int m, int n, int k, float alpha,
//                                  float *matA, int lda, float *matB, int ldb,
//                                  float beta, float *matC, int ldc) {


      // void symm(cublasSideMode_t side, cublasFillMode_t uplo, T_data alpha,
      //           carma_obj<T_data> *matA, int lda, carma_obj<T_data> *matB, int ldb,
      //           T_data beta, int ldc);
      .def("symm", &Class::symm)
      // void syrk(cublasFillMode_t uplo, char transa, T_data alpha,
      //           carma_obj<T_data> *matA, int lda, T_data beta, int ldc);
      .def("syrk", &Class::syrk)
      // void syrkx(cublasFillMode_t uplo, char transa, T_data alpha,
      //            carma_obj<T_data> *matA, int lda, carma_obj<T_data> *matB, int ldb,
      //            T_data beta, int ldc);
      .def("syrkx", &Class::syrkx)
      // void geam(char transa, char transb, T_data alpha, carma_obj<T_data> *matA,
      //           int lda, T_data beta, carma_obj<T_data> *matB, int ldb, int ldc);
      .def("geam", &Class::geam)
      // void dgmm(cublasSideMode_t side, carma_obj<T_data> *matA, int lda,
      //           carma_obj<T_data> *vectx, int incx, int ldc);
      .def("dgmm", &Class::dgmm)

      // /**< Curand */
      .def("is_rng_init", &Class::is_rng_init)
      // int init_prng();
      .def("init_prng", (int (Class::*)()) & Class::init_prng)
      // int init_prng(long seed);
      .def("init_prng", (int (Class::*)(long)) & Class::init_prng)
      // int destroy_prng();
      .def("destroy_prng", &Class::destroy_prng)
      // int prng(T_data *output, char gtype, float alpha, float beta);
      .def("prng", (int (Class::*)(T *, char, float, float)) & Class::prng)
      // int prng(T_data *output, char gtype, float alpha);
      .def("prng", (int (Class::*)(T *, char, float)) & Class::prng)
      // int prng(char gtype, float alpha, float beta);
      .def("prng", (int (Class::*)(char, float, float)) & Class::prng)
      // int prng(char gtype, float alpha);
      .def("prng", (int (Class::*)(char, float)) & Class::prng)
      // int prng(char gtype);
      .def("prng", (int (Class::*)(char)) & Class::prng)

      .def("random", [](Class &data, int seed, char gtype) {
        data.init_prng(seed);
        data.prng(gtype);
      },
           py::arg("seed") = 1234, py::arg("j") = 'U')

      // int prng_montagn( float init_montagn );
      .def("prng_montagn", &Class::prng_montagn)

      // int init_prng_host(int seed);
      .def("init_prng_host", (int (Class::*)(int)) & Class::init_prng_host)
      // int prng_host(char gtype);
      .def("prng_host", (int (Class::*)(char)) & Class::prng_host)
      // int prng_host(char gtype, T_data stddev);
      .def("prng_host", (int (Class::*)(char, T)) & Class::prng_host)
      // int prng_host(char gtype, T_data stddev, T_data alpha);
      .def("prng_host", (int (Class::*)(char, T, T)) & Class::prng_host)
      // int destroy_prng_host();
      .def("destroy_prng_host", &Class::destroy_prng_host)

      ;
  // CU functions clip
  // template<class T_data>
  // void clip_array(T_data *d_data, T_data min, T_data max, int N, carma_device *device);

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
  //                   curandState *state, char gtype, int n, float alpha, float beta);
  // template<class T>
  // int carma_curand_montagn(curandState *state, T *d_odata, int N, carma_device *device);

  // CU functions fft
  // template<class T_in, class T_out>
  // cufftType carma_select_plan();
  // template<class T_in, class T_out>
  // void carma_initfft(const long *dims_data, cufftHandle *plan, cufftType tPlan);
  // template<class T_in, class T_out>
  // int carma_fft(T_in *input, T_out *output, int dir, cufftHandle plan);

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
  // int fillarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
  //                 carma_device *device);
  // template<class T>
  // int fillarray2d2(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
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
  // int carma_initfftconv(caObjS *data_in, caObjS *kernel_in, caObjS *padded_data,
  // caObjC *padded_spectrum, int kernelY, int kernelX);

  // CPP functions fftconv
  // int carma_fftconv(caObjS *data_out, caObjS *padded_data,
  //                   caObjC *padded_spectrum, int kernelY, int kernelX);

  // MAGMA functions

  // template<class T>
  // int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
  //               carma_obj<T> *mod2act, carma_obj<T> *mes2mod);
  mod.def(("naga_svd_" + suffix).c_str(), &carma_svd<T>);

  // TODO after carma_host_obj
  // template<class T>
  // int carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals);
  // mod.def( ("naga_syevd_" + suffix).c_str(), &carma_syevd<T>);

  // template<class T, int method>
  // int carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals);
  // template<class T>
  // int carma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
  // template<class T>
  // int carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
  //                   carma_host_obj<T> *eigenvals);
  // template<class T>
  // int carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
  //                   carma_host_obj<T> *eigenvals, carma_host_obj<T> *U);
  // template<class T>
  // int carma_getri(carma_obj<T> *d_iA);
  mod.def(("naga_getri_" + suffix).c_str(), &carma_getri<T>);

  // template<class T>
  // int carma_potri(carma_obj<T> *d_iA);
  mod.def(("naga_potri_" + suffix).c_str(), &carma_potri<T>);

  // TODO after carma_host_obj
  // template<class T>
  // int carma_potri_m(long num_gpus, carma_host_obj<T> *h_A, carma_obj<T> *d_iA);

  // MAGMA functions (direct access)
  // template<class T>
  // int carma_syevd(char jobz, long N, T *mat, T *eigenvals);
  // template<class T, int method>
  // int carma_syevd(char jobz, long N, T *mat, T *eigenvals);
  // template<class T>
  // int carma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
  // template<class T>
  // int carma_potri_m(long num_gpus, long N, T *h_A, T *d_iA);

  // CULA functions
  // template<class T>
  // int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
  //                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod);

  // int snapTransformSize(unsigned int dataSize);
}
