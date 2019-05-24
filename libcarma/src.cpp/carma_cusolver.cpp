#include <carma_cusolver.h>

#include <type_list.hpp>

#ifdef CAN_DO_HALF
using TypeListObj = GenericTypeList<uint16_t, int, unsigned int, float, double,
                                    half, cuFloatComplex,
                                    cuDoubleComplex>;  // , tuple_t<float>>;
#else
using TypeListObj =
    GenericTypeList<uint16_t, int, unsigned int, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;
#endif

#define CHECK_CUSOLVER(fct, dealloc)                                    \
  do {                                                                  \
    cusolverStatus_t info = fct;                                        \
    if (info != CUSOLVER_STATUS_SUCCESS) {                              \
      printf("%s@%d cuSolver returned error %d.\n", __FILE__, __LINE__, \
             (int)info);                                                \
      dealloc;                                                          \
      return EXIT_FAILURE;                                              \
    }                                                                   \
  } while (0)

template <class T, typename Fn_bufferSize, typename Fn>
int carma_syevd_gen(Fn_bufferSize const &ptr_syevd_gpu_bufferSize,
                    Fn const &ptr_syevd_gpu, cusolverEigMode_t jobz, int N,
                    int lda, T *d_mat, T *d_eigenvals) {
  cusolverDnHandle_t cusolverH = NULL;

  int *devInfo = NULL;
  T *d_work = NULL;
  int lwork = 0;

  int info_gpu = 0;

  // step 1: create cusolver/cublas handle
  CHECK_CUSOLVER(cusolverDnCreate(&cusolverH), nullptr);

  // step 2: query working space of syevd
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  CHECK_CUSOLVER(ptr_syevd_gpu_bufferSize(cusolverH, jobz, uplo, N, d_mat, lda,
                                          d_eigenvals, &lwork),
                 nullptr);
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int));

  // step 3: compute spectrum
  CHECK_CUSOLVER(ptr_syevd_gpu(cusolverH, jobz, uplo, N, d_mat, lda,
                               d_eigenvals, d_work, lwork, devInfo),
                 nullptr);

  if (devInfo) cudaFree(devInfo);
  if (d_work) cudaFree(d_work);
  cudaDeviceSynchronize();
  return EXIT_SUCCESS;
}

template <class T>
int carma_syevd(cusolverEigMode_t jobz, long N, T *mat, T *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(cusolverEigMode_t jobz, long N, float *mat,
                       float *eigenvals) {
  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat, eigenvals);
}

template <>
int carma_syevd<double>(cusolverEigMode_t jobz, long N, double *mat,
                        double *eigenvals) {
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat, eigenvals);
}

template <class T>
int carma_syevd(cusolverEigMode_t jobz, carma_obj<T> *mat,
                carma_obj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(cusolverEigMode_t jobz, caObjS *mat, caObjS *eigenvals) {
  long N = mat->getDims(1);

  if (N != mat->getDims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat->getData(), eigenvals->getData());
}

template <>
int carma_syevd<double>(cusolverEigMode_t jobz, caObjD *mat,
                        caObjD *eigenvals) {
  long N = mat->getDims(1);

  if (N != mat->getDims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat->getData(), eigenvals->getData());
}

struct CarmaCusolverInterfacer {
  template <typename T_data>
  static void call() {
    force_keep((int (*)(cusolverEigMode_t, long, T_data *, T_data *)) &
               carma_syevd<T_data>);
    force_keep(
        (int (*)(cusolverEigMode_t, carma_obj<T_data> *, carma_obj<T_data> *)) &
        carma_syevd<T_data>);
  }
};

void declare_carma() { apply<CarmaCusolverInterfacer, TypeListObj>(); }
