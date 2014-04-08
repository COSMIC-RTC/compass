#include <carma_cublas.h>
#include <carma_obj.h>

#define carma_checkCublasStatus(status) carma_checkCublasStatus_v2(status, __LINE__, __FILE__)

cublasStatus_t carma_checkCublasStatus_v2(cublasStatus_t status, int line,
    string file)
    /**< Generic CUBLAS check status routine */
    {
  switch (status) {
  case CUBLAS_STATUS_SUCCESS:
    return status;
  case CUBLAS_STATUS_NOT_INITIALIZED:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_NOT_INITIALIZED !!!!!\n";
    break;
  case CUBLAS_STATUS_ALLOC_FAILED:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_ALLOC_FAILED !!!!!\n";
    break;
  case CUBLAS_STATUS_INVALID_VALUE:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_ALLOC_FAILED !!!!!\n";
    break;
  case CUBLAS_STATUS_ARCH_MISMATCH:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_ARCH_MISMATCH !!!!!\n";
    break;
  case CUBLAS_STATUS_MAPPING_ERROR:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_MAPPING_ERROR !!!!!\n";
    break;
  case CUBLAS_STATUS_EXECUTION_FAILED:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_EXECUTION_FAILED !!!!!\n";
    break;
  case CUBLAS_STATUS_INTERNAL_ERROR:
    cerr << "!!!!! Cublas error : CUBLAS_STATUS_INTERNAL_ERROR !!!!!\n";
    break;
  }
  cerr << "!!!!! Cublas error in " << file << "@" << line << " !!!!!\n";
  return status;
}

cublasStatus_t carma_initCublas(cublasHandle_t *cublas_handle)
/**< Generic CUBLAS init routine */
{
  return carma_checkCublasStatus(cublasCreate(cublas_handle));
}

cublasStatus_t carma_shutdownCublas(cublasHandle_t cublas_handle)
/**< Generic CUBLAS shutdown routine */
{
  return carma_checkCublasStatus(cublasDestroy(cublas_handle));
}

cublasOperation_t carma_char2cublasOperation(char operation) {
  switch (operation) {
  case 't':
    return CUBLAS_OP_T;
  case 'c':
    return CUBLAS_OP_C;
  default:
    return CUBLAS_OP_N;
  }
}

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____  
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T_data>
int
carma_wheremax(cublasHandle_t cublas_handle, int n, T_data *vect, int incx);
/**< Generic template for Iamax executable selection */
template<>
int carma_wheremax<float>(cublasHandle_t cublas_handle, int n, float *vect,
    int incx) {
  int result = 0;
  carma_checkCublasStatus(cublasIsamax(cublas_handle,n, vect, incx,&result));
  return result;
}
template<>
int carma_wheremax<double>(cublasHandle_t cublas_handle, int n, double *vect,
    int incx) {
  int result = 0;
  carma_checkCublasStatus(cublasIdamax(cublas_handle,n, vect, incx,&result));
  return result;
}

/** These templates are used to select the proper Iamin executable from T_data*/
template<class T_data>
int
carma_wheremin(cublasHandle_t cublas_handle, int n, T_data *vect, int incx);
/**< Generic template for Iamin executable selection */
template<>
int carma_wheremin<float>(cublasHandle_t cublas_handle, int n, float *vect,
    int incx) {
  int result = 0;
  carma_checkCublasStatus(cublasIsamin(cublas_handle,n, vect, incx,&result));
  return result;
}
template<>
int carma_wheremin<double>(cublasHandle_t cublas_handle, int n, double *vect,
    int incx) {
  int result = 0;
  carma_checkCublasStatus(cublasIdamin(cublas_handle,n, vect, incx,&result));
  return result;
}

/** These templates are used to select the proper asum executable from T_data*/
template<class T_data>
T_data
carma_getasum(cublasHandle_t cublas_handle, int n, T_data *vect, int incx);
/**< Generic template for asum executable selection */
template<>
float carma_getasum<float>(cublasHandle_t cublas_handle, int n, float *vect,
    int incx) {
  float result = 0;
  carma_checkCublasStatus(cublasSasum(cublas_handle,n, vect, incx,&result));
  return result;
}
template<>
double carma_getasum<double>(cublasHandle_t cublas_handle, int n, double *vect,
    int incx) {
  double result = 0;
  carma_checkCublasStatus(cublasDasum(cublas_handle,n, vect, incx,&result));
  return result;
}

/** These templates are used to select the proper axpy executable from T_data*/
template<class T_data>
cublasStatus_t
carma_axpy(cublasHandle_t cublas_handle, int n, T_data alpha, T_data *vectx,
    int incx, T_data *vecty, int incy);
/**< Generic template for axpy executable selection */
template<>
cublasStatus_t carma_axpy<float>(cublasHandle_t cublas_handle, int n,
    float alpha, float *vectx, int incx, float *vecty, int incy) {
  return carma_checkCublasStatus(cublasSaxpy(cublas_handle,n, &alpha, vectx,incx,vecty,incy));
}
template<>
cublasStatus_t carma_axpy<double>(cublasHandle_t cublas_handle, int n,
    double alpha, double *vectx, int incx, double *vecty, int incy) {
  return carma_checkCublasStatus(cublasDaxpy(cublas_handle,n, &alpha, vectx,incx,vecty,incy));
}

/** These templates are used to select the proper copy executable from T_data*/
template<class T_data>
cublasStatus_t
carma_copy(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
    T_data *vecty, int incy);
/**< Generic template for copy executable selection */
template<>
cublasStatus_t carma_copy<float>(cublasHandle_t cublas_handle, int n,
    float *vectx, int incx, float *vecty, int incy) {
  return carma_checkCublasStatus(cublasScopy(cublas_handle,n, vectx,incx,vecty,incy));
}
template<>
cublasStatus_t carma_copy<double>(cublasHandle_t cublas_handle, int n,
    double *vectx, int incx, double *vecty, int incy) {
  return carma_checkCublasStatus(cublasDcopy(cublas_handle,n,vectx,incx,vecty,incy));
}

/** These templates are used to select the proper dot executable from T_data*/
template<class T_data>
T_data
carma_dot(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
    T_data *vecty, int incy);
/**< Generic template for dot executable selection */
template<>
float carma_dot<float>(cublasHandle_t cublas_handle, int n, float *vectx,
    int incx, float *vecty, int incy) {
  float result = 0;
  carma_checkCublasStatus(
      cublasSdot(cublas_handle,n, vectx,incx,vecty,incy,&result));
  return result;
}
template<>
double carma_dot<double>(cublasHandle_t cublas_handle, int n, double *vectx,
    int incx, double *vecty, int incy) {
  double result = 0;
  carma_checkCublasStatus(
      cublasDdot(cublas_handle,n, vectx,incx,vecty,incy,&result));
  return result;
}

/** These templates are used to select the proper nrm2 executable from T_data*/
template<class T_data>
T_data
carma_nrm2(cublasHandle_t cublas_handle, int n, T_data *vect, int incx);
/**< Generic template for nrm2 executable selection */
template<>
float carma_nrm2<float>(cublasHandle_t cublas_handle, int n, float *vect,
    int incx) {
  float result = 0;
  carma_checkCublasStatus(cublasSnrm2(cublas_handle,n, vect,incx,&result));
  return result;
}
template<>
double carma_nrm2<double>(cublasHandle_t cublas_handle, int n, double *vect,
    int incx) {
  double result = 0;
  carma_checkCublasStatus(cublasDnrm2(cublas_handle,n, vect,incx,&result));
  return result;
}

/** These templates are used to select the proper rot executable from T_data*/
template<class T_data>
cublasStatus_t
carma_rot(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
    T_data *vecty, int incy, T_data sc, T_data ss);
/**< Generic template for rot executable selection */
template<>
cublasStatus_t carma_rot<float>(cublasHandle_t cublas_handle, int n,
    float *vectx, int incx, float *vecty, int incy, float sc, float ss) {
  return carma_checkCublasStatus(cublasSrot(cublas_handle,n,vectx,incx,vecty,incy,&sc,&ss));
}
template<>
cublasStatus_t carma_rot<double>(cublasHandle_t cublas_handle, int n,
    double *vectx, int incx, double *vecty, int incy, double sc, double ss) {
  return carma_checkCublasStatus(cublasDrot(cublas_handle,n,vectx,incx,vecty,incy,&sc,&ss));
}

/** These templates are used to select the proper scal executable from T_data*/
template<class T_data>
cublasStatus_t
carma_scal(cublasHandle_t cublas_handle, int n, T_data alpha, T_data *vectx,
    int incx);
/**< Generic template for scal executable selection */
template<>
cublasStatus_t carma_scal<float>(cublasHandle_t cublas_handle, int n,
    float alpha, float *vectx, int incx) {
  return carma_checkCublasStatus(cublasSscal(cublas_handle,n,&alpha,vectx,incx));
}
template<>
cublasStatus_t carma_scal<double>(cublasHandle_t cublas_handle, int n,
    double alpha, double *vectx, int incx) {
  return carma_checkCublasStatus(cublasDscal(cublas_handle,n,&alpha,vectx,incx));
}

/** These templates are used to select the proper swap executable from T_data*/
template<class T_data>
cublasStatus_t
carma_swap(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
    T_data *vecty, int incy);
/**< Generic template for swap executable selection */
template<>
cublasStatus_t carma_swap<float>(cublasHandle_t cublas_handle, int n,
    float *vectx, int incx, float *vecty, int incy) {
  return carma_checkCublasStatus(cublasSswap(cublas_handle,n,vectx,incx,vecty,incy));
}
template<>
cublasStatus_t carma_swap<double>(cublasHandle_t cublas_handle, int n,
    double *vectx, int incx, double *vecty, int incy) {
  return carma_checkCublasStatus(cublasDswap(cublas_handle,n,vectx,incx,vecty,incy));
}

/** These templates are used to select the proper gemv executable from T_data*/
template<class T_data>
cublasStatus_t
carma_gemv(cublasHandle_t cublas_handle, char trans, int m, int n, T_data alpha,
    T_data *matA, int lda, T_data *vectx, int incx, T_data beta, T_data *vecty,
    int incy);
/**< Generic template for gemv executable selection */
template<>
cublasStatus_t carma_gemv<float>(cublasHandle_t cublas_handle, char trans,
    int m, int n, float alpha, float *matA, int lda, float *vectx, int incx,
    float beta, float *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublasOperation(trans);
  return carma_checkCublasStatus(cublasSgemv(cublas_handle,trans2,m,n,&alpha,matA,lda,vectx,incx,&beta,vecty,incy));
}
template<>
cublasStatus_t carma_gemv<double>(cublasHandle_t cublas_handle, char trans,
    int m, int n, double alpha, double *matA, int lda, double *vectx, int incx,
    double beta, double *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublasOperation(trans);
  return carma_checkCublasStatus(cublasDgemv(cublas_handle,trans2,m,n,&alpha,matA,lda,vectx,incx,&beta,vecty,incy));
}

/** These templates are used to select the proper ger executable from T_data*/
template<class T_data>
cublasStatus_t
carma_ger(cublasHandle_t cublas_handle, int m, int n, T_data alpha,
    T_data *vectx, int incx, T_data *vecty, int incy, T_data *matA, int lda);
/**< Generic template for ger executable selection */
template<>
cublasStatus_t carma_ger<float>(cublasHandle_t cublas_handle, int m, int n,
    float alpha, float *vectx, int incx, float *vecty, int incy, float *matA,
    int lda) {
  return carma_checkCublasStatus(cublasSger(cublas_handle,m,n,&alpha,vectx,incx,vecty,incy,matA,lda));
}
template<>
cublasStatus_t carma_ger<double>(cublasHandle_t cublas_handle, int m, int n,
    double alpha, double *vectx, int incx, double *vecty, int incy,
    double *matA, int lda) {
  return carma_checkCublasStatus(cublasDger(cublas_handle,m,n,&alpha,vectx,incx,vecty,incy,matA,lda));
}

/** These templates are used to select the proper symv executable from T_data*/
template<class T_data>
cublasStatus_t
carma_symv(cublasHandle_t cublas_handle, cublasFillMode_t uplo, int n,
    T_data alpha, T_data *matA, int lda, T_data *vectx, int incx, T_data beta,
    T_data *vecty, int incy);
/**< Generic template for symv executable selection */
template<>
cublasStatus_t carma_symv<float>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, int n, float alpha, float *matA, int lda,
    float *vectx, int incx, float beta, float *vecty, int incy) {
  return carma_checkCublasStatus(cublasSsymv(cublas_handle,uplo,n,&alpha,matA,lda,vectx,incx,&beta,vecty,incy));
}
template<>
cublasStatus_t carma_symv<double>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, int n, double alpha, double *matA, int lda,
    double *vectx, int incx, double beta, double *vecty, int incy) {
  return carma_checkCublasStatus(cublasDsymv(cublas_handle,uplo,n,&alpha,matA,lda,vectx,incx,&beta,vecty,incy));
}

/** These templates are used to select the proper gemm executable from T_data*/
template<class T_data>
cublasStatus_t
carma_gemm(cublasHandle_t cublas_handle, char transa, char transb, int m, int n,
    int k, T_data alpha, T_data *matA, int lda, T_data *matB, int ldb,
    T_data beta, T_data *matC, int ldc);
/**< Generic template for gemm executable selection */
template<>
cublasStatus_t carma_gemm<float>(cublasHandle_t cublas_handle, char transa,
    char transb, int m, int n, int k, float alpha, float *matA, int lda,
    float *matB, int ldb, float beta, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  cublasOperation_t transb2 = carma_char2cublasOperation(transb);
  return carma_checkCublasStatus(cublasSgemm(cublas_handle,transa2,transb2,m,n,k,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}
template<>
cublasStatus_t carma_gemm<double>(cublasHandle_t cublas_handle, char transa,
    char transb, int m, int n, int k, double alpha, double *matA, int lda,
    double *matB, int ldb, double beta, double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  cublasOperation_t transb2 = carma_char2cublasOperation(transb);
  return carma_checkCublasStatus(cublasDgemm(cublas_handle,transa2,transb2,m,n,k,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}

/** These templates are used to select the proper symm executable from T_data*/
template<class T_data>
cublasStatus_t
carma_symm(cublasHandle_t cublas_handle, cublasSideMode_t side,
    cublasFillMode_t uplo, int m, int n, T_data alpha, T_data *matA, int lda,
    T_data *matB, int ldb, T_data beta, T_data *matC, int ldc);
/**< Generic template for symm executable selection */
template<>
cublasStatus_t carma_symm<float>(cublasHandle_t cublas_handle,
    cublasSideMode_t side, cublasFillMode_t uplo, int m, int n, float alpha,
    float *matA, int lda, float *matB, int ldb, float beta, float *matC,
    int ldc) {
  return carma_checkCublasStatus(cublasSsymm(cublas_handle,side,uplo,m,n,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}
template<>
cublasStatus_t carma_symm<double>(cublasHandle_t cublas_handle,
    cublasSideMode_t side, cublasFillMode_t uplo, int m, int n, double alpha,
    double *matA, int lda, double *matB, int ldb, double beta, double *matC,
    int ldc) {
  return carma_checkCublasStatus(cublasDsymm(cublas_handle,side,uplo,m,n,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}

/** These templates are used to select the proper syrk executable from T_data*/
template<class T_data>
cublasStatus_t
carma_syrk(cublasHandle_t cublas_handle, cublasFillMode_t uplo, char transa,
    int n, int k, T_data alpha, T_data *matA, int lda, T_data beta,
    T_data *matC, int ldc);
/**< Generic template for syrk executable selection */
template<>
cublasStatus_t carma_syrk<float>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, char transa, int n, int k, float alpha, float *matA,
    int lda, float beta, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  return carma_checkCublasStatus(cublasSsyrk(cublas_handle,uplo,transa2,n,k,&alpha,matA,lda,&beta,matC,ldc));
}
template<>
cublasStatus_t carma_syrk<double>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, char transa, int n, int k, double alpha,
    double *matA, int lda, double beta, double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  return carma_checkCublasStatus(cublasDsyrk(cublas_handle,uplo,transa2,n,k,&alpha,matA,lda,&beta,matC,ldc));
}

/** These templates are used to select the proper syrkx executable from T_data*/
template<class T_data>
cublasStatus_t
carma_syrkx(cublasHandle_t cublas_handle, cublasFillMode_t uplo, char transa,
    int n, int k, T_data alpha, T_data *matA, int lda, T_data *matB, int ldb,
    T_data beta, T_data *matC, int ldc);
/**< Generic template for syrkx executable selection */
template<>
cublasStatus_t carma_syrkx<float>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, char transa, int n, int k, float alpha, float *matA,
    int lda, float *matB, int ldb, float beta, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  return carma_checkCublasStatus(cublasSsyrkx(cublas_handle,uplo,transa2,n,k,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}
template<>
cublasStatus_t carma_syrkx<double>(cublasHandle_t cublas_handle,
    cublasFillMode_t uplo, char transa, int n, int k, double alpha,
    double *matA, int lda, double *matB, int ldb, double beta, double *matC,
    int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  return carma_checkCublasStatus(cublasDsyrkx(cublas_handle,uplo,transa2,n,k,&alpha,matA,lda,matB,ldb,&beta,matC,ldc));
}

/** These templates are used to select the proper geam executable from T_data*/
template<class T_data>
cublasStatus_t
carma_geam(cublasHandle_t cublas_handle, char transa, char transb, int m, int n,
    T_data alpha, T_data *matA, int lda, T_data beta, T_data *matB, int ldb,
    T_data *matC, int ldc);
/**< Generic template for geam executable selection */
template<>
cublasStatus_t carma_geam<float>(cublasHandle_t cublas_handle, char transa,
    char transb, int m, int n, float alpha, float *matA, int lda, float beta,
    float *matB, int ldb, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  cublasOperation_t transb2 = carma_char2cublasOperation(transb);
  return carma_checkCublasStatus(cublasSgeam(cublas_handle,transa2,transb2,m,n,&alpha,matA,lda,&beta,matB,ldb,matC,ldc));
}
template<>
cublasStatus_t carma_geam<double>(cublasHandle_t cublas_handle, char transa,
    char transb, int m, int n, double alpha, double *matA, int lda, double beta,
    double *matB, int ldb, double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublasOperation(transa);
  cublasOperation_t transb2 = carma_char2cublasOperation(transb);
  return carma_checkCublasStatus(cublasDgeam(cublas_handle,transa2,transb2,m,n,&alpha,matA,lda,&beta,matB,ldb,matC,ldc));
}

/** These templates are used to select the proper dgmm executable from T_data*/
template<class T_data>
cublasStatus_t
carma_dgmm(cublasHandle_t cublas_handle, cublasSideMode_t side, int m, int n,
    T_data *matA, int lda, T_data *vectx, int incx, T_data *matC, int ldc);
/**< Generic template for geam executable selection */
template<>
cublasStatus_t carma_dgmm<float>(cublasHandle_t cublas_handle,
    cublasSideMode_t side, int m, int n, float *matA, int lda, float *vectx,
    int incx, float *matC, int ldc) {
  return carma_checkCublasStatus(cublasSdgmm(cublas_handle,side,m,n,matA,lda,vectx,incx,matC,ldc));
}
template<>
cublasStatus_t carma_dgmm<double>(cublasHandle_t cublas_handle,
    cublasSideMode_t side, int m, int n, double *matA, int lda, double *vectx,
    int incx, double *matC, int ldc) {
  return carma_checkCublasStatus(cublasDdgmm(cublas_handle,side,m,n,matA,lda,vectx,incx,matC,ldc));
}

/*
 ____ _        _    ____ ____        _       __ _       _ _   _
 / ___| |      / \  / ___/ ___|    __| | ___ / _(_)_ __ (_) |_(_) ___  _ __  
 | |   | |     / _ \ \___ \___ \   / _` |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \
 | |___| |___ / ___ \ ___) |__) | | (_| |  __/  _| | | | | | |_| | (_) | | | |
 \____|_____/_/   \_\____/____/   \__,_|\___|_| |_|_| |_|_|\__|_|\___/|_| |_|

 */

template<class T_data>
int carma_obj<T_data>::host2deviceVect(T_data *data, int incx, int incy) {
  /** \brief host2device generic data transfer method for a blas object vector.
   * \param data : input data (x)
   * \param incx : increment in x
   * \param incy : increment in y
   *
   * this method fills a vector with the input data
   */

  carma_checkCublasStatus(
      cublasSetVector(this->nb_elem, sizeof(T_data), data, incx, this->d_data,
          incy));
  return EXIT_SUCCESS;
}
template int
caObjS::host2deviceVect(float *data, int incx, int incy);
template int
caObjD::host2deviceVect(double *data, int incx, int incy);

template<class T_data>
int carma_obj<T_data>::device2hostVect(T_data *data, int incx, int incy) {
  /** \brief device2host generic data transfer method for a blas object vector.
   * \param data : output data (y)
   * \param incx : increment in x
   * \param incy : increment in y
   *
   * this method fills output with the vector data
   */

  carma_checkCublasStatus(
      cublasGetVector(this->nb_elem, sizeof(T_data), this->d_data, incx, data,
          incy));
  return EXIT_SUCCESS;
}

template<class T_data>
int carma_obj<T_data>::host2deviceMat(T_data *data, int lda, int ldb) {
  /** \brief host2device generic data transfer method.
   * \param data : input data  (A)
   * \param mat : matrix to fill(B)
   * \param lda : leading dim of A (# of rows)
   * \param ldb : leading dim of B (# of rows)
   *
   * this method fills mat with the input data
   */
  carma_checkCublasStatus(
      cublasSetMatrix(this->dims_data[1], this->dims_data[2], sizeof(T_data),
          data, lda, this->d_data, ldb));
  return EXIT_SUCCESS;
}

template int
caObjS::host2deviceMat(float *data, int lda, int ldb);
template int
caObjD::host2deviceMat(double *data, int lda, int ldb);

template<class T_data>
int carma_obj<T_data>::device2hostMat(T_data *data, int lda, int ldb) {
  /** \brief device2host generic data transfer.
   * \param data : output data  (A)
   * \param mat : matrix to copy(B)
   * \param lda : leading dim of A (# of rows)
   * \param ldb : leading dim of B (# of rows)
   *
   * this method fills output with the mat data
   */
  carma_checkCublasStatus(
      cublasGetMatrix(this->dims_data[1], this->dims_data[2], sizeof(T_data),
          this->d_data, lda, data, ldb));
  return EXIT_SUCCESS;
}

/*
 ____  _        _    ____  _ 
 | __ )| |      / \  / ___|/ |
 |  _ \| |     / _ \ \___ \| |
 | |_) | |___ / ___ \ ___) | |
 |____/|_____/_/   \_\____/|_|

 */

template<class T_data>
int carma_obj<T_data>::imax(int incx) {
  /** \brief imax method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the maximum magnitude element in vect
   */
  return carma_wheremax(current_context->get_cublasHandle(), this->nb_elem,
      this->d_data, incx);;
}
template int
carma_obj<float>::imax(int incx);
template int
carma_obj<double>::imax(int incx);

template<class T_data>
int carma_obj<T_data>::imin(int incx) {
  /** \brief imin method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the minimum magnitude element in vect
   */
  return carma_wheremin(current_context->get_cublasHandle(), this->nb_elem,
      this->d_data, incx);;
}
template int
carma_obj<float>::imin(int incx);
template int
carma_obj<double>::imin(int incx);

template<class T_data>
T_data carma_obj<T_data>::asum(int incx) {
  /** \brief asum method
   * \param incx : increment in x
   *
   * this method computes the sum of the absolute values of the elements
   */
  return carma_getasum(current_context->get_cublasHandle(), this->nb_elem,
      this->d_data, incx);;
}
template float
carma_obj<float>::asum(int incx);
template double
carma_obj<double>::asum(int incx);

template<class T_data>
T_data carma_obj<T_data>::nrm2(int incx) {
  /** \brief getNrm2 method
   * \param n    : vect size
   * \param vect : vector to sum (x)
   * \param incx : increment in x
   *
   * this method computes the Euclidean norm of vect
   */
  return carma_nrm2(current_context->get_cublasHandle(), this->nb_elem,
      this->d_data, incx);;
}
template float
carma_obj<float>::nrm2(int incx);
template double
carma_obj<double>::nrm2(int incx);

template<class T_data>
void carma_obj<T_data>::scale(T_data alpha, int incx) {
  /** \brief vectScale method
   * \param n    : vect size
   * \param alpha : scale factor
   * \param vect : vector to scale (x)
   * \param incx : increment in x
   *
   * this method replaces vector x with alpha * x

   */
  carma_scal(current_context->get_cublasHandle(), this->nb_elem, alpha,
      this->d_data, incx);
}
template void
carma_obj<float>::scale(float alpha, int incx);
template void
carma_obj<double>::scale(double alpha, int incx);

template<class T_data>
void carma_obj<T_data>::copy(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief vectCopy method
   * \param n    : vect size
   * \param vectx : vector to copy (x)
   * \param incx : increment in x
   * \param vecty : vector to fill (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_copy(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
      incx, this->d_data, incy);
}
template void
carma_obj<float>::copy(caObjS *, int incx, int incy);
template void
carma_obj<double>::copy(caObjD *, int incx, int incy);

template<class T_data>
void carma_obj<T_data>::swap(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief vectSwap method
   * \param n    : vect size
   * \param vectx : first vector to swap (x)
   * \param incx : increment in x
   * \param vecty : second vector to swap (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_swap(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
      incx, this->d_data, incy);
}
template void
carma_obj<float>::swap(caObjS *, int incx, int incy);
template void
carma_obj<double>::swap(caObjD *, int incx, int incy);

template<class T_data>
void carma_obj<T_data>::axpy(T_data alpha, carma_obj<T_data> *source, int incx,
    int incy) {
  /** \brief computeAXPY method
   * \param alpha : scale factor
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incy : increment in y
   *
   * this method multiplies vector x by scalar alpha and adds the result to vector y
   */
  carma_axpy(current_context->get_cublasHandle(), this->nb_elem, alpha,
      source->d_data, incx, this->d_data, incy);
}
template void
carma_obj<float>::axpy(float alpha, caObjS *source, int incx, int incy);
template void
carma_obj<double>::axpy(double alpha, caObjD *source, int incx, int incy);

template<class T_data>
T_data carma_obj<T_data>::dot(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief ccomputeDot method
   * \param n    : vect size
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incy : increment in y
   *
   * this method computes the dot product of two vectors
   */
  return carma_dot(current_context->get_cublasHandle(), this->nb_elem,
      source->d_data, incx, this->d_data, incy);
}
template float
carma_obj<float>::dot(caObjS *source, int incx, int incy);
template double
carma_obj<double>::dot(caObjD *source, int incx, int incy);

template<class T_data>
void carma_obj<T_data>::rot(carma_obj<T_data> *source, int incx, int incy,
    T_data sc, T_data ss) {
  /** \brief ccomputeRot method
   * \param n    : vect size
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incx : increment in y
   * \param sc : cosinus of rotation angle
   * \param ss : sinus of rotation angle
   *
   * this method computes the dot product of two vectors
   */
  carma_rot(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
      incx, this->d_data, incy, sc, ss);
}
template void
carma_obj<float>::rot(caObjS *source, int incx, int incy, float, float);
template void
carma_obj<double>::rot(caObjD *source, int incx, int incy, double sc,
    double ss);

/*
 ____  _        _    ____ ____  
 | __ )| |      / \  / ___|___ \
 |  _ \| |     / _ \ \___ \ __) |
 | |_) | |___ / ___ \ ___) / __/
 |____/|_____/_/   \_\____/_____|

 */

template<class T_data>
void carma_obj<T_data>::gemv(char trans, T_data alpha, carma_obj<T_data> *matA,
    int lda, carma_obj<T_data> *vectx, int incx, T_data beta, int incy) {
  /** \brief gemv method.
   * \param trans : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment for x
   * \param incy : increment for y
   *
   * this method performs one of the matrix‐vector operations y = alpha * op(A) * x + beta * y
   */
  carma_gemv(current_context->get_cublasHandle(), trans, matA->dims_data[1],
      matA->dims_data[2], alpha, matA->d_data, lda, vectx->d_data, incx, beta,
      this->d_data, incy);
}
template void
carma_obj<float>::gemv(char, float, caObjS *, int, caObjS*, int, float, int);
template void
carma_obj<double>::gemv(char, double, caObjD *, int, caObjD *, int, double,
    int);

template<class T_data>
void carma_obj<T_data>::ger(T_data alpha, carma_obj<T_data> *vectx, int incx,
    carma_obj<T_data> *vecty, int incy, int lda) {
  /** \brief ger method.
   * \param alpha : alpha
   * \param vectx : m-element vector
   * \param incx : increment for x
   * \param vecty : vector y
   * \param incy : increment for y
   * \param lda : leading dim of A (# of rows)
   *
   * this method performs the symmetric rank 1 operation A = alpha * x * y T + A
   */
  carma_ger(current_context->get_cublasHandle(), this->dims_data[1],
      this->dims_data[2], alpha, vectx->d_data, incx, vecty->d_data, incy,
      this->d_data, lda);
}
template void
carma_obj<float>::ger(float, caObjS *, int, caObjS*, int, int);
template void
carma_obj<double>::ger(double, caObjD *, int, caObjD *, int, int);

template<class T_data>
void carma_obj<T_data>::symv(cublasFillMode_t uplo, T_data alpha,
    carma_obj<T_data> *matA, int lda, carma_obj<T_data> *vectx, int incx,
    T_data beta, int incy) {
  /** \brief symv method.
   * \param uplo : upper or lower fill mode
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment for x
   * \param incy : increment for y
   *
   * this method performs one of the matrix‐vector operations y = alpha * op(A) * x + beta * y
   */
  carma_symv(current_context->get_cublasHandle(), uplo, matA->dims_data[1],
      alpha, matA->d_data, lda, vectx->d_data, incx, beta, this->d_data, incy);
}
template void
carma_obj<float>::symv(cublasFillMode_t, float, caObjS *, int, caObjS*, int,
    float, int);
template void
carma_obj<double>::symv(cublasFillMode_t, double, caObjD *, int, caObjD *, int,
    double, int);

/*
 ____  _        _    ____ _____ 
 | __ )| |      / \  / ___|___ /
 |  _ \| |     / _ \ \___ \ |_ \
 | |_) | |___ / ___ \ ___) |__) |
 |____/|_____/_/   \_\____/____/

 */

template<class T_data>
void carma_obj<T_data>::gemm(char transa, char transb, T_data alpha,
    carma_obj<T_data> *matA, int lda, carma_obj<T_data> *matB, int ldb,
    T_data beta, int ldc) {
  /** \brief generic gemm method.
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param transb : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the matrix‐matrix operations:
   * C = alpha * op(A) * op(B) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  int k = (
      ((transa == 'N') || (transa == 'n')) ?
          matA->dims_data[2] : matA->dims_data[1]);
  carma_gemm(current_context->get_cublasHandle(), transa, transb,
      this->dims_data[1], this->dims_data[2], k, alpha, matA->d_data, lda,
      matB->d_data, ldb, beta, this->d_data, ldc);
}
template void
carma_obj<float>::gemm(char, char, float, caObjS *, int, caObjS*, int, float,
    int);
template void
carma_obj<double>::gemm(char, char, double, caObjD *, int, caObjD *, int,
    double, int);

template<class T_data>
void carma_obj<T_data>::symm(cublasSideMode_t side, cublasFillMode_t uplo,
    T_data alpha, carma_obj<T_data> *matA, int lda, carma_obj<T_data> *matB,
    int ldb, T_data beta, int ldc) {
  /** \brief generic symm method.
   * \param side : which side of the equation is symmetric matrix A
   * \param uplo : fill mode of matrix A
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * A * B + beta * C ,
   * or
   * C = alpha * B * A + beta * C ,
   * where A is a symmetric matrix
   */
  carma_symm(current_context->get_cublasHandle(), side, uplo,
      this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
      matB->d_data, ldb, beta, this->d_data, ldc);
}
template void
carma_obj<float>::symm(cublasSideMode_t, cublasFillMode_t, float, caObjS *, int,
    caObjS*, int, float, int);
template void
carma_obj<double>::symm(cublasSideMode_t, cublasFillMode_t, double, caObjD *,
    int, caObjD *, int, double, int);

template<class T_data>
void carma_obj<T_data>::syrk(cublasFillMode_t uplo, char transa, T_data alpha,
    carma_obj<T_data> *matA, int lda, T_data beta, int ldc) {
  /** \brief generic syrk method.
   * \param uplo : fill mode of matrix A
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * op(A) * transpose(op(A)) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  carma_syrk(current_context->get_cublasHandle(), uplo, transa,
      this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda, beta,
      this->d_data, ldc);
}
template void
carma_obj<float>::syrk(cublasFillMode_t, char, float, caObjS *, int, float,
    int);
template void
carma_obj<double>::syrk(cublasFillMode_t, char, double, caObjD *, int, double,
    int);

template<class T_data>
void carma_obj<T_data>::syrkx(cublasFillMode_t uplo, char transa, T_data alpha,
    carma_obj<T_data> *matA, int lda, carma_obj<T_data> *matB, int ldb,
    T_data beta, int ldc) {
  /** \brief generic syrkx method.
   * \param uplo : fill mode of matrix A
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * op(A) * transpose(op(B)) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  carma_syrkx(current_context->get_cublasHandle(), uplo, transa,
      this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
      matB->d_data, ldb, beta, this->d_data, ldc);
}
template void
carma_obj<float>::syrkx(cublasFillMode_t, char, float, caObjS *, int, caObjS*,
    int, float, int);
template void
carma_obj<double>::syrkx(cublasFillMode_t, char, double, caObjD *, int,
    caObjD *, int, double, int);

template<class T_data>
void carma_obj<T_data>::geam(char transa, char transb, T_data alpha,
    carma_obj<T_data> *matA, int lda, T_data beta, carma_obj<T_data> *matB,
    int ldb, int ldc) {
  /** \brief generic geam method.
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param transb : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param beta : beta
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param ldc : leading dim of C (# of rows)
   *
   * C = alpha * op(A) + beta * op(B),
   * where  op(X) = X  or  op(X) = X^T
   */

  carma_geam(current_context->get_cublasHandle(), transa, transb,
      this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda, beta,
      matB->d_data, ldb, this->d_data, ldc);
}
template void
carma_obj<float>::geam(char, char, float, caObjS *, int, float, caObjS*, int,
    int);
template void
carma_obj<double>::geam(char, char, double, caObjD *, int, double, caObjD *,
    int, int);

template<class T_data>
void carma_obj<T_data>::dgmm(cublasSideMode_t side, carma_obj<T_data> *matA,
    int lda, carma_obj<T_data> *vectx, int incx, int ldc) {
  /** \brief generic dgmm method.
   * \param side : side of equation for matrix A
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment on x
   * \param ldc : leading dim of C (# of rows)
   *
   * C = A * diag(X) or C = diag(X) * A
   */

  carma_dgmm(current_context->get_cublasHandle(), side, this->dims_data[1],
      this->dims_data[2], matA->d_data, lda, vectx->d_data, incx, this->d_data,
      ldc);
}
template void
carma_obj<float>::dgmm(cublasSideMode_t, caObjS *, int, caObjS*, int, int);
template void
carma_obj<double>::dgmm(cublasSideMode_t, caObjD *, int, caObjD *, int, int);

