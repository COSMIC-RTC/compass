#include "yoga_template.h"
#include <yoga_api.h>

#include "magma.h"
#include "magma_lapack.h"
#include <cublas_v2.h>

#define magma_dsyevd_gpu_t	magma_dsyevd_gpu	    // slow dsymv from cublas legacy api
//#define magma_dsyevd_gpu_t	magma_dsyevd_gpu_magmablas  // fast dsymv from magma
//#define magma_dsyevd_gpu_t	magma_dsyevd_gpu_kblas	    // fast dsymv from kblas

#define magma_ssyevd_gpu_t	magma_ssyevd_gpu	 

extern "C" {

  void
  Y_yoga_magma_syevd(int argc)
  {
    magma_init();  
                                                         
    cublasHandle_t handle;
    if( CUBLAS_STATUS_SUCCESS != cublasCreate(&handle) ) {                          
        fprintf(stderr, "ERROR: cublasInit failed\n");                     
        magma_finalize();                                                  
        exit(-1);                                                          
    }    
                                                                      
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    //magma_print_devices();
    magma_int_t *iwork;
    magma_int_t N, n2, info, lwork, liwork, lda, ldda, aux_iwork[1];
    void *h_A;

    if (!yarg_subroutine()) {
      // retreive matrix and parameters from yorick session
      char jobz = ygets_c(argc-1);
      char uplo = ygets_c(argc-2);
      void *A_handler = ygeta_any(argc - 3, &ntot, dims, &yType);
      yObj_struct *d_R_handler = yoga_getyObj(argc, 4);

      if (yType == Y_OPAQUE) {
	if (A_handler->type != d_R_handler->type) 
	  y_error("Please provide an array and a yoga_obj with same precision\n");
      } else {
	h_A = (void *)A_handler;
	if (yType != d_R_handler->type) 
	  y_error("Please provide an array and a yoga_obj with same precision\n");
      }

      void *h_R, *h_work, *w1;
      void *d_A;

      if (d_R_handler->type == Y_DOUBLE) {
	if (yType == Y_OPAQUE) {
	  (caObjD *)d_A = (caObjD *)(A_handler->carma_object);
	  N    = (caObjD *)d_A->getDims()[1];
	} else N = dims[1];

	n2   = N*N;
	lda  = N;
	ldda = N;//((N + 31)/32)*32;
	long dims_eigen[2];
	dims_eigen[0] = 1; dims_eigen[1] = N; 
	// create magma workspace
	double aux_work[1];
	
	caObjD *d_R = (caObjD *)(d_R_handler->carma_object);
	if (N * ldda != d_R->getNbElem()) 
	  y_error("Please provide arrays of proper size\n");
	
	/* Query for workspace sizes */      
	magma_dsyevd_gpu_t( jobz, uplo,
			    N, (double *)d_R->getData(), 
			    ldda, (double *)w1,
			    (double *)h_R, lda,
			    aux_work,  -1,
			    aux_iwork, -1,
			    &info );
	
	lwork  = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];
	posix_memalign( (void**) &iwork, 32, (liwork)*sizeof(int) );
	
	cudaMallocHost( (void**) &h_R, (N*lda)*sizeof(double) );
	cudaMallocHost( (void**) &h_work, (lwork)*sizeof(double) );
	/* pushing eigen values on the yorick stack */
	w1 = (double *)ypush_d(dims_eigen);
	
	/* Initialize the matrix */
	//fprintf(stderr, "%s: %s@%d\n", __FILE__, __FUNCTION__, __LINE__);
	if (yType == Y_OPAQUE) {
	  d_R->copy(d_A, 1, 1);
	} else {
	  magma_dsetmatrix( N, N, (double *)h_A, lda, (double *)d_R->getData(), ldda );
	}
	magma_dsyevd_gpu_t( jobz, uplo,
			    N, (double *)d_R->getData(), 
			    ldda, (double *)w1,
			    (double *)h_R, lda,
			    (double *)h_work, lwork,
			    iwork, liwork,
			    &info );
      }
      if (d_R_handler->type == Y_FLOAT) {
	if (yType == Y_OPAQUE) {
	  (caObjS *)d_A = (caObjS *)(A_handler->carma_object);
	  N    = (caObjS *)d_A->getDims()[1];
	} else N = dims[1];

	n2   = N*N;
	lda  = N;
	ldda = N;//((N + 31)/32)*32;
	long dims_eigen[2];
	dims_eigen[0] = 1; dims_eigen[1] = N; 
	// create magma workspace
	float aux_work[1];
	
	caObjS *d_R = (caObjS *)(d_R_handler->carma_object);
	if (N * ldda != d_R->getNbElem()) 
	  y_error("Please provide arrays of proper size\n");
	
	/* Query for workspace sizes */      
	magma_ssyevd_gpu_t( jobz, uplo,
			    N, (float *)d_R->getData(), 
			    ldda, (float *)w1,
			    (float *)h_R, lda,
			    aux_work,  -1,
			    aux_iwork, -1,
			    &info );
	
	lwork  = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];
	posix_memalign( (void**) &iwork, 32, (liwork)*sizeof(int) );
	
	cudaMallocHost( (void**) &h_R, (N*lda)*sizeof(float) );
	cudaMallocHost( (void**) &h_work, (lwork)*sizeof(float) );
	/* pushing eigen values on the yorick stack */
	w1 = (float *)ypush_f(dims_eigen);
	
	/* Initialize the matrix */
	//fprintf(stderr, "%s: %s@%d\n", __FILE__, __FUNCTION__, __LINE__);
	if (yType == Y_OPAQUE) {
	  d_R->copy(d_A, 1, 1);
	} else {
	  magma_ssetmatrix( N, N, (float *)h_A, lda, (float *)d_R->getData(), ldda );
	}
	
	magma_ssyevd_gpu_t( jobz, uplo,
			    N, (float *)d_R->getData(), 
			    ldda, (float *)w1,
			    (float *)h_R, lda,
			    (float *)h_work, lwork,
			    iwork, liwork,
			    &info );
      }

      cudaFreeHost(h_R);
      cudaFreeHost(h_work);
      free(iwork);
      magma_finalize();                                                      
      cublasDestroy(handle);

    } else y_error("can only be called as a subroutine");
  }

  void
  Y_yoga_magma_getri(int argc)
  {
    magma_init();  
                                                         
    cublasHandle_t handle;
    if( CUBLAS_STATUS_SUCCESS != cublasCreate(&handle) ) {                          
        fprintf(stderr, "ERROR: cublasInit failed\n");                     
        magma_finalize();                                                  
        exit(-1);                                                          
    }    
                                                                      
    //magma_print_devices();
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    if (yarg_subroutine()) {
      // retreive matrix and parameters from yorick session
      void *dwork;
      magma_int_t N, n2, lda, ldda, info, lwork, ldwork;
      magma_int_t *ipiv;

      void *h_A = ygeta_any(argc - 1, &ntot, dims, &yType);
      N = dims[1];
      n2   = N*N;
      if (n2 != ntot) 
	y_error("Please provide a square matrix\n");

      lda  = N;
      ldda = ((N + 31)/32)*32;
      ldwork = N * magma_get_dgetri_nb( N );

      posix_memalign( (void**) &ipiv, 32, (N)*sizeof(int) );

      yObj_struct *d_A_handler = yoga_getyObj(argc, 2);

      if (yType != d_A_handler->type) 
	y_error("Please provide an array and a yoga_obj with same precision\n");

      if (yType == Y_DOUBLE) {
	// query for workspace size
	lwork = -1;
	double tmp;
	lapackf77_dgetri( &N, NULL, &lda, NULL, &tmp, &lwork, &info );
	lwork = int(tmp);

	// create magma workspace
	caObjD *d_A = (caObjD *)(d_A_handler->carma_object);
	if (N * ldda != d_A->getNbElem()) 
	  y_error("Please provide arrays of proper size\n");
	
	cudaMalloc( (void**) &dwork, (ldwork)*sizeof(double) );

	/* Factor the matrix */
	magma_dsetmatrix( N, N, (double *)h_A, lda, d_A->getData(), ldda );
	magma_dgetrf_gpu( N, N, d_A->getData(), ldda, ipiv, &info );
	magma_dgetmatrix( N, N, d_A->getData(), ldda, (double *)h_A, lda );

	magma_dgetri_gpu( N,    d_A->getData(), ldda, ipiv, (double *)dwork, ldwork, &info );

      }
      if (yType == Y_FLOAT) {
	// query for workspace size
	lwork = -1;
	float tmp;
	lapackf77_sgetri( &N, NULL, &lda, NULL, &tmp, &lwork, &info );
	lwork = int(tmp);

	// create magma workspace
	caObjS *d_A = (caObjS *)(d_A_handler->carma_object);
	if (N * ldda != d_A->getNbElem()) 
	  y_error("Please provide arrays of proper size\n");
	
	cudaMalloc( (void**) &dwork, (ldwork)*sizeof(float) );

	/* Factor the matrix */
	magma_ssetmatrix( N, N, (float *)h_A, lda, d_A->getData(), ldda );
	magma_sgetrf_gpu( N, N, d_A->getData(), ldda, ipiv, &info );
	magma_sgetmatrix( N, N, d_A->getData(), ldda, (float *)h_A, lda );

	magma_sgetri_gpu( N,    d_A->getData(), ldda, ipiv, (float *)dwork, ldwork, &info );
      }

      cudaFree(dwork);
      free(ipiv);
      magma_finalize();                                                      
      cublasDestroy(handle);

    } else y_error("cannot be called as a subroutine");
  }
}
