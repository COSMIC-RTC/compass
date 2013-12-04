#include "yoga_template.h"
#include <yoga_api.h>

#include "magma.h"
#include "magma_lapack.h"
//#include "testings.h"

#define magma_dsyevd_gpu_t	magma_dsyevd_gpu		// slow dsymv from cublas legacy api
//#define magma_dsyevd_gpu_t	magma_dsyevd_gpu_magmablas	// fast dsymv from magma
//#define magma_dsyevd_gpu_t	magma_dsyevd_gpu_kblas		// fast dsymv from kblas

extern "C" {

  void
  Y_yoga_immult(int argc)
  {
    if (yarg_subroutine()) {
      yObj_struct *handle_out = (yObj_struct *) yget_obj(argc - 1, &yObj);
      yObj_struct *handle_in = (yObj_struct *) yget_obj(argc - 2, &yObj);
      
      carma_context *context_handle = _getCurrentContext();
      context_handle->set_activeDevice(handle_in->device);
      if (handle_in->type == Y_FLOAT) {
	caObjS *carma_out_handler = (caObjS *)(handle_out->carma_object);
	caObjS *carma_in_handler = (caObjS *)(handle_in->carma_object);
	multim((float *)carma_out_handler->getData(),(float *)carma_in_handler->getData(),carma_in_handler->getNbElem());
      } else y_error("double not implemented yet");
    } else y_error("can only be called as a subroutine");
  }

  void
  Y_yoga_magma_evd(int argc)
  {
    int yType;
    long ntot;
    long dims[Y_DIMSIZE];
    if (yarg_subroutine()) {
      // retreive matrix and parameters from yorick session
      char jobz = ygets_c(argc-1);
      char uplo = ygets_c(argc-2);
      void *h_A = ygeta_any(argc - 3, &ntot, dims,&yType);
      if (yType != Y_DOUBLE) y_error("Please provide a double precision array\n");

      // d_R becomes a carma object retreived from Yorick 
      // !!!!!!!!! FIX ME !!!!!!!!!!!!!!!!!!!!!
      // need to test size and type of provided carma object
      //magma_malloc( (void**) &d_R, (N*lda)*sizeof(double) );
      yObj_struct *d_R_handler = (yObj_struct *) yget_obj(argc - 4, &yObj);
      caObjD *d_R = (caObjD *)(d_R_handler->carma_object);

      // create magma workspace
      double *h_R, *h_work;
      magma_int_t *iwork;

      magma_int_t N, n2, info, lwork, liwork, lda, ldda, aux_iwork[1];
      magma_int_t ione     = 1;
      magma_int_t ISEED[4] = {0,0,0,1};
      double      aux_work[1];

      N = dims[1];
      n2   = N*N;
      lda  = N;
      ldda = ((N + 31)/32)*32;

      // eigenvalues are pushed on the interpreter stack
      //magma_malloc_cpu( (void**) &w1, (N)*sizeof(double) );
      long dims_eigen[2];
      dims_eigen[0] = 1; dims_eigen[1] = N; 
      double *w1 = ypush_d(dims_eigen);

       /* Query for workspace sizes */
      magma_dsyevd_gpu_t( jobz, uplo,
			  N, d_R->getData(), ldda, w1,
			  h_R, lda,
			  aux_work,  -1,
			  aux_iwork, -1,
			  &info );
      lwork  = (magma_int_t) aux_work[0];
      liwork = aux_iwork[0];

      cudaMallocHost( (void**) &h_R, (N*lda)*sizeof(double) );
      cudaMallocHost( (void**) &h_work, (lwork)*sizeof(double) );
      posix_memalign( (void**) &iwork, 32, (liwork)*sizeof(int) );
      //magma_malloc_cpu( (void**) &iwork, (liwork)*sizeof(magma_int_t) );

      /* Initialize the matrix */
      lapackf77_dlarnv( &ione, ISEED, &n2, (double *)h_A );
      magma_dsetmatrix( N, N, (double *)h_A, lda, d_R->getData(), ldda );

      magma_dsyevd_gpu_t( jobz, uplo,
			  N, d_R->getData(), ldda, w1,
			  h_R, lda,
			  h_work, lwork,
			  iwork, liwork,
			  &info );

      cudaFreeHost(h_R);
      cudaFreeHost(h_work);
      free(iwork);

    } else y_error("can only be called as a subroutine");
  }

}
