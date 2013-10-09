/*
 * yoga_cublas.h
 *
 *  Created on: Apr 23, 2012
 *      Author: sevin
 */

#ifndef YOGA_CUBLAS_H_
#define YOGA_CUBLAS_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cublas_v2.h>

using namespace std;

cublasStatus_t yoga_checkCublasStatus(cublasStatus_t status);

void yoga_initCublas();
void yoga_shutdownCublas();

cublasOperation_t yoga_char2cublasOperation(char operation);

/*
 _____ _____ __  __ ____  _        _  _____ _____ ____
|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/

 */

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T_data> int yoga_wheremax(int n, T_data *vect, int incx);
/**< Generic template for Iamax executable selection */

/** These templates are used to select the proper Iamin executable from T_data*/
template<class T_data> int yoga_wheremin(int n, T_data *vect, int incx);
/**< Generic template for Iamin executable selection */

/** These templates are used to select the proper asum executable from T_data*/
template<class T_data> T_data yoga_getasum(int n, T_data *vect, int incx);
/**< Generic template for asum executable selection */

/** These templates are used to select the proper axpy executable from T_data*/
template<class T_data> void yoga_axpy(int n, T_data alpha, T_data *vectx, int incx, T_data *vecty, int incy);
/**< Generic template for axpy executable selection */

/** These templates are used to select the proper dot executable from T_data*/
template<class T_data> T_data yoga_dot(int n, T_data *vectx, int incx, T_data *vecty, int incy);
/**< Generic template for dot executable selection */

/** These templates are used to select the proper nrm2 executable from T_data*/
template<class T_data> T_data yoga_nrm2(int n, T_data *vect, int incx);
/**< Generic template for nrm2 executable selection */

/** These templates are used to select the proper rot executable from T_data*/
template<class T_data> void yoga_rot(int n, T_data *vectx, int incx, T_data *vecty, int incy, T_data sc, T_data ss);
/**< Generic template for rot executable selection */

/** These templates are used to select the proper scal executable from T_data*/
template<class T_data> void yoga_scal(int n, T_data alpha, T_data *vectx, int incx);
/**< Generic template for scal executable selection */

/** These templates are used to select the proper swap executable from T_data*/
template<class T_data> void yoga_swap(int n, T_data *vectx, int incx, T_data *vecty, int incy);
/**< Generic template for swap executable selection */

/** These templates are used to select the proper gemv executable from T_data*/
template<class T_data> void yoga_gemv(char trans, int m, int n, T_data alpha, T_data *matA, int lda, T_data *vectx, int incx, T_data beta, T_data *vecty, int incy);
/**< Generic template for gemv executable selection */

/** These templates are used to select the proper ger executable from T_data*/
template<class T_data> void yoga_ger(int m, int n, T_data alpha, T_data *vectx, int incx, T_data *vecty, int incy, T_data *matA, int lda);
/**< Generic template for ger executable selection */

/** These templates are used to select the proper gemm executable from T_data*/
template<class T_data> void yoga_gemm(char transa, char transb, int m, int n, int k, T_data alpha, T_data *matA, int lda, T_data *matB, int ldb, T_data beta, T_data *matC, int ldc);
/**< Generic template for gemm executable selection */


/*
__        ______      _    ____  ____  _____ ____  ____
\ \      / /  _ \    / \  |  _ \|  _ \| ____|  _ \/ ___|
 \ \ /\ / /| |_) |  / _ \ | |_) | |_) |  _| | |_) \___ \
  \ V  V / |  _ <  / ___ \|  __/|  __/| |___|  _ < ___) |
   \_/\_/  |_| \_\/_/   \_\_|   |_|   |_____|_| \_\____/

 */

/** The following are wrappers to be called from Yorick */

/*
 BLAS 1
*/

int _yogaInitCublas();

int _yogaShutdownCublas();



#endif /* YOGA_CUBLAS_H_ */
