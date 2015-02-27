/*
 * carma_sparse.h
 *
 *  Created on: Apr 23, 2012
 *      Author: sevin
 */

#ifndef CARMA_SPARSE_H_
#define CARMA_SPARSE_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cusparse_v2.h>

using namespace std;

cusparseStatus_t carma_checkCusparseStatus(cusparseStatus_t status);

cusparseStatus_t carma_initCusparse(cusparseHandle_t *cusparse_handle);
cusparseStatus_t carma_shutdownCusparse(cusparseHandle_t cusparse_handle);

cusparseOperation_t carma_char2cusparseOperation(char operation);

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

#endif /* CARMA_SPARSE_H_ */
