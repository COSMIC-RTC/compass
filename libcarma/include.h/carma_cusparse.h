// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_cusparse.h
//! \ingroup   libcarma
//! \brief     this file provides the cusparse features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#ifndef CARMA_SPARSE_H_
#define CARMA_SPARSE_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cusparse_v2.h>
#include <string>

#define carma_check_cusparse_status(status) \
  carma_check_cusparse_status_v2(status, __LINE__, __FILE__)

cusparseStatus_t carma_check_cusparse_status_v2(cusparseStatus_t status, int line,
                                              std::string file);
cusparseStatus_t carma_init_cusparse(cusparseHandle_t *cusparse_handle);
cusparseStatus_t carma_shutdown_cusparse(cusparseHandle_t cusparse_handle);

cusparseOperation_t carma_char2cusparse_operation(char operation);

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

#endif /* CARMA_SPARSE_H_ */
