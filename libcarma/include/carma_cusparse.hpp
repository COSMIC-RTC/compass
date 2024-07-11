// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_cusparse.hpp
//! \ingroup   libcarma
//! \brief     this file provides the cusparse features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef CARMA_SPARSE_H_
#define CARMA_SPARSE_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cusparse_v2.h>
#include <string>

#define carma_check_cusparse_status(status) \
  carma_check_cusparse_status_v2(status, __LINE__, __FILE__)

cusparseStatus_t carma_check_cusparse_status_v2(cusparseStatus_t status, int32_t line,
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
