// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_ksparse.cpp
//! \ingroup   libcarma
//! \brief     this file provides the ksparse features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <carma_obj.h>
#include <carma_sparse_obj.h>

#ifdef USE_KSPARSE
#include <ksparse.h>
#define TEST_USE_KSPARSE(...) __VA_ARGS__
#else
#define TEST_USE_KSPARSE(...)                          \
  DEBUG_TRACE("!!!!!! KSPARSE not compiled !!!!!!\n"); \
  return EXIT_FAILURE;
#endif

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

template <class T_data, int (*ksparse_bsrmv)(
                            int matrix_dim, int block_dim, T_data alpha,
                            T_data* A, int* bsr_row_ptr, int* bsr_col_ind,
                            const T_data* __restrict x, T_data beta, T_data* y)>
int carma_kgemv_gen(int matrix_dim, int block_dim, T_data alpha, T_data* A,
                    int* bsr_row_ptr, int* bsr_col_ind,
                    const T_data* __restrict x, T_data beta, T_data* y) {
  TEST_USE_KSPARSE(return ksparse_bsrmv(matrix_dim, block_dim, alpha, A,
                                        bsr_row_ptr, bsr_col_ind, x, beta, y));
}

template <>
int carma_kgemv<float>(CarmaSparseObj<float>* A, float alpha,
                       const float* __restrict x, float beta, float* y) {
  if (A->format != "BSR") {
    DEBUG_TRACE("carma_kgemv needs a BSR matrix as input");
  }
  TEST_USE_KSPARSE(return carma_kgemv_gen<float, ksparse_sbsrmv>(
      A->dims_data[1], A->block_dim, alpha, A->d_data, A->d_rowind, A->d_colind,
      x, beta, y));
}
template <>
int carma_kgemv<double>(CarmaSparseObj<double>* A, double alpha,
                        const double* __restrict x, double beta, double* y) {
  if (A->format != "BSR") {
    DEBUG_TRACE("carma_kgemv needs a BSR matrix as input");
  }
  TEST_USE_KSPARSE(return carma_kgemv_gen<double, ksparse_dbsrmv>(
      A->dims_data[1], A->block_dim, alpha, A->d_data, A->d_rowind, A->d_colind,
      x, beta, y));
}
