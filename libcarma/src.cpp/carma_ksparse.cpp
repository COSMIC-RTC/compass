// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_ksparse.cpp
//! \ingroup   libcarma
//! \brief     this file provides the ksparse features to carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
int carma_kgemv<float>(carma_sparse_obj<float>* A, float alpha,
                       const float* __restrict x, float beta, float* y) {
  if (A->format != "BSR") {
    DEBUG_TRACE("carma_kgemv needs a BSR matrix as input");
  }
  TEST_USE_KSPARSE(return carma_kgemv_gen<float, ksparse_sbsrmv>(
      A->dims_data[1], A->blockDim, alpha, A->d_data, A->d_rowind, A->d_colind,
      x, beta, y));
}
template <>
int carma_kgemv<double>(carma_sparse_obj<double>* A, double alpha,
                        const double* __restrict x, double beta, double* y) {
  if (A->format != "BSR") {
    DEBUG_TRACE("carma_kgemv needs a BSR matrix as input");
  }
  TEST_USE_KSPARSE(return carma_kgemv_gen<double, ksparse_dbsrmv>(
      A->dims_data[1], A->blockDim, alpha, A->d_data, A->d_rowind, A->d_colind,
      x, beta, y));
}
