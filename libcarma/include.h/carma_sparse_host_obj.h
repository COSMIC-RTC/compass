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

//! \file      carma_sparse_host_obj.h
//! \ingroup   libcarma
//! \class     carma_sparse_host_obj
//! \brief     this class provides wrappers to the generic carma sparse host object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef CARMASPARSEHOSTOBJ_H_
#define CARMASPARSEHOSTOBJ_H_
#include "carma_host_obj.h"

template <class T_data>
class carma_sparse_obj;

template <class T_data>
class carma_sparse_host_obj {
 public:
  carma_sparse_host_obj();
  carma_sparse_host_obj(carma_sparse_obj<T_data>& sm);
  carma_sparse_host_obj(carma_sparse_host_obj<T_data>& sm);
  carma_sparse_host_obj(const long* dims, T_data* M, char order);
  virtual ~carma_sparse_host_obj();

  // delete all arrays and create new for nnz=new_nnz
  void resize(int new_nnz, int dim1_, int dim2_);

  void operator=(carma_sparse_obj<T_data>& M);
  void operator=(carma_sparse_host_obj<T_data>& M);

  void init_from_matrix(const long* dims, T_data* M, char majorDim);
  void copy_into_matrix(T_data* M, char majorDim);
  void resize2rowMajor();
  void resize2colMajor();
  char get_majorDim() const { return majorDim; }

  /**< General Utilities */
  operator T_data*() { return h_data; }
  T_data* operator[](int index) { return &h_data[index]; }
  T_data* getData() { return h_data; }
  T_data* getData(int index) { return &h_data[index]; }
  const long* getDims() { return dims_data; }
  long getDims(int i) { return dims_data[i]; }
  int getNzElem() { return nz_elem; }

 private:
  void _create(int nnz_, int dim1_, int dim2_);  // create new arrays
  void _clear();                                 // clear arrays

  char majorDim;  // U - undefined
  // R - row major
  // C - col major

  MemAlloc mallocType;  ///< type of host alloc

 public:
  long dims_data[3];  ///< dimensions of the array
  int nz_elem;        ///< number of elements in the array

  // see Sparse Matrix Storage Formats (coordinate format)
  // ONE-BASED!!!!
  T_data* h_data;
  int* rowind;
  int* colind;
};

// Multiply sparce matrix by vector
// y := alpha*A*x + betta * y
template <class T_data>
void carma_gemv(T_data alpha, carma_sparse_host_obj<T_data>* A,
                carma_host_obj<T_data>* x, T_data betta,
                carma_host_obj<T_data>* y,
                void (*ptr_coomv)(char* transa, long* m, long* k, T_data* alpha,
                                  char* matdescra, T_data* val, int* rowind,
                                  int* colind, int* nnz, T_data* x,
                                  T_data* beta, T_data* y));

// Multiply sparce matrix by dense matrix
// C := alpha*op(A)*B + betta * y
// op_A could be N (do nothing) or T (transpose)
template <class T_data>
void carma_gemm(char op_A, T_data alpha, carma_sparse_host_obj<T_data>* A,
                carma_host_obj<T_data>* B, T_data betta,
                carma_host_obj<T_data>* C);

// Multiply sparce matrix by dense matrix
// C := alpha*A*B + betta * y
template <class T_data>
inline void carma_gemm(T_data alpha, carma_sparse_host_obj<T_data>* A,
                       carma_host_obj<T_data>* B, T_data betta,
                       carma_host_obj<T_data>* C) {
  kp_gemm('N', alpha, A, B, betta, C);
}

#endif /* CARMASPARSEHOSTOBJ_H_ */
