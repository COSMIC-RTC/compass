
#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <unistd.h>
#include "kp_KFPP.h"


template<typename T> void kernel_add(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_sub(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_mult(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_div(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_add_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_sub_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_mult_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_div_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_sqrt(T* d_cu1, int length);
template<typename T> void kernel_inv(T* d_cu1, int length);
template<typename T, typename U> void kernel_memcpy(T* d_cu_dst, const U* d_cu_src, int length);
template<typename T> void kernel_memset(T* d_cu_dst, T valeur, int length);
template<typename T> void kernel_A1_2vec(T* Atur_d_cu, T* Btur_d_cu, int* colind, int* rowind, T* values, int nb_az, int nnz);
template<typename T> void kernel_memcpy_diag(T *dev_dst, const T *dev_src, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);
template<typename T> void kernel_memset_diag(T *dev_dst, T val, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);
template<typename T> void kernel_sparse2full(T *dev_dst, int *dev_src_rowind, int *dev_src_colind, T* dev_src_values, int nnz, int src_dim1, int src_dim2);
template<typename T> void kernel_add_diag_const(T* d_cu1, T val, int dim1);
template<typename T> void kernel_get_diag(T* dst_diag, T* src_M, int dim1);
template<typename T> void kernel_diag_mult(T* mat_dst, T* mat_src, T* vec_src, int dim1, int nb_elements);
template<typename T> void kernel_diag_mult2(T* mat_dst, T* mat_src, T* vec_src, int dim1, int nb_elements);
template<typename T> void kernel_diag_mult3(T* mat_dst1, T* mat_dst2, T* mat_src, T* vec_src, int dim1, int nb_elements);
template<typename T> void kernel_set_submatrix(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_col);
template<typename T> void kernel_set_Linf(T* mat_Linf, T* mat_Hinf, T* atur, T* btur, int nb_n, int nb_p, int ordreAR);
//template<typename T> void kernel_set_Linf_AR2(T* mat_Linf, T* mat_Hinf, T* atur, T* btur, int nb_n, int nb_p);
template<typename T> void kernel_add_A1(T* mat, T* atur, T* btur, int nb_az, int ordreAR);



#endif
