
#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <unistd.h>
#include <iostream>
#include "kp_KFPP.h"

// d_cu1[i] += d_cu2[i] (pour i de 0 a length-1)
template <typename T>
void kernel_add(T* d_cu1, const T* d_cu2, int length);
// d_cu1[i] -= d_cu2[i] (pour i de 0 a length-1)
template <typename T>
void kernel_sub(T* d_cu1, const T* d_cu2, int length);
// d_cu1[i] *= d_cu2[i] (pour i de 0 a length-1)
template <typename T>
void kernel_mult(T* d_cu1, const T* d_cu2, int length);
// d_cu1[i] /= d_cu2[i] (pour i de 0 a length-1)
template <typename T>
void kernel_div(T* d_cu1, const T* d_cu2, int length);
// d_cu1[i] += valeur (pour i de 0 a length-1)
template <typename T>
void kernel_add_const(T* d_cu1, T valeur, int length);
// d_cu1[i] -= valeur (pour i de 0 a length-1)
template <typename T>
void kernel_sub_const(T* d_cu1, T valeur, int length);
// d_cu1[i] *= valeur (pour i de 0 a length-1)
template <typename T>
void kernel_mult_const(T* d_cu1, T valeur, int length);
// d_cu1[i] /= valeur (pour i de 0 a length-1)
template <typename T>
void kernel_div_const(T* d_cu1, T valeur, int length);

// d_cu1[i] = sqrt(d_cu1[i]) (pour i de 0 a length-1)
template <typename T>
void kernel_sqrt(T* d_cu1, int length);
// d_cu1[i] = 1/d_cu1[i]     (pour i de 0 a length-1)
template <typename T>
void kernel_inv(T* d_cu1, int length);

// d_cu_dst[i] = d_cu_src[i] (pour i de 0 a length-1)
template <typename T, typename U>
void kernel_memcpy(T* d_cu_dst, const U* d_cu_src, int length);
// d_cu_dst[i] = valeur (pour i de 0 a length-1)
template <typename T>
void kernel_memset(T* d_cu_dst, T valeur, int length);

// Atur_d_cu = diag(A(0:nb_az-1,0:nb_az-1)) et Btur_d_cu =
// diag(A(0:nb_az-1,nb_az:2*nb_az-1)) (avec A matrice creuse)
template <typename T>
void kernel_A1_2vec(T* Atur_d_cu, T* Btur_d_cu, int* colind, int* rowind,
                    T* values, int nb_az, int nnz);

// kernel_memcpy_diag(A, b, i, j, N, D). diag(A(i:i+N,j:j+N)) = b
// Equivalent Matlab : for k=1:N,A((i+1)+(k-1),(j+1)+(k-1))=b(k);end
// i et j ZERO-BASED INDEXING
// A matrice column-major, doit avoir au moins i+N lignes et au moins j+N
// colonnes b doit etre un vecteur de taille N N : nombre de lignes (et de
// colonnes) de la sous matrice de A D : nombre de lignes de de A
template <typename T>
void kernel_memcpy_diag(T* dev_dst, const T* dev_src, int rowDebut,
                        int colDebut, int numElements, int dev_dest_dim1);

// kernel_memset_diag(A, k, i, j, N, D). diag(A(i:i+N,j:j+N)) = k (avec k reel)
// Equivalent Matlab : for k=1:N,A((i+1)+(k-1),(j+1)+(k-1))=b(k);end
// i et j ZERO-BASED INDEXING
// A matrice column-major, doit avoir au moins i+N lignes et au moins j+N
// colonnes k doit etre real N : nombre de lignes (et de colonnes) de la sous
// matrice de A D : nombre de lignes de de A
template <typename T>
void kernel_memset_diag(T* dev_dst, T val, int rowDebut, int colDebut,
                        int numElements, int dev_dest_dim1);

// A_full = full(A); A_full est la matrice dense correspondant a la matrice
// creuse A
template <typename T>
void kernel_sparse2full(T* dev_dst, int* dev_src_rowind, int* dev_src_colind,
                        T* dev_src_values, int nnz, int src_dim1, int src_dim2);
// A(i,i) += valeur (pour i de 0 a length-1)
template <typename T>
void kernel_add_diag_const(T* d_cu1, T val, int dim1);
// v_diag = diag(M);
template <typename T>
void kernel_get_diag(T* dst_diag, T* src_M, int dim1);
// A(i,j) = v(i) * B(i,j) (pour i de 0 a B_dim1 et j de 0 a
// B_nb_elements/B_dim1)
template <typename T>
void kernel_diag_mult(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                      int nb_elements);
// A(i,j) = B(i,j) * v(j) (pour i de 0 a B_dim1 et j de 0 a
// B_nb_elements/B_dim1)
template <typename T>
void kernel_diag_mult2(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                       int nb_elements);

//    A1(i,j) = B(i,j) * v(i)
// et A2(i,j) = B(i,j) * v(i) * v(j)
// (pour i de 0 a B_dim1 et j de 0 a B_nb_elements/B_dim1)
template <typename T>
void kernel_diag_mult3(T* mat_dst1, T* mat_dst2, T* mat_src, T* vec_src,
                       int dim1, int nb_elements);
// A(i,j) = B(i+r1,j+c1) (pour i de 0 a nb_rows-1 et j de 0 a nb_col-1)
template <typename T>
void kernel_set_submatrix(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1,
                          int nb_rows, int nb_col);

// Linf(i,j) = Hinf(i-nb_n/2,j) (pour i de nb_n/2 a nb_n-1)
// si AR1 :
//      Linf(i,j) = atur(i)*Hinf(i,j) (pour i de 0 a nb_n/2-1)
// si AR2 :
//      Linf(i,j) = atur(i)*Hinf(i,j) + btur(i)*Hinf(i+nbn/2,j)   (pour i de 0 a
//      nb_n/2-1)
template <typename T>
void kernel_set_Linf(T* mat_Linf, T* mat_Hinf, T* atur, T* btur, int nb_n,
                     int nb_p, int ordreAR);

// mat += A1 (avec A1, diagonale par bloc : A1=[Atur Btur ; eye(nb_az)
// zeros(nb_az)] avec Atur (resp. Btur), matrice diagonale de dimension nb_az
// dont la diagonale est le vecteur atur (resp. btur))
template <typename T>
void kernel_add_A1(T* mat, T* atur, T* btur, int nb_az, int ordreAR);

#endif
