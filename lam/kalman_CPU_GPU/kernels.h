
#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <unistd.h>
#include "kp_real.h"

// d_cu1[i] += d_cu2[i] (pour i de 0 a length-1)
extern "C" void kernel_add_int(int* d_cu1, const int* d_cu2, int length);  
extern "C" void kernel_add_real(real* d_cu1, const real* d_cu2, int length);  
// d_cu1[i] -= d_cu2[i] (pour i de 0 a length-1)
extern "C" void kernel_sub_int(int* d_cu1, const int* d_cu2, int length);  
extern "C" void kernel_sub_real(real* d_cu1, const real* d_cu2, int length); 
// d_cu1[i] *= d_cu2[i] (pour i de 0 a length-1)
extern "C" void kernel_mult_int(int* d_cu1, const int* d_cu2, int length); 
extern "C" void kernel_mult_real(real* d_cu1, const real* d_cu2, int length); 
// d_cu1[i] /= d_cu2[i] (pour i de 0 a length-1)
extern "C" void kernel_div(real* d_cu1, const real* d_cu2, int length);  
// d_cu1[i] += valeur (pour i de 0 a length-1)
extern "C" void kernel_add_const_int(int* d_cu1, int valeur, int length);  
extern "C" void kernel_add_const_real(real* d_cu1, real valeur, int length); 
// d_cu1[i] -= valeur (pour i de 0 a length-1)
extern "C" void kernel_sub_const_int(int* d_cu1, int valeur, int length);  
extern "C" void kernel_sub_const_real(real* d_cu1, real valeur, int length); 
// d_cu1[i] *= valeur (pour i de 0 a length-1)
extern "C" void kernel_mult_const_int(int* d_cu1, int valeur, int length); 
extern "C" void kernel_mult_const_real(real* d_cu1, real valeur, int length); 
// d_cu1[i] /= valeur (pour i de 0 a length-1)
extern "C" void kernel_div_const(real* d_cu1, real valeur, int length);  

// d_cu1[i] = sqrt(d_cu1[i]) (pour i de 0 a length-1)
extern "C" void kernel_sqrt(real* d_cu1, int length); 
// d_cu1[i] = 1/d_cu1[i]     (pour i de 0 a length-1)
extern "C" void kernel_inv(real* d_cu1, int length);  



//d_cu_dst[i] = d_cu_src[i] (pour i de 0 a length-1)
extern "C" void kernel_memcpy_real(real* d_cu_dst, const real* d_cu_src, int length); 
//d_cu_dst[i] = d_cu_src[i] (pour i de 0 a length-1)
extern "C" void kernel_memcpy_int(int* d_cu_dst, const int* d_cu_src, int length);



//d_cu_dst[i] = valeur (pour i de 0 a length-1)
extern "C" void kernel_memset_real(real* d_cu_dst, real valeur, int length);
//d_cu_dst[i] = valeur (pour i de 0 a length-1)
extern "C" void kernel_memset_int(int* d_cu_dst, int valeur, int length); 

// Atur_d_cu = diag(A(0:nb_az-1,0:nb_az-1)) et Btur_d_cu = diag(A(0:nb_az-1,nb_az:2*nb_az-1)) (avec A matrice creuse)
extern "C" void kernel_A1_2vec(real* Atur_d_cu, real* Btur_d_cu, int* A_colind, int* A_rowind, real* A_values, int nb_az, int nnz);

// kernel_memcpy_diag(A, b, i, j, N, D). diag(A(i:i+N,j:j+N)) = b
// Equivalent Matlab : for k=1:N,A((i+1)+(k-1),(j+1)+(k-1))=b(k);end
// i et j ZERO-BASED INDEXING
// A matrice column-major, doit avoir au moins i+N lignes et au moins j+N colonnes
// b doit etre un vecteur de taille N
// N : nombre de lignes (et de colonnes) de la sous matrice de A
// D : nombre de lignes de de A
extern "C" void kernel_memcpy_diag(real *dev_dst, const real *dev_src, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);

// kernel_memset_diag(A, k, i, j, N, D). diag(A(i:i+N,j:j+N)) = k (avec k reel)
// Equivalent Matlab : for k=1:N,A((i+1)+(k-1),(j+1)+(k-1))=b(k);end
// i et j ZERO-BASED INDEXING
// A matrice column-major, doit avoir au moins i+N lignes et au moins j+N colonnes
// k doit etre real
// N : nombre de lignes (et de colonnes) de la sous matrice de A
// D : nombre de lignes de de A
extern "C" void kernel_memset_diag(real *dev_dst, real val, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);

// A_full = full(A); A_full est la matrice dense correspondant a la matrice creuse A
extern "C" void kernel_sparse2full(real* A_full, int* A_rowind, int* A_colind, real* A_values, int A_nnz, int A_dim1, int A_dim2);

// A(i,i) += valeur (pour i de 0 a length-1)
extern "C" void kernel_add_diag_const(real* A, real valeur, int A_dim1);

// v_diag = diag(M);
extern "C" void kernel_get_diag(real* v_diag, real* M_d_cu, int M_dim1);

// A(i,j) = v(i) * B(i,j) (pour i de 0 a B_dim1 et j de 0 a B_nb_elements/B_dim1)
extern "C" void kernel_diag_mult(real* A_d_cu, real* B_d_cu, real* v_d_cu, int B_dim1, int B_nb_elements);

// A(i,j) = B(i,j) * v(j) (pour i de 0 a B_dim1 et j de 0 a B_nb_elements/B_dim1)
extern "C" void kernel_diag_mult2(real* A_d_cu, real* B_d_cu, real* v_d_cu, int B_dim1, int B_nb_elements);
//    A1(i,j) = B(i,j) * v(i) 
// et A2(i,j) = B(i,j) * v(i) * v(j)
// (pour i de 0 a B_dim1 et j de 0 a B_nb_elements/B_dim1)
extern "C" void kernel_diag_mult3(real* A1_d_cu, real* A2_d_cu, real* B_d_cu, real* v_d_cu, int B_dim1, int B_nb_elements);
// A(i,j) = B(i+r1,j+c1) (pour i de 0 a nb_rows-1 et j de 0 a nb_col-1)
extern "C" void kernel_set_submatrix(real* A, real* B, int B_dim1, int r1, int c1, int nb_rows, int nb_col);

// Linf(i,j) = Hinf(i-nb_n/2,j) (pour i de nb_n/2 a nb_n-1)
// si AR1 :
//      Linf(i,j) = atur(i)*Hinf(i,j) (pour i de 0 a nb_n/2-1)
// si AR2 : 
//      Linf(i,j) = atur(i)*Hinf(i,j) + btur(i)*Hinf(i+nbn/2,j)   (pour i de 0 a nb_n/2-1)
extern "C" void kernel_set_Linf(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_p, int ordreAR);
//extern "C" void kernel_set_Linf_AR2(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_p);

// mat += A1 (avec A1, diagonale par bloc (atur, btur, Id))
extern "C" void kernel_add_A1(real* mat, real* atur, real* btur, int nb_az, int ordreAR);


#endif
