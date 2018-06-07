#include "kernels.h"

// prototypes
template <typename T>
__global__ void Add_lam(T* A, const T* B, int numElements);
template <typename T>
__global__ void Sub_lam(T* A, const T* B, int numElements);
template <typename T>
__global__ void Mult_lam(T* A, const T* B, int numElements);
template <typename T>
__global__ void Div_lam(T* A, const T* B, int numElements);

template <typename T>
__global__ void AddConst_lam(T* A, T val, int numElements);
template <typename T>
__global__ void SubConst_lam(T* A, T val, int numElements);
template <typename T>
__global__ void MultConst_lam(T* A, T val, int numElements);
template <typename T>
__global__ void DivConst_lam(T* A, T val, int numElements);

template <typename T>
__global__ void Sqrt_lam(T* A, int numElements);
template <typename T>
__global__ void Inv_lam(T* A, int numElements);

template <typename T, typename U>
__global__ void Memcpy_lam(T* d_cu_dst, const U* d_cu_src, int length);

template <typename T>
__global__ void Memset_lam(T* dev_dst, T valeur, int numElements);

template <typename T>
__global__ void A1_2vec_lam(T* Atur_d_cu, T* Btur_d_cu, int* colind,
                            int* rowind, T* values, int nb_az, int nnz);

template <typename T>
__global__ void MemcpyDiag_lam(T* dev_dst, const T* dev_src, int rowDebut,
                               int colDebut, int numElements,
                               int dev_dest_dim1);
template <typename T>
__global__ void MemsetDiag_lam(T* dev_dst, T val, int rowDebut, int colDebut,
                               int numElements, int dev_dest_dim1);

template <typename T>
__global__ void Sparse2full_lam(T* dev_dst, int* dev_src_rowind,
                                int* dev_src_colind, T* dev_src_values, int nnz,
                                int src_dim1);

template <typename T>
__global__ void AddDiagConst_lam(T* dev_dst, T val, int dim1);

template <typename T>
__global__ void GetDiag_lam(T* dst_diag, T* src_M, int dim1);

template <typename T>
__global__ void DiagMult_lam(T* mat_dst, T* mat_src, T* vec_src,
                             int mat_src_dim1, int nb_elements);
template <typename T>
__global__ void DiagMult2_lam(T* mat_dst, T* mat_src, T* vec_src,
                              int mat_src_dim1, int nb_elements);
template <typename T>
__global__ void DiagMult3_lam(T* mat_dst1, T* mat_dst2, T* mat_src, T* vec_src,
                              int mat_src_dim1, int nb_elements);
template <typename T>
__global__ void SetSubmatrix_lam(T* mat_dst, T* mat_src, int src_dim1, int r1,
                                 int c1, int nb_rows, int nb_elements);
template <typename T>
__global__ void SetLinf_lam(T* mat_Linf, T* mat_Hinf, T* atur, T* btur,
                            int nb_n, int nb_elements, int ordreAR);
// template<typename T> __global__ void SetLinfAR2(T* mat_Linf, T* mat_Hinf, T*
// atur, T* btur, int nb_n, int nb_elements);
template <typename T>
__global__ void AddA1_lam(T* mat, T* atur, T* btur, int nb_az, int nb_elements,
                          int ordreAR);

// Kernel wrappers
template <typename T>
void kernel_add(T* d_cu1, const T* d_cu2, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Add_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, d_cu2, length);
}
template <typename T>
void kernel_sub(T* d_cu1, const T* d_cu2, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Sub_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, d_cu2, length);
}
template <typename T>
void kernel_mult(T* d_cu1, const T* d_cu2, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Mult_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, d_cu2, length);
}
template <typename T>
void kernel_div(T* d_cu1, const T* d_cu2, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Div_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, d_cu2, length);
}
template <typename T>
void kernel_add_const(T* d_cu1, T valeur, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  AddConst_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, valeur, length);
}
template <typename T>
void kernel_sub_const(T* d_cu1, T valeur, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  SubConst_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, valeur, length);
}
template <typename T>
void kernel_mult_const(T* d_cu1, T valeur, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  MultConst_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, valeur, length);
}
template <typename T>
void kernel_div_const(T* d_cu1, T valeur, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  DivConst_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, valeur, length);
}
template <typename T>
void kernel_sqrt(T* d_cu1, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Sqrt_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, length);
}

template <typename T>
void kernel_inv(T* d_cu1, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Inv_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, length);
}

template <typename T, typename U>
void kernel_memcpy(T* d_cu_dst, const U* d_cu_src, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Memcpy_lam<T, U>
      <<<blocksPerGrid, threadsPerBlock>>>(d_cu_dst, d_cu_src, length);
}

template <typename T>
void kernel_memset(T* d_cu_dst, T valeur, int length) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
  Memset_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu_dst, valeur, length);
}
template <typename T>
void kernel_A1_2vec(T* Atur_d_cu, T* Btur_d_cu, int* colind, int* rowind,
                    T* values, int nb_az, int nnz) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
  A1_2vec_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      Atur_d_cu, Btur_d_cu, colind, rowind, values, nb_az, nnz);
}
template <typename T>
void kernel_memcpy_diag(T* dev_dst, const T* dev_src, int rowDebut,
                        int colDebut, int numElements, int dev_dest_dim1) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
  MemcpyDiag_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      dev_dst, dev_src, rowDebut, colDebut, numElements, dev_dest_dim1);
}
template <typename T>
void kernel_memset_diag(T* dev_dst, T val, int rowDebut, int colDebut,
                        int numElements, int dev_dest_dim1) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
  MemsetDiag_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      dev_dst, val, rowDebut, colDebut, numElements, dev_dest_dim1);
}
template <typename T>
void kernel_sparse2full(T* dev_dst, int* dev_src_rowind, int* dev_src_colind,
                        T* dev_src_values, int nnz, int src_dim1,
                        int src_dim2) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
  T zero = (T)0.0;
  kernel_memset<T>(dev_dst, zero, src_dim1 * src_dim2);
  Sparse2full_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      dev_dst, dev_src_rowind, dev_src_colind, dev_src_values, nnz, src_dim1);
}
template <typename T>
void kernel_add_diag_const(T* d_cu1, T val, int dim1) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
  AddDiagConst_lam<T><<<blocksPerGrid, threadsPerBlock>>>(d_cu1, val, dim1);
}
template <typename T>
void kernel_get_diag(T* dst_diag, T* src_M, int dim1) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
  GetDiag_lam<T><<<blocksPerGrid, threadsPerBlock>>>(dst_diag, src_M, dim1);
}
template <typename T>
void kernel_diag_mult(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                      int nb_elements) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  DiagMult_lam<T><<<blocksPerGrid, threadsPerBlock>>>(mat_dst, mat_src, vec_src,
                                                      dim1, nb_elements);
}
template <typename T>
void kernel_diag_mult2(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                       int nb_elements) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  DiagMult2_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      mat_dst, mat_src, vec_src, dim1, nb_elements);
}
template <typename T>
void kernel_diag_mult3(T* mat_dst1, T* mat_dst2, T* mat_src, T* vec_src,
                       int dim1, int nb_elements) {
  int threadsPerBlock = 256;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  DiagMult3_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      mat_dst1, mat_dst2, mat_src, vec_src, dim1, nb_elements);
}
template <typename T>
void kernel_set_submatrix(T* mat_dst, T* mat_src, int src_dim1, int r1, int c1,
                          int nb_rows, int nb_col) {
  int threadsPerBlock = 256;
  int nb_elements = nb_rows * nb_col;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  SetSubmatrix_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      mat_dst, mat_src, src_dim1, r1, c1, nb_rows, nb_elements);
}

template <typename T>
void kernel_set_Linf(T* mat_Linf, T* mat_Hinf, T* atur, T* btur, int nb_n,
                     int nb_p, int ordreAR) {
  int threadsPerBlock = 256;
  int nb_elements = nb_n * nb_p;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  SetLinf_lam<T><<<blocksPerGrid, threadsPerBlock>>>(
      mat_Linf, mat_Hinf, atur, btur, nb_n, nb_elements, ordreAR);
}

/*
extern "C" void kernel_set_Linf_AR2_KFPP(KFPP* mat_Linf, KFPP* mat_Hinf, KFPP*
atur, KFPP* btur, int nb_n, int nb_p) {kernel_set_Linf_AR2_tpl<KFPP>(mat_Linf,
mat_Hinf, atur, btur, nb_n, nb_p);} extern "C" void
kernel_set_Linf_AR2_double(double* mat_Linf, double* mat_Hinf, double* atur,
double* btur, int nb_n, int nb_p) {kernel_set_Linf_AR2_tpl<double>(mat_Linf,
mat_Hinf, atur, btur, nb_n, nb_p);} template<typename T> void
kernel_set_Linf_AR2(T* mat_Linf, T* mat_Hinf, T* atur, T* btur, int nb_n, int
nb_p)
{
        int threadsPerBlock = 256;
        int nb_elements = nb_n*nb_p;
        int blocksPerGrid = (nb_elements + threadsPerBlock - 1) /
threadsPerBlock; SetLinfAR2<T><<<blocksPerGrid,threadsPerBlock>>>(mat_Linf,
mat_Hinf, atur, btur, nb_n, nb_elements);

}*/
template <typename T>
void kernel_add_A1(T* mat, T* atur, T* btur, int nb_az, int ordreAR) {
  int threadsPerBlock = 256;
  int nb_elements = 4 * nb_az * nb_az;
  int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
  AddA1_lam<T><<<blocksPerGrid, threadsPerBlock>>>(mat, atur, btur, nb_az,
                                                   nb_elements, ordreAR);
}

// Kernels

template <typename T>
__global__ void Add_lam(T* A, const T* B, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] + B[i];
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void Sub_lam(T* A, const T* B, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] - B[i];
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void Mult_lam(T* A, const T* B, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] * B[i];
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void Div_lam(T* A, const T* B, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] / B[i];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void AddConst_lam(T* A, T val, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] + val;
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void SubConst_lam(T* A, T val, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] - val;
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void MultConst_lam(T* A, T val, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] * val;
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void DivConst_lam(T* A, T val, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = A[i] / val;
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void Sqrt_lam(T* A, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = sqrt(A[i]);
    i += gridDim.x * blockDim.x;
  }
}
template <typename T>
__global__ void Inv_lam(T* A, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    A[i] = 1 / A[i];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T, typename U>
__global__ void Memcpy_lam(T* dev_dst, const U* dev_src, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    dev_dst[i] = (T)dev_src[i];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void Memset_lam(T* dev_dst, T valeur, int numElements) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    dev_dst[i] = valeur;
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void A1_2vec_lam(T* Atur_d_cu, T* Btur_d_cu, int* colind,
                            int* rowind, T* values, int nb_az, int nnz) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < nnz) {
    int col = colind[i];
    int row = rowind[i];
    KFPP val = values[i];

    // si col<=nb_az && row<=nb_az
    if (col <= nb_az && row <= nb_az) {
      if (col == row)
        Atur_d_cu[col - 1] = val;
      else
        Atur_d_cu[row - 1] = 0;

    }
    // si col>nb_az && row<=nb_az
    else if (col > nb_az && row <= nb_az) {
      if ((col - nb_az) == row)
        Btur_d_cu[row - 1] = val;
      else
        Btur_d_cu[row - 1] = 0;
    }
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void MemcpyDiag_lam(T* dev_dst, const T* dev_src, int rowDebut,
                               int colDebut, int numElements,
                               int dev_dest_dim1) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    int j = (i + colDebut) * dev_dest_dim1 + (i + rowDebut);
    dev_dst[j] = dev_src[i];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void MemsetDiag_lam(T* dev_dst, T val, int rowDebut, int colDebut,
                               int numElements, int dev_dest_dim1) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < numElements) {
    int j = (i + colDebut) * dev_dest_dim1 + (i + rowDebut);
    dev_dst[j] = val;
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void Sparse2full_lam(T* dev_dst, int* dev_src_rowind,
                                int* dev_src_colind, T* dev_src_values, int nnz,
                                int src_dim1) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < nnz) {
    int j = (dev_src_colind[i] - 1) * src_dim1 + (dev_src_rowind[i] - 1);
    dev_dst[j] = dev_src_values[i];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void AddDiagConst_lam(T* A, T val, int dim1) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < dim1) {
    int j = i * (dim1 + 1);
    A[j] = A[j] + val;
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void GetDiag_lam(T* dst_diag, T* src_M, int dim1) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;

  while (i < dim1) {
    int j = i * (dim1 + 1);
    dst_diag[i] = src_M[j];
    i += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void DiagMult_lam(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                             int nb_elements) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int col = (int)(ind / dim1);
    int ligne = ind - dim1 * col;

    mat_dst[ind] = mat_src[ind] * vec_src[ligne];
    ind += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void DiagMult2_lam(T* mat_dst, T* mat_src, T* vec_src, int dim1,
                              int nb_elements) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int col = (int)(ind / dim1);
    int ligne = ind - dim1 * col;

    mat_dst[ind] = mat_src[ind] * vec_src[ligne] * vec_src[col];
    ind += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void DiagMult3_lam(T* mat_dst1, T* mat_dst2, T* mat_src, T* vec_src,
                              int dim1, int nb_elements) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int col = (int)(ind / dim1);
    int ligne = ind - dim1 * col;
    KFPP val = mat_src[ind] * vec_src[ligne];
    mat_dst1[ind] = val;
    mat_dst2[ind] = val * vec_src[col];
    ind += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void SetSubmatrix_lam(T* mat_dst, T* mat_src, int src_dim1, int r1,
                                 int c1, int nb_rows, int nb_elements) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int col = (int)(ind / nb_rows);
    int ligne = ind - nb_rows * col;
    mat_dst[ind] = mat_src[(col + c1) * src_dim1 + (ligne + r1)];
    ind += gridDim.x * blockDim.x;
  }
}

/*template<typename T> __global__ void SetLinfAR1(T* mat_Linf, T* mat_Hinf, T*
atur, int nb_n, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    while(ind<nb_elements)
    {
        int col = (int)(ind/nb_n) ;
        int ligne = ind - nb_n*col;
        int nb_az = nb/2;
        if(ligne < nb_az)
        {
             mat_Linf[ind] = atur[ligne] * mat_Hinf[ind];
        }
        else
        {
             mat_Linf[ind] = mat_Hinf[col*nb_n + ligne-nb_az];
        }
        ind += gridDim.x * blockDim.x;
    }
}*/

template <typename T>
__global__ void SetLinf_lam(T* mat_Linf, T* mat_Hinf, T* atur, T* btur,
                            int nb_n, int nb_elements, int ordreAR) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int col = (int)(ind / nb_n);
    int ligne = ind - nb_n * col;
    int nb_az = nb_n / 2;
    if (ligne < nb_az) {
      mat_Linf[ind] = atur[ligne] * mat_Hinf[ind];
      if (ordreAR == 2)
        mat_Linf[ind] += btur[ligne] * mat_Hinf[col * nb_n + ligne + nb_az];
    } else {
      mat_Linf[ind] = mat_Hinf[col * nb_n + ligne - nb_az];
    }
    ind += gridDim.x * blockDim.x;
  }
}

template <typename T>
__global__ void AddA1_lam(T* mat, T* atur, T* btur, int nb_az, int nb_elements,
                          int ordreAR) {
  int ind = blockDim.x * blockIdx.x + threadIdx.x;
  while (ind < nb_elements) {
    int nb_n = 2 * nb_az;
    int col = (int)(ind / nb_n);
    int ligne = ind - nb_n * col;
    if (ligne < nb_az) {
      if (col == ligne)
        mat[ind] += atur[ligne];
      else if ((ordreAR == 2) && (col == (ligne + nb_az)))
        mat[ind] += btur[ligne];
    } else if (ligne == col + nb_az)
      mat[ind] += 1.0;
    ind += gridDim.x * blockDim.x;
  }
}

#include "kernel_def.hu"
