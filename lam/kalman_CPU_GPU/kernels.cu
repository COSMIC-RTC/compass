#include "kernels.h"

//prototypes
template<typename T> __global__ void Add(T *A, const T *B, int numElements);
template<typename T> __global__ void Sub(T *A, const T *B, int numElements);
template<typename T> __global__ void Mult(T *A, const T *B, int numElements);
__global__ void Div(real *A, const real *B, int numElements);

template<typename T> __global__ void AddConst(T *A, T val, int numElements);
template<typename T> __global__ void SubConst(T *A, T val, int numElements);
template<typename T> __global__ void MultConst(T *A, T val, int numElements);
__global__ void DivConst(real *A, real val, int numElements);

__global__ void Sqrt(real *A, int numElements);
__global__ void Inv(real *A, int numElements);

template<typename T> __global__ void Memcpy(T* d_cu_dst, const T* d_cu_src, int length);


template<typename T> __global__ void Memset(T* dev_dst, T valeur, int numElements);

__global__ void A1_2vec(real* Atur_d_cu, real* Btur_d_cu, int* colind, int* rowind, real* values, int nb_az, int nnz);

__global__ void MemcpyDiag(real *dev_dst, const real *dev_src, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);
__global__ void MemsetDiag(real *dev_dst, real val, int rowDebut, int colDebut, int numElements, int dev_dest_dim1);

__global__ void Sparse2full(real *dev_dst, int *dev_src_rowind, int *dev_src_colind, real* dev_src_values, int nnz, int src_dim1);

__global__ void AddDiagConst(real* dev_dst, real val, int dim1);

__global__ void GetDiag(real* dst_diag, real* src_M,  int dim1);

__global__ void DiagMult(real* mat_dst, real* mat_src, real* vec_src, int mat_src_dim1, int nb_elements);
__global__ void DiagMult2(real* mat_dst, real* mat_src, real* vec_src, int mat_src_dim1, int nb_elements);
__global__ void DiagMult3(real* mat_dst1, real* mat_dst2, real* mat_src, real* vec_src, int mat_src_dim1, int nb_elements);
__global__ void SetSubmatrix(real* mat_dst, real* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_elements);
__global__ void SetLinf(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_elements, int ordreAR);
//__global__ void SetLinfAR2(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_elements);
__global__ void AddA1(real* mat, real* atur, real* btur, int nb_az, int nb_elements, int ordreAR);


template<typename T> void kernel_add(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_sub(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_mult(T* d_cu1, const T* d_cu2, int length);
template<typename T> void kernel_add_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_sub_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_mult_const(T* d_cu1, T valeur, int length);
template<typename T> void kernel_memcpy(T* d_cu_dst, const T* d_cu_src, int length);
template<typename T> void kernel_memset(T* d_cu_dst, T valeur, int length);





// Kernel wrappers
extern "C" void kernel_add_int(int* d_cu1, const int* d_cu2, int length)
{kernel_add<int>(d_cu1, d_cu2, length);}
extern "C" void kernel_add_real(real* d_cu1, const real* d_cu2, int length)
{kernel_add<real>(d_cu1, d_cu2, length);}
template<typename T> void kernel_add(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Add<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
}

extern "C" void kernel_sub_int(int* d_cu1, const int* d_cu2, int length)
{kernel_sub<int>(d_cu1, d_cu2, length);}
extern "C" void kernel_sub_real(real* d_cu1, const real* d_cu2, int length)
{kernel_sub<real>(d_cu1, d_cu2, length);}
template<typename T> void kernel_sub(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Sub<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
}

extern "C" void kernel_mult_int(int* d_cu1, const int* d_cu2, int length)
{kernel_mult<int>(d_cu1, d_cu2, length);}
extern "C" void kernel_mult_real(real* d_cu1, const real* d_cu2, int length)
{kernel_mult<real>(d_cu1, d_cu2, length);}
template<typename T> void kernel_mult(T* d_cu1, const T* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Mult<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
}

extern "C" void kernel_div(real* d_cu1, const real* d_cu2, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Div<<<blocksPerGrid,threadsPerBlock>>>(d_cu1, d_cu2, length);
}

extern "C" void kernel_add_const_int(int* d_cu1, int valeur, int length)
{kernel_add_const<int>(d_cu1, valeur, length);}
extern "C" void kernel_add_const_real(real* d_cu1, real valeur, int length)
{kernel_add_const<real>(d_cu1, valeur, length);}
template<typename T> void kernel_add_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	AddConst<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
}

extern "C" void kernel_sub_const_int(int* d_cu1, int valeur, int length)
{kernel_sub_const<int>(d_cu1, valeur, length);}
extern "C" void kernel_sub_const_real(real* d_cu1, real valeur, int length)
{kernel_sub_const<real>(d_cu1, valeur, length);}
template<typename T> void kernel_sub_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	SubConst<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
}

extern "C" void kernel_mult_const_int(int* d_cu1, int valeur, int length)
{kernel_mult_const<int>(d_cu1, valeur, length);}
extern "C" void kernel_mult_const_real(real* d_cu1, real valeur, int length)
{kernel_mult_const<real>(d_cu1, valeur, length);}
template<typename T> void kernel_mult_const(T* d_cu1, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	MultConst<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
}

extern "C" void kernel_div_const(real* d_cu1, real valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	DivConst<<<blocksPerGrid,threadsPerBlock>>>(d_cu1, valeur, length);
}

extern "C" void kernel_sqrt(real* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Sqrt<<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
}
extern "C" void kernel_inv(real* d_cu1, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Inv<<<blocksPerGrid,threadsPerBlock>>>(d_cu1, length);
}





extern "C" void kernel_memcpy_real(real* d_cu_dst, const real* d_cu_src, int length)
{kernel_memcpy<real>(d_cu_dst, d_cu_src, length);}
extern "C" void kernel_memcpy_int(int* d_cu_dst, const int* d_cu_src, int length)
{kernel_memcpy<int>(d_cu_dst, d_cu_src, length);}
template<typename T> void kernel_memcpy(T* d_cu_dst, const T* d_cu_src, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Memcpy<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, d_cu_src, length);
}


extern "C" void kernel_memset_int(int* d_cu_dst, const int valeur, int length)
{kernel_memset<int>(d_cu_dst, valeur, length);}
extern "C" void kernel_memset_real(real* d_cu_dst, const real valeur, int length)
{kernel_memset<real>(d_cu_dst, valeur, length);}
template<typename T> void kernel_memset(T* d_cu_dst, T valeur, int length)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
	Memset<T><<<blocksPerGrid,threadsPerBlock>>>(d_cu_dst, valeur, length);
}	

extern "C" void kernel_A1_2vec(real* Atur_d_cu, real* Btur_d_cu, int* colind, int* rowind, real* values, int nb_az, int nnz)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
	A1_2vec<<<blocksPerGrid,threadsPerBlock>>>(Atur_d_cu, Btur_d_cu, colind, rowind, values, nb_az, nnz);
}	
extern "C" void kernel_memcpy_diag(real *dev_dst, const real *dev_src, int rowDebut, int colDebut, int numElements, int dev_dest_dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
	MemcpyDiag<<<blocksPerGrid,threadsPerBlock>>>(dev_dst, dev_src, rowDebut, colDebut, numElements, dev_dest_dim1);
}

extern "C" void kernel_memset_diag(real *dev_dst, real val, int rowDebut, int colDebut, int numElements, int dev_dest_dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (numElements + threadsPerBlock - 1) / threadsPerBlock;
	MemsetDiag<<<blocksPerGrid,threadsPerBlock>>>(dev_dst, val, rowDebut, colDebut, numElements, dev_dest_dim1);
}

extern "C" void kernel_sparse2full(real *dev_dst, int *dev_src_rowind, int *dev_src_colind, real* dev_src_values, int nnz, int src_dim1, int src_dim2)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nnz + threadsPerBlock - 1) / threadsPerBlock;
	kernel_memset_real(dev_dst, 0.0, src_dim1*src_dim2);
	Sparse2full<<<blocksPerGrid,threadsPerBlock>>>(dev_dst, dev_src_rowind, dev_src_colind, dev_src_values, nnz, src_dim1);
}

extern "C" void kernel_add_diag_const(real* d_cu1, real val, int dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
	AddDiagConst<<<blocksPerGrid,threadsPerBlock>>>(d_cu1, val, dim1);
}

extern "C" void kernel_get_diag(real* dst_diag, real* src_M, int dim1)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (dim1 + threadsPerBlock - 1) / threadsPerBlock;
	GetDiag<<<blocksPerGrid,threadsPerBlock>>>(dst_diag, src_M, dim1);
}

extern "C" void kernel_diag_mult(real* mat_dst, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	DiagMult<<<blocksPerGrid,threadsPerBlock>>>(mat_dst, mat_src, vec_src, dim1, nb_elements);
}

extern "C" void kernel_diag_mult2(real* mat_dst, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	DiagMult2<<<blocksPerGrid,threadsPerBlock>>>(mat_dst, mat_src, vec_src, dim1, nb_elements);
}

extern "C" void kernel_diag_mult3(real* mat_dst1, real* mat_dst2, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	DiagMult3<<<blocksPerGrid,threadsPerBlock>>>(mat_dst1, mat_dst2, mat_src, vec_src, dim1, nb_elements);
}

extern "C" void kernel_set_submatrix(real* mat_dst, real* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_col)
{
	int threadsPerBlock = 256;
	int nb_elements = nb_rows*nb_col;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	SetSubmatrix<<<blocksPerGrid,threadsPerBlock>>>(mat_dst, mat_src, src_dim1, r1, c1, nb_rows, nb_elements);
}


extern "C" void kernel_set_Linf(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_p, int ordreAR)
{
	int threadsPerBlock = 256;
	int nb_elements = nb_n*nb_p;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	SetLinf<<<blocksPerGrid,threadsPerBlock>>>(mat_Linf, mat_Hinf, atur, btur, nb_n, nb_elements, ordreAR);

}

/*extern "C" void kernel_set_Linf_AR2(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_p)
{
	int threadsPerBlock = 256;
	int nb_elements = nb_n*nb_p;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	SetLinfAR2<<<blocksPerGrid,threadsPerBlock>>>(mat_Linf, mat_Hinf, atur, btur, nb_n, nb_elements);

}*/

extern "C" void kernel_add_A1(real* mat, real* atur, real* btur, int nb_az, int ordreAR)
{
	int threadsPerBlock = 256;
	int nb_elements = 4*nb_az*nb_az;
	int blocksPerGrid = (nb_elements + threadsPerBlock - 1) / threadsPerBlock;
	AddA1<<<blocksPerGrid,threadsPerBlock>>>(mat, atur, btur, nb_az, nb_elements, ordreAR);
}


// Kernels

template<typename T> __global__ void Add(T* A, const T* B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] + B[i];
    }
}
template<typename T> __global__ void Sub(T* A, const T* B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] - B[i];
    }
}
template<typename T> __global__ void Mult(T *A, const T *B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] * B[i];
    }
}
__global__ void Div(real *A, const real *B, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] / B[i];
    }
}

template<typename T> __global__ void AddConst(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] + val;
    }
}
template<typename T> __global__ void SubConst(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] - val;
    }
}
template<typename T> __global__ void MultConst(T *A, T val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] * val;
    }
}
__global__ void DivConst(real* A, real val, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = A[i] / val;
    }
}

__global__ void Sqrt(real *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = sqrt(A[i]) ;
    }
}
__global__ void Inv(real *A, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        A[i] = 1/A[i] ;
    }
}


template<typename T> __global__ void Memcpy(T* dev_dst, const T* dev_src, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        dev_dst[i] = dev_src[i] ;
    }
}

template<typename T> __global__ void Memset(T* dev_dst, T valeur, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        dev_dst[i] = valeur ;
    }
}


__global__ void A1_2vec(real* Atur_d_cu, real* Btur_d_cu, int* colind, int* rowind, real* values, int nb_az, int nnz)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    int col = colind[i];
    int row = rowind[i];
    real val = values[i];
    


    if(i<nnz)
    {	
	// si col<=nb_az && row<=nb_az
	if (col<=nb_az && row<=nb_az)
	{
	    if(col==row) Atur_d_cu[col-1]=val;
	    else Atur_d_cu[row-1]=0;

	}
	//si col>nb_az && row<=nb_az
	else if(col>nb_az  && row<=nb_az)
        {
	    if((col-nb_az)==row) Btur_d_cu[row-1]=val;
            else Btur_d_cu[row-1]=0;

	}

    }
}


__global__ void MemcpyDiag(real *dev_dst, const real *dev_src, int rowDebut, int colDebut, int numElements, int dev_dest_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = (i+colDebut) * dev_dest_dim1 + (i+rowDebut);

    if (i < numElements)
    {
        dev_dst[j] = dev_src[i] ;
    }
}

__global__ void MemsetDiag(real *dev_dst, real val, int rowDebut, int colDebut, int numElements, int dev_dest_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = (i+colDebut) * dev_dest_dim1 + (i+rowDebut);

    if (i < numElements)
    {
        dev_dst[j] = val ;
    }
}

__global__ void Sparse2full(real *dev_dst, int *dev_src_rowind, int *dev_src_colind, real* dev_src_values, int nnz, int src_dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < nnz)
    {
    int j = (dev_src_colind[i]-1) * src_dim1 + (dev_src_rowind[i]-1);
        dev_dst[j] =  dev_src_values[i];
    }
}

__global__ void AddDiagConst(real *A, real val, int dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = i*(dim1+1);

    if (i < dim1)
    {
        A[j] = A[j] + val;
    }
}

__global__ void GetDiag(real* dst_diag, real* src_M, int dim1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = i*(dim1+1);

    if (i < dim1)
    {
        dst_diag[i] = src_M[j];
    }
}

__global__ void DiagMult(real* mat_dst, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
        int col = (int)(ind/dim1) ;
        int ligne = ind - dim1*col;

        mat_dst[ind] = mat_src[ind] * vec_src[ligne];
    }
}

__global__ void DiagMult2(real* mat_dst, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
        int col = (int)(ind/dim1) ;
        int ligne = ind - dim1*col;

        mat_dst[ind] = mat_src[ind] * vec_src[ligne] * vec_src[col];
    }
}

__global__ void DiagMult3(real* mat_dst1, real* mat_dst2, real* mat_src, real* vec_src, int dim1, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
        int col = (int)(ind/dim1) ;
        int ligne = ind - dim1*col;
	real val = mat_src[ind] * vec_src[ligne];
        mat_dst1[ind] = val;
	mat_dst2[ind] = val * vec_src[col];
    }
}

__global__ void SetSubmatrix(real* mat_dst, real* mat_src, int src_dim1, int r1, int c1, int nb_rows, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
	    int col = (int)(ind/nb_rows) ;
            int ligne = ind - nb_rows*col;
	    mat_dst[ind] = mat_src[(col+c1)*src_dim1+(ligne+r1)];        
    }
}

/*__global__ void SetLinfAR1(real* mat_Linf, real* mat_Hinf, real* atur, int nb_n, int nb_elements)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
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
    }
}*/

__global__ void SetLinf(real* mat_Linf, real* mat_Hinf, real* atur, real* btur, int nb_n, int nb_elements, int ordreAR)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
	int col = (int)(ind/nb_n) ;
        int ligne = ind - nb_n*col;
	int nb_az = nb_n/2;
	if(ligne < nb_az)
	{
             mat_Linf[ind] = atur[ligne] * mat_Hinf[ind];
	     if (ordreAR==2) 
		     mat_Linf[ind] += btur[ligne] * mat_Hinf[col*nb_n + ligne+nb_az];
	}
	else
	{
	     mat_Linf[ind] = mat_Hinf[col*nb_n + ligne-nb_az];
	}

    }
}


__global__ void AddA1(real* mat, real* atur, real* btur, int nb_az, int nb_elements, int ordreAR)
{
    int ind = blockDim.x*blockIdx.x+threadIdx.x;
    if(ind<nb_elements)
    {
	int nb_n = 2*nb_az;
	int col = (int)(ind/nb_n) ;
        int ligne = ind - nb_n*col;
	if (ligne < nb_az)
	{
		if (col == ligne) 
			mat[ind] += atur[ligne];
		else if ((ordreAR==2) && (col == (ligne+nb_az))) 
			mat[ind] += btur[ligne];
	}
	else if(ligne == col+nb_az) 
		mat[ind] += 1.0;

    }

}



