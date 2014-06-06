//kp_cu_matrix.cpp


#include "kp_cu_matrix.h"
#include "kp_smatrix.h"
#include "mkl_blas.h"
#include "mkl.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;



kp_cu_matrix::kp_cu_matrix()
{
	_create(0,0);
}

kp_cu_matrix::kp_cu_matrix(int dim1_, int dim2_)
{
	_create(dim1_, dim2_);
}

kp_cu_matrix::kp_cu_matrix(const kp_cu_matrix& m)
{
	_create(m.dim1, m.dim2);
	//kp_cu_cudaMemcpy(d_cu, m.d_cu, sizeof(real) * dim1 * dim2, cudaMemcpyDeviceToDevice); 
    	kernel_memcpy_real(d_cu, m.d_cu, dim1*dim2);	
}

kp_cu_matrix::kp_cu_matrix(const kp_matrix& m)
{
	_create(m.dim1, m.dim2);
	kp_cu_cudaMemcpy(d_cu, m.d, sizeof(real) * dim1 * dim2, cudaMemcpyHostToDevice);     
}

kp_cu_matrix::kp_cu_matrix(const vector< vector<real> >& m )
{
	size_t d1,d2;
	if (m.size() == 0)
	{
		d1 = 0;
		d2 = 0;
	}
	else
	{
		d1 = m.size();
		d2 = m[0].size();
		for (size_t i = 0 ; i < m.size(); i++)
		if (d2 != m[i].size())
		{
			cout<<"error | kp_cu_matrix::kp_cu_matrix | vector< vector<double> > is not matrix"<<endl;
			exit(EXIT_FAILURE);
		}	
	}
	_create(d1,d2);

	//real* d_tmp = new real[d1*d2];
	for (size_t i = 0 ; i < m.size() ; i++)
		for (size_t j = 0 ; j < m[i].size() ; j++)
			//d_tmp[j*d1+i] = m[i][j];
			kp_cu_cudaMemcpy(d_cu+(j*d1+i), &(m[i][j]), sizeof(real), cudaMemcpyHostToDevice);

	//kp_cu_cudaMemcpy(d_cu, d_tmp, sizeof(real) * dim1 * dim2, cudaMemcpyHostToDevice );     

	//delete[] d_tmp;
}

kp_cu_matrix::~kp_cu_matrix()
{
   _clear();
}


void kp_cu_matrix::resize(int dim1_, int dim2_)
{

   if (dim1 * dim2 != dim1_ * dim2_)
     {
	_clear();
	_create(dim1_, dim2_);
     }
   else
     {
	dim1 = dim1_;
	dim2 = dim2_;
     }
}


// kp_cu_matrix::at()


void kp_cu_matrix::operator=(const kp_cu_matrix&m)
{
   resize(m.dim1, m.dim2);
   //memcpy(d, m.d, sizeof(real) * dim1 * dim2);     
   kernel_memcpy_real(d_cu, m.d_cu, dim1*dim2);
   
}

void kp_cu_matrix::operator=(const kp_matrix& m)
{
   resize(m.dim1, m.dim2);
   kp_cu_cudaMemcpy(d_cu, m.d, dim1*dim2*sizeof(real),cudaMemcpyHostToDevice);
   
}


void kp_cu_matrix::zeros()
{
   kp_cu_cudaMemset(d_cu, 0, sizeof(real) * dim1 * dim2);
   
}

void kp_cu_matrix::operator/=(real val)
{
  kernel_div_const(d_cu, val, dim1*dim2);

}
void kp_cu_matrix::operator*=(real val)
{
  kernel_mult_const_real(d_cu, val, dim1*dim2);

}
void kp_cu_matrix::operator+=(real val)
{
  kernel_add_const_real(d_cu, val, dim1*dim2);

}
void kp_cu_matrix::operator-=(real val)
{
  kernel_sub_const_real(d_cu, val, dim1*dim2);

}


void kp_cu_matrix::operator+= (const kp_cu_matrix& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_cu_matrix::operator+= | dimension problems"<<endl;
     }
   kernel_add_real(d_cu, M.d_cu, dim1*dim2);

}
void kp_cu_matrix::operator-= (const kp_cu_matrix& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_cu_matrix::operator-= | dimension problems"<<endl;
     }
   kernel_sub_real(d_cu, M.d_cu, dim1*dim2);

}












/*// calcule op(A)+op(B) (possibilite d'avoir A = *this)
void kp_cu_matrix::sum(char opA, char opB, kp_cu_matrix& A, const kp_cu_matrix& B)
{
	real alpha = 1;
	real beta = 1;
	cublasStatus_t stat;

	real* A_data;

	cublasOperation_t transa, transb;


	if (opA != 'N' && opA != 'T' && opB !='N' && opB != 'T')
	{
		cerr<<"error | kp_cu_matrix::sum | opA ou opB non reconnu"<<endl;
		exit(EXIT_FAILURE);
	}

	if (opA=='N') transa = CUBLAS_OP_N ; else transa = CUBLAS_OP_T ;
	if (opB=='N') transb = CUBLAS_OP_N ; else transb = CUBLAS_OP_T ;

	if (opA == opB)
	{
		if (dim1 != M.dim1 || dim2 != M.dim2)
		{
			cerr<<"error | kp_cu_matrix::sum | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if (dim1 != M.dim2 || dim2 != M.dim1)
		{
			cerr<<"error | kp_cu_matrix::sum | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}

	if (this == &A)
	{
		cudaMalloc((void** )&A_data, dim1*dim2*sizeof(A.d_cu[0]));
		cudaMemcpy(A_data, A.d_cu, sizeof(A.d_cu[0]) * A.dim1 * A.dim2, cudaMemcpyDeviceToDevice );			
	}
	else
	{
		A_data = A.d_cu;
	}


	if (opA=='N') resize(A.dim1, A.dim2);
	else resize(A.dim2, A.dim1);




	#ifdef KP_SINGLE
		stat = cublasSgeam(handle, transa, transb, dim1, dim2,
			&alpha, A_data, A.dim1,
			&beta, B.d_cu, B.dim1,
			d_cu, dim1);

	#else
		stat = cublasDgeam(handle, transa, transb, dim1, dim2,
			&alpha, A_data, A.dim1,
			&beta, B.d_cu, A.dim1,
			d_cu, dim1);

	#endif
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::sum | erreur lors de l'addition sur GPU" << endl;
		exit(EXIT_FAILURE);
	}

	if (this == &A) cudaFree(A_data);


}*/





// M[r1:r2-1 : c1:c2-1] avec  0 <= r1,r2 <= DIM1-1 ; 0 <= c1,c2 <= DIM2-1 ; r1<r2 ; c1<c2 
void kp_cu_matrix::init_from_matrix(const kp_matrix& M, int r1, int r2, int c1, int c2)
{
	if (r1 < 0 || r1 > r2 || r2 > M.dim1 || 
			c1 < 0 || c1 > c2 || c2 > M.dim2)
	{
		cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
		cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
		cerr<<"Error | kp_cu_matrix::init_from_matrix | index problems"<<endl;
		exit(EXIT_FAILURE);
	}
	resize(r2 - r1, c2 - c1);
	//real* d_tmp = new real[r2-r1, c2-c1];

	for (int j = 0 ; j < c2 - c1 ; j++)
		//d_tmp[j * dim1 + i] = M(i + r1, j + c1); 	
		kp_cu_cudaMemcpy(d_cu+(j*dim1), &(M(0, j+c1)), (r2-r1)*sizeof(real), cudaMemcpyHostToDevice);


	//delete[] d_tmp;

}

void kp_cu_matrix::init_from_matrix(const kp_cu_matrix& M, int r1, int r2, int c1, int c2)
{
	if (r1 < 0 || r1 > r2 || r2 > M.dim1 || 
			c1 < 0 || c1 > c2 || c2 > M.dim2)
	{
		cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
		cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
		cerr<<"Error | kp_cu_matrix::init_from_matrix | index problems"<<endl;
		exit(EXIT_FAILURE);
	}
	resize(r2 - r1, c2 - c1);
	//real* d_tmp = new real[r2-r1, c2-c1];

	kernel_set_submatrix(d_cu, M.d_cu, M.dim1, r1, c1, r2-r1, c2-c1);


	//delete[] d_tmp;

}







void kp_cu_matrix::init_from_transpose(const kp_matrix& M)
{
	resize(M.dim2, M.dim1);
	real* d_tmp = new real[dim1*dim2];
	for (int i = 0 ; i < dim1; i++)
		for (int j = 0 ; j < dim2 ; j++)
			d_tmp[j*dim1+i] = M(j,i);

	kp_cu_cudaMemcpy(d_cu, d_tmp, sizeof(real)*dim1*dim2, cudaMemcpyHostToDevice);

	delete[] d_tmp;

	
}

void kp_cu_matrix::init_from_transpose(cublasHandle_t handle, const kp_cu_matrix& M)
{
	resize(M.dim2, M.dim1);
	
	real alpha = 1;
	real beta = 0;


	cublasStatus_t stat;

	#ifdef KP_SINGLE
		stat = cublasSgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, M.dim2, M.dim1,
			&alpha, M.d_cu, M.dim1,
			&beta, M.d_cu, M.dim1,
			d_cu, M.dim2);

	#else
		stat = cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, M.dim2, M.dim1,
			&alpha, M.d_cu, M.dim1,
			&beta, M.d_cu, M.dim1,
			d_cu, M.dim2);

	#endif
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::init_from_transpose | erreur lors de la transposition sur GPU" << endl;
		throw KP_CUBLAS_GEAM;
		//exit(EXIT_FAILURE);
	}
}



void kp_cu_matrix::set_from_submatrix(const kp_matrix& subM, int r, int c)
{
	if (subM.dim1 + r > dim1 || subM.dim2 + c > dim2 || r < 0 || c < 0)
	{
		cerr<<"Error | kp_cu_matrix::set_from_submatrix | dimension problem"<<endl;
		cerr<<dim1<<"x"<<dim2<<" "<<subM.dim1<<"x"<<subM.dim2<<" "<<r<<" "<<c<<endl;
		exit(EXIT_FAILURE);
	}
	for (int j = 0 ; j < subM.dim2 ; j++)
		kp_cu_cudaMemcpy(d_cu+((j+c)*dim1+r), &(subM(0,j)), sizeof(real)*subM.dim1, cudaMemcpyHostToDevice);
	
   
}


void kp_cu_matrix::set_from_submatrix(const kp_cu_matrix& subM, int r, int c)
{
	if (subM.dim1 + r > dim1 || subM.dim2 + c > dim2 || r < 0 || c < 0)
	{
		cerr<<"Error | kp_cu_matrix::set_from_submatrix | dimension problem"<<endl;
		cerr<<dim1<<"x"<<dim2<<" "<<subM.dim1<<"x"<<subM.dim2<<" "<<r<<" "<<c<<endl;
		exit(EXIT_FAILURE);
	}
	for (int j = 0 ; j < subM.dim2 ; j++)
		//kp_cu_cudaMemcpy(d_cu+((j+c)*dim1+r), subM.d_cu+(j*subM.dim1), sizeof(real)*subM.dim1, cudaMemcpyDeviceToDevice);
		kernel_memcpy_real(d_cu+((j+c)*dim1+r), subM.d_cu+(j*subM.dim1), subM.dim1);
	
   
}


/*
void kp_cu_matrix::kp_cu_alloc()
{
	cudaError_t cudaStat;
	cudaStat = cudaMalloc ((void** )&d_cu, dim1*dim2*sizeof(d[0]));
	if (cudaStat != cudaSuccess) 
	{
		 cerr << "device memory allocation failed"<<endl;
		 exit(EXIT_FAILURE);
	}
}


void kp_cu_matrix::kp_cu_set_matrix()
{
	if ( sizeof(d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}
	cublasStatus_t stat;
	stat = cublasSetMatrix (dim1, dim2, sizeof(d[0]), d, dim1, d_cu, dim1);
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | erreur de copie de la matrice vers le GPU" << endl;
		exit(EXIT_FAILURE);
	}
}
void kp_cu_matrix::kp_cu_set_matrix(const int nb_row, const int nb_col)
{
	if ( sizeof(d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}

	cublasStatus_t stat;
	stat = cublasSetMatrix (nb_row, nb_col, sizeof(d[0]), d, dim1, d_cu, dim1);
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | erreur de copie de la matrice vers le GPU" << endl;
		exit(EXIT_FAILURE);
	}
}
// copie de la (sous-)matrice l'objet courant (*this) du CPU vers le GPU dans la (sous-)matrice A 
void kp_cu_matrix::kp_cu_set_matrix(const kp_cu_matrix& A, const int nb_row, const int nb_col)
{
	if ( sizeof(A.d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | A.d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}

	cublasStatus_t stat;
	stat = cublasSetMatrix (nb_row, nb_col, sizeof(d[0]), d, dim1, A.d_cu, A.dim1);
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_set_matrix | erreur de copie de la matrice vers le GPU" << endl;
		exit(EXIT_FAILURE);
	}
}




void kp_cu_matrix::kp_cu_get_matrix()
{
	if ( sizeof(d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}
	cublasStatus_t stat;
	stat = cublasGetMatrix (dim1, dim2, sizeof(d[0]), d_cu, dim1, d, dim1);
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | erreur de copie de la matrice vers le CPU" << endl;
		exit(EXIT_FAILURE);
	}
}
void kp_cu_matrix::kp_cu_get_matrix(const int nb_row, const int nb_col)
{
	if ( sizeof(d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}

	cublasStatus_t stat;
	stat = cublasGetMatrix (nb_row, nb_col, sizeof(d[0]), d_cu, dim1, d, dim1);
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | erreur de copie de la matrice vers le CPU" << endl;
		exit(EXIT_FAILURE);
	}
}
// copie de la (sous-)matrice de A du GPU vers le CPU dans la (sous-)matrice de l'objet courant (*this) 
void kp_cu_matrix::kp_cu_get_matrix(const kp_cu_matrix& A, const int nb_row, const int nb_col)
{
	if ( sizeof(A.d_cu[0]) != sizeof(d[0]) )
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | A.d_cu et d sont de type differents" << endl;
		exit(EXIT_FAILURE);
	}

	cublasStatus_t stat;
	stat = cublasGetMatrix (nb_row, nb_col, sizeof(d[0]), A.d_cu, A.dim1, d, dim1);
	if (stat != CUBLAS_STATUSi_SUCCESS) 
	{
		cerr << "error | kp_cu_matrix::kp_cu_get_matrix | erreur de copie de la matrice vers le CPU" << endl;
		exit(EXIT_FAILURE);
	}
}*/


void kp_cu_matrix::inverse()
{
	//kp_cu_timer temps_getrf,temps_getri; 
	if (dim1 != dim2)
	{
		cerr<<"Error | kp_cu_matrix::inverse | matrix is not square"<<endl;
		exit(EXIT_FAILURE);
	}

	int status;
	magma_int_t info;
	int* ipiv = new int[dim1];
//temps_getrf.start();
	#ifdef KP_SINGLE
		status = magma_sgetrf_gpu(dim1, dim2, d_cu, dim1, ipiv, &info);
	#else
		status = magma_dgetrf_gpu(dim1, dim2, d_cu, dim1, ipiv, &info);

	#endif
//temps_getrf.pause();
//cout << "temps getrf"<<temps_getrf.rez()<<endl;

		
	if (status != 0)
	{	
		cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_<type>getrf_gpu"<<endl;
		exit(EXIT_FAILURE);
	}


	#ifdef KP_SINGLE
		magma_int_t lwork = magma_get_sgetri_nb(dim1)*dim1;
	#else
		magma_int_t lwork = magma_get_dgetri_nb(dim1)*dim1;	
	#endif


	real* dwork;
	kp_cu_cudaMalloc((void**)&dwork, sizeof(real) * lwork);
//temps_getri.start();
	#ifdef KP_SINGLE
		status = magma_sgetri_gpu(dim1, d_cu, dim1, ipiv, dwork, lwork, &info);
		if (status != 0)
		{
			cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_sgetri_gpu"<<endl;
			exit(EXIT_FAILURE);
		}
	#else
		status = magma_dgetri_gpu(dim1, d_cu, dim1, ipiv, dwork, lwork, &info);
		if (status != 0)
		{
			cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_dgetri_gpu"<<endl;
			exit(EXIT_FAILURE);
		}
	#endif

	//temps_getri.pause();
//cout << "temps getri"<<temps_getri.rez()<<endl<<endl;


	delete[] ipiv;
	kp_cu_cudaFree(dwork);

}


void kp_cu_matrix::init_from_smatrix(const kp_cu_smatrix& M)
{
	resize(M.dim1, M.dim2);
	//zeros();
	kernel_sparse2full(d_cu, M.rowind_cu, M.colind_cu, M.values_cu, M.nnz, M.dim1, M.dim2);	
}


void kp_cu_matrix::_create(int dim1_, int dim2_)
{
	dim1 = dim1_;
	dim2 = dim2_;
	
	if(dim1*dim2>0)
	{
		try
		{
			kp_cu_cudaMalloc((void**)&d_cu, sizeof(real)*dim1*dim2);
		}
		catch (int err)
		{
			cerr<<"Erreur allocation CUDA"<<endl;
		}
	}
	
}

void kp_cu_matrix::_clear()  
{
     if(dim1*dim2>0) 
     {
   	if (d_cu == NULL)
     	{
		cerr<<"Error | kp_cu_matrix::_clear| d_cu == NULL ; dim1="<<dim1<<", dim2="<<dim2<<endl;
		exit(EXIT_FAILURE);
	}
	     kp_cu_cudaFree(d_cu);
     }
}







void kp_cu_geam(cublasHandle_t handle, char opA, char opB, real alpha, real beta, const kp_cu_matrix& A, const kp_cu_matrix& B, kp_cu_matrix& C )
{
	cublasStatus_t stat;
	cublasOperation_t transa, transb;

	real* A_data, *B_data;

	if( &C == &B || &C == &B)
	{
		cerr<<"error | kp_cu_geam | C doit etre different de A et B"<<endl;
		exit(EXIT_FAILURE);

	}


	if (opA != 'N' && opA != 'T' && opB !='N' && opB != 'T')
	{
		cerr<<"error | kp_cu_geam | opA ou opB non reconnu"<<endl;
		exit(EXIT_FAILURE);
	}

	if (opA=='N') transa = CUBLAS_OP_N ; else transa = CUBLAS_OP_T ;
	if (opB=='N') transb = CUBLAS_OP_N ; else transb = CUBLAS_OP_T ;

	if (opA == opB)
	{
		if (A.dim1 != B.dim1 || A.dim2 != B.dim2)
		{
			cerr<<"error | kp_cu_geam | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if (A.dim1 != B.dim2 || A.dim2 != B.dim1)
		{
			cerr<<"error | kp_cu_geam | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}
	

	if (&C == &A)
	{
		kp_cu_cudaMalloc((void** )&A_data, A.dim1*A.dim2*sizeof(A.d_cu[0]));
		//kp_cu_cudaMemcpy(A_data, A.d_cu, sizeof(A.d_cu[0]) * A.dim1 * A.dim2, cudaMemcpyDeviceToDevice);
		kernel_memcpy_real(A_data, A.d_cu, A.dim1 * A.dim2);
				
	}
	else
	{
		A_data = A.d_cu;
	}
	if (&C == &B)
	{		
		kp_cu_cudaMalloc((void** )&B_data, B.dim1*B.dim2*sizeof(B.d_cu[0]));
		//kp_cu_cudaMemcpy(B_data, B.d_cu, sizeof(B.d_cu[0]) * B.dim1 * B.dim2, cudaMemcpyDeviceToDevice);
		kernel_memcpy_real(B_data, B.d_cu, B.dim1 * B.dim2);

	}
	else
	{
		B_data = B.d_cu;
	}

	if (opA=='N') C.resize(A.dim1, A.dim2);
	else C.resize(A.dim2, A.dim1);




	#ifdef KP_SINGLE
		stat = cublasSgeam(handle, transa, transb, C.dim1, C.dim2,
			&alpha, A_data, A.dim1,
			&beta, B_data, B.dim1,
			C.d_cu, C.dim1);

	#else
		stat = cublasDgeam(handle, transa, transb, C.dim1, C.dim2,
			&alpha, A_data, A.dim1,
			&beta, B_data, A.dim1,
			C.d_cu, C.dim1);

	#endif
	if (stat != CUBLAS_STATUS_SUCCESS) 
	{
		cerr << "error | kp_cu_geam | erreur lors de l'addition sur GPU" << endl;
		throw KP_CUBLAS_GEAM;
		//exit(EXIT_FAILURE);
	}

	if (&C == &A)
	{
		kp_cu_cudaFree(A_data);
	}
	if (&C == &B)
	{
		kp_cu_cudaFree(B_data);
	}




}



void kp_cu_gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C)
{
	int opA_dim1, opA_dim2;	
	int opB_dim1, opB_dim2;

	cublasOperation_t transa,transb;
	cublasStatus_t stat1;

	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &transa);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2, &transb);
	

	if (C.dim1 != opA_dim1 || opA_dim2 != opB_dim1 || C.dim2 != opB_dim2)
	{
		cerr<<"Error | kp_cu_gemm | diminsion problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.dim1 == 0 || B.dim1 == 0 || C.dim1 == 0)
	{
		cerr<<"A.dim1="<<A.dim1<<" B.dim1="<<B.dim1<<" C.dim1="<<C.dim1<<endl;
		cerr<<"Error | kp_cu_gemm| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}


	#ifdef KP_SINGLE
	
	/*if (A.precision == 'D')
	{
	double *dAd, *dBd, *dCd;

	float* hAf = new float[A.dim1*A.dim2];
	float* hBf = new float[B.dim1*B.dim2];
	float* hCf = new float[C.dim1*C.dim2];

	double* hAd = new double[A.dim1*A.dim2];
	double* hBd = new double[B.dim1*B.dim2];
	double* hCd = new double[C.dim1*C.dim2];

	kp_cu_cudaMalloc((void**)&dAd, A.dim1*A.dim2*sizeof(double));
	kp_cu_cudaMalloc((void**)&dBd, B.dim1*B.dim2*sizeof(double));
	kp_cu_cudaMalloc((void**)&dCd, C.dim1*C.dim2*sizeof(double));

	kp_cu_cudaMemcpy(hAf, A.d_cu, A.dim1*A.dim2*sizeof(float), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(hBf, B.d_cu, B.dim1*B.dim2*sizeof(float), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(hCf, C.d_cu, C.dim1*C.dim2*sizeof(float), cudaMemcpyDeviceToHost);

	for (int i=0;i<A.dim1*A.dim2;i++) hAd[i] = (double) hAf[i];
	for (int i=0;i<B.dim1*B.dim2;i++) hBd[i] = (double) hBf[i];
	for (int i=0;i<C.dim1*C.dim2;i++) hCd[i] = (double) hCf[i];

	kp_cu_cudaMemcpy(dAd, hAd, A.dim1*A.dim2*sizeof(double), cudaMemcpyHostToDevice);
	kp_cu_cudaMemcpy(dBd, hBd, B.dim1*B.dim2*sizeof(double), cudaMemcpyHostToDevice);
	kp_cu_cudaMemcpy(dCd, hCd, C.dim1*C.dim2*sizeof(double), cudaMemcpyHostToDevice);


	double alpha_d = (double) alpha;
 	double beta_d  = (double) beta;




		stat1= cublasDgemm(handle, transa, transb,\
			opA_dim1, opB_dim2, opA_dim2, &alpha_d, dAd, A.dim1, dBd, B.dim1, &beta_d, dCd, C.dim1);

	cudaDeviceSynchronize();
	cudaMemcpy(hCd, dCd, C.dim1*C.dim2*sizeof(double), cudaMemcpyDeviceToHost);
	
	for (int i=0;i<C.dim1*C.dim2;i++) hCf[i] = (float) hCd[i];

	cudaMemcpy(C.d_cu, hCf, C.dim1*C.dim2*sizeof(float), cudaMemcpyHostToDevice);

	delete[] hAf ; delete[] hBf ; delete[] hCf ;
	delete[] hAd ; delete[] hBd ; delete[] hCd ;
	kp_cu_cudaFree(dAd) ; kp_cu_cudaFree(dBd) ; kp_cu_cudaFree(dCd);
	}

	else*/
	{
	stat1= cublasSgemm(handle, transa, transb,\
		       	opA_dim1, opB_dim2, opA_dim2, &alpha, A.d_cu, A.dim1, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
	}



	#else
		stat1= cublasDgemm(handle, transa, transb,\
			opA_dim1, opB_dim2, opA_dim2, &alpha, A.d_cu, A.dim1, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
	#endif

	if (stat1 != CUBLAS_STATUS_SUCCESS)
	{
		cerr << "error | kp_cu_matrix | kp_cu_gemm failed : "<< stat1 << endl;
		throw KP_CUBLAS_GEMM;
	}
	

	
}

void kp_cu_gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix& A, const kp_cu_vector& x, real beta, kp_cu_vector& y)
{
	cublasStatus_t cublasStat;
	cublasOperation_t transa;

	int opA_dim1, opA_dim2;
	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &transa);

	if (opA_dim2 != x.size() || opA_dim1 != y.size())
	{
		cerr<<"Error | kp_cu_gemv | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.dim1 == 0 || x.d_cu == 0 || y.d_cu == 0)
	{
		cerr<<"A.dim1="<<A.dim1<<"x.d_cu="<<x.d_cu<<" y.d_cu="<<y.d_cu<<endl;
		cerr<<"Error | kp_cu_gemv| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}

	int stride = 1;
	#ifdef KP_SINGLE
		cublasStat = cublasSgemv(handle, transa, A.dim1, A.dim2, &alpha, A.d_cu, A.dim1, x.d_cu, stride, &beta, y.d_cu, stride);

	#else
		cublasStat = cublasDgemv(handle, transa, A.dim1, A.dim2, &alpha, A.d_cu, A.dim1, x.d_cu, stride, &beta, y.d_cu, stride);
	#endif

	if (cublasStat != CUBLAS_STATUS_SUCCESS)
	{
		cerr << "error | kp_cu_matrix | kp_cu_gemv failed !" << endl;
		throw KP_CUBLAS_GEMV;
	}

   
}





void kp_cu_check_op_set_dim(int op, const kp_cu_matrix&M, int& dim1, int& dim2, cublasOperation_t* trans)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
	*trans = CUBLAS_OP_N;
     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
	*trans = CUBLAS_OP_T;
     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim | op should ge either N either T"<<endl;
	exit(EXIT_FAILURE);
   
     }
}
void kp_cu_vertcat(kp_cu_matrix& rez, const vector<kp_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim1 += ms[i]->dim1;
	if (ms[i]->dim2 != ms[0]->dim2)
	  {
	     cerr<<"Error | kp_cu_vertcat | inconsistant second dimension: "<<ms[i]->dim2<<" "<<ms[0]->dim2<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(sum_dim1, ms[0]->dim2);
   int dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), dim1, 0);
	dim1 += ms[i]->dim1;
     }
}

void kp_cu_vertcat(kp_cu_matrix& rez, const vector<kp_cu_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim1 += ms[i]->dim1;
	if (ms[i]->dim2 != ms[0]->dim2)
	  {
	     cerr<<"Error | kp_cu_vertcat | inconsistant second dimension: "<<ms[i]->dim2<<" "<<ms[0]->dim2<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(sum_dim1, ms[0]->dim2);
   int dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), dim1, 0);
	dim1 += ms[i]->dim1;
     }
}

void kp_cu_horizcat(kp_cu_matrix& rez, const vector<kp_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim2 += ms[i]->dim2;
	if (ms[i]->dim1 != ms[0]->dim1)
	  {
	     cerr<<"Error | kp_cu_horizcat | inconsistant first dimension: "<<ms[i]->dim1<<" "<<ms[0]->dim1<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(ms[0]->dim1, sum_dim2);
   int dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), 0, dim2);
	dim2 += ms[i]->dim2;
     }
}

void kp_cu_horizcat(kp_cu_matrix& rez, const vector<kp_cu_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim2 += ms[i]->dim2;
	if (ms[i]->dim1 != ms[0]->dim1)
	  {
	     cerr<<"Error | kp_cu_horizcat | inconsistant first dimension: "<<ms[i]->dim1<<" "<<ms[0]->dim1<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(ms[0]->dim1, sum_dim2);
   int dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), 0, dim2);
	dim2 += ms[i]->dim2;
     }
}





