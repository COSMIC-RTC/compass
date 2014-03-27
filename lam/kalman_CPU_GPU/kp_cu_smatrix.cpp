//kp_cu_smatrix.cpp

#include "kp_cu_smatrix.h"
#include <mkl_spblas.h>
#include <iostream>
#include <algorithm>
#include <string.h>
#include "kp_cu2kp.h"
#include <fstream>
using namespace std;

#ifdef __WITH_FLOPS_CALCULATOR__
#include "kp_flopscalc.h"
#endif


kp_cu_smatrix::kp_cu_smatrix()
{
	_create(0,0,0);
}

kp_cu_smatrix::kp_cu_smatrix(const kp_cu_smatrix& M)
{
	_create(0,0,0);
	resize(M.nnz, M.dim1, M.dim2);

	kernel_memcpy_real(values_cu, M.values_cu, nnz);
   	kernel_memcpy_int(rowind_cu, M.rowind_cu, nnz);   
   	kernel_memcpy_int(colind_cu, M.colind_cu, nnz);
	

	majorDim = M.majorDim;

	cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
	cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
	cusparseSetMatIndexBase(descr,cusparseGetMatIndexBase(M.descr));
	cusparseSetMatType(descr, cusparseGetMatType(M.descr));




#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)nnz * 3);
#endif

}

void kp_cu_smatrix::operator=(const kp_cu_smatrix& M)
{
	resize(M.nnz, M.dim1, M.dim2);
	kernel_memcpy_real(values_cu, M.values_cu, nnz);
   	kernel_memcpy_int(rowind_cu, M.rowind_cu, nnz);   
   	kernel_memcpy_int(colind_cu, M.colind_cu, nnz);
	

	majorDim = M.majorDim;

	cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
	cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
	cusparseSetMatIndexBase(descr,cusparseGetMatIndexBase(M.descr));
	cusparseSetMatType(descr, cusparseGetMatType(M.descr));
#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)nnz * 3);
#endif
}

void kp_cu_smatrix::operator=(const kp_smatrix& M)
{
	resize(M.nnz, M.dim1, M.dim2);
	kp_cu_cudaMemcpy(values_cu, M.values, nnz*sizeof(real), cudaMemcpyHostToDevice);
   	kp_cu_cudaMemcpy(rowind_cu, M.rowind, nnz*sizeof(int), cudaMemcpyHostToDevice);   
   	kp_cu_cudaMemcpy(colind_cu, M.colind, nnz*sizeof(int), cudaMemcpyHostToDevice);

	majorDim = M.get_majorDim();

#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)nnz * 3);
#endif

}


kp_cu_smatrix::kp_cu_smatrix(const kp_smatrix& M)
{
	_create(0,0,0);;
	resize(M.nnz, M.dim1, M.dim2);

	kp_cu_cudaMemcpy(values_cu, M.values, nnz*sizeof(real), cudaMemcpyHostToDevice);
   	kp_cu_cudaMemcpy(rowind_cu, M.rowind, nnz*sizeof(int), cudaMemcpyHostToDevice);   
   	kp_cu_cudaMemcpy(colind_cu, M.colind, nnz*sizeof(int), cudaMemcpyHostToDevice);

	majorDim = M.get_majorDim();

#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)nnz * 3);
#endif

}

void kp_cu_smatrix::resize(int nnz_, int dim1_, int dim2_)
{
   if (nnz != nnz_)
     {
	_clear();
	_create(nnz_, dim1_, dim2_);
     }
   else
     {
	dim1 = dim1_;
	dim2 = dim2_;
	majorDim = 'U';
     }
   isCSRconverted = false;
   isCSRconvertedT = false;

}

void kp_cu_smatrix::_create(int nnz_, int dim1_, int dim2_)
{
   cusparseStatus_t status;
   nnz  = nnz_;
   dim1 = dim1_;
   dim2 = dim2_;


   if(nnz>0)
   {
   	kp_cu_cudaMalloc((void**)&values_cu, nnz*sizeof(real));
	kp_cu_cudaMalloc((void**)&rowind_cu, nnz*sizeof(int));
	kp_cu_cudaMalloc((void**)&colind_cu, nnz*sizeof(int));
   	kp_cu_cudaMalloc((void**)&csrRowPtr, (dim1+1)*sizeof(int));
   	kp_cu_cudaMalloc((void**)&csrRowPtrT, (dim2+1)*sizeof(int));	
   }



   majorDim = 'U' ;
   isCSRconverted=false;
   isCSRconvertedT=false;


   status = cusparseCreateMatDescr(&descr); 
   if (status != CUSPARSE_STATUS_SUCCESS) 
	{
		cerr<<"Error | kp_cu_smatrix::_create | Matrix descriptor initialization failed"<<endl;
		throw KP_CUSPARSE_MAT_DESCR1;
		//exit(EXIT_FAILURE);
		
	}
   status = cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);

   if (status != CUSPARSE_STATUS_SUCCESS) 
	{
		cerr<<"Error | kp_cu_smatrix::_create | set IndexBase failed"<<endl;
		throw KP_CUSPARSE_SET_MAT_BASE;
		//exit(EXIT_FAILURE);
		
	}


}

void kp_cu_smatrix::_clear()
{
   cusparseStatus_t status;



   if (nnz>0)
   {
	if (csrRowPtr == NULL || csrRowPtrT==NULL || values_cu==NULL || rowind_cu==NULL || colind_cu==NULL )
	{
		cerr<<"Error | kp_cu_smatrix::_clear | double clear"<<endl;
		exit(EXIT_FAILURE);
	}
   	kp_cu_cudaFree(values_cu);
   	kp_cu_cudaFree(rowind_cu);
   	kp_cu_cudaFree(colind_cu);
   	kp_cu_cudaFree(csrRowPtr);
   	kp_cu_cudaFree(csrRowPtrT);
   }
   values_cu   = NULL;
   rowind_cu   = NULL;
   colind_cu   = NULL;
   csrRowPtr = NULL;
   csrRowPtrT= NULL;
   nnz      = 0;
   dim1     = 0;
   dim2     = 0;
   majorDim =  'U' ;
   isCSRconverted = false;
   isCSRconvertedT = false;

   status = cusparseDestroyMatDescr(descr);

   descr = 0;

    if (status != CUSPARSE_STATUS_SUCCESS) {
        cerr<<"Error | kp_cu_smatrix::_clear | Matrix descriptor destruction failed"<<endl;
	throw KP_CUSPARSE_MAT_DESCR2;
	//exit(EXIT_FAILURE);
    }   
}



void kp_cu_smatrix::init_from_transpose(const kp_cu_smatrix& M)
{
   resize(M.nnz, M.dim1, M.dim2);

   kernel_memcpy_real(values_cu, M.values_cu, sizeof(real) * nnz);
   kernel_memcpy_int(colind_cu, M.rowind_cu, sizeof(int) * nnz);
   kernel_memcpy_int(rowind_cu, M.colind_cu, sizeof(int) * nnz);
   

   if (M.majorDim == 'C') majorDim = 'R';
  else if (M.majorDim == 'R') majorDim = 'C';
  else majorDim = 'U';
     
   
#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)nnz * 3);
#endif
}

bool kp_cu_smatrix::isColumnMajor()
{
    bool colMajor = true;
    int colm1=0;
    int rowm1=0;

    kp_smatrix A_tmp;
    kp_cu2kp_smatrix(A_tmp,*this);


    for (int i=0;i<nnz;i++)
    {

        if (A_tmp.colind[i]==colm1 )
        {
            colMajor = colMajor &&(A_tmp.rowind[i]>rowm1);
            rowm1=A_tmp.rowind[i];
        }
        else
        {
            rowm1=A_tmp.rowind[i];
            colMajor = colMajor && (A_tmp.colind[i]>colm1);
            colm1=A_tmp.colind[i];
        }
    }
    return colMajor;
}



kp_cu_smatrix::~kp_cu_smatrix()
{
   _clear();
}




// y = alpha * op_A(A) * x + beta * y
void kp_cu_gemv(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_vector& x, real beta, kp_cu_vector& y)
{
	int opA_dim1, opA_dim2;
	cusparseOperation_t trans;
	cusparseStatus_t status;
       
	kp_cu_timer t1=kp_cu_timer();
	kp_cu_timer t2=kp_cu_timer();
	kp_cu_timer t3=kp_cu_timer(); 
	kp_cu_timer t4=kp_cu_timer(); 
	kp_cu_timer t5=kp_cu_timer(); 


	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &trans);

	if (opA_dim1 != y.size() || opA_dim2 != x.size() )
	{
		cerr<<"Error | kp_cu_gemv (sparse) | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}

	if (A.majorDim == 'R') //si row-major
	{	

		if (!A.isCSRconverted)
		{


			status= cusparseXcoo2csr(handle, A.rowind_cu, A.nnz, A.dim1, A.csrRowPtr, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemv (sparse) | Conversion from COO to CSR format failed"<<endl;
				throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR"<<endl;
			A.isCSRconverted = true;
		}
	

		//t5.start();
		#ifdef KP_SINGLE
		status= cusparseScsrmv(handle, trans, A.dim1, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, x.d_cu, &beta, y.d_cu);		
		#else
		status= cusparseDcsrmv(handle, trans, A.dim1, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, x.d_cu, &beta, y.d_cu);
		#endif
		//t5.pause();
		//cout << "temps de multiplication R = " << t5.rez()<<endl;

	}


	else if (A.majorDim == 'C') //si column-major
	{
		//t1.start();

		//t2.start();
		if (!A.isCSRconvertedT)
		{
			status= cusparseXcoo2csr(handle, A.colind_cu, A.nnz, A.dim2, A.csrRowPtrT, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemv (sparse) | Conversion from COO to CSR format failed"<<endl;
				throw KP_CUSPARSE_COO2CSR;				
				//exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR T"<<endl;
			A.isCSRconvertedT = true;

		}
		//t2.pause();
		//cout << "temps de conversion = " << t2.rez()<<endl;

		if(op_A == 'N')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A
			// En mettant transA=CUSPARSE_OPERATION_TRANSPOSE en parametre de cusparseDcsrmm,
			// on effectue en fait, C = alpha * ((A)T)T * B  + beta * C

			
			//t3.start();
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmv(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.d_cu, &beta, y.d_cu);
			#else
			status= cusparseDcsrmv(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.d_cu, &beta, y.d_cu);
			#endif
			//t3.pause();
			//cout << "temps multiplication CN = "<<t3.rez()<<endl;
		}

		else if (op_A == 'T')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A


			//t4.start();			
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.d_cu, &beta, y.d_cu);

			#else
			status= cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.d_cu, &beta, y.d_cu);
			#endif
			//t4.pause();
			//cout << "temps multiplication CT = "<<t4.rez()<<endl;

		}
		else
		{
			cerr<<"Error | kp_cu_gemv (sparse) | op_A transpose character is unknown."<<endl;
			exit(EXIT_FAILURE);

		}
		//t1.pause();
		//cout << "temps total = "<< t1.rez()<<endl;
	}

	else //si indeterminate
	{
		cerr<<"Error | kp_cu_gemv (sparse) | Sparse matrix has been defined neither as a row-major nor as a column-major"<<endl;
		exit(EXIT_FAILURE);
	}
	
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{
        	cerr<<"Error | kp_cu_gemv (sparse) | Matrix-vector multiplication failed"<<endl;
		throw KP_CUSPARSE_GEMV;
		//exit(EXIT_FAILURE);
	}

}

// alpha * op_A(A) * B + beta * C 
void kp_cu_gemm(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C)
{
	ofstream fichier;

	int opA_dim1, opA_dim2;
	cusparseOperation_t trans;
	cusparseStatus_t status;

	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &trans);
			
	if (C.dim1 != opA_dim1 || opA_dim2 != B.dim1 || C.dim2 != B.dim2)
	{
		cerr<<"Error | kp_cu_gemm (sparse) | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}

	if (A.majorDim == 'R') //si row-major
	{				
		if (!A.isCSRconverted)
		{
			status= cusparseXcoo2csr(handle, A.rowind_cu, A.nnz, A.dim1, A.csrRowPtr, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR"<<endl;
			A.isCSRconverted = true;
		}
			

	
		#ifdef KP_SINGLE


	/*double *dAd, *dBd, *dCd;

	float* hAf = new float[A.nnz];
	float* hBf = new float[B.dim1*B.dim2];
	float* hCf = new float[C.dim1*C.dim2];

	double* hAd = new double[A.nnz];
	double* hBd = new double[B.dim1*B.dim2];
	double* hCd = new double[C.dim1*C.dim2];

	kp_cu_cudaMalloc((void**)&dAd, A.nnz*sizeof(double));
	kp_cu_cudaMalloc((void**)&dBd, B.dim1*B.dim2*sizeof(double));
	kp_cu_cudaMalloc((void**)&dCd, C.dim1*C.dim2*sizeof(double));

	kp_cu_cudaMemcpy(hAf, A.values_cu, A.nnz*sizeof(float), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(hBf, B.d_cu, B.dim1*B.dim2*sizeof(float), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(hCf, C.d_cu, C.dim1*C.dim2*sizeof(float), cudaMemcpyDeviceToHost);
	for (int i=0;i<A.nnz;i++) hAd[i] = (double) hAf[i];
	for (int i=0;i<B.dim1*B.dim2;i++) hBd[i] = (double) hBf[i];
	for (int i=0;i<C.dim1*C.dim2;i++) hCd[i] = (double) hCf[i];

	kp_cu_cudaMemcpy(dAd, hAd, A.nnz*sizeof(double), cudaMemcpyHostToDevice);
	kp_cu_cudaMemcpy(dBd, hBd, B.dim1*B.dim2*sizeof(double), cudaMemcpyHostToDevice);
	kp_cu_cudaMemcpy(dCd, hCd, C.dim1*C.dim2*sizeof(double), cudaMemcpyHostToDevice);


	double alpha_d = (double) alpha;
 	double beta_d  = (double) beta;




	status= cusparseDcsrmm(handle, trans, A.dim1, B.dim2, A.dim2, A.nnz, &alpha_d, A.descr, dAd,\
			       	A.csrRowPtr, A.colind_cu, dBd, B.dim1, &beta_d, dCd, C.dim1);	
	cudaDeviceSynchronize();
	cudaMemcpy(hCd, dCd, C.dim1*C.dim2*sizeof(double), cudaMemcpyDeviceToHost);
	
	for (int i=0;i<C.dim1*C.dim2;i++) hCf[i] = (float) hCd[i];

	cudaMemcpy(C.d_cu, hCf, C.dim1*C.dim2*sizeof(float), cudaMemcpyHostToDevice);

	delete[] hAf ; delete[] hBf ; delete[] hCf ;
	delete[] hAd ; delete[] hBd ; delete[] hCd ;
	kp_cu_cudaFree(dAd) ; kp_cu_cudaFree(dBd) ; kp_cu_cudaFree(dCd);*/







		
		status= cusparseScsrmm(handle, trans, A.dim1, B.dim2, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);		
		#else
		status= cusparseDcsrmm(handle, trans, A.dim1, B.dim2, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
		#endif



	}

	else if (A.majorDim == 'C') //si column-major
	{
		
		if (!A.isCSRconvertedT)
		{
			status= cusparseXcoo2csr(handle, A.colind_cu, A.nnz, A.dim2, A.csrRowPtrT, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR T"<<endl;
			A.isCSRconvertedT = true;

		}

		

		if(op_A == 'N')
		{


			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A
			// En mettant transA=CUSPARSE_OPERATION_TRANSPOSE en parametre de cusparseDcsrmm,
			// on effectue en fait, C = alpha * ((A)T)T * B  + beta * C
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmm(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, B.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
			#else
			status= cusparseDcsrmm(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, B.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
			#endif


		}

		else if (op_A == 'T')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A

	
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, B.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);

			#else
			status= cusparseDcsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, B.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.d_cu, B.dim1, &beta, C.d_cu, C.dim1);
			#endif
		}
		else
		{
			cerr<<"Error | kp_cu_gemm (sparse) | op_A transpose character is unknown."<<endl;
			exit(EXIT_FAILURE);

		}
		


	}
	else //si indeterminate
	{
		cerr<<"Error | kp_cu_gemm (sparse) | Sparse matrix has been defined neither as a row-major nor as a column-major"<<endl;
		exit(EXIT_FAILURE);
	}

	

	
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{
        	cerr<<"Error | kp_cu_gemm (sparse) | Matrix-matrix multiplication failed"<<endl;
		throw KP_CUSPARSE_GEMM;
		//exit(EXIT_FAILURE);
	}

}


/*void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublashandle,  char op_A, char op_B, real alpha, kp_cu_smatrix& A, kp_cu_matrix& B, real beta, kp_cu_matrix& C)
{

	int opA_dim1, opA_dim2;
	int opB_dim1, opB_dim2;

	cusparseOperation_t transA;
	cusparseOperation_t transB;
	cusparseOperation_t transB2;

	cusparseStatus_t status;

	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &transA);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2, &transB);

	int ldb=0;
	kp_cu_matrix B2;
	B2.init_from_transpose(cublashandle,B);
	
	if (op_B=='T') 
		transB2 = CUSPARSE_OPERATION_NON_TRANSPOSE;
	else if (op_B=='N')
		transB2 = CUSPARSE_OPERATION_TRANSPOSE;


			
	if (C.dim1 != opA_dim1 || opA_dim2 != opB_dim1 || C.dim2 != opB_dim2)
	{
		cerr<<"Error | kp_cu_gemm (sparse) | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}

	if (A.majorDim == 'R') //si row-major
	{			
		if (!A.isCSRconverted)
		{
			status= cusparseXcoo2csr(handle, A.rowind_cu, A.nnz, A.dim1, A.csrRowPtr, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
					exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR"<<endl;
			A.isCSRconverted = true;
		}
		



	
		#ifdef KP_SINGLE
		status= cusparseScsrmm2(handle, transA, transB2, A.dim1, C.dim2, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, B2.d_cu, B2.dim1, &beta, C.d_cu, C.dim1);		
		#else
		status= cusparseDcsrmm2(handle, transA, transB2, A.dim1, C.dim2, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, B2.d_cu, B2.dim1, &beta, C.d_cu, C.dim1);
		#endif



	}

	else if (A.majorDim == 'C') //si column-major
	{
		if (!A.isCSRconvertedT)
		{
			status= cusparseXcoo2csr(handle, A.colind_cu, A.nnz, A.dim2, A.csrRowPtrT, CUSPARSE_INDEX_BASE_ONE); 
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				exit(EXIT_FAILURE);
			}
			//cout << "conversion CSR T"<<endl;
			A.isCSRconvertedT = true;
		}

		if(op_A == 'N')
		{


			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A
			// En mettant transA=CUSPARSE_OPERATION_TRANSPOSE en parametre de cusparseDcsrmm,
			// on effectue en fait, C = alpha * ((A)T)T * B  + beta * C
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmm2(handle, CUSPARSE_OPERATION_TRANSPOSE, transB2, A.dim2, C.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B2.d_cu, B2.dim1, &beta, C.d_cu, C.dim1);
			#else
			status= cusparseDcsrmm2(handle, CUSPARSE_OPERATION_TRANSPOSE, transB2, A.dim2, C.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B2.d_cu, B2.dim1, &beta, C.d_cu, C.dim1);
			#endif

		}

		else if (op_A == 'T')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A

	
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			#ifdef KP_SINGLE
			status= cusparseScsrmm2(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, transB2, A.dim2, C.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B2.d_cu, B2.dim, &beta, C.d_cu, C.dim1);

			#else
			status= cusparseDcsrmm2(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, transB2, A.dim2, C.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B2.d_cu, B2.dim1, &beta, C.d_cu, C.dim1);
			#endif
		}
		else
		{
			cerr<<"Error | kp_cu_gemm (sparse) | op_A transpose character is unknown."<<endl;
			exit(EXIT_FAILURE);

		}
		


	}
	else //si indeterminate
	{
		cerr<<"Error | kp_cu_gemm (sparse) | Sparse matrix has been defined neither as a row-major nor as a column-major"<<endl;
		exit(EXIT_FAILURE);
	}

	

	
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{
		if(status == CUSPARSE_STATUS_NOT_INITIALIZED) cout<<"err1"<<endl;
		if(status == CUSPARSE_STATUS_ALLOC_FAILED) cout<<"err2"<<endl;
		if(status == CUSPARSE_STATUS_INVALID_VALUE) cout<<"err3"<<endl;
		if(status == CUSPARSE_STATUS_ARCH_MISMATCH) cout<<"err4"<<endl;
		if(status == CUSPARSE_STATUS_EXECUTION_FAILED) cout<<"err5"<<endl;
		if(status == CUSPARSE_STATUS_INTERNAL_ERROR) cout<<"err6"<<endl;
		if(status == CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED) cout<<"err7"<<endl;



        	cerr<<"Error | kp_cu_gemm (sparse) | Matrix-matrix multiplication failed"<<" "<<status<<endl;
		exit(EXIT_FAILURE);
	}
	cout<<"mult OK"<<endl;
}*/


void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix&M, int& dim1, int& dim2, cusparseOperation_t* trans)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
	*trans = CUSPARSE_OPERATION_NON_TRANSPOSE ;
     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
	*trans = CUSPARSE_OPERATION_TRANSPOSE ;
     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim (sparse) | op should be either N or T"<<endl;
	exit(EXIT_FAILURE);
     }
}
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix&M, int& dim1, int& dim2, cusparseOperation_t* trans)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
	*trans = CUSPARSE_OPERATION_NON_TRANSPOSE ;
     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
	*trans = CUSPARSE_OPERATION_TRANSPOSE ;
     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim | op should be either N or T"<<endl;
	exit(EXIT_FAILURE);
     }
}








