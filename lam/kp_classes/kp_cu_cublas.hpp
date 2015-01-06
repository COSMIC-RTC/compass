//kp_cu_cublas.cpp


#ifndef __SEGER__KP_CU_CUBLAS_HPP__
#define __SEGER__KP_CU_CUBLAS_HPP__

//FULL MATRIX (CUBLAS)



template <typename real>
void kp_cu_geam(cublasHandle_t handle, char opA, char opB, real alpha, real beta, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, kp_cu_matrix<real>& C )
{
	cublasOperation_t transa, transb;


	if (opA != 'N' && opA != 'T' && opB !='N' && opB != 'T')
	{
		cerr<<"error | kp_cu_geam | opA ou opB non reconnu"<<endl;
		exit(EXIT_FAILURE);
	}

	if (opA=='N') transa = CUBLAS_OP_N ; else transa = CUBLAS_OP_T ;
	if (opB=='N') transb = CUBLAS_OP_N ; else transb = CUBLAS_OP_T ;

	if (opA == opB)
	{
		if (A.getDim1() != B.getDim1() || A.getDim2() != B.getDim2())
		{
			cerr<<"error | kp_cu_geam | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if (A.getDim1() != B.getDim2() || A.getDim2() != B.getDim1())
		{
			cerr<<"error | kp_cu_geam | dimension problems"<<endl;
			exit(EXIT_FAILURE);
		}
	}
	

	if (opA=='N') C.resize(A.getDim1(), A.getDim2());
	else C.resize(A.getDim2(), A.getDim1());


	kp_cu_geam_core(handle, transa, transb, C.getDim1(), C.getDim2(),
			&alpha, A.getData(), A.getDim1(),
			&beta, B.getData(), B.getDim1(),
			C.getData(), C.getDim1());
}

template <typename real>
void kp_cu_gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, real beta, kp_cu_matrix<real>& C)
{

	int opA_dim1, opA_dim2;	
	int opB_dim1, opB_dim2;


	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2);
	
	if (C.getDim1() != opA_dim1 || opA_dim2 != opB_dim1 || C.getDim2() != opB_dim2)
	{
		cerr<<"Error | kp_cu_gemm | diminsion problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.getDim1() < 1 || B.getDim1() < 1 || C.getDim1() < 1)
	{
		cerr<<"A.getDim1()="<<A.getDim1()<<" B.getDim1()="<<B.getDim1()<<" C.getDim1()="<<C.getDim1()<<endl;
		cerr<<"Error | kp_cu_gemm| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}

	cublasOperation_t transa,transb;

	kp_cu_op2cublastrans(op_A, &transa);
	kp_cu_op2cublastrans(op_B, &transb);


	kp_cu_gemm_core(handle, transa, transb,
		       	opA_dim1, opB_dim2, opA_dim2, 
			&alpha, A.getData(), A.getDim1(), 
			B.getData(), B.getDim1(), &beta,
		       	C.getData(), C.getDim1());

}


template <typename real>
void kp_cu_gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x, real beta, kp_cu_vector<real>& y)
{

	int opA_dim1, opA_dim2;
	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);

	if (opA_dim2 != x.size() || opA_dim1 != y.size())
	{
		cerr<<"Error | kp_cu_gemv | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.getDim1() < 1 || x.size() < 1 || y.size() < 1)
	{
		cerr<<"A.getDim1()="<<A.getDim1()<<"x.size="<<x.size()<<" y.size()="<<y.size()<<endl;
		cerr<<"Error | kp_cu_gemv| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}

	cublasOperation_t transa;
	kp_cu_op2cublastrans(op_A, &transa);

	int stride = 1;

	kp_cu_gemv_core(handle, transa, A.getDim1(), A.getDim2(), &alpha, A.getData(), A.getDim1(), x.getData(), stride, &beta, y.getData(), stride);
}


template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix<real>&M, int& dim1, int& dim2)
{
   if (op == 'N')
     {
	dim1 = M.getDim1();
	dim2 = M.getDim2();
     }
   else if (op == 'T')
     {
	dim1 = M.getDim2();
	dim2 = M.getDim1();
     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim | op should ge either N or T"<<endl;
	exit(EXIT_FAILURE);
   
     }
}







// SPRASE MATRIX (CUSPARSE)



// y = alpha * op_A(A) * x + beta * y
template <typename real>
void kp_cu_gemv(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix<real>& A, const kp_cu_vector<real>& x, real beta, kp_cu_vector<real>& y)
{
	int opA_dim1, opA_dim2;
	cusparseOperation_t trans;
       
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
                A.convert2csr(handle);

		kp_cu_gemv_core(handle, trans, A.dim1, A.dim2, A.nnz, &alpha, A.descr, A.values_cu,\
			       	A.csrRowPtr, A.colind_cu, x.getData(), &beta, y.getData());

	}


	else if (A.majorDim == 'C') //si column-major
	{
		A.convert2csrT(handle);

		if(op_A == 'N')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A
			// En mettant transA=CUSPARSE_OPERATION_TRANSPOSE en parametre de cusparseDcsrmm,
			// on effectue en fait, C = alpha * ((A)T)T * B  + beta * C

			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			kp_cu_gemv_core(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.getData(), &beta, y.getData());
		}

		else if (op_A == 'T')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A

			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			
			kp_cu_gemv_core(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, A.dim1, A.nnz, &alpha,\
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, x.getData(), &beta, y.getData());

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

}

// alpha * op_A(A) * B + beta * C 
template <typename real>
void kp_cu_gemm(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix<real>& A, const kp_cu_matrix<real>& B, real beta, kp_cu_matrix<real>& C)
{

	int opA_dim1, opA_dim2;
	cusparseOperation_t trans;

	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &trans);
			
	if (C.getDim1() != opA_dim1 || opA_dim2 != B.getDim1() || C.getDim2() != B.getDim2())
	{
		cerr<<"Error | kp_cu_gemm (sparse) | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}

	if (A.majorDim == 'R') //si row-major
	{				
		A.convert2csr(handle);
		kp_cu_gemm_core(handle, trans, A.dim1, B.getDim2(), A.dim2, A.nnz, &alpha, A.descr, A.values_cu,
			       	A.csrRowPtr, A.colind_cu, B.getData(), B.getDim1(), &beta, C.getData(), C.getDim1());
	}

	else if (A.majorDim == 'C') //si column-major
	{
	        A.convert2csrT(handle);	

		if(op_A == 'N')
		{

			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A
			// En mettant transA=CUSPARSE_OPERATION_TRANSPOSE en parametre de cusparseDcsrmm,
			// on effectue en fait, C = alpha * ((A)T)T * B  + beta * C
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE

			kp_cu_gemm_core(handle, CUSPARSE_OPERATION_TRANSPOSE, A.dim2, B.getDim2(), A.dim1, A.nnz, &alpha,
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.getData(), B.getDim1(), &beta, C.getData(), C.getDim1());
		}

		else if (op_A == 'T')
		{
			// La matrice creuse en entree de cusparseDcsrmm doit etre row-major
			// Pour cela on met en entree de cusparseDcsrmm ce qui equivaut a (A)T
			// Comme A est col-major, en inversant colind_cu et row_ind_cu on obtient
			// une matrice row-major, qui est la transposee de A

	
			// les parametres en entree correspondent a (A)T avec transA=CUSPARSE_OPERATION_NON_TRANSPOSE
			
			kp_cu_gemm_core(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, A.dim2, B.getDim2(), A.dim1, A.nnz, &alpha,
				       A.descr, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.getData(), B.getDim1(), &beta, C.getData(), C.getDim1());
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

}


/*template <typename real>
void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublashandle,  char op_A, char op_B, real alpha, kp_cu_smatrix<real>& A, kp_cu_matrix<real>& B, real beta, kp_cu_matrix<real>& C)
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
                A.convert2csr(handle);	

	
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
	        A.convert2csrT(handle);

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


// C = op_A(A) * op_B(B) 
template <typename real>
void kp_cu_sgemm(cusparseHandle_t cusparseHandle, char op_A, char op_B, kp_cu_smatrix<real>& A, kp_cu_smatrix<real>& B, kp_cu_smatrix<real>& C)
{
	int opA_dim1, opA_dim2;
	int opB_dim1, opB_dim2;

	int nnzC=-1, baseC=-1;
	int *nnzTotalDevHostPtr = &nnzC;

	cusparseOperation_t transA, transB;
	cusparseOperation_t transAinv, transBinv;

	
	cusparseStatus_t status;

	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2, &transA, &transAinv);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2, &transB, &transBinv);
	
	if (C.dim1 != opA_dim1 || opA_dim2 != opB_dim1 || C.dim2 != opB_dim2)
	{
		cerr<<"Error | kp_cu_sgemm (sparse) | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}

	C.resize(1, C.dim1,C.dim2);
	if (A.get_majorDim() == 'R') //si row-major
	{	
		A.convert2csr(cusparseHandle);		

		if(B.get_majorDim() == 'R')	
		{
			B.convert2csr(cusparseHandle);

			status = cusparseXcsrgemmNnz(cusparseHandle, transA, transB, C.dim1, C.dim2, opA_dim2, 
        A.descr, A.nnz, A.csrRowPtr, A.colind_cu,
        B.descr, B.nnz, B.csrRowPtr, B.colind_cu,
        C.descr, C.csrRowPtr, nnzTotalDevHostPtr);

			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				//throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			

	
			if (nnzTotalDevHostPtr != NULL)
				nnzC = *nnzTotalDevHostPtr;
			else
			{
				kp_cu_cudaMemcpy(&nnzC , C.csrRowPtr+C.dim1, sizeof(int), cudaMemcpyDeviceToHost);
				kp_cu_cudaMemcpy(&baseC, C.csrRowPtr  , sizeof(int), cudaMemcpyDeviceToHost);
				nnzC -= baseC;
			}
			C.resize(nnzC, C.dim1,C.dim2);
			
			kp_cu_sgemm_core(cusparseHandle, transA, transB, C.dim1, C.dim2, opA_dim2, A.descr, A.nnz, A.values_cu, A.csrRowPtr, A.colind_cu, B.descr, B.nnz, B.values_cu, B.csrRowPtr, B.colind_cu, C.descr, C.values_cu, C.csrRowPtr, C.colind_cu);	
		}

		else if(B.get_majorDim() == 'C')
		{
			B.convert2csrT(cusparseHandle);	

			status = cusparseXcsrgemmNnz(cusparseHandle, transA, transBinv, C.dim1, C.dim2, opA_dim2, 
        A.descr, A.nnz, A.csrRowPtr, A.colind_cu,
        B.descr, B.nnz, B.csrRowPtrT, B.rowind_cu,
        C.descr, C.csrRowPtr, nnzTotalDevHostPtr);
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				//throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			

				
			if (nnzTotalDevHostPtr != NULL)
				nnzC = *nnzTotalDevHostPtr;
			else
			{
				kp_cu_cudaMemcpy(&nnzC , C.csrRowPtr+C.dim1, sizeof(int), cudaMemcpyDeviceToHost);
				kp_cu_cudaMemcpy(&baseC, C.csrRowPtr  , sizeof(int), cudaMemcpyDeviceToHost);
				nnzC -= baseC;
			}
			C.resize(nnzC, C.dim1,C.dim2);
		
			kp_cu_sgemm_core(cusparseHandle, transA, transBinv, C.dim1, C.dim2, opA_dim2, A.descr, A.nnz, A.values_cu, A.csrRowPtr, A.colind_cu, B.descr, B.nnz, B.values_cu, B.csrRowPtrT, B.rowind_cu, C.descr, C.values_cu, C.csrRowPtr, C.colind_cu);

		}
		else
		{
			cerr<<"Error | kp_cu_sgemm (sparse) | Sparse matrix B has been defined neither as a row-major nor as a column-major"<<endl;
			exit(EXIT_FAILURE);
		}


	}

	else if (A.get_majorDim() == 'C') //si column-major
	{
		A.convert2csrT(cusparseHandle);	

		if(B.get_majorDim() == 'R')	
		{
			B.convert2csr(cusparseHandle);

			status = cusparseXcsrgemmNnz(cusparseHandle, transAinv, transB, C.dim1, C.dim2, opA_dim2, 
        A.descr, A.nnz, A.csrRowPtrT, A.rowind_cu,
        B.descr, B.nnz, B.csrRowPtr, B.colind_cu,
        C.descr, C.csrRowPtr, nnzTotalDevHostPtr);
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				//throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			
				
			if (nnzTotalDevHostPtr != NULL)
				nnzC = *nnzTotalDevHostPtr;
			else
			{
				kp_cu_cudaMemcpy(&nnzC , C.csrRowPtr+C.dim1, sizeof(int), cudaMemcpyDeviceToHost);
				kp_cu_cudaMemcpy(&baseC, C.csrRowPtr  , sizeof(int), cudaMemcpyDeviceToHost);
				nnzC -= baseC;
			}
			C.resize(nnzC, C.dim1,C.dim2);

			kp_cu_sgemm_core(cusparseHandle, transAinv, transB, C.dim1, C.dim2, opA_dim2, A.descr, A.nnz, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.descr, B.nnz, B.values_cu, B.csrRowPtr, B.colind_cu, C.descr, C.values_cu, C.csrRowPtr, C.colind_cu);
		}


		else if(B.get_majorDim() == 'C')
		{
			B.convert2csrT(cusparseHandle);	
			
			status = cusparseXcsrgemmNnz(cusparseHandle, transAinv, transBinv, C.dim1, C.dim2, opA_dim2, 
        A.descr, A.nnz, A.csrRowPtrT, A.rowind_cu,
        B.descr, B.nnz, B.csrRowPtrT, B.rowind_cu,
        C.descr, C.csrRowPtr, nnzTotalDevHostPtr);
			if (status != CUSPARSE_STATUS_SUCCESS) 
			{
				cerr<<"Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format failed"<<endl;
				//throw KP_CUSPARSE_COO2CSR;
				//exit(EXIT_FAILURE);
			}
			

				
			if (nnzTotalDevHostPtr != NULL)
				nnzC = *nnzTotalDevHostPtr;
			else
			{
				kp_cu_cudaMemcpy(&nnzC , C.csrRowPtr+C.dim1, sizeof(int), cudaMemcpyDeviceToHost);
				kp_cu_cudaMemcpy(&baseC, C.csrRowPtr  , sizeof(int), cudaMemcpyDeviceToHost);
				nnzC -= baseC;
			}
			C.resize(nnzC, C.dim1,C.dim2);
		
			kp_cu_sgemm_core(cusparseHandle, transAinv, transBinv, C.dim1, C.dim2, opA_dim2, A.descr, A.nnz, A.values_cu, A.csrRowPtrT, A.rowind_cu, B.descr, B.nnz, B.values_cu, B.csrRowPtrT, B.rowind_cu, C.descr, C.values_cu, C.csrRowPtr, C.colind_cu);

		}
		else
		{
			cerr<<"Error | kp_cu_sgemm (sparse) | Sparse matrix B has been defined neither as a row-major nor as a column-major"<<endl;
			exit(EXIT_FAILURE);
		}


	}
	else //si indeterminate
	{
		cerr<<"Error | kp_cu_sgemm (sparse) | Sparse matrix A has been defined neither as a row-major nor as a column-major"<<endl;
		exit(EXIT_FAILURE);
	}


	status=cusparseXcsr2coo(cusparseHandle, C.csrRowPtr, C.nnz, C.dim1, C.rowind_cu, cusparseGetMatIndexBase(C.descr));

	if (status != CUSPARSE_STATUS_SUCCESS) 
	{
        	cerr<<"Error | kp_cu_sgemm (sparse) | Conversion from CSR to COO failed"<<endl;
		//throw KP_CUSPARSE_GEMM;
		//exit(EXIT_FAILURE);
	}
	C.isCSRconverted = true;

}


template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans)
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

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans)
{
   if (op == 'N')
     {
	dim1 = M.getDim1();
	dim2 = M.getDim2();
	*trans = CUSPARSE_OPERATION_NON_TRANSPOSE ;
     }
   else if (op == 'T')
     {
	dim1 = M.getDim2();
	dim2 = M.getDim1();
	*trans = CUSPARSE_OPERATION_TRANSPOSE ;
     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim | op should be either N or T"<<endl;
	exit(EXIT_FAILURE);
     }
}

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans, cusparseOperation_t* transinv)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
	*trans = CUSPARSE_OPERATION_NON_TRANSPOSE ;
	*transinv = CUSPARSE_OPERATION_TRANSPOSE ;

     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
	*trans = CUSPARSE_OPERATION_TRANSPOSE ;
	*transinv = CUSPARSE_OPERATION_NON_TRANSPOSE ;

     }
   else
     {
	cerr<<"Error | kp_cu_check_op_set_dim (sparse) | op should be either N or T"<<endl;
	exit(EXIT_FAILURE);
     }
}



#endif
