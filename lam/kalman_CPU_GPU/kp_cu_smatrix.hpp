//kp_cu_smatrix.cpp

#ifndef __SEGER__KP_CU_SMATRIX_HPP__
#define __SEGER__KP_CU_SMATRIX_HPP__
template <typename real>
kp_cu_smatrix<real>::kp_cu_smatrix()
{       
        values_cu=NULL;
        rowind_cu=NULL;
        colind_cu=NULL;
        csrRowPtr=NULL ; 
        csrRowPtrT=NULL;
	_create(0,0,0);
}

template <typename real>
template <typename T>
kp_cu_smatrix<real>::kp_cu_smatrix(const kp_cu_smatrix<T>& M)
{ 
        values_cu=NULL;
        rowind_cu=NULL;
        colind_cu=NULL;
        csrRowPtr=NULL ; 
        csrRowPtrT=NULL;
	_create(0,0,0);
	resize(M.nnz, M.dim1, M.dim2);
	this->operator=(M);
}
template <typename real>
kp_cu_smatrix<real>::kp_cu_smatrix(const kp_cu_smatrix<real>& M)
{ 
        values_cu=NULL;
        rowind_cu=NULL;
        colind_cu=NULL;
        csrRowPtr=NULL ; 
        csrRowPtrT=NULL;
	_create(0,0,0);
	resize(M.nnz, M.dim1, M.dim2);
	this->operator=(M);
}



template <typename real>
template <typename T> 
void kp_cu_smatrix<real>::operator=(const kp_cu_smatrix<T>& M)
{
	resize(M.nnz, M.dim1, M.dim2);
   	//kernel_memcpy(rowind_cu, M.rowind_cu, nnz);   
   	//kernel_memcpy(colind_cu, M.colind_cu, nnz);
	kp_cu_cudaMemcpy(rowind_cu, M.rowind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
	kp_cu_cudaMemcpy(colind_cu, M.colind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
	kernel_memcpy(values_cu, M.values_cu, nnz);
        
        if (M.isCSRconverted)
	{
   	        //kernel_memcpy(csrRowPtr, M.csrRowPtr, dim1+1);  
		kp_cu_cudaMemcpy(csrRowPtr, M.csrRowPtr, sizeof(int) * (dim1+1), cudaMemcpyDeviceToDevice); 
		isCSRconverted = true;
	}
        if(M.isCSRconvertedT)
	{	
   	        //kernel_memcpy(csrRowPtrT, M.rowind_cu, dim2+1);  
		kp_cu_cudaMemcpy(csrRowPtrT, M.csrRowPtrT, sizeof(int) * (dim2+1), cudaMemcpyDeviceToDevice); 
		isCSRconvertedT = true;
	}
	

	majorDim = M.get_majorDim();

	cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
	cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
	cusparseSetMatIndexBase(descr,cusparseGetMatIndexBase(M.descr));
	cusparseSetMatType(descr, cusparseGetMatType(M.descr));
}

template <typename real>
void kp_cu_smatrix<real>::operator=(const kp_cu_smatrix<real>& M)
{
	resize(M.nnz, M.dim1, M.dim2);
	//kernel_memcpy_real(values_cu, M.values_cu, nnz);
   	//kernel_memcpy_int(rowind_cu, M.rowind_cu, nnz);   
   	//kernel_memcpy_int(colind_cu, M.colind_cu, nnz);
	kp_cu_cudaMemcpy(values_cu, M.values_cu, sizeof(real) * nnz, cudaMemcpyDeviceToDevice); 
	kp_cu_cudaMemcpy(rowind_cu, M.rowind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
	kp_cu_cudaMemcpy(colind_cu, M.colind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
        
        if (M.isCSRconverted)
   	        //kernel_memcpy_int(csrRowPtr, M.csrRowPtr, dim1+1);  
		kp_cu_cudaMemcpy(csrRowPtr, M.csrRowPtr, sizeof(int) * (dim1+1), cudaMemcpyDeviceToDevice); 
        if(M.isCSRconvertedT) 
   	        //kernel_memcpy_int(csrRowPtrT, M.rowind_cu, dim2+1);  
		kp_cu_cudaMemcpy(csrRowPtrT, M.csrRowPtrT, sizeof(int) * (dim2+1), cudaMemcpyDeviceToDevice); 
	

	majorDim = M.majorDim;

	cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
	cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
	cusparseSetMatIndexBase(descr,cusparseGetMatIndexBase(M.descr));
	cusparseSetMatType(descr, cusparseGetMatType(M.descr));
}




template <typename real>
template <typename T>
void kp_cu_smatrix<real>::operator=(const kp_smatrix<T>& M)
{
	resize(M.nnz, M.dim1, M.dim2);
	if (nnz>0)
	{
   	   kp_cu_cudaMemcpy(rowind_cu, M.rowind, nnz*sizeof(int), cudaMemcpyHostToDevice);   
   	   kp_cu_cudaMemcpy(colind_cu, M.colind, nnz*sizeof(int), cudaMemcpyHostToDevice);
           if(typeid(real)==typeid(T))
	      kp_cu_cudaMemcpy(values_cu, M.values, nnz*sizeof(real), cudaMemcpyHostToDevice);
	   else
	   {
	      real* data_tmp = new real[nnz];
	      for (int i=0 ; i<nnz ; i++)
	         data_tmp[i] = (real) M.values[i];  
	      kp_cu_cudaMemcpy(values_cu, data_tmp, nnz*sizeof(real), cudaMemcpyHostToDevice);
 	      delete[] data_tmp ; data_tmp=NULL;
	   }
	}

	majorDim = M.get_majorDim();


}


template <typename real>
template <typename T>
kp_cu_smatrix<real>::kp_cu_smatrix(const kp_smatrix<T>& M)
{
	_create(0,0,0);;
	this->operator=(M);
}

template <typename real>
void kp_cu_smatrix<real>::resize(int nnz_, int dim1_, int dim2_)
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

template <typename real>
void kp_cu_smatrix<real>::_create(int nnz_, int dim1_, int dim2_)
{
   cusparseStatus_t status;
   nnz  = nnz_;
   dim1 = dim1_;
   dim2 = dim2_;


   if(nnz>0)
   { 
        if(values_cu != NULL)
            cerr<<"Error | kp_cu_smatrix::_create | values_cu has already been allocated";
        else
	{
   	    kp_cu_cudaMalloc((void**)&values_cu, nnz*sizeof(real));
            kernel_memset(values_cu, (real)0.0, nnz);
	}

        if(rowind_cu != NULL)
            cerr<<"Error | kp_cu_smatrix::_create | rowind_cu has already been allocated";
        else
	{
	    kp_cu_cudaMalloc((void**)&rowind_cu, nnz*sizeof(int));
            kernel_memset(rowind_cu, 0, nnz);
	}

        if(colind_cu != NULL)
            cerr<<"Error | kp_cu_smatrix::_create | colind_cu has already been allocated";
        else
	{
	    kp_cu_cudaMalloc((void**)&colind_cu, nnz*sizeof(int));
            kernel_memset(colind_cu, 0, nnz);
	}

        if(csrRowPtr != NULL)
            cerr<<"Error | kp_cu_smatrix::_create | csrRowPtr has already been allocated";
        else
	{
   	    kp_cu_cudaMalloc((void**)&csrRowPtr, (dim1+1)*sizeof(int));
            kernel_memset(csrRowPtr, 0, (dim1+1));
	}

        if(csrRowPtrT != NULL)
            cerr<<"Error | kp_cu_smatrix::_create | csrRowPtrT has already been allocated";
        else
	{
   	    kp_cu_cudaMalloc((void**)&csrRowPtrT, (dim2+1)*sizeof(int));
            kernel_memset(csrRowPtrT, 0, (dim2+1));
	}	    
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

template <typename real>
void kp_cu_smatrix<real>::_clear()
{
   cusparseStatus_t status;



   if (nnz>0)
   {
	if (csrRowPtr == NULL || csrRowPtrT==NULL || values_cu==NULL || rowind_cu==NULL || colind_cu==NULL )
	{
		cerr<<"Error | kp_cu_smatrix::_clear | double clear"<<endl;
		exit(EXIT_FAILURE);
	}
        /*kernel_memset(rowind_cu, 0, nnz);
        kernel_memset(colind_cu, 0, nnz);
        kernel_memset(values_cu, (real)0.0, nnz);
        kernel_memset(csrRowPtr, 0, (dim1+1));
        kernel_memset(csrRowPtrT, 0, (dim2+1));*/
   	
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



template <typename real>
void kp_cu_smatrix<real>::init_from_transpose(const kp_cu_smatrix<real>& M)
{
   resize(M.nnz, M.dim2, M.dim1);

   //kernel_memcpy(values_cu, M.values_cu, sizeof(real) * nnz);
   //kernel_memcpy(colind_cu, M.rowind_cu, sizeof(int) * nnz);
   //kernel_memcpy(rowind_cu, M.colind_cu, sizeof(int) * nnz);
	kp_cu_cudaMemcpy(values_cu, M.values_cu, sizeof(real) * nnz, cudaMemcpyDeviceToDevice); 
	kp_cu_cudaMemcpy(rowind_cu, M.rowind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
	kp_cu_cudaMemcpy(colind_cu, M.colind_cu, sizeof(int) * nnz, cudaMemcpyDeviceToDevice); 
   

   if (M.majorDim == 'C') majorDim = 'R';
  else if (M.majorDim == 'R') majorDim = 'C';
  else majorDim = 'U';
     
   
}

template <typename real>
bool kp_cu_smatrix<real>::isColumnMajor()
{
    bool colMajor = true;
    int colm1=0;
    int rowm1=0;

    kp_smatrix<real> A_tmp;
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


template <typename real>
void kp_cu_smatrix<real>::convert2csr(cusparseHandle_t handle)
{
   cusparseStatus_t status;
   if (!isCSRconverted)
   {

      status= cusparseXcoo2csr(handle, rowind_cu, nnz, dim1, csrRowPtr, CUSPARSE_INDEX_BASE_ONE); 

      if (status != CUSPARSE_STATUS_SUCCESS) 
      {
         cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
         throw KP_CUSPARSE_COO2CSR;
         //exit(EXIT_FAILURE);
      }
      isCSRconverted = true;
   }
}

template <typename real>
void kp_cu_smatrix<real>::convert2csrT(cusparseHandle_t handle)
{
   if (!isCSRconvertedT)
   {

   cusparseStatus_t status;
   status= cusparseXcoo2csr(handle, colind_cu, nnz, dim2, csrRowPtrT, CUSPARSE_INDEX_BASE_ONE); 
   if (status != CUSPARSE_STATUS_SUCCESS) 
   {
      cerr<<"Error | kp_cu_gemm (sparse) | Conversion from COO to CSR format failed"<<endl;
      throw KP_CUSPARSE_COO2CSR;
      //exit(EXIT_FAILURE);
   }
   isCSRconvertedT = true;
   }
}


template <typename real>
kp_cu_smatrix<real>::~kp_cu_smatrix()
{
   _clear();
}



#endif
