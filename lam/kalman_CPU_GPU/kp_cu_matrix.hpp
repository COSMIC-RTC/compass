//kp_cu_matrix.cpp

#ifndef __SEGER__KP_CU_MATRIX_HPP__
#define __SEGER__KP_CU_MATRIX_HPP__

#ifndef KP_WITH_CARMA
template <typename real>
kp_cu_matrix<real>::kp_cu_matrix()
{
        d_cu =NULL;
	_create(0,0);
}


template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(int dim1_, int dim2_)
{
        d_cu =NULL;
	_create(dim1_, dim2_);
}


template <typename real>
template <typename T>
kp_cu_matrix<real>::kp_cu_matrix(const kp_cu_matrix<T>& m)
{
        d_cu =NULL;
	_create(m.getDim1(), m.getDim2());
	this->operator=(m);
}
template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(const kp_cu_matrix<real>& m)
{
        d_cu =NULL;
	_create(m.dim1, m.dim2);
	this->operator=(m);
}

template <typename real>
template <typename T>
kp_cu_matrix<real>::kp_cu_matrix(const kp_matrix<T>& m)
{
        d_cu =NULL;
	_create(m.getDim1(), m.getDim2());
	this->operator=(m);
}

template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(const vector< vector<real> >& m )
{
        d_cu =NULL;
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
        if(d1*d2>0)
        {
	        //real* d_tmp = new real[d1*d2];
	        for (size_t i = 0 ; i < m.size() ; i++)
	        	for (size_t j = 0 ; j < m[i].size() ; j++)
		        	//d_tmp[j*d1+i] = m[i][j];
			        kp_cu_cudaMemcpy(d_cu+(j*d1+i), &(m[i][j]), sizeof(real), cudaMemcpyHostToDevice);
        }
        else
                _clear();

	//kp_cu_cudaMemcpy(d_cu, d_tmp, sizeof(real) * dim1 * dim2, cudaMemcpyHostToDevice );     

	//delete[] d_tmp;
}




template <typename real>
void kp_cu_matrix<real>::resize(int dim1_, int dim2_)
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


template <typename real> template <typename T>
void kp_cu_matrix<real>::operator=(const kp_cu_matrix<T>&m)
{
   resize(m.getDim1(), m.getDim2());
   //memcpy(d, m.d, sizeof(real) * dim1 * dim2);
   if (dim1*dim2>0) 
      kernel_memcpy(d_cu, m.getData(), dim1*dim2);
   else
      _clear();
}

template <typename real>
void kp_cu_matrix<real>::operator=(const kp_cu_matrix<real>&m)
{
   resize(m.dim1, m.dim2);
   //memcpy(d, m.d, sizeof(real) * dim1 * dim2);
   if (dim1*dim2>0) 
      kernel_memcpy(d_cu, m.d_cu, dim1*dim2);
   else
      _clear();
}



template <typename real>
template <typename T>
void kp_cu_matrix<real>::operator=(const kp_matrix<T>& m)
{
   resize(m.getDim1(), m.getDim2());
   if (dim1*dim2>0)   
   {  
      if(typeid(real) == typeid(T))
      {
         kp_cu_cudaMemcpy(d_cu, m.getData(), dim1*dim2*sizeof(real),cudaMemcpyHostToDevice);
      }
      else
      {
         const T* data_src = m.getData();
         real* data_tmp = new real[dim1*dim2];
         for (int i=0 ; i<dim1*dim2 ; i++)
            data_tmp[i] = (real) data_src[i];
         kp_cu_cudaMemcpy(d_cu, data_tmp, dim1*dim2*sizeof(real),cudaMemcpyHostToDevice);
	 delete[] data_tmp ; data_tmp=NULL;
      }
   }
   else
        _clear();
   
}

template <typename real>
void kp_cu_matrix<real>::init_from_transpose(cublasHandle_t handle, const kp_cu_matrix<real>& M)
{
   resize(M.dim2, M.dim1);
   if (dim1*dim2>0)
   { 
	
	real alpha = 1;
	real beta = 0;

	cublasOperation_t trans = CUBLAS_OP_T;
	cublasOperation_t not_trans = CUBLAS_OP_N;

	int d1 = M.getDim1();
	int d2 = M.getDim2();

	kp_cu_geam_core(handle, trans, not_trans, d2, d1,
			&alpha, M.getData(), d1,
			&beta, getData(), d2,
			getData(), d2);
   }
   else
        _clear();
}



/*// calcule op(A)+op(B) (possibilite d'avoir A = *this)
template <typename real>
void kp_cu_matrix::sum(char opA, char opB, kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B)
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




	kp_cu_geam_core(handle, transa, transb, dim1, dim2,
			&alpha, A_data, A.dim1,
			&beta, B.d_cu, B.dim1,
			d_cu, dim1);

	if (this == &A) cudaFree(A_data);


}*/

template <typename real>
void kp_cu_matrix<real>::inverse()
{
	//kp_cu_timer temps_getrf,temps_getri; 
	if (dim1 != dim2)
	{
		cerr<<"Error | kp_cu_matrix::inverse | matrix is not square"<<endl;
		exit(EXIT_FAILURE);
	}

	int info;
	int* ipiv = new int[dim1];

	kp_cu_getrf_core(dim1, dim2, d_cu, dim1, ipiv, &info);

	magma_int_t lwork = (kp_cu_get_getri_nb_core<real>(dim1))*dim1;

	real* dwork;
	kp_cu_cudaMalloc((void**)&dwork, sizeof(real) * lwork);

	kp_cu_getri_core(dim1, d_cu, dim1, ipiv, dwork, lwork, &info);

	delete[] ipiv;
	kp_cu_cudaFree(dwork);
        dwork=NULL;

}


template <typename real>
void kp_cu_matrix<real>::init_from_smatrix(const kp_cu_smatrix<real>& M)
{
	resize(M.getDim1(), M.getDim2());
	//zeros();
   if (dim1*dim2>0)
	kernel_sparse2full(d_cu, M.rowind_cu, M.colind_cu, M.values_cu, M.nnz, M.getDim1(), M.getDim2());
   else
        _clear();   
}


template <typename real>
void kp_cu_matrix<real>::_create(int dim1_, int dim2_)
{
	dim1 = dim1_;
	dim2 = dim2_;
	
	if(dim1*dim2>0)
	{
                if(d_cu != NULL)
                        cerr<<"Error | kp_cu_matrix::_create | memory has already been allocated";
		kp_cu_cudaMalloc((void**)&d_cu, sizeof(real)*dim1*dim2);
        	kernel_memset(d_cu, (real)0.0, dim1*dim2);
	}
        else
                d_cu=NULL;	
}


template <typename real>
void kp_cu_matrix<real>::_clear()  
{
     if(dim1*dim2>0) 
     {
   	if (d_cu == NULL)
     	{
		cerr<<"Error | kp_cu_matrix::_clear| d_cu == NULL ; dim1="<<dim1<<", dim2="<<dim2<<endl;
		exit(EXIT_FAILURE);
	}
        //kernel_memset(d_cu, (real)0.0, dim1 * dim2);
	kp_cu_cudaFree(d_cu);
        d_cu=NULL;
     }
     else
        d_cu = NULL;
}

template <typename real>
void kp_cu_matrix<real>::gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, real beta)
{
	int opA_dim1, opA_dim2;	
	int opB_dim1, opB_dim2;
	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2);
	if (dim1 != opA_dim1 || opA_dim2 != opB_dim1 || dim2 != opB_dim2)
	{
		cerr<<"Error | kp_cu_matrix::gemm | diminsion problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.dim1 < 1 || B.dim1 < 1 || dim1 < 1)
	{
		cerr<<"A.dim1="<<A.dim1<<" B.dim1="<<B.dim1<<" C.dim1="<<dim1<<endl;
		cerr<<"Error | kp_cu_matrix::gemm| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}
   
	kp_cu_gemm(handle, op_A, op_B, alpha, A, B, beta, *this);
}



#else


template <typename real>
kp_cu_matrix<real>::kp_cu_matrix()
{
	d_cu=NULL;
	dim1 = 0;
	dim2=0;
	carma_object=NULL;

}
template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(carma_context* context_)
{
        d_cu =NULL;
	carma_object=NULL;
	_create(context_,0,0);
}

template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(carma_context* context_, int dim1_, int dim2_)
{
        d_cu =NULL;
	carma_object=NULL;
	_create(context_, dim1_, dim2_);
}

template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(const kp_cu_matrix<real>& m)
{
        d_cu =NULL;
	carma_object=NULL;
	carma_object = new carma_obj(&m.carma_object);
	dim1 = (int) carma_object->getDims(1);
	dim2 = (int) carma_object->getDims(2);
	d_cu = carma_object->getData();
}

template <typename real>
kp_cu_matrix<real>::kp_cu_matrix(carma_context* context_, const kp_matrix<real>& m)
{
        d_cu =NULL;
	_create(context_, m.getDim1(), m.getDim2());
  	carma_object->host2device(m.getData(), m.carma_object->getNbElem());

}


template <typename real>
void kp_cu_matrix<real>::resize(int dim1_, int dim2_)
{
	_clear();
	_create(carma_object->getContext(), dim1_, dim2_);

}

template <typename real>
void kp_cu_matrix<real>::operator=(const kp_cu_matrix<real>&m)
{
   resize(m.dim1, m.dim2);
   if (dim1*dim2>0)     
	carma_object->copyFrom(m.getData(), m.carma_object->getNbElem());		
   else
        _clear();
   
}

template <typename real>
void kp_cu_matrix<real>::operator=(const kp_matrix<real>& m)
{
   resize(m.getDim1(), m.getDim2());
   if (dim1*dim2>0)     
  	carma_object->host2device(m.getData(), m.carma_object->getNbElem());
   else
        _clear();
   
}


template <typename real>
void kp_cu_matrix<real>::init_from_transpose(const kp_cu_matrix<real>& M)
{
   resize(M.dim2, M.dim1);
   if (dim1*dim2>0)
   { 
	carma_object->transpose(M.carma_object);		
   }
   else
           _clear();
}


template <typename real>
void kp_cu_matrix<real>::gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, real beta)
{
	int opA_dim1, opA_dim2;	
	int opB_dim1, opB_dim2;


	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
	kp_cu_check_op_set_dim(op_B, B, opB_dim1, opB_dim2);
	

	if (dim1 != opA_dim1 || opA_dim2 != opB_dim1 || dim2 != opB_dim2)
	{
		cerr<<"Error | kp_cu_matrix::gemm | diminsion problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.dim1 < 1 || B.dim1 < 1 || dim1 < 1)
	{
		cerr<<"A.dim1="<<A.dim1<<" B.dim1="<<B.dim1<<" C.dim1="<<dim1<<endl;
		cerr<<"Error | kp_cu_matrix::gemm| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}

    
    carma_object->gemm(op_A, op_B, alpha, A.carma_object, A.dim1, B.carma_object, B.dim1, beta, dim1);
}

template <typename real>
void kp_cu_matrix<real>::_create(carma_context* context_, int dim1_, int dim2_)
{
	if(dim1_*dim2_>0)
	{
                if(d_cu!=NULL || carma_object!=NULL)
		{
                        cerr<<"Error | kp_cu_matrix::_create | memory has already been allocated";
			exit(EXIT_FAILURE);
		}
		long dims_data[3]={2, dim1_, dim2_};
		carma_object = new carma_obj(context_, dims_data);
		dim1 = (int) carma_object->getDims(1);
		dim2 = (int) carma_object->getDims(2);
		d_cu = carma_object->getData();
	}
        else
	{
		dim1 =0 ;
		dim2 = 0;
                d_cu=NULL;
		carma_object=NULL;
	}
}

template <typename real>
void kp_cu_matrix<real>::_clear()  
{
	if(dim1*dim2>0) 
	{
   		if (d_cu==NULL || carma_object==NULL)
     		{
			cerr<<"Error | kp_cu_matrix::_clear| d_cu == NULL ; dim1="<<dim1<<", dim2="<<dim2<<endl;
			exit(EXIT_FAILURE);
		}
		delete carma_object;
	}
	carma_object=NULL;
        d_cu=NULL;
}

#endif


template <typename real>
kp_cu_matrix<real>::~kp_cu_matrix()
{
   _clear();
}

template <typename real>
void kp_cu_matrix<real>::zeros()
{
   if (dim1*dim2>0)     
        kernel_memset(d_cu, (real)0.0, dim1 * dim2);
   
}

template <typename real>
void kp_cu_matrix<real>::init_from_transpose(const kp_matrix<real>& M)
{
	resize(M.getDim2(), M.getDim1());
   if (dim1*dim2>0)
   { 
	real* d_tmp = new real[dim1*dim2];
	for (int i = 0 ; i < dim1; i++)
		for (int j = 0 ; j < dim2 ; j++)
			d_tmp[j*dim1+i] = M(j,i);

	kp_cu_cudaMemcpy(d_cu, d_tmp, sizeof(real)*dim1*dim2, cudaMemcpyHostToDevice);

	delete[] d_tmp;
   }
   else
           _clear();

	
}

// M[r1:r2-1 : c1:c2-1] avec  0 <= r1,r2 <= DIM1-1 ; 0 <= c1,c2 <= DIM2-1 ; r1<r2 ; c1<c2 
template <typename real>
void kp_cu_matrix<real>::init_from_matrix(const kp_matrix<real>& M, int r1, int r2, int c1, int c2)
{
	if (r1 < 0 || r1 > r2 || r2 > M.getDim1() || 
			c1 < 0 || c1 > c2 || c2 > M.getDim2())
	{
		//cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
		//cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
		cerr<<"Error | kp_cu_matrix::init_from_matrix | index problems"<<endl;
		exit(EXIT_FAILURE);
	}
	resize(r2 - r1, c2 - c1);
	//real* d_tmp = new real[r2-r1, c2-c1];

   if (dim1*dim2>0) 
   {    
	for (int j = 0 ; j < c2 - c1 ; j++)
		//d_tmp[j * dim1 + i] = M(i + r1, j + c1); 	
		kp_cu_cudaMemcpy(getData(j*dim1), M.getData((j+c1)*M.getDim1()), (r2-r1)*sizeof(real), cudaMemcpyHostToDevice);
   }
   else
           _clear();


	//delete[] d_tmp;

}

template <typename real>
void kp_cu_matrix<real>::init_from_matrix(const kp_cu_matrix<real>& M, int r1, int r2, int c1, int c2)
{
	if (r1 < 0 || r1 > r2 || r2 > M.dim1 || 
			c1 < 0 || c1 > c2 || c2 > M.dim2)
	{
		//cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
		//cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
		cerr<<"Error | kp_cu_matrix::init_from_matrix | index problems"<<endl;
		exit(EXIT_FAILURE);
	}
	resize(r2 - r1, c2 - c1);
	//real* d_tmp = new real[r2-r1, c2-c1];

   if (dim1*dim2>0) 
	kernel_set_submatrix(d_cu, M.d_cu, M.dim1, r1, c1, r2-r1, c2-c1);
   else
           _clear();


	//delete[] d_tmp;

}

template <typename real>
void kp_cu_matrix<real>::operator/=(real val)
{
   if (dim1*dim2>0)     
        kernel_div_const(d_cu, val, dim1*dim2);

}
template <typename real>
void kp_cu_matrix<real>::operator*=(real val)
{
   if (dim1*dim2>0)     
  kernel_mult_const(d_cu, val, dim1*dim2);

}
template <typename real>
void kp_cu_matrix<real>::operator+=(real val)
{
   if (dim1*dim2>0)     
  kernel_add_const(d_cu, val, dim1*dim2);

}
template <typename real>
void kp_cu_matrix<real>::operator-=(real val)
{
   if (dim1*dim2>0)     
  kernel_sub_const(d_cu, val, dim1*dim2);

}

template <typename real>
void kp_cu_matrix<real>::operator+= (const kp_cu_matrix<real>& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_cu_matrix::operator+= | dimension problems"<<endl;
     }
   if (dim1*dim2>0)     
   kernel_add(d_cu, M.d_cu, dim1*dim2);

}
template <typename real>
void kp_cu_matrix<real>::operator-= (const kp_cu_matrix<real>& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_cu_matrix::operator-= | dimension problems"<<endl;
     }
   if (dim1*dim2>0)     
   kernel_sub(d_cu, M.d_cu, dim1*dim2);

}

template <typename real>
void kp_cu_matrix<real>::set_from_submatrix(const kp_matrix<real>& subM, int r, int c)
{
	if (subM.getDim1() + r > dim1 || subM.getDim2() + c > dim2 || r < 0 || c < 0)
	{
		cerr<<"Error | kp_cu_matrix::set_from_submatrix | dimension problem"<<endl;
		cerr<<dim1<<"x"<<dim2<<" "<<subM.getDim1()<<"x"<<subM.getDim2()<<" "<<r<<" "<<c<<endl;
		exit(EXIT_FAILURE);
	}
	for (int j = 0 ; j < subM.getDim2() ; j++)
		kp_cu_cudaMemcpy(getData((j+c)*dim1+r), subM.getData(j*subM.getDim1()), sizeof(real)*subM.getDim1(), cudaMemcpyHostToDevice);
	
   
}


template <typename real>
void kp_cu_matrix<real>::set_from_submatrix(const kp_cu_matrix<real>& subM, int r, int c)
{
	if (subM.dim1 + r > dim1 || subM.dim2 + c > dim2 || r < 0 || c < 0)
	{
		cerr<<"Error | kp_cu_matrix::set_from_submatrix | dimension problem"<<endl;
		cerr<<dim1<<"x"<<dim2<<" "<<subM.getDim1()<<"x"<<subM.getDim2()<<" "<<r<<" "<<c<<endl;
		exit(EXIT_FAILURE);
	}
	for (int j = 0 ; j < subM.getDim2() ; j++)
		kp_cu_cudaMemcpy(getData((j+c)*dim1+r), subM.getData(j*subM.getDim1()), sizeof(real)*subM.getDim1(), cudaMemcpyDeviceToDevice);
		//kernel_memcpy(d_cu+((j+c)*dim1+r), subM.d_cu+(j*subM.dim1), subM.dim1);
	
   
}








template <typename real>
void kp_cu_vertcat(kp_cu_matrix<real>& rez, const vector<kp_matrix<real>*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim1 += ms[i]->getDim1();
	if (ms[i]->getDim2() != ms[0]->getDim2())
	  {
	     cerr<<"Error | kp_cu_vertcat | inconsistant second dimension: "<<ms[i]->getDim2()<<" "<<ms[0]->getDim2()<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(sum_dim1, ms[0]->getDim2());
   int dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), dim1, 0);
	dim1 += ms[i]->getDim1();
     }
}

template <typename real>
void kp_cu_vertcat(kp_cu_matrix<real>& rez, const vector<kp_cu_matrix<real>*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim1 += ms[i]->getDim1();
	if (ms[i]->getDim2() != ms[0]->getDim2())
	  {
	     cerr<<"Error | kp_cu_vertcat | inconsistant second dimension: "<<ms[i]->getDim2()<<" "<<ms[0]->getDim2()<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(sum_dim1, ms[0]->getDim2());
   int dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), dim1, 0);
	dim1 += ms[i]->getDim1();
     }
}

template <typename real>
void kp_cu_horizcat(kp_cu_matrix<real>& rez, const vector<kp_matrix<real>*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim2 += ms[i]->getDim2();
	if (ms[i]->getDim1() != ms[0]->getDim1())
	  {
	     cerr<<"Error | kp_cu_horizcat | inconsistant first dimension: "<<ms[i]->getDim1()<<" "<<ms[0]->getDim1()<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(ms[0]->getDim1(), sum_dim2);
   int dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), 0, dim2);
	dim2 += ms[i]->getDim2();
     }
}

template <typename real>
void kp_cu_horizcat(kp_cu_matrix<real>& rez, const vector<kp_cu_matrix<real>*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim2 += ms[i]->getDim2();
	if (ms[i]->getDim1() != ms[0]->getDim1())
	  {
	     cerr<<"Error | kp_cu_horizcat | inconsistant first dimension: "<<ms[i]->getDim1()<<" "<<ms[0]->getDim1()<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(ms[0]->getDim1(), sum_dim2);
   int dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), 0, dim2);
	dim2 += ms[i]->getDim2();
     }
}




#endif
