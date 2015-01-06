//kp_cu_vector.cpp

#ifndef __SEGER__KP_CU_VECTOR_HPP__
#define __SEGER__KP_CU_VECTOR_HPP__

#ifndef KP_WITH_CARMA
template <typename real>
kp_cu_vector<real>::kp_cu_vector(int s_)  
{
       d_cu =NULL;
       _create(s_);
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector()        
{
        d_cu =NULL;
        _create(0);
}


template <typename real>
template <typename T>
kp_cu_vector<real>::kp_cu_vector(const kp_cu_vector<T>& v)
{
        d_cu =NULL;
	_create(v.size());
	this->operator=(v);
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector(const kp_cu_vector<real>& v)
{
        d_cu =NULL;
	_create(v.size());
	this->operator=(v);
}



template <typename real> template <typename T>
kp_cu_vector<real>::kp_cu_vector(const kp_vector<T>& v)
{
        d_cu =NULL;
	_create(v.size());
	this->operator=(v);
}


template <typename real>
kp_cu_vector<real>::kp_cu_vector(const int vec_size, real val)
{
   d_cu =NULL;
   _create(vec_size);
   if (s>0)
      kernel_memset(d_cu, val, s);
}

template<typename real> template<typename T>
void kp_cu_vector<real>::operator=( const kp_vector<T>& v)
{
   resize(v.size());

   if (s>0)
   {
      if(typeid(T) == typeid(real))
         kp_cu_cudaMemcpy(d_cu, v.getData(), sizeof(real) * v.size(), cudaMemcpyHostToDevice);
      else
      {
         real* data_tmp = new real[v.size()];
	 for (int i=0 ; i< v.size() ; i++)
            data_tmp[i] = (real) v[i];
         kp_cu_cudaMemcpy(d_cu, data_tmp, sizeof(real) * v.size(), cudaMemcpyHostToDevice);
	 delete[] data_tmp ; data_tmp=NULL;
      }
   }
}




template <typename real> template <typename T>
void kp_cu_vector<real>::operator=( const kp_cu_vector<T>& v)
{
   resize(v.size());
   if (s>0)
      kernel_memcpy(d_cu, v.getData(), v.size());
}

template <typename real> 
void kp_cu_vector<real>::operator=( const kp_cu_vector<real>& v)
{
   resize(v.size());
   if (s>0)
      kernel_memcpy(d_cu, v.d_cu, v.size());
}




template <typename real>
void kp_cu_vector<real>::init_from_matrix_column(const kp_cu_matrix<real>& M, int col)
{
   resize(M.getDim1()); //number of rows in matrix (because we set from column)
	if (s>0)
   //kernel_memcpy(d_cu, M.d_cu+col*M.dim1, M.dim1);
   kp_cu_cudaMemcpy(d_cu, M.getData(col*M.getDim1()), M.getDim1()*sizeof(real), cudaMemcpyDeviceToDevice);
}


template <typename real>
void kp_cu_vector<real>::init_from_matrix_column(const kp_matrix<real>& M, int col)
{
   resize(M.getDim1()); //number of rows in matrix (because we set from column)
	if (s>0)
   kp_cu_cudaMemcpy(d_cu, M.getData(col*M.getDim1()), M.getDim1()*sizeof(real), cudaMemcpyHostToDevice);
}

template <typename real>
void kp_cu_vector<real>::init_from_vector(const kp_vector<real>& v, int ind_begin, int ind_end)
{
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
	if (s>0)
   kp_cu_cudaMemcpy(d_cu, v.getData(ind_begin), (ind_end - ind_begin)*sizeof(real), cudaMemcpyHostToDevice); 

}

template <typename real>
void kp_cu_vector<real>::init_from_vector(const kp_cu_vector<real>& v, int ind_begin, int ind_end)
{
   if (this == &v)
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
	if (s>0)
   kp_cu_cudaMemcpy(d_cu, v.getData(ind_begin), (ind_end - ind_begin)*sizeof(real), cudaMemcpyDeviceToDevice); 
   //kernel_memcpy(d_cu, v.d_cu+ind_begin, ind_end - ind_begin);
}




template <typename real>
void kp_cu_vector<real>::resize(int s_)
{
   if (s != s_)
     {
	_clear();
	_create(s_);
     }
}

template <typename real>
void kp_cu_vector<real>::_create(int s_)
{
        s = s_; 
        if(s>0) 
        {
                if(d_cu != NULL)
                    cerr<<"Error | kp_cu_vector::_create | d_cu has already been allocated";
                else
		{
                    kp_cu_cudaMalloc((void**)&d_cu, sizeof(real)*s);
            	    kernel_memset(d_cu, (real)0.0, s);
		}
        }
        else 
                d_cu=NULL;
}


template <typename real>
void kp_cu_vector<real>::_clear()        
{
	if (d_cu == NULL && s>0)
	{
		cerr<<"Error | kp_cu_vector::_clear | d_cu == NULL"<<endl;
		exit(EXIT_FAILURE);
	}
	if(s>0)
        { 
                kernel_memset(d_cu, (real)0.0, s);
                kp_cu_cudaFree(d_cu);
        }
	d_cu=NULL;
}
 



template <typename real>
void kp_cu_vector<real>::gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x, real beta)
{
	int opA_dim1, opA_dim2;
	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);

	if (opA_dim2 != x.size() || opA_dim1 != size())
	{
		cerr<<"Error | kp_cu_gemv | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.getDim1() < 1 || x.size() < 1 || size() < 1)
	{
		cerr<<"A.dim1="<<A.getDim1()<<"x.size="<<x.size()<<" y.size()="<<size()<<endl;
		cerr<<"Error | kp_cu_gemv| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}


	kp_cu_gemv(handle, op_A, alpha, A, x, beta, *this);

}





#else
template <typename real>
kp_cu_vector<real>::kp_cu_vector()        
{
        d_cu = NULL;
	carma_object = NULL;
        s = 0;
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector(carma_context* context_)
{
        d_cu =NULL;
	carma_object=NULL;
	_create(context_,0);
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector(const kp_cu_vector<real>& v)
{
        d_cu =NULL;
	carma_object=NULL;
	carma_object = new carma_obj(&v.carma_object);
	s = (int) carma_object->getDims(1);
	d_cu = carma_object->getData();
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector(carma_context* context_, const kp_vector<real>& v)
{
        d_cu =NULL;
	carma_object=NULL;
	_create(context_, v.size());
  	carma_object->host2device(v.getData(), v.getNbElem());
}

template <typename real>
kp_cu_vector<real>::kp_cu_vector(carma_context* context_, const int vec_size, real val)
{
        d_cu =NULL;
	carma_object=NULL;
   	_create(context_, vec_size);
	if (s>0)
   		kernel_memset(carma_object->getData(), val, s);
}


template <typename real>
void kp_cu_vector<real>::operator=(const kp_vector<real>& v)
{
	resize(v.size());
	if (s>0)
  		carma_object->host2device(v.getData(), v.carma_object->getNbElem());
	else
        	_clear();
		
}

template <typename real>
void kp_cu_vector<real>::operator=( const kp_cu_vector<real>& v)
{
   	resize(v.size());
	if (s>0)
		carma_object->copyFrom(v.getData(), v.carma_object->getNbElem());		
	else
		_clear();
}



template <typename real>
void kp_cu_vector<real>::init_from_matrix_column(const kp_cu_matrix<real>& M, int col)
{
   resize(M.getDim1()); //number of rows in matrix (because we set from column)
	if (s>0)
		carma_object->copyFrom(M.getData(col*M.getDim1()), M.getDim1());		
}


template <typename real>
void kp_cu_vector<real>::init_from_matrix_column(const kp_matrix<real>& M, int col)
{
   resize(M.getDim1()); //number of rows in matrix (because we set from column)
	if (s>0)
  		carma_object->host2device(M.getData(col*M.getDim1()),M.getDim1());
}

template <typename real>
void kp_cu_vector<real>::init_from_vector(const kp_vector<real>& v, int ind_begin, int ind_end)
{
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
	if (s>0)
  		carma_object->host2device(v.getData(ind_begin), ind_end-ind_begin);

}

template <typename real>
void kp_cu_vector<real>::init_from_vector(const kp_cu_vector<real>& v, int ind_begin, int ind_end)
{
   if (this == &v)
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
	if (s>0)
		carma_object->copyFrom(v.getData(ind_begin), ind_end-ind_begin);		
}




template <typename real>
void kp_cu_vector<real>::gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x, real beta)
{
	int opA_dim1, opA_dim2;
	kp_cu_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);

	if (opA_dim2 != x.size() || opA_dim1 != size())
	{
		cerr<<"Error | kp_cu_gemv | dimension problem"<<endl;
		exit(EXIT_FAILURE);
	}
	if (A.dim1 < 1 || x.size() < 1 || size() < 1)
	{
		cerr<<"A.dim1="<<A.dim1<<"x.size()="<<x.size()<<" y.size()="<<size()<<endl;
		cerr<<"Error | kp_cu_gemv| leading dimension should be > 0"<<endl;
		exit(EXIT_FAILURE);
	}

	carma_object->gemv(op_A, alpha, A.carma_object, A.dim1, x.carma_object, 1, beta, 1);

}





template <typename real>
void kp_cu_vector<real>::_create(carma_context* context_, int s_)
{
        s = s_; 
        if(s>0) 
        {
                if(d_cu != NULL)
		{
                    cerr<<"Error | kp_cu_vector::_create | d_cu has already been allocated";
		    exit(EXIT_FAILURE);
		}
		long dims_data[2]={1, s_};
		carma_object = new carma_obj(context_, dims_data);
		s = (int) carma_object->getDims(1);
		d_cu = carma_object->getData();
        }
        else 
                d_cu=NULL;
}


template <typename real>
void kp_cu_vector<real>::_clear()        
{
	if ((d_cu == NULL || carma_object==NULL) && s>0)
	{
		cerr<<"Error | kp_cu_vector::_clear | d_cu == NULL"<<endl;
		exit(EXIT_FAILURE);
	}
	if(s>0)
        { 
		delete carma_object;
        }
	carma_object=NULL;
	d_cu=NULL;
}


#endif


template <typename real>
kp_cu_vector<real>::~kp_cu_vector()
{
   _clear();
}



template <typename real>
void kp_cu_vector<real>::zeros()
{
	if (s>0)
   kernel_memset(d_cu, (real)0.0, s);
}


template <typename real>
void kp_cu_vector<real>::operator+=(const kp_cu_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator+= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
	if (s>0)
   kernel_add(d_cu, v.d_cu, s);
}
template <typename real>
void kp_cu_vector<real>::operator-=(const kp_cu_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator-= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
	if (s>0)
   kernel_sub(d_cu, v.d_cu, s);
}
template <typename real>
void kp_cu_vector<real>::operator*=(const kp_cu_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator*= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
	if (s>0)
   kernel_mult(d_cu, v.d_cu, s);
}
template <typename real>
void kp_cu_vector<real>::operator/=(const kp_cu_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator/= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
	if (s>0)
   kernel_div(d_cu, v.d_cu, s);
}



template <typename real>
void kp_cu_vector<real>::operator+=(real val)
{
	if (s>0)
   kernel_add_const(d_cu, val, s);
}
template <typename real>
void kp_cu_vector<real>::operator-=(real val)
{
	if (s>0)
   kernel_sub_const(d_cu, val, s);
}
template <typename real>
void kp_cu_vector<real>::operator*=(real val)
{
	if (s>0)
   kernel_mult_const(d_cu, val, s);
}
template <typename real>
void kp_cu_vector<real>::operator/=(real val)
{
	if (s>0)
   kernel_div_const(d_cu, val, s);
}



template <typename real>
void kp_cu_vector<real>::init_from_idx(const kp_vector<real>& v, const vector<int>& idx)
{
   resize(idx.size());
	if (s>0)
        {
   for (size_t i = 0 ; i < idx.size(); i++)
     {
	if (idx[i] < 0 || idx[i] >= v.size() )
	  {
	     cerr<<"Error | kp_vector::init_from_idx | Indexing error"<<endl;
	     cerr<<"i="<<i<<" idx[i]="<<idx[i]<<" v.size()="<<v.size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	kp_cu_cudaMemcpy(getData(i), v.getData(idx[i]), sizeof(real), cudaMemcpyHostToDevice);
     }
        }
}
template <typename real>
void kp_cu_vector<real>::init_from_idx(const kp_cu_vector<real>& v, const vector<int>& idx)
{
   resize(idx.size());
	if (s>0)
        {
   for (size_t i = 0 ; i < idx.size(); i++)
     {
	if (idx[i] < 0 || idx[i] >= v.size() )
	  {
	     cerr<<"Error | kp_vector::init_from_idx | Indexing error"<<endl;
	     cerr<<"i="<<i<<" idx[i]="<<idx[i]<<" v.size()="<<v.size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	kp_cu_cudaMemcpy(getData(i), v.getData(idx[i]), sizeof(real), cudaMemcpyDeviceToDevice);
	//kernel_memcpy(d_cu+i, v.d_cu+idx[i], 1);
     }
        }
}

template <typename real>
void kp_cu_vector<real>::set_from_subvector(const kp_vector<real>& subv, const vector<int>& idx)
{
   if ((unsigned int)subv.size() != idx.size())
     {
	cerr<<"Error | kp_cu_vector::set_from_subvector | subv.size() != idx.size()"<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < idx.size() ; i++)
     {
	if (idx[i] < 0 || idx[i] >= this->size() )
	  {
	     cerr<<"Error | kp_cu_vector::set_from_subvector | Indexing error"<<endl;
	     cerr<<"idx[i]="<<idx[i]<<" size"<<size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	kp_cu_cudaMemcpy(getData(idx[i]), subv.getData(i), sizeof(real), cudaMemcpyHostToDevice);
     }
}
template <typename real>
void kp_cu_vector<real>::set_from_subvector(const kp_cu_vector<real>& subv, const vector<int>& idx)
{
   if (this == &subv)
     {
	cerr<<"Error | kp_cu_vector::set_from_subvector | the same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if ((unsigned int)subv.size() != idx.size())
     {
	cerr<<"Error | kp_cu_vector::set_from_subvector | subv.size() != idx.size()"<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < idx.size() ; i++)
     {
	if (idx[i] < 0 || idx[i] >= this->size() )
	  {
	     cerr<<"Error | kp_cu_vector::set_from_subvector | Indexing error"<<endl;
	     cerr<<"idx[i]="<<idx[i]<<" size"<<size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	kp_cu_cudaMemcpy(getData(idx[i]), subv.getData(i), sizeof(real), cudaMemcpyDeviceToDevice);
	//kernel_memcpy(d_cu+idx[i], subv.d_cu+i, 1);
     }
}


template <typename real>
void kp_cu_inverse(kp_cu_vector<real>& v)
{
	if (v.size()>0)
   kernel_inv(v.getData(), v.size());
   
}

template <typename real>
void sqrt(kp_cu_vector<real>& v)
{
	if (v.size()>0)
   kernel_sqrt(v.getData(), v.size());
   
}

#endif
