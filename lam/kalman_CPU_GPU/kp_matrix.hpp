//kp_matrix.cpp

#ifndef __SEGER__KP_MATRIX_HPP__
#define __SEGER__KP_MATRIX_HPP__

#ifndef KP_WITH_CARMA
//        
template<typename real>
template<typename T>
kp_matrix<real>::kp_matrix(const kp_matrix<T>& m)
{
   _create(m.getDim1(), m.getDim2());
   for (int i=0 ; i<dim1*dim2 ; i++)
      d[i] = (real) m.getData()[i];     
}

template<typename real>
kp_matrix<real>::kp_matrix(const kp_matrix<real>& m)
{
   _create(m.dim1, m.dim2);
   for (int i=0 ; i<dim1*dim2 ; i++)
      d[i] = (real) m.d[i];     
}


//                                                                                           
template<typename real> 
kp_matrix<real>::kp_matrix(const vector< vector<double> >& m)
{
   size_t d1, d2;   
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
	       cout<<"error | kp_matrix::kp_matrix | vector< vector<double> > is not matrix"<<endl;
	       exit(EXIT_FAILURE);
	    }	
     }
   _create(d1,d2);
   for (size_t i = 0 ; i < m.size() ; i++)
     for (size_t j = 0 ; j < m[i].size() ; j++)
       el(i,j) = m[i][j];
}
//                                                                                           
template<typename real>
void kp_matrix<real>::resize(int dim1_, int dim2_)
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
//                                                                                      
/*real& kp_matrix::at(int i, int j)
{
   if (i < 0 || i >= dim1 || j < 0 || j >= dim2)
     {
	cerr<<"error | kp_matrix::at | index problem"<<endl;
	exit(EXIT_FAILURE);
     }
   return this->operator()(i,j);
}*/
//                                                                                           
template<typename real> 
template<typename T>                                                                                   
void kp_matrix<real>::operator=(const kp_matrix<T>&m)
{
   resize(m.getDim1(), m.getDim2());
   for (int i=0 ; i<dim1*dim2 ; i++)
      d[i] = (real) m.getData()[i];     
   
}
template<typename real> 
void kp_matrix<real>::operator=(const kp_matrix<real>&m)
{
   resize(m.getDim1(), m.getDim2());
   memcpy(d, m.d, sizeof(real) * dim1 * dim2);     
   
}

//                                                                                           
template<typename real> 
void kp_matrix<real>::zeros()
{
   memset(d, 0, sizeof(real) * dim1 * dim2);
   
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::operator/=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] /= val;
   
}
//                                                                                           
template<typename real>
void kp_matrix<real>::operator*=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] *= val;
   
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::operator+=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] += val;
   
}
//                                                                                           
template<typename real>
void kp_matrix<real>::operator-=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] -= val;

}
//                                                                                           
template<typename real>
void kp_matrix<real>::operator+= (const kp_matrix<real>& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_matrix::operator+= | dimension problems"<<endl;
     }
   for (int i = 0 ; i < dim1 * dim2 ; i++)
     d[i] += M.d[i];

}
//  

template<typename real>
void kp_matrix<real>::inverse()
{  
   
   if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::inverse | matrix is not square"<<endl;
	exit(EXIT_FAILURE);
     }
   int info;
   int* ipiv = new int[max(1,min(dim1,dim2))];

   kp_getrf_core(&dim1, &dim2, d, ipiv, &info);
   
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetrf "<<info<<endl;
	exit(EXIT_FAILURE);
     }

      
   real wkopt;
   int lwork = -1;
   
   kp_getri_core(&dim1, d, ipiv, &wkopt, &lwork, &info);
   
   lwork = (int)wkopt;
   real *work = new real[lwork];
   
   kp_getri_core(&dim1, d, ipiv, work, &lwork, &info);
   
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetri "<<info<<endl;
	exit(EXIT_FAILURE);
     }
   delete [] ipiv;
   delete [] work;
}
//                                                                                           
template<typename real>
void kp_matrix<real>::inverse_symmetric()
{   
   inverse();
/*      if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | matrix is not square"<<endl;
	exit(EXIT_FAILURE);
     }
   int* ipiv = new int[max(1,dim1)];
   double wkopt;
   int lwork = -1;
   int info;
   dsytrf( "L", &dim1, d, &dim1, ipiv, &wkopt, &lwork, &info );
   
   lwork = max((int)wkopt, 2 * dim1);
   double *work = new double[lwork];
   
   dsytrf( "L", &dim1, d, &dim1, ipiv, work, &lwork, &info );
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | some error in dsytrf "<<info<<endl;
	exit(EXIT_FAILURE);
     }
   dsytri("L", &dim1, d, &dim1, ipiv, work, &info);
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | some error in dsytri "<<info<<endl;
	exit(EXIT_FAILURE);
     }   
   delete [] ipiv;
   delete [] work;
*/
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::_create(int dim1_, int dim2_)
{
   dim1 = dim1_;
   dim2 = dim2_;
   d = new real[dim1 * dim2];
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::_clear()  
{
   if (d == NULL)
     {
	cerr<<"Error | kp_matrix::_clear| d == NULL"<<endl;
	exit(EXIT_FAILURE);
     }
   if (d != NULL)
     {
	delete [] d;	
     }
   d = NULL;
}



#else

template<typename real> 
kp_matrix<real>::kp_matrix(const kp_matrix<real>& m)
{
   carma_host_object = new carma_host_obj(&m.carma_host_object);
   dim1 = (int) carma_host_object->getDims(1);
   dim2 = (int) carma_host_object->getDims(2);
   d = carma_host_object->getData();
}

template<typename real> 
void kp_matrix<real>::_create(int dim1_, int dim2_)
{
	dim1 = dim1_;
	dim2 = dim2_;
	long dims_data[3] = {2, dim1_, dim2_};
	carma_host_object = new carma_host_obj(dims_data);
	d = carma_host_object->getData();
}
                                                                                          
template<typename real> 
void kp_matrix<real>::_clear()  
{
	if (carma_host_object==NULL || d==NULL)
	{
		cerr<<"Error | kp_matrix::_clear| d == NULL"<<endl;
		exit(EXIT_FAILURE);
	}
	delete  carma_host_object;
	carma_host_object = NULL;
	d = NULL;
}
template<typename real>
void kp_matrix<real>::resize(int dim1_, int dim2_)
{
	_clear();
	_create(dim1_, dim2_);
}

template<typename real> 
void kp_matrix<real>::operator=(const kp_matrix<real>&m)
{
   _clear();
   kp_matrix(m);     
}
                                                                                           
template<typename real>
void kp_matrix<real>::zeros()
{
   dim1_avant = dim1;
   dim2_avant = dim2
   _clear();
   kp_matrix(dim1_avant, dim2_avant);
}

template<typename real> 
void kp_matrix<real>::operator/=(real val)
{
   real* data = carma_host_object->getData();
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     data[i] /= val;
}
//                                                                                           
template<typename real>
void kp_matrix<real>::operator*=(real val)
{
   real* data = carma_host_object->getData();
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     data[i] *= val;
}
//                                                                                             
template<typename real> 
void kp_matrix<real>::operator+=(real val)
{
   real* data = carma_host_object->getData();
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     data[i] += val;
}
//                                                                                           
template<typename real>
void kp_matrix<real>::operator-=(real val)
{
   real* data = carma_host_object->getData();
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     data[i] -= val;
}

template<typename real>
void kp_matrix<real>::operator+= (const kp_matrix<real>& M)
{
   real* data = carma_host_object->getData();
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_matrix::operator+= | dimension problems"<<endl;
     }
   for (int i = 0 ; i < dim1 * dim2 ; i++)
     data[i] += *(M.getData(i));
}

template<typename real> 
void kp_matrix<real>::inverse()
{  
   
   if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::inverse | matrix is not square"<<endl;
	exit(EXIT_FAILURE);
     }
   int info;
   int* ipiv = new int[max(1,min(dim1,dim2))];
   
   kp_getrf_core(&dim1, &dim2, carma_host_object->getData(), ipiv, &info);
   
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetrf "<<info<<endl;
	exit(EXIT_FAILURE);
     }

      
   real wkopt;
   int lwork = -1;
   
   kp_getri_core(&dim1, carma_host_object->getData(), ipiv, &wkopt, &lwork, &info);
   
   lwork = (int)wkopt;
   real *work = new real[lwork];
   
   kp_getri_core(&dim1, carma_host_object->getData(), ipiv, work, &lwork, &info);
   
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetri "<<info<<endl;
	exit(EXIT_FAILURE);
     }
   delete [] ipiv;
   delete [] work;
}



#endif




template<typename real>
kp_matrix<real>::kp_matrix()
{
   _create(0,0);
}
//                                                                                            
template<typename real>
kp_matrix<real>::kp_matrix(int dim1_, int dim2_)
{
   _create(dim1_, dim2_);
}

template<typename real>
kp_matrix<real>::~kp_matrix()
{
   _clear();
}
//                                                                                           


template<typename real> 
void kp_matrix<real>::cols_mean(kp_vector<real>& rez)
{
   rez.resize(dim2);
   rez.zeros();
   if (dim1 == 0 || dim2 == 0)
     return;
   for (size_t j = 0 ; j < (unsigned int)dim2 ; j++) 
     for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
       {
	  rez[j] += el(i, j);
       }
   rez /= (real) dim1;
}
//                                                                                           
template<typename real>
void kp_matrix<real>::rows_mean(kp_vector<real>& rez)
{
   rez.resize(dim1);
   rez.zeros();
   if (dim1 == 0 || dim2 == 0)
     return;
   for (size_t j = 0 ; j < (unsigned int)dim2 ; j++) 
     for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
       {
	  rez[i] += el(i, j);
       }
   rez /= (real) dim2;
}
//                                                                                           
template<typename real>
double kp_matrix<real>::trace()
{
   if (dim1 != dim2)
     {
	cout<<"error | kp_matrix::trace | dim1 != dim2 (matrix must be square)"<<endl;
	exit(EXIT_FAILURE);
     }
   double sum = 0;
   for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
     sum += el(i,i);
   return sum;
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::init_from_matrix(const kp_matrix<real>& M, int r1, int r2, int c1, int c2)
{
   if (r1 < 0 || r1 > r2 || r2 > M.dim1 || 
       c1 < 0 || c1 > c2 || c2 > M.dim2)
     {
	cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
	cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
	cerr<<"Error | kp_matrix::init_from_matrix | index problems"<<endl;
	exit(EXIT_FAILURE);
     }
   resize(r2 - r1, c2 - c1);
   
   for (int i = 0 ; i < r2 - r1 ; i++)
     for (int j = 0 ; j < c2 - c1 ; j++)
       el(i, j) = M(i + r1, j + c1); 
}
//                                                                                           
template<typename real>
void kp_matrix<real>::init_from_rowidx(const kp_matrix<real>& M, const vector<int>& rowidx)
{
   if (this == &M)
     {
	cerr<<"Error | kp_matrix::init_from_rowidx | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }
   resize(rowidx.size(), M.dim2);
   for (unsigned int i = 0 ; i < rowidx.size(); i++)
     {
	if (rowidx[i] < 0 || rowidx[i] >= M.dim1)
	  {
	     cerr<<"error | kp_matrix::init_from_rowidx | index error "<<endl;
	  }
     }
   for (int j = 0 ; j < M.dim2 ; j++)
     for (unsigned int i = 0 ; i < rowidx.size(); i++)
       {
	  el(i,j) = M(rowidx[i], j);
       }
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::set_from_rowsubmatrix(const kp_matrix<real>& subM, const vector<int>& rowidx)
{
   if (this == &subM)
     {
	cerr<<"Error | kp_matrix::kp_matrix::set_from_rowsubmatrix | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }
   if ((unsigned int)subM.dim1 != rowidx.size() || subM.dim2 != dim2)
     {
	cerr<<"Error | kp_matrix::set_from_rowsubmatrix | dimension probles"<<endl;
	exit(EXIT_FAILURE);
     }
   for (unsigned int i = 0 ; i < rowidx.size(); i++)
     {
	if (rowidx[i] < 0 || rowidx[i] >= dim1)
	  {
	     cerr<<"error | kp_matrix::set_from_rowsubmatrix | index error "<<endl;
	  }
     }
   for (int j = 0 ; j < dim2 ; j++)
     for (unsigned int i = 0 ; i < rowidx.size() ; i++)
       {
	  el(rowidx[i], j) = subM(i,j);	  
       }
}
//                                                                                           
template<typename real>
void kp_matrix<real>::make_rowidx(int col, vector<int>& rowidx)
{
   if (col < 0 || col >= dim2)
     {
	cerr<<"Error | kp_matrix::kp_make_rowidx | index problem"<<endl;
	exit(EXIT_FAILURE);
     }
   rowidx.clear();
   for (int i = 0; i < dim1 ; i++)
     {
	if (el(i, col) != 0)
	  {
	     if (fabs(el(i, col)) < 0.1)
	       {
		  cerr<<"Error | kp_matrix::kp_make_rowidx | looks like strange"<<endl;
		  exit(EXIT_FAILURE);
	       }
	     rowidx.push_back(i);
	  }
     }
   
}
//                                                                                           
template<typename real>
void kp_matrix<real>::init_from_transpose(const kp_matrix<real>& M)
{
   resize(M.dim2, M.dim1);
   for (int i = 0 ; i < dim1; i++)
     for (int j = 0 ; j < dim2 ; j++)
       el(i,j) = M(j,i);
}
//                                                                                           
template<typename real>
void kp_matrix<real>::mult_each_row(const kp_vector<real>& v)
{
   if (v.size() != dim1)
     {
	cerr<<"Error | kp_matrix::mult_each_row | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < dim2 ; j++)
     for (int i = 0 ; i < dim1 ; i++)
       el(i,j) *= v[i];
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::set_from_submatrix(const kp_matrix<real>& subM, int r, int c)
{
   if (subM.dim1 + r > dim1 || subM.dim2 + c > dim2 || r < 0 || c < 0)
     {
	cerr<<"Error | kp_matrix::set_from_submatrix | dimension problem"<<endl;
	cerr<<dim1<<"x"<<dim2<<" "<<subM.dim1<<"x"<<subM.dim2<<" "<<r<<" "<<c<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < subM.dim2 ; j++)
     for (int i = 0 ; i < subM.dim1 ; i++)
       {
	  el(i + r, j + c) = subM(i,j);
       }
}
//                                                                                           
template<typename real>
void kp_matrix<real>::init_from_smatrix(const kp_smatrix<real>&m)
{
   resize(m.dim1, m.dim2);
   zeros();
   for (int i = 0 ; i < m.nnz ; i++)
     {
	if (m.rowind[i] < 1 || m.rowind[i] > dim1 || m.colind[i] < 1 || m.colind[i] > dim2)
	  {
	     cerr<<"error | kp_matrix::init_from_smatrix | smatrix problem"<<endl;
	     exit(EXIT_FAILURE);
	  }
	el(m.rowind[i] - 1 , m.colind[i] - 1) = m.values[i];
     }   
}
//                                                                                           
template<typename real>
void kp_matrix<real>::each_column_substract_mult(const kp_vector<real> &v, double val)
{
   kp_matrix& M = *this;
   
   if (v.size() != M.dim1)
     {
	cerr<<"Error | kp_matrix::each_column_substract_mult | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < M.dim2 ; j++)
     for (int i = 0 ; i < M.dim1 ; i++)
       {
	  M(i,j) -= v[i];
	  M(i,j) *= val;
       }
}

template<typename real>
void kp_matrix<real>::each_column_add(const kp_vector<real> &v)
{
   if (v.size() != dim1)
     {
	cerr<<"Error | kp_matrix::each_column_add | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < dim2 ; j++)
     for (int i = 0 ; i < dim1 ; i++)
       {
	  el(i,j) += v[i];
       }
}
//                                                                                           
template<typename real> 
void kp_matrix<real>::add_unit_to_diagonal()
{
   if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::add_unit_to_diagonal | matrix isn't square"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < dim1 ; i++)
     el(i,i) += 1;
}


                                                                                          
template<typename real> 
void kp_gemv(char op_A, real alpha, const kp_matrix<real>& A, const kp_vector<real>& x, real beta, kp_vector<real>& y)
{
   int opA_dim1, opA_dim2;
   kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
//   cout<<"gemv "<<op_A<<" "<<opA_dim1<<"x"<<opA_dim2<<endl;
   if (opA_dim2 != x.size() || opA_dim1 != y.size())
     {
	cerr<<"Error | kp_gemv | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   if (A.getDim1() == 0 || x.size() == 0 || y.size() == 0)
     {
	cerr<<"A.dim1="<<A.getDim1()<<"x.size()="<<x.size()<<" y.size()="<<y.size()<<endl;
	cerr<<"Error | kp_gemv| leading dimension should be > 0"<<endl;
	exit(EXIT_FAILURE);
     }

   const int A_dim1 = A.getDim1();
   const int A_dim2 = A.getDim2();

   int one = 1;

   kp_gemv_core(&op_A, &A_dim1, &A_dim2, &alpha, A.getData(), x.getData(), &one, &beta, y.getData());

}
//                                                                                           
template<typename real>
void kp_gemm(char op_A, char op_B, real alpha, const kp_matrix<real>& A, const kp_matrix<real>& B, real beta, kp_matrix<real>& C)
{
   int opA_dim1, opA_dim2;
   kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
   
   int opB_dim1, opB_dim2;
   kp_check_op_set_dim(op_B, B, opB_dim1, opB_dim2);


   const int A_dim1 = A.getDim1();
   const int B_dim1 = B.getDim1();
   const int C_dim1 = C.getDim1();
   const int C_dim2 = C.getDim2();

   
   if (C.getDim1() != opA_dim1 || opA_dim2 != opB_dim1 || C.getDim2() != opB_dim2)
     {
	cerr<<"Error | kp_gemm | diminsion problem"<<endl;
	exit(EXIT_FAILURE);
     }
   if (A.getDim1() == 0 || B.getDim1() == 0 || C.getDim1() == 0)
     {
	cerr<<"A.dim1="<<A.getDim1()<<" B.dim1="<<B.getDim1()<<" C.dim1="<<C.getDim1()<<endl;
	cerr<<"Error | kp_gemm| leading dimension should be > 0"<<endl;
	exit(EXIT_FAILURE);
     }
   
   kp_gemm_core(&op_A, &op_B, &opA_dim1, &C_dim2, &opA_dim2, &alpha, A.getData(), &A_dim1, B.getData(), &B_dim1, &beta, C.getData(), &C_dim1);
}
//   
template<typename real> 
void kp_syevd(kp_matrix<real>& Q, kp_vector<real>& DQ)
{
   if (Q.getDim1() != Q.getDim2() || Q.getDim1() != DQ.size())
     {
	cerr<<"Error | kp_syevd | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   int n = Q.getDim1();
   
   int lwork   = 2 * n * n +  6 * n + 1;
   real * work = new real[lwork];
   int liwork  = 5 * n + 3;
   int* iwork  = new int[liwork];
   int info;
  
   kp_syevd_core("V", "U", &n, Q.getData(), &n, DQ.getData(), work, &lwork, iwork, &liwork, &info);
   
   delete [] work;
   delete [] iwork;
   

}
//                                                                                            
template<typename real> 
void kp_check_op_set_dim(int op, const kp_matrix<real>&M, int& dim1, int& dim2)
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
	cerr<<"Error | kp_check_op_set_dim | op should ge either N either T"<<endl;
	exit(EXIT_FAILURE);
     }
}

//                                                                                            
template<typename real> 
double kp_calc_diff(const kp_matrix<real>& m1, const kp_matrix<real>& m2)
{
   if (m1.getDim1() != m2.getDim1() || m1.getDim2() != m2.getDim2())
     {
	cerr<<"error | kp_calc_diff | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   double rez = 0;
   for (int j = 0 ; j < m1.getDim2(); j++)
     for (int i = 0 ; i < m1.getDim1() ; i++)
       rez  += fabs(m1(i,j) - m2(i,j));
   return rez;
}
//                                                                                            
template<typename real>
void kp_vertcat(kp_matrix<real>& rez, const vector<kp_matrix<real>*>& ms)
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
	     cerr<<"Error | kp_vertcat | inconsistant second dimension: "<<ms[i]->getDim2()<<" "<<ms[0]->getDim2()<<endl;
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
//                                                                                  
template<typename real> 
void kp_horizcat(kp_matrix<real>& rez, const vector<kp_matrix<real>*>& ms)
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
	     cerr<<"Error | kp_horizcat | inconsistant first dimension: "<<ms[i]->getDim1()<<" "<<ms[0]->getDim1()<<endl;
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
//                                                                                           
template<typename real>
ostream& operator<<(ostream& out, kp_matrix<real>& m)
{
   for (int i = 0 ; i < m.getDim1() ; i++)
     {
	for (int j = 0 ; j < m.getDim2() ; j++)
	  {
	     if (j)
	       out<<" ";
	     out<<m(i,j);
	  }
	out<<endl;
     }
   return out;
}

#endif
