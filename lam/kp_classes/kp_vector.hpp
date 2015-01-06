//kp_vector.cpp

#ifndef __SEGER__KP_VECTOR_HPP__
#define __SEGER__KP_VECTOR_HPP__

#ifndef KP_WITH_CARMA
template<typename real>
template<typename T>
kp_vector<real>::kp_vector(const kp_vector<T>&v)
{
   _create(v.size());
   this->operator=(v); 
}
template<typename real>
kp_vector<real>::kp_vector(const kp_vector<real>&v)
{
   _create(v.size());
   this->operator=(v); 
}
	
template<typename real>
void kp_vector<real>::resize(int s_)
{
   if (s != s_)
     {
	_clear();
	_create(s_);
     }
}
//                                                                                           
template<typename real> 
template<typename T>  
void kp_vector<real>::operator=( const kp_vector<T>& v)
{
   resize(v.size());
   if (typeid(real) == typeid(T))
      memcpy(d, v.getData(), sizeof(real) * v.size());
   else
      for (int i=0 ; i<v.size() ; i++)
         d[i] = (real) v[i];
}
template<typename real> 
void kp_vector<real>::operator=( const kp_vector<real>& v)
{
   resize(v.size());
   memcpy(d, v.d, sizeof(real) * v.size());
}

//                                                                                            
template<typename real>  
void kp_vector<real>::zeros()
{
   memset(d, 0, sizeof(real) * s);
   
}


#else

template<typename real> 
kp_vector<real>::kp_vector(const kp_vector<real>&v)
{
   carma_host_object = new carma_host_obj(&v.carma_host_object);
   s = (int) carma_host_object->getDims(1);
   d = carma_host_object->getData();
}
//                                                                                         


template<typename real>
void kp_vector<real>::resize(int s_)
{
	_clear();
	_create(s_);
}

template<typename real>
void kp_vector<real>::operator=( const kp_vector<real>& v)
{
   _clear();
   kp_vector(v);     
}


template<typename real>
void kp_vector<real>::zeros()
{
   s_avant = s;
   _clear();
   kp_vector(s_avant);
}

template<typename real>
void kp_vector<real>::_create(int s_)
{
	s = s_; 
	long dims_data[2] = {1, s_};
	carma_host_object = new carma_host_obj(dims_data);
	d = carma_host_object->getData();
}	


//                                                                                            
template<typename real> 
void kp_vector<real>::_clear()        
{
	if (d==NULL || carma_host_object==NULL)
	{
		cerr<<"Error | kp_vector::_clear | d == NULL"<<endl;
		exit(EXIT_FAILURE);
	}
	delete  carma_host_object;
	carma_host_object = NULL;
	d = NULL;
}

#endif




template<typename real>
kp_vector<real>::kp_vector(){_create(0);}

template<typename real> 
kp_vector<real>::kp_vector(int s_){_create(s_);}

template<typename real>                                                                                   
kp_vector<real>::kp_vector(const vector<double>& v)
{
   _create(v.size());
   for (size_t i = 0 ; i < v.size() ; i++)
     d[i] = v[i];
}
                                                                                         
template<typename real> 
kp_vector<real>::kp_vector(const int vec_size, const double val)
{
   _create(vec_size);
   for (int i = 0 ; i < vec_size ; i++)
     d[i] = val;
}

template<typename real> 
kp_vector<real>::~kp_vector(){_clear();}


template<typename real> 
void kp_vector<real>::operator/=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] /= val;
}
//                                                                                            
template<typename real> 
void kp_vector<real>::operator*=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] *= val;
}
//                                                                                            
template<typename real>
void kp_vector<real>::operator+=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] += val;
}
//                                                                                            
template<typename real>
void kp_vector<real>::operator-=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] -= val;
}
//                                                                                            
template<typename real>
void kp_vector<real>::operator+=(const kp_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator+= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] += v.d[i];
}
//                                                                                            
template<typename real>
void kp_vector<real>::operator-=(const kp_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator-= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] -= v.d[i];
}
//                                                                                            
template<typename real>
void kp_vector<real>::operator*=(const kp_vector<real>& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator*= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] *= v.d[i];
}

template<typename real>
void kp_vector<real>::init_from_matrix_column(const kp_matrix<real>& M, int col)
{
   resize(M.getDim1()); //number of rows in matrix (because we set from column)
   for (int i = 0 ; i < M.getDim1() ; i++)
     el(i) = M(i, col); 
}
//                                                                                            
template<typename real>
void kp_vector<real>::init_from_vector(const kp_vector<real>& v, int ind_begin, int ind_end)
{
   if (this == &v)
     {
	cerr<<"Error | kp_vector::init_from_vector | same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
   for (int i = 0 ; i < ind_end - ind_begin ; i++)
     {
	d[i] = v[ind_begin + i]; 
     }
}
//                                                                                            
template<typename real>
void kp_vector<real>::init_from_idx(const kp_vector<real>& v, const vector<int>& idx)
{
   resize(idx.size());
   for (size_t i = 0 ; i < idx.size(); i++)
     {
	if (idx[i] < 0 || idx[i] >= v.size() )
	  {
	     cerr<<"Error | kp_vector::init_from_idx | Indexing error"<<endl;
	     cerr<<"i="<<i<<" idx[i]="<<idx[i]<<" v.size()="<<v.size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	d[i] = v[idx[i]];
     }
}
//                                                                                            
template<typename real> 
void kp_vector<real>::set_from_subvector(const kp_vector<real>& subv, const vector<int>& idx)
{
   if (this == &subv)
     {
	cerr<<"Error | kp_vector::set_from_subvector | the same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if ((unsigned int)subv.size() != idx.size())
     {
	cerr<<"Error | kp_vector::set_from_subvector | subv.size() != idx.size()"<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < idx.size() ; i++)
     {
	if (idx[i] < 0 || idx[i] >= this->size() )
	  {
	     cerr<<"Error | kp_vector::set_from_subvector | Indexing error"<<endl;
	     cerr<<"idx[i]="<<idx[i]<<" size"<<size()<<endl;
	     exit(EXIT_FAILURE);
	  }
	d[idx[i]] = subv[i];
     }
}
//                                                                                            
template<typename real>
void kp_vector<real>::set_from_subvector(const kp_vector<real>& subv, int start)
{
   if (size() < (subv.size() + start))
     {
	cerr<<"Error | kp_vector::set_from_subvector | size problem"<<endl;
	cerr<<"size()="<<size()<<" start="<<start<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < (unsigned int)subv.size() ; i++)
     d[i + start] = subv[i];
}
//                                                                                            
template<typename real> 
real kp_vector<real>::var()
{
   long double sum   = 0;
   long double sum_2 = 0;
   for (int i = 0 ; i < s ; i++)
     {
	sum   += d[i];
	sum_2 += d[i] * d[i];
     }
   
   long double rez = (sum_2 - sum * sum / (long double)s) / (long double)(s - 1); 
   return rez;
}
//                                                                                            
template<typename real> 
real kp_vector<real>::mean()
{
   long double sum   = 0;
   for (int i = 0 ; i < s ; i++)
     {
	sum   += d[i];
     }   
   
   return sum / (long double)s; 
}
//                                                                                            
template<typename real> 
real kp_vector<real>::sum_sqr()
{
   real rez = 0;
   for (size_t i = 0 ; i < (unsigned int)size() ; i++)
     {
	rez += el(i) * el(i);
     }
   return rez;
}

template<typename real>
void kp_vector<real>::_create(int s_)
{
	s = s_; 
	d = new real[s];
}	


//                                                                                            
template<typename real>
void kp_vector<real>::_clear()        
{
   if (d == NULL)
     {
	cerr<<"Error | kp_vector::_clear | d == NULL"<<endl;
	exit(EXIT_FAILURE);
     }
   delete [] d;
   d = NULL;
}



//                                                                                            
template<typename real> 
void kp_inverse(kp_vector<real>& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     *(v.getData(i)) = 1 / v[i];
}
//                                                                                           
template<typename real>
void sqrt(kp_vector<real>& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     *(v.getData(i))  = sqrt(v[i]);
}
//                                                                                           
template<typename real>
double kp_calc_diff(const kp_vector<real>& v1,const kp_vector<real>& v2)
{
   if (v1.size() != v2.size())
     {
	cerr<<"Error kp_calc_diff different sizes"<<endl;
	cout<<v1.size()<<" "<<v2.size()<<endl;
	exit(EXIT_FAILURE);
     }
     
   double rez = 0;
   for (int i = 0 ; i < v1.size() ; i++)
     rez += fabs(v1[i] - v2[i]);
   
   return rez;
}
//                                                                                           
template<typename real>
void kp_vertcat(kp_vector<real>& rez, vector<kp_vector<real>*> vs)
{
   int sum_s = 0;
   for (size_t i = 0 ; i < vs.size() ; i++)
     sum_s += vs[i]->size();
   
   rez.resize(sum_s);
   int s = 0;
   for (size_t i = 0 ; i < vs.size() ; i++)
     {
	rez.set_from_subvector(*(vs[i]), s);
	s += vs[i]->size();
     }
}
//                                                                                           
template<typename real>
ostream& operator<<(ostream& out, kp_vector<real>& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     {
	if (i)
	  cout<<" ";
	cout<<v[i];
     }
   return out;
}
#endif
