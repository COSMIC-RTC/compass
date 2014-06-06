//kp_vector.cpp
#include "kp_vector.h"
#include "kp_matrix.h"
#include <math.h>
#include <iostream>
#include <string.h>
using namespace std;

#ifdef __WITH_FLOPS_CALCULATOR__
#include "kp_flopscalc.h"
#endif

kp_vector::kp_vector(const kp_vector&v)
{
   _create(v.size());
   memcpy(d, v.d, sizeof(real) * v.size());
}
//                                                                                         
kp_vector::kp_vector(const vector<double>& v)
{
   _create(v.size());
   for (size_t i = 0 ; i < v.size() ; i++)
     d[i] = v[i];
}
//                                                                                         
kp_vector::kp_vector(const int vec_size, const double val)
{
   _create(vec_size);
   for (int i = 0 ; i < vec_size ; i++)
     d[i] = val;
}
//                                                                                         
void kp_vector::resize(int s_)
{
   if (s != s_)
     {
	_clear();
	_create(s_);
     }
}
//                                                                                           
void kp_vector::operator=( const kp_vector& v)
{
   resize(v.size());
   memcpy(d, v.d, sizeof(real) * v.size());
  
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
//                                                                                            
void kp_vector::zeros()
{
   memset(d, 0, sizeof(real) * s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
//                                                                                            
void kp_vector::operator/=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] /= val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::operator*=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] *= val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::operator+=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] += val;
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::operator-=(real val)
{
   for (int i = 0 ; i < s ; i++)
     d[i] -= val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::operator+=(const kp_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator+= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] += v.d[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
//                                                                                            
void kp_vector::operator-=(const kp_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator-= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] -= v.d[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
//                                                                                            
void kp_vector::operator*=(const kp_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_vector::operator*= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < size() ; i++)
     d[i] *= v.d[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
//                                                                                            
void kp_vector::init_from_matrix_column(const kp_matrix& M, int col)
{
   resize(M.dim1); //number of rows in matrix (because we set from column)
   for (int i = 0 ; i < M.dim1 ; i++)
     el(i) = M(i, col); 
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::init_from_vector(const kp_vector& v, int ind_begin, int ind_end)
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
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::init_from_idx(const kp_vector& v, const vector<int>& idx)
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
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}
//                                                                                            
void kp_vector::set_from_subvector(const kp_vector& subv, const vector<int>& idx)
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
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(idx.size());
   #endif
}
//                                                                                            
void kp_vector::set_from_subvector(const kp_vector& subv, int start)
{
   if (size() < (subv.size() + start))
     {
	cerr<<"Error | kp_vector::set_from_subvector | size problem"<<endl;
	cerr<<"size()="<<size()<<" start="<<start<<endl;
	exit(EXIT_FAILURE);
     }
   for (size_t i = 0 ; i < (unsigned int)subv.size() ; i++)
     d[i + start] = subv[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(subv.size());
   #endif

}
//                                                                                            
real kp_vector::var()
{
   long double sum   = 0;
   long double sum_2 = 0;
   for (int i = 0 ; i < s ; i++)
     {
	sum   += d[i];
	sum_2 += d[i] * d[i];
     }
   
   long double rez = (sum_2 - sum * sum / (long double)s) / (long double)(s - 1); 
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)size() * 3);
   #endif
   
   return rez;
}
//                                                                                            
real kp_vector::mean()
{
   long double sum   = 0;
   for (int i = 0 ; i < s ; i++)
     {
	sum   += d[i];
     }   
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
   
   return sum / (long double)s; 
}
//                                                                                            
real kp_vector::sum_sqr()
{
   real rez = 0;
   for (size_t i = 0 ; i < (unsigned int)size() ; i++)
     {
	rez += el(i) * el(i);
     }
   return rez;
}
//                                                                                            
void kp_vector::_clear()        
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
void kp_inverse(kp_vector& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     v[i] = 1 / v[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
//                                                                                           
void sqrt(kp_vector& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     v[i] = sqrt(v[i]);

   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
//                                                                                           
double kp_calc_diff(const kp_vector& v1,const kp_vector& v2)
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
void kp_vertcat(kp_vector& rez, vector<kp_vector*> vs)
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
ostream& operator<<(ostream& out, kp_vector& v)
{
   for (int i = 0 ; i < v.size() ; i++)
     {
	if (i)
	  cout<<" ";
	cout<<v[i];
     }
   return out;
}
