//kp_cu_vector.cpp
#include "kp_cu_vector.h"
#include "kp_cu_matrix.h"
#include <math.h>
#include <iostream>
#include <string.h>
using namespace std;

#ifdef __WITH_FLOPS_CALCULATOR__
#include "kp_flopscalc.h"
#endif







kp_cu_vector::kp_cu_vector(const kp_cu_vector& v)
{
	_create(v.s);
	//kp_cu_cudaMemcpy(d_cu, v.d_cu, sizeof(real) * s, cudaMemcpyDeviceToDevice);
	kernel_memcpy_real(d_cu, v.d_cu, s);
}

kp_cu_vector::kp_cu_vector(const kp_vector& v)
{
	_create(v.s);
	kp_cu_cudaMemcpy(d_cu, v.d, sizeof(real) * s, cudaMemcpyHostToDevice);
}

kp_cu_vector::kp_cu_vector(const int vec_size, real val)
{
   _create(vec_size);
   kernel_memset_real(d_cu, val, s);
}

void kp_cu_vector::operator=( const kp_vector& v)
{
   resize(v.size());
   kp_cu_cudaMemcpy(d_cu, v.d, sizeof(real) * v.size(), cudaMemcpyHostToDevice);
  
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
void kp_cu_vector::operator=( const kp_cu_vector& v)
{
   resize(v.size());
   //kp_cu_cudaMemcpy(d_cu, v.d_cu, sizeof(real) * v.size(), cudaMemcpyDeviceToDevice);
   kernel_memcpy_real(d_cu, v.d_cu, v.size());
   
#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}


void kp_cu_vector::operator+=(const kp_cu_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator+= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_add_real(d_cu, v.d_cu, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
void kp_cu_vector::operator-=(const kp_cu_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator-= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_sub_real(d_cu, v.d_cu, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
void kp_cu_vector::operator*=(const kp_cu_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator*= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_mult_real(d_cu, v.d_cu, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
void kp_cu_vector::operator/=(const kp_cu_vector& v)
{
   if (size() != v.size())
     {
	cerr<<"error | kp_cu_vector::operator/= | different vector sizes"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_div(d_cu, v.d_cu, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}



void kp_cu_vector::operator+=(real val)
{
   kernel_add_const_real(d_cu, val, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
void kp_cu_vector::operator-=(real val)
{
   kernel_sub_const_real(d_cu, val, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
void kp_cu_vector::operator*=(real val)
{
   kernel_mult_const_real(d_cu, val, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
void kp_cu_vector::operator/=(real val)
{
   kernel_div_const(d_cu, val, s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}



void kp_cu_vector::init_from_matrix_column(const kp_cu_matrix& M, int col)
{
   resize(M.dim1); //number of rows in matrix (because we set from column)
   //kp_cu_cudaMemcpy(d_cu, M.d_cu+j*M.dim1, M.dim1*sizeof(real), cudaMemcpyDeviceToDevice);
   kernel_memcpy_real(d_cu, M.d_cu+col*M.dim1, M.dim1);

   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}


void kp_cu_vector::init_from_matrix_column(const kp_matrix& M, int col)
{
   resize(M.dim1); //number of rows in matrix (because we set from column)
   kp_cu_cudaMemcpy(d_cu, M.d+col*M.dim1, M.dim1*sizeof(real), cudaMemcpyHostToDevice);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}

void kp_cu_vector::init_from_vector(const kp_vector& v, int ind_begin, int ind_end)
{
   if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
     {
	cerr<<"Error | kp_cu_vector::init_from_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
	exit(EXIT_FAILURE);
     }
   resize(ind_end - ind_begin);
   kp_cu_cudaMemcpy(d_cu, &v[ind_begin], (ind_end - ind_begin)*sizeof(real), cudaMemcpyHostToDevice); 
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif

}

void kp_cu_vector::init_from_vector(const kp_cu_vector& v, int ind_begin, int ind_end)
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
   //kp_cu_cudaMemcpy(d_cu, &v.d_cu[ind_begin], (ind_end - ind_begin)*sizeof(real), cudaMemcpyDeviceToDevice); 
   kernel_memcpy_real(d_cu, v.d_cu+ind_begin, ind_end - ind_begin);

   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}

void kp_cu_vector::init_from_idx(const kp_vector& v, const vector<int>& idx)
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
	kp_cu_cudaMemcpy(d_cu+i, v.d+idx[i], sizeof(real), cudaMemcpyHostToDevice);
     }
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}
void kp_cu_vector::init_from_idx(const kp_cu_vector& v, const vector<int>& idx)
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
	//kp_cu_cudaMemcpy(d_cu+i, v.d_cu+idx[i], sizeof(real), cudaMemcpyDeviceToDevice);
	kernel_memcpy_real(d_cu+i, v.d_cu+idx[i], 1);
     }
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}

void kp_cu_vector::set_from_subvector(const kp_vector& subv, const vector<int>& idx)
{
   if (subv.size() != idx.size())
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
	kp_cu_cudaMemcpy(d_cu+idx[i], subv.d+i, sizeof(real), cudaMemcpyHostToDevice);
     }
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(idx.size());
   #endif
}
void kp_cu_vector::set_from_subvector(const kp_cu_vector& subv, const vector<int>& idx)
{
   if (this == &subv)
     {
	cerr<<"Error | kp_cu_vector::set_from_subvector | the same vector"<<endl;
	exit(EXIT_FAILURE);
     }
   if (subv.size() != idx.size())
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
	//kp_cudaMemcpy(d_cu+idx[i], subv.d_cu+i, sizeof(real), cudaMemcpyDeviceToDevice);
	kernel_memcpy_real(d_cu+idx[i], subv.d_cu+i, 1);
     }
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(idx.size());
   #endif
}


void kp_cu_vector::zeros()
{
   kp_cu_cudaMemset(d_cu, 0, sizeof(real) * s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(size());
   #endif
}



void kp_cu_vector::resize(int s_)
{
   if (s != s_)
     {
	_clear();
	_create(s_);
     }
}

void kp_cu_vector::_clear()        
{
	if (d_cu == NULL && s>0)
	{
		cerr<<"Error | kp_cu_vector::_clear | d_cu == NULL"<<endl;
		exit(EXIT_FAILURE);
	}
	if(s>0) kp_cu_cudaFree(d_cu);
	d_cu=NULL;
}
  
void kp_cu_inverse(kp_cu_vector& v)
{
   kernel_inv(v.d_cu, v.s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}

void sqrt(kp_cu_vector& v)
{
   kernel_sqrt(v.d_cu, v.s);
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other(v.size());
   #endif
}
