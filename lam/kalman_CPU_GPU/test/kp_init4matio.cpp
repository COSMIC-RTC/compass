//kp_init_from_matio.cpp

#include "kp_init4matio.h"
#include <iostream>
#include <string.h>
#include <math.h>

void kp_memcpy_double2real(real* in, double*from, size_t s)
{
#ifdef KP_SINGLE
   for (size_t i = 0 ; i < s ; i++)
     in[i] = from[i];
#else
     memcpy(in, from, s * sizeof(double));   
#endif
}

void kp_init4matio_smatrix(kp_smatrix& M, mat_t* mat, string name)
{
   matvar_t* var = kp_matio_readvar(mat, name);
   mat_sparse_t *s = (mat_sparse_t*)var->data;
   if (var->class_type != MAT_C_SPARSE || var->data_type != MAT_T_DOUBLE || var->rank != 2 ||
      var->isComplex)
     {
	cerr<<"Error | kp_init4matio_smatrix | variable "<<name<<" is not sparce real matrix"<<endl;
	exit(EXIT_FAILURE);
     }

   if (var->rank != 2 || var->dims[1] + 1 != s->njc || s->nir != s->ndata || s->jc[var->dims[1]] != s->ndata)
     {
	cerr<<"Error | kp_init4matio_smatrix | problems with dimensions"<<endl;
	exit(EXIT_FAILURE);
     }
   M.resize(s->ndata, var->dims[0], var->dims[1]);   
   
   kp_memcpy_double2real(M.values, (double*)s->data, s->ndata);   
   int count = 0;
   for (int i = 0 ; i < M.dim2 ; i++)
     for (int j = s->jc[i] ; j < s->jc[i + 1] ; j++)
       {
	  M.colind[count] = i + 1;
	  M.rowind[count] = s->ir[count] + 1;
	  count++;
       }
   if (count != M.nnz)
     {
	cerr<<"Error | kp_init4matio_smatrix | count != M.nnz"<<endl;
	exit(EXIT_FAILURE);
     }   
   Mat_VarFree(var);
}
//                                                                                            
void kp_init4matio_matrix(kp_matrix&M, mat_t* mat, string name)
{
   matvar_t* var = kp_matio_readvar(mat, name);
   if (var->class_type != MAT_C_DOUBLE || var->data_type != MAT_T_DOUBLE || var->rank != 2 ||
       var->isComplex || var->data_size != sizeof(double))
     {
	cerr<<"Error | kp_init4matio_matrix | variable "<<name<<" is not double real matrix"<<endl;
	cerr<<var->class_type<<" "<<var->data_type<<" "<<var->rank<<endl;
	exit(EXIT_FAILURE);
     }
   M.resize(var->dims[0], var->dims[1]);
   kp_memcpy_double2real(M.d, (double*)var->data, var->dims[0] * var->dims[1]);
   Mat_VarFree(var);   
}
//                                                                                            
void kp_init4matio_vector(kp_vector&V, mat_t* mat, string name)
{
   matvar_t* var = kp_matio_readvar(mat, name);
   if (var->class_type != MAT_C_DOUBLE || var->data_type != MAT_T_DOUBLE || var->rank != 2 || 
      var->isComplex || var->data_size != sizeof(double) || (var->dims[1] != 1 && var->dims[0] != 1))
     {
	cerr<<"Error | kp_init4matio_vector | variable "<<name<<" is not double real vector"<<endl;
	cerr<<"class="<<var->class_type<<" data_type"<<var->data_type<<" rank="<<var->rank<<endl;
	cerr<<"data_size="<<var->data_size<<endl;
	cerr<<"dims[0]"<<var->dims[0]<<endl;
	cerr<<"dims[1]"<<var->dims[1]<<endl;
	exit(EXIT_FAILURE);
     }
   int s = 0;
   if (var->dims[0] != 0 && var->dims[1] != 0)
     s = max(var->dims[0], var->dims[1]);
   V.resize(s);
   kp_memcpy_double2real(V.d, (double*)var->data, s);
   Mat_VarFree(var);
}
//                                                                                            
void kp_init4matio_vector(vector<double>&v, mat_t* mat, string name)
{
   kp_vector kpv;
   kp_init4matio_vector(kpv,mat,name);
   v.resize(kpv.size());
   for (size_t i = 0 ; i < v.size() ; i++)
     v[i] = kpv[i];
}
//                                                                                            
void kp_init4matio_vector(vector<int>&v, mat_t* mat, string name)
{   
   
   kp_vector kpv;
   kp_init4matio_vector(kpv,mat,name);
   v.resize(kpv.size());
   for (size_t i = 0 ; i < v.size() ; i++)
     {
	v[i] = kpv[i];
	if (fabs((double)v[i] - kpv[i]) > 1e-10)
	  {
	     cerr<<"Error | kp_init4matio_vector{int} | "<<kpv[i]<<"not look like integer"<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
}
//                                                                                            
double kp_init4matio_double(mat_t* mat, string name)
{
   double rez;
   matvar_t* var = kp_matio_readvar(mat, name);
   if (var->class_type != MAT_C_DOUBLE || var->data_type != MAT_T_DOUBLE || var->rank != 2 || 
      var->isComplex || var->data_size != sizeof(double) || var->dims[1] != 1 || var->dims[1] != 1)
     {
	cerr<<"Error | kp_init4matio_vector | variable "<<name<<" is not double value"<<endl;
	exit(EXIT_FAILURE);
     }
   
   rez = ((double*)var->data)[0];
   Mat_VarFree(var);
   return rez;
}
//                                                                                           
int kp_init4matio_int(mat_t* mat, string name)
{
   double rezd = kp_init4matio_double(mat,name);
   int rezi = (int) rint(rezd);
   if (fabs((double)rezi - rezd) > 1e-10)
     {
	cerr<<"Error | kp_init4matio_int | "<<rezd<<"not look like integer"<<endl;
	exit(EXIT_FAILURE);
     }
   return rezi;
}
//                                                                                            
matvar_t * kp_matio_readvar(mat_t* mat, string name)
{
   matvar_t* var = Mat_VarRead(mat, name.c_str());
   if (var == NULL)
     {
	Mat_Rewind(mat);
	var = Mat_VarRead(mat, name.c_str());
     }
   if (var == NULL)
     {
	cerr<<"error | kp_matio_readvar | cannot read "<<name<<endl;
	exit(EXIT_FAILURE);
     
     }
   return var;
}
//                                                                                            
void   kp_init4matio_index_vector(vector<int>&v,   mat_t* mat, string name)
{
   kp_init4matio_vector(v, mat, name);
   for (size_t i =0 ; i < v.size() ; i++)
     {
	v[i] -= 1;
	if (v[i] < 0)
	  {
	     cerr<<"Error | kp_init4matio_index_vector | not a index "<<v[i]<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
}
