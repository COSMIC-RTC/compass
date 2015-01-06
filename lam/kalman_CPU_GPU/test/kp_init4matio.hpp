template<typename T> void kp_init4matio_smatrix(kp_smatrix<T>& M, mat_t* mat, string name)
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
   
   for (size_t i = 0 ; i < s->ndata ; i++)
     M.values[i] = (T) (((double*)(s->data))[i]);
   
     
   int count = 0;
   for (int i = 0 ; i < M.getDim2() ; i++)
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
template<typename T> void kp_init4matio_matrix(kp_matrix<T>&M, mat_t* mat, string name)
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
   for (size_t i = 0 ; i < var->dims[0] * var->dims[1] ; i++)
     (M.getData())[i] = (T) (((double*)(var->data))[i]);
   Mat_VarFree(var);   
}
//                                                                                            
template<typename T> void kp_init4matio_vector(kp_vector<T>&V, mat_t* mat, string name)
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
   for (size_t i = 0 ; i < s ; i++)
     (V.getData())[i] = (T) (((double*)(var->data))[i]);
   Mat_VarFree(var);
}
