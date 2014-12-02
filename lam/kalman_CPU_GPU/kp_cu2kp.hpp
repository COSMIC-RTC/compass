
#ifndef __SEGER__KP_CU2KP_HPP__
#define __SEGER__KP_CU2KP_HPP__
template<typename T, typename U>
void kp_cu2kp_smatrix(kp_smatrix<T>& A, const kp_cu_smatrix<U>& M)
{
	A.resize(M.nnz, M.dim1, M.dim2);
	kp_cu_cudaMemcpy(A.rowind, M.rowind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(A.colind, M.colind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);


	if (typeid(T) != typeid(U))
	{
	   U* data_tmp = new U[A.nnz];
	   kp_cu_cudaMemcpy(data_tmp, M.values_cu, A.nnz*sizeof(U), cudaMemcpyDeviceToHost);
	   for (int i=0 ; i<A.nnz ;i++)
              A.values[i] = (T) data_tmp[i];
	   delete [] data_tmp ; data_tmp=NULL;
	}
	else
	{
	   kp_cu_cudaMemcpy(A.values, M.values_cu, A.nnz*sizeof(T), cudaMemcpyDeviceToHost);
	}	
}

template<typename T, typename U>
void kp_cu2kp_matrix(kp_matrix<T>& A, const kp_cu_matrix<U>& M)
{
	A.resize(M.getDim1(), M.getDim2());

	if (typeid(T) != typeid(U))
	{
	   int nb_el = A.getDim1()*A.getDim2();
   	   U* data_tmp = new U[nb_el];
	   kp_cu_cudaMemcpy(data_tmp, M.getData(), sizeof(U) * nb_el, cudaMemcpyDeviceToHost);
	   for (int i=0 ; i<nb_el ;i++)
           *(A.getData(i)) = (T) data_tmp[i];
	   delete [] data_tmp ; data_tmp=NULL;
	}
	else
	   kp_cu_cudaMemcpy(A.getData(), M.getData(), sizeof(T) * M.getDim1() * M.getDim2(), cudaMemcpyDeviceToHost);
}

template<typename T, typename U>
void kp_cu2kp_vector(kp_vector<T>& u, const kp_cu_vector<U>& v)
{
	u.resize(v.size());
	if (typeid(T) != typeid(U))
	{
	   U* data_tmp = new U[v.size()];
	   kp_cu_cudaMemcpy(data_tmp, v.getData(), sizeof(U) * v.size(), cudaMemcpyDeviceToHost);
	   for (int i=0 ; i<u.size() ;i++)
              *(u.getData(i)) = (T) data_tmp[i];
	   delete [] data_tmp ; data_tmp=NULL;	
	}
	else
	   kp_cu_cudaMemcpy(u.getData(), v.getData(), sizeof(T) * v.size(), cudaMemcpyDeviceToHost);
}


template<typename real>
void init_from_kp_cu_vector2kp_vector(kp_vector<real>& u, const kp_cu_vector<real>& v, int ind_begin, int ind_end)
{
	if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
	{
		cerr<<"Error | init_from_kp_cu_vector2kp_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
		exit(EXIT_FAILURE);
	}
	u.resize(ind_end - ind_begin);
	kp_cu_cudaMemcpy(u.getData(), v.getData() + ind_begin, (ind_end - ind_begin)*sizeof(real), cudaMemcpyDeviceToHost );
}
#endif
