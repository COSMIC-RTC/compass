#include "kp_cu2kp.h"


void kp_cu2kp_smatrix(kp_smatrix& A, const kp_cu_smatrix& M)
{
	A.resize(M.nnz, M.dim1, M.dim2);
	kp_cu_cudaMemcpy(A.rowind, M.rowind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(A.colind, M.colind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(A.values, M.values_cu, A.nnz*sizeof(real), cudaMemcpyDeviceToHost);

	kp_cu_cudaMemcpy(A.rowind, M.rowind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(A.colind, M.colind_cu, A.nnz*sizeof(int), cudaMemcpyDeviceToHost);
	kp_cu_cudaMemcpy(A.values, M.values_cu, A.nnz*sizeof(real), cudaMemcpyDeviceToHost);

}


void kp_cu2kp_matrix(kp_matrix& A, const kp_cu_matrix& M)
{
	A.resize(M.dim1, M.dim2);
	kp_cu_cudaMemcpy(A.d, M.d_cu, sizeof(real) * M.dim1 * M.dim2, cudaMemcpyDeviceToHost);
	
}


void kp_cu2kp_vector(kp_vector& u, const kp_cu_vector& v)
{
	u.resize(v.size());
	kp_cu_cudaMemcpy(u.d, v.d_cu, sizeof(real) * v.s, cudaMemcpyDeviceToHost);
}

void init_from_kp_cu_vector2kp_vector(kp_vector& u, const kp_cu_vector& v, int ind_begin, int ind_end)
{
	if (ind_begin < 0 || ind_begin > ind_end || ind_end > v.size())
	{
		cerr<<"Error | init_from_kp_cu_vector2kp_vector | v.size() = "<<v.size()<<" ind_begin="<<ind_begin<<" ind_end="<<ind_end<<endl;
		exit(EXIT_FAILURE);
	}
	u.resize(ind_end - ind_begin);
	kp_cu_cudaMemcpy(u.d, v.d_cu + ind_begin, (ind_end - ind_begin)*sizeof(real), cudaMemcpyDeviceToHost );

}


