#include "kp_carma_tools.h"

int main()
{
   carma_context contex;
   //kp_carma_obj_to_kp_matrix
     {
	 long dims[3];
	dims[0] = 2;
	dims[1] = 2;
	dims[2] = 2;
	
	carma_obj<float> obj(&contex, dims);
	carma_host_obj<float> hobj(dims);
	
	for (size_t i = 0; i < 4; i++)
	  *(hobj[i]) = i;
	
	hobj.cpy_obj(&obj, cudaMemcpyHostToDevice);
	
	kp_matrix k;
	kp_carma_obj_to_kp_matrix(obj, k);
	
	cout<<"dim1="<<kp_carma_obj_get_dim1(obj)<<endl;
	for (size_t i = 0 ;i < k.dim1 ; i++)
	  for (size_t j = 0 ;j < k.dim2 ; j++)
	    cout<<i<<" "<<j<<" "<<k(i,j)<<endl;
     }
   cout<<endl;
   //kp_carma_obj_to_kp_vector
     {
	 long dims[2];
	dims[0] = 1;
	dims[1] = 5;
	
	carma_obj<float> obj(&contex, dims);
	carma_host_obj<float> hobj(dims);
	
	for (size_t i = 0; i < 5; i++)
	  *(hobj[i]) = i;
	
	hobj.cpy_obj(&obj, cudaMemcpyHostToDevice);
	
	kp_vector v;
	kp_carma_obj_to_kp_vector(obj, v);
	
	for (size_t i = 0 ;i < v.size() ; i++)
	    cout<<i<<" "<<v[i]<<endl;
     }
   cout<<endl;
   //kp_kp_vector_to_carma_obj
     {
	int N = 10;
	kp_vector v(N);
	for (size_t i = 0 ; i < N ; i++)
	  v[i] = N - i;
	
	long dims[2];
	dims[0] = 1;
	dims[1] = N;
	
	carma_obj<float> obj(&contex, dims);
	carma_host_obj<float> hobj(dims);
	kp_kp_vector_to_carma_obj(v, obj);
	hobj.cpy_obj(&obj, cudaMemcpyDeviceToHost);
	
	for (size_t i = 0; i < N; i++)
	  cout<< *(hobj[i]) <<endl;
	
     }
   
}
