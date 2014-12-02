//kp_carma_tools.hpp
//some tools to convert carma_obj <-> kp_matrix, kp_vector



template<typename T, typename U>
void kp_carma_host_obj_to_kp_matrix(carma_host_obj<T>& ch, kp_matrix<U>& k)
{
   if (ch.getDims(0) != 2)
     {
	cerr<<"Error | kp_carma_host_obj_to_kp_matrix | carma_host_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }

   k.resize(ch.getDims(1), ch.getDims(2));
   
   size_t s = ch.getDims(1) * ch.getDims(2);
   
   for (size_t i = 0 ; i < s ; i++)
     *(k.getData(i)) = (U)(ch.getData()[i]);
}

template<typename T, typename U>
void kp_carma_obj_to_kp_matrix(carma_obj<T>& c, kp_matrix<U>& k)
{
   if (c.getDims(0) != 2)
     {
	cerr<<"Error | kp_carma_obj_to_kp_matrix | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<T> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims(1), ch.getDims(2));
   
   size_t s = c.getDims(1) * c.getDims(2);
   
   for (size_t i = 0 ; i < s ; i++)
     *(k.getData(i)) = (U) (ch.getData()[i]);
}
/*void kp_carma_obj_to_kp_matrix(carma_obj<int>& c, kp_matrix& k)
{
   if (c.getDims(0) != 2)
     {
	cerr<<"Error | kp_carma_obj_to_kp_matrix | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<int> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims(1), ch.getDims(2));
   
   size_t s = c.getDims(1) * c.getDims(2);
   
   for (size_t i = 0 ; i < s ; i++)
     *(k.getData(i)) = (real) (ch.getData()[i]);
}*/

template<typename T, typename U>
void kp_carma_host_obj_to_kp_vector(carma_host_obj<T>& ch, kp_vector<U>& k)
{
   if (ch.getDims(0) != 1)
     {
	cerr<<"Error | kp_carma_host_obj_to_kp_vector | carma_host_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   k.resize(ch.getDims(1));  
   for (size_t i = 0 ; i < (size_t)k.size() ; i++)
     *(k.getData(i)) = (U)(ch.getData()[i]);
}
//                                                                                     
template<typename T, typename U>
void kp_carma_obj_to_kp_vector(carma_obj<T>& c, kp_vector<U>& k)
{
   if (c.getDims(0) != 1)
     {
	cerr<<"Error | kp_carma_obj_to_kp_vector | carma_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<T> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims(1));  
   for (size_t i = 0 ; i < (unsigned int)k.size() ; i++)
     *(k.getData(i)) = (U) (ch.getData()[i]);
}
//                                                                                     
/*void kp_carma_obj_to_kp_vector(carma_obj<int>& c, kp_vector& k)
{
   if (c.getDims(0) != 1)
     {
	cerr<<"Error | kp_carma_obj_to_kp_vector | carma_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<int> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims(0));  
   for (size_t i = 0 ; i < (unsigned int)k.size() ; i++)
      *(k.getData(i)) = (real) (ch.getData()[i]);
}*/

//       
template<typename T>                                                                              
size_t kp_carma_obj_get_dim1(carma_obj<T>& cM)
{
   if (cM.getDims(0) != 2)
     {
	cerr<<"Error | kp_carma_obj_get_dim1 | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   return cM.getDims(1);
}
//                                                                                     
template<typename T>                                                                              
size_t kp_carma_obj_get_dim2(carma_obj<T>& cM)
{
   if (cM.getDims(0) != 2)
     {
	cerr<<"Error | kp_carma_obj_get_dim2 | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   return cM.getDims(0);
}
//                                                                                     
template<typename T, typename U>
void kp_kp_vector_to_carma_obj(const kp_vector<T>& k, carma_obj<U>& c)
{
   if (c.getDims(0) != 1)
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | carma_obj is not a vector(1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   if (k.size() != c.getDims(1))
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | sizes are not equal"<<endl;
	exit(EXIT_FAILURE);
     }
   U* f = new U[k.size()];
   for (size_t i = 0 ;i < (unsigned int)k.size() ; i++)
     f[i] = (U) (*(k.getData(i)));
   
   c.host2device(f);   
   delete[] f; f=NULL;
}

	
template<typename T, typename U>
void kp_kp_cu_vector_to_carma_obj(const kp_cu_vector<T>& cu_k, carma_obj<U>& c)
{
   if (c.getDims(0) != 1)
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | carma_obj is not a vector(1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   if (cu_k.size() != c.getDims(1))
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | sizes are not equal"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_memcpy(c.getData(), cu_k.getData(), cu_k.size());
}	
template<typename T, typename U>
void kp_carma_obj_to_kp_cu_vector(const carma_obj<T>& c, kp_cu_vector<U>& cu_k)
{
   if (c.getDims(0) != 1)
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | carma_obj is not a vector(1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   if (cu_k.size() != c.getDims(1))
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | sizes are not equal"<<endl;
	exit(EXIT_FAILURE);
     }
   kernel_memcpy(cu_k.getData(), c.getData(), cu_k.size());
}
