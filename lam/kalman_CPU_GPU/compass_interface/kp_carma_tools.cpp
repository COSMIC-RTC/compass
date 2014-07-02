//kp_carma_tools.cpp
//some tools to convert carma_obj <-> kp_matrix, kp_vector

#include "kp_carma_tools.h"
#include <carma_host_obj.h>
#include "kp_real.h"

//
void kp_carma_host_obj_to_kp_matrix(carma_host_obj<float>& ch, kp_matrix& k)
{
   if (ch.getDims()[0] != 2)
     {
	cerr<<"Error | kp_carma_host_obj_to_kp_matrix | carma_host_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }

   k.resize(ch.getDims()[1], ch.getDims()[2]);
   
   size_t s = ch.getDims()[1] * ch.getDims()[2];
   
   for (size_t i = 0 ; i < s ; i++)
     k.d[i] = (real)(ch.getData()[i]);
}
void kp_carma_obj_to_kp_matrix(carma_obj<float>& c, kp_matrix& k)
{
   if (c.getDims()[0] != 2)
     {
	cerr<<"Error | kp_carma_obj_to_kp_matrix | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<float> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims()[1], ch.getDims()[2]);
   
   size_t s = c.getDims()[1] * c.getDims()[2];
   
   for (size_t i = 0 ; i < s ; i++)
     k.d[i] = ch.getData()[i];
}
void kp_carma_obj_to_kp_matrix(carma_obj<int>& c, kp_matrix& k)
{
   if (c.getDims()[0] != 2)
     {
	cerr<<"Error | kp_carma_obj_to_kp_matrix | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<int> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims()[1], ch.getDims()[2]);
   
   size_t s = c.getDims()[1] * c.getDims()[2];
   
   for (size_t i = 0 ; i < s ; i++)
     k.d[i] = ch.getData()[i];
}
//
void kp_carma_host_obj_to_kp_vector(carma_host_obj<float>& ch, kp_vector& k)
{
   if (ch.getDims()[0] != 1)
     {
	cerr<<"Error | kp_carma_obj_to_kp_vector | carma_host_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   k.resize(ch.getDims()[1]);  
   for (size_t i = 0 ; i < (size_t)k.size() ; i++)
     k.d[i] = (real)(ch.getData()[i]);
}
//                                                                                     
void kp_carma_obj_to_kp_vector(carma_obj<float>& c, kp_vector& k)
{
   if (c.getDims()[0] != 1)
     {
	cerr<<"Error | kp_carma_obj_to_kp_vector | carma_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<float> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims()[1]);  
   for (size_t i = 0 ; i < (unsigned int)k.size() ; i++)
     k.d[i] = ch.getData()[i];
}
//                                                                                     
void kp_carma_obj_to_kp_vector(carma_obj<int>& c, kp_vector& k)
{
   if (c.getDims()[0] != 1)
     {
	cerr<<"Error | kp_carma_obj_to_kp_vector | carma_obj is not a vector (1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   carma_host_obj<int> ch(c.getDims());
   ch.cpy_obj(&c, cudaMemcpyDeviceToHost);
   k.resize(ch.getDims()[1]);  
   for (size_t i = 0 ; i < (unsigned int)k.size() ; i++)
     k.d[i] = ch.getData()[i];
}

//                                                                                     
size_t kp_carma_obj_get_dim1(carma_obj<float>& cM)
{
   if (cM.getDims()[0] != 2)
     {
	cerr<<"Error | kp_carma_obj_get_dim1 | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   return cM.getDims()[1];
}
//                                                                                     
size_t kp_carma_obj_get_dim2(carma_obj<float>& cM)
{
   if (cM.getDims()[0] != 2)
     {
	cerr<<"Error | kp_carma_obj_get_dim2 | carma_obj is not a matrix (2D)"<<endl;
	exit(EXIT_FAILURE);
     }
   return cM.getDims()[2];
}
//                                                                                     
void kp_kp_vector_to_carma_obj(const kp_vector& k, carma_obj<float>& c)
{
   if (c.getDims()[0] != 1)
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | carma_obj is not a vector(1D)"<<endl;
	exit(EXIT_FAILURE);
     }
   if (k.size() != c.getDims()[1])
     {
	cerr<<"Error | kp_kp_vector_to_carma_obj | sizes are not equal"<<endl;
	exit(EXIT_FAILURE);
     }
   float *f = new float[k.size()];
   for (size_t i = 0 ;i < (unsigned int)k.size() ; i++)
     f[i] = k[i];
   
   c.host2device(f);   
   delete f;
}
