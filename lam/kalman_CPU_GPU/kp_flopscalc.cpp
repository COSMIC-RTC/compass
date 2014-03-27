//kp_flopscalc.cpp


#include "kp_flopscalc.h"
#include <omp.h>
#include <cstdlib>
#include <iostream>
using namespace std;

kp_fpc_data kp_fpc_global;

//                                                                            
long kp_fpc_data::sum()
{
   long rez = 0;
   for (int i = 0 ; i < A_SIZE ; i++)
     rez += a[i];
   return rez;
}
//                                                                            
void kp_fpc_data::reset()
{
   for (int i = 0 ; i < A_SIZE ; i++)
     a[i] = 0;
}
//                                                                            
void kp_fpc_data::check()
{
   if (omp_in_parallel())
     {
	cerr<<"Error | kp_fpc_data::check_activ() | in omp parallel region"<<endl;
	exit(EXIT_FAILURE);
     }
}
//                                                                            
void kp_fpc_data::operator+=(const kp_fpc_data& d)
{
   for (int i = 0 ; i < A_SIZE ; i++)
     a[i] += d.a[i];
}

