//kp_randn.cpp
#include "kp_randn.h"
#include <gsl/gsl_randist.h>

void kp_randn(kp_vector& v, gsl_rng* rng)
{
   for (int i = 0 ; i < v.size() ; i++)
     v[i] = gsl_ran_ugaussian(rng);
}
//                                                                                            
void kp_randn(kp_matrix& m, gsl_rng* rng)
{
   for (int j = 0 ; j < m.dim2 ; j++)
     for (int i = 0 ; i < m.dim1 ; i++)
     m(i,j) = gsl_ran_ugaussian(rng);
}

