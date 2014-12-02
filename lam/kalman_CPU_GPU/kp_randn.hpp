//kp_randn.cpp

#ifndef __SEGER__KP_RANDN_HPP__
#define __SEGER__KP_RANDN_HPP__
template<typename real>
void kp_randn(kp_vector<real>& v, gsl_rng* rng)
{
   for (int i = 0 ; i < v.size() ; i++)
     v[i] = (real) gsl_ran_ugaussian(rng);
}
//                                                                                            
template<typename real>
void kp_randn(kp_matrix<real>& m, gsl_rng* rng)
{
   for (int j = 0 ; j < m.getDim2() ; j++)
     for (int i = 0 ; i < m.getDim1() ; i++)
     m(i,j) = (real) gsl_ran_ugaussian(rng);
}

#endif
