//kp_randn.h
//init vector and vectors as Normally distributed pseudorandom numbers

#ifndef __SEGER__KP_RANDN_H__
#define __SEGER__KP_RANDN_H__
#include <gsl/gsl_rng.h>
#include "kp_matrix.h"


void kp_randn(kp_vector& v, gsl_rng* rng);
void kp_randn(kp_matrix& m, gsl_rng* rng);


#endif
