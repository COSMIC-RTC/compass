//kp_randn.h
//init vector and vectors as Normally distributed pseudorandom numbers

#ifndef __SEGER__KP_RANDN_H__
#define __SEGER__KP_RANDN_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "kp_vector.h"
#include "kp_matrix.h"


template<typename real>
void kp_randn(kp_vector<real>& v, gsl_rng* rng);

template<typename real>
void kp_randn(kp_matrix<real>& m, gsl_rng* rng);

#include "kp_randn.hpp"

#endif
