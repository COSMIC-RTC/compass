#ifndef __A1_2VEC_H__
#define __A1_2VEC_H__

#include "kp_smatrix.h"
#include "kp_vector.h"

template <typename real>
void A1_2vec(kp_vector<real>& A1_00, kp_vector<real>& A1_01,
             const kp_smatrix<real>& A1);

#include "A1_2vec.hpp"

#endif
