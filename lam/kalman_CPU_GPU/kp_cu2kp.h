
#ifndef __KP_CU2KP_H__
#define __KP_CU2KP_H__



#include "kp_smatrix.h"
#include "kp_cu_smatrix.h"
#include "kp_vector.h"
#include "cusparse_v2.h"

void kp_cu2kp_smatrix(kp_smatrix& A, const kp_cu_smatrix& M);
void kp_cu2kp_matrix(kp_matrix& A, const kp_cu_matrix& M);
void kp_cu2kp_vector(kp_vector& u, const kp_cu_vector& v);
void init_from_kp_cu_vector2kp_vector(kp_vector& u, const kp_cu_vector& v, int ind_begin, int ind_end);

#endif




