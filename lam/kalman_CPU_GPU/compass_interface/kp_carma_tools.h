//kp_carma_tools.h
//some tools to convert carma_obj <-> kp_matrix, kp_vector


#ifndef __SEGER__KP_CARMA_TOOLS_H__
#define __SEGER__KP_CARMA_TOOLS_H__

#include <carma_obj.h>
#include <kp_matrix.h>

void kp_carma_obj_to_kp_matrix(carma_obj<float>& c, kp_matrix& k);
void kp_carma_obj_to_kp_vector(carma_obj<float>& c, kp_vector& k);
size_t kp_carma_obj_get_dim1(carma_obj<float>& M);
void kp_kp_vector_to_carma_obj(const kp_vector& k, carma_obj<float>& c);
#endif
