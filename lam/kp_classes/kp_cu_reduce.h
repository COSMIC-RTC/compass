
#ifndef __LAM_KP_CU_REDUCE_H__
#define __LAM_KP_CU_REDUCE_H__

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "kp_KFPP.h"

template <typename real>
class kp_cu_vector;

template <typename real>
real kp_cu_reduce(const kp_cu_vector<real>& v);

#endif
