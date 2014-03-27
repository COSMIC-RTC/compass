#include "kp_cu_reduce.h"



real kp_cu_reduce(const kp_cu_vector& v)
{
	thrust::device_ptr<real> dev_ptr(v.d_cu);
	real somme = thrust::reduce(dev_ptr, dev_ptr+v.size());
	return somme;
}
