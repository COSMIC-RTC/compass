#include "kp_cu_reduce.h"
#include "kp_cu_vector.h"

template <>
double kp_cu_reduce<double>(const kp_cu_vector<double>& v) {
  thrust::device_ptr<double> dev_ptr(v.getData());
  double somme = thrust::reduce(dev_ptr, dev_ptr + v.size());
  return somme;
}

template <>
float kp_cu_reduce<float>(const kp_cu_vector<float>& v) {
  thrust::device_ptr<float> dev_ptr(v.getData());
  float somme = thrust::reduce(dev_ptr, dev_ptr + v.size());
  return somme;
}
