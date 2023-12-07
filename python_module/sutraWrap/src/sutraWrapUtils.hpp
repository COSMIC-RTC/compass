#ifndef SUTRA_WRAP_UTILS_HPP_
#define SUTRA_WRAP_UTILS_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

template <typename T>
using ArrayFStyle = py::array_t<T, py::array::f_style | py::array::forcecast>;

#endif  // SUTRA_WRAP_UTILS_HPP_