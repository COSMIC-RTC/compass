#ifndef _DECLARE_NAME_H_
#define _DECLARE_NAME_H_

#include <carma.h>
#include <wyrm>

namespace py = pybind11;

using complex64 = cuFloatComplex;
static_assert(sizeof(complex64) == 8, "Bad size");

using complex128 = cuDoubleComplex;
static_assert(sizeof(complex128) == 16, "Bad size");

namespace pybind11 { namespace detail {

// Similar to enums in `pybind11/numpy.h`. Determined by doing:
// python3 -c 'import numpy as np; print(np.dtype(np.float16).num)'
constexpr int NPY_COMPLEX64 = 14;
constexpr int NPY_COMPLEX128 = 15;

// Kinda following: https://github.com/pybind/pybind11/blob/9bb3313162c0b856125e481ceece9d8faa567716/include/pybind11/numpy.h#L1000
template <>
struct npy_format_descriptor<complex64> {
  static pybind11::dtype dtype() {
    handle ptr = npy_api::get().PyArray_DescrFromType_(NPY_COMPLEX64);
    return reinterpret_borrow<pybind11::dtype>(ptr);
  }
  static std::string format() {
    // following: https://docs.python.org/3/library/struct.html#format-characters
    return "Zf";
  }
  static constexpr auto name() {
    return _("complex64");
  }
};

// Kinda following: https://github.com/pybind/pybind11/blob/9bb3313162c0b856125e481ceece9d8faa567716/include/pybind11/numpy.h#L1000
template <>
struct npy_format_descriptor<complex128> {
  static pybind11::dtype dtype() {
    handle ptr = npy_api::get().PyArray_DescrFromType_(NPY_COMPLEX128);
    return reinterpret_borrow<pybind11::dtype>(ptr);
  }
  static std::string format() {
    // following: https://docs.python.org/3/library/struct.html#format-characters
    return "Zd";
  }
  static constexpr auto name() {
    return _("complex128");
  }
};

template <typename T>
struct npy_scalar_caster {
  PYBIND11_TYPE_CASTER(T, _("PleaseOverride"));
  using Array = wy::array_f<T>;

  bool load(handle src, bool convert) {
    // Taken from Eigen casters. Permits either scalar dtype or scalar array.
    handle type = dtype::of<T>().attr("type");  // Could make more efficient.
    if (!convert && !isinstance<Array>(src) && !isinstance(src, type))
      return false;
    Array tmp = Array::ensure(src);
    if (tmp && tmp.size() == 1 && tmp.ndim() == 0) {
      this->value = *tmp.data();
      return true;
    }
    return false;
  }

  static handle cast(T src, return_value_policy, handle) {
    Array tmp({1});
    tmp.mutable_at(0) = src;
    tmp.resize({});
    // You could also just return the array if you want a scalar array.
    object scalar = tmp[tuple()];
    return scalar.release();
  }
};

template <>
struct type_caster<complex64> : npy_scalar_caster<complex64> {
  static constexpr auto name = _("complex64");
};

template <>
struct type_caster<complex128> : npy_scalar_caster<complex128> {
  static constexpr auto name = _("complex128");
};

}}  // namespace pybind11::detail

#ifdef CAN_DO_HALF

using float16 = __half;
static_assert(sizeof(float16) == 2, "Bad size");

namespace pybind11 { namespace detail {

// Similar to enums in `pybind11/numpy.h`. Determined by doing:
// python3 -c 'import numpy as np; print(np.dtype(np.float16).num)'
constexpr int NPY_FLOAT16 = 23;

// Kinda following: https://github.com/pybind/pybind11/blob/9bb3313162c0b856125e481ceece9d8faa567716/include/pybind11/numpy.h#L1000
template <>
struct npy_format_descriptor<float16> {
  static pybind11::dtype dtype() {
    handle ptr = npy_api::get().PyArray_DescrFromType_(NPY_FLOAT16);
    return reinterpret_borrow<pybind11::dtype>(ptr);
  }
  static std::string format() {
    // following: https://docs.python.org/3/library/struct.html#format-characters
    return "e";
  }
  static constexpr auto name() {
    return _("float16");
  }
};

template <>
struct type_caster<float16> : npy_scalar_caster<float16> {
  static constexpr auto name = _("float16");
};

}}  // namespace pybind11::detail

#endif

template <typename T>
constexpr static char const* explicit_name();

#define DECLARE_NAME_CARMA(Type, Name)       \
  template <>                          \
  constexpr char const* explicit_name<Type>() { \
    return #Name;                      \
  }

// DECLARE_NAME_CARMA(u8, Uint8)
// DECLARE_NAME_CARMA(i8, Int8)

// DECLARE_NAME_CARMA(u16, Uint16)
// DECLARE_NAME_CARMA(i16, Int16)

// DECLARE_NAME_CARMA(u32, Uint32)
// DECLARE_NAME_CARMA(i32, Int32)

// DECLARE_NAME_CARMA(u64, Uint64)
// DECLARE_NAME_CARMA(i64, Int64)

DECLARE_NAME_CARMA(int, int);
DECLARE_NAME_CARMA(unsigned int, uint);
DECLARE_NAME_CARMA(uint16_t, uint16);

#ifdef CAN_DO_HALF
DECLARE_NAME_CARMA(half, half);
#endif
DECLARE_NAME_CARMA(float, float);
DECLARE_NAME_CARMA(double, double);

DECLARE_NAME_CARMA(cuFloatComplex, float_complex);
DECLARE_NAME_CARMA(cuDoubleComplex, double_complex);
// DECLARE_NAME_CARMA(tuple_t<float>, tuple_float);

template <typename T>
std::string appendName(std::string str) {
  return str + explicit_name<T>();
}
#endif
