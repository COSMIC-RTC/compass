#ifndef _DECLARE_NAME_H_
#define _DECLARE_NAME_H_

#include <carma.h>


template <typename T>
constexpr static char const* explicit_name();

#define DECLARE_NAME(Type, Name)       \
  template <>                          \
  constexpr char const* explicit_name<Type>() { \
    return #Name;                      \
  }

// DECLARE_NAME(u8, Uint8)
// DECLARE_NAME(i8, Int8)

// DECLARE_NAME(u16, Uint16)
// DECLARE_NAME(i16, Int16)

// DECLARE_NAME(u32, Uint32)
// DECLARE_NAME(i32, Int32)

// DECLARE_NAME(u64, Uint64)
// DECLARE_NAME(i64, Int64)

DECLARE_NAME(int, int);
DECLARE_NAME(unsigned int, uint);

#ifdef CAN_DO_HALF
DECLARE_NAME(half, half);
#endif
DECLARE_NAME(float, float);
DECLARE_NAME(double, double);

DECLARE_NAME(cuFloatComplex, float_complex);
DECLARE_NAME(cuDoubleComplex, double_complex);
// DECLARE_NAME(tuple_t<float>, tuple_float);

template <typename T>
std::string appendName(std::string str) {
  return str + explicit_name<T>();
}
#endif
