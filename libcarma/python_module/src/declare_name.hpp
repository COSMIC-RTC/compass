#ifndef _DECLARE_NAME_H_
#define _DECLARE_NAME_H_

#include <carma.h>


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
