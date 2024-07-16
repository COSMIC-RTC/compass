// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      declare_name.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for libcarma
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _DECLARE_NAME_H_
#define _DECLARE_NAME_H_

#include <pybind11/pybind11.h>
#include "pybind11/complex.h"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <cstdint>
#include <cuComplex.h>

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

DECLARE_NAME_CARMA(int32_t, int32_t);
DECLARE_NAME_CARMA(uint32_t, uint);
DECLARE_NAME_CARMA(uint16_t, uint16);
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
