// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      declare_name.hpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for libcarma
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _DECLARE_NAME_H_
#define _DECLARE_NAME_H_

#include <string>
#include <cstdint>
#include <cuComplex.h>

#ifdef CAN_DO_HALF
#include <cuda_fp16.h>
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
