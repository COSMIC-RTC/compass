// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_fft.hpp
//! \ingroup   libcarma
//! \class     CarmaFFT
//! \brief     this class provides the fft features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#ifndef _CARMA_FFT_H_
#define _CARMA_FFT_H_

#include <carma_obj.hpp>
#include <cufft.h>

template <class T_in, class T_out>
class CarmaFFT {
 protected:
  CarmaObj<T_in> *d_input;    ///< Input data
  CarmaObj<T_out> *d_output;  ///< Output data
  cufftHandle plan;            ///< FFT plan
  cufftType type_plan;             ///< FFT plan type
  int32_t inplace;  ///< flag to select inplace transform or not (1 or 0)

 public:
  CarmaFFT(int64_t *dims_data, int32_t inplace);
  ~CarmaFFT();

  int32_t host2device(T_in *data);
  int32_t device2host(T_out *data);
  /**< Memory transfers both ways */

  int32_t compute(int32_t dir);
  /**< Compute on the class members */
  int32_t compute(T_in *input, T_out *output, int32_t dir);
  /**< Compute on any array */
};

typedef CarmaFFT<cuFloatComplex, cuFloatComplex> caFFT_C2C;
typedef CarmaFFT<cufftReal, cuFloatComplex> caFFT_R2C;
typedef CarmaFFT<cuFloatComplex, cufftReal> caFFT_C2R;

typedef CarmaFFT<cuDoubleComplex, cuDoubleComplex> caFFT_Z2Z;
typedef CarmaFFT<cufftDoubleReal, cuDoubleComplex> caFFT_D2Z;
typedef CarmaFFT<cuDoubleComplex, cufftDoubleReal> caFFT_Z2D;

#endif  // _CARMA_FFT_H_
