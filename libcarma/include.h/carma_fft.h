// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_fft.h
//! \ingroup   libcarma
//! \class     CarmaFFT
//! \brief     this class provides the fft features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24


#ifndef _CARMA_FFT_H_
#define _CARMA_FFT_H_

#include <carma_obj.h>
#include <cufft.h>

template <class T_in, class T_out>
class CarmaFFT {
 protected:
  CarmaObj<T_in> *d_input;    ///< Input data
  CarmaObj<T_out> *d_output;  ///< Output data
  cufftHandle plan;            ///< FFT plan
  cufftType type_plan;             ///< FFT plan type
  int inplace;  ///< flag to select inplace transform or not (1 or 0)

 public:
  CarmaFFT(long *dims_data, int inplace);
  ~CarmaFFT();

  int host2device(T_in *data);
  int device2host(T_out *data);
  /**< Memory transfers both ways */

  int compute(int dir);
  /**< Compute on the class members */
  int compute(T_in *input, T_out *output, int dir);
  /**< Compute on any array */
};

typedef CarmaFFT<cuFloatComplex, cuFloatComplex> caFFT_C2C;
typedef CarmaFFT<cufftReal, cuFloatComplex> caFFT_R2C;
typedef CarmaFFT<cuFloatComplex, cufftReal> caFFT_C2R;

typedef CarmaFFT<cuDoubleComplex, cuDoubleComplex> caFFT_Z2Z;
typedef CarmaFFT<cufftDoubleReal, cuDoubleComplex> caFFT_D2Z;
typedef CarmaFFT<cuDoubleComplex, cufftDoubleReal> caFFT_Z2D;

#endif  // _CARMA_FFT_H_
