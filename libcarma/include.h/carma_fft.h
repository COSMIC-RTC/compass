// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_fft.h
//! \ingroup   libcarma
//! \class     CarmaFFT
//! \brief     this class provides the fft features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License


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
