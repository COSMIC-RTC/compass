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

//! \file      carma_fft.h
//! \ingroup   libcarma
//! \class     carma_fft
//! \brief     this class provides the fft features to carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef _CARMA_FFT_H_
#define _CARMA_FFT_H_

#include <carma_obj.h>
#include <cufft.h>

template <class T_in, class T_out>
class carma_fft {
 protected:
  carma_obj<T_in> *d_input;    ///< Input data
  carma_obj<T_out> *d_output;  ///< Output data
  cufftHandle plan;            ///< FFT plan
  cufftType tPlan;             ///< FFT plan type
  int inplace;  ///< flag to select inplace transform or not (1 or 0)

 public:
  carma_fft(long *dims_data, int inplace);
  ~carma_fft();

  int host2device(T_in *data);
  int device2host(T_out *data);
  /**< Memory transfers both ways */

  int compute(int dir);
  /**< Compute on the class members */
  int compute(T_in *input, T_out *output, int dir);
  /**< Compute on any array */
};

typedef carma_fft<cuFloatComplex, cuFloatComplex> caFFT_C2C;
typedef carma_fft<cufftReal, cuFloatComplex> caFFT_R2C;
typedef carma_fft<cuFloatComplex, cufftReal> caFFT_C2R;

typedef carma_fft<cuDoubleComplex, cuDoubleComplex> caFFT_Z2Z;
typedef carma_fft<cufftDoubleReal, cuDoubleComplex> caFFT_D2Z;
typedef carma_fft<cuDoubleComplex, cufftDoubleReal> caFFT_Z2D;

extern "C" {

/** This is a collection of wrappers for yorick. */
int _fftCUInitC2C(caFFT_C2C **handle, long *dims_data, int inplace);
int _fftCUInitZ2Z(caFFT_Z2Z **handle, long *dims_data, int inplace);
int _fftCUInitR2C(caFFT_R2C **handle, long *dims_data);
int _fftCUInitD2Z(caFFT_D2Z **handle, long *dims_data);
int _fftCUInitC2R(caFFT_C2R **handle, long *dims_data);
int _fftCUInitZ2D(caFFT_Z2D **handle, long *dims_data);
/**< Init wrappers */

int _fftCUFreeC2C(caFFT_C2C **handle);
int _fftCUFreeZ2Z(caFFT_Z2Z **handle);
int _fftCUFreeR2C(caFFT_R2C **handle);
int _fftCUFreeD2Z(caFFT_D2Z **handle);
int _fftCUFreeC2R(caFFT_C2R **handle);
int _fftCUFreeZ2D(caFFT_Z2D **handle);
/**< Free wrappers */

int _fftCUhost2deviceC2C(caFFT_C2C *handle, cuFloatComplex *data);
int _fftCUhost2deviceZ2Z(caFFT_Z2Z *handle, cuDoubleComplex *data);
int _fftCUhost2deviceR2C(caFFT_R2C *handle, cufftReal *data);
int _fftCUhost2deviceD2Z(caFFT_D2Z *handle, cufftDoubleReal *data);
int _fftCUhost2deviceC2R(caFFT_C2R *handle, cuFloatComplex *data);
int _fftCUhost2deviceZ2D(caFFT_Z2D *handle, cuDoubleComplex *data);
/**< host2device wrappers */

int _fftCUdevice2hostC2C(caFFT_C2C *handle, cuFloatComplex *data);
int _fftCUdevice2hostZ2Z(caFFT_Z2Z *handle, cuDoubleComplex *data);
int _fftCUdevice2hostR2C(caFFT_R2C *handle, cuFloatComplex *data);
int _fftCUdevice2hostD2Z(caFFT_D2Z *handle, cuDoubleComplex *data);
int _fftCUdevice2hostC2R(caFFT_C2R *handle, cufftReal *data);
int _fftCUdevice2hostZ2D(caFFT_Z2D *handle, cufftDoubleReal *data);
/**< device2host wrappers */

int _fftCUcomputeC2C(caFFT_C2C *handle, int dir);
int _fftCUcomputeZ2Z(caFFT_Z2Z *handle, int dir);
int _fftCUcomputeR2C(caFFT_R2C *handle, int dir);
int _fftCUcomputeD2Z(caFFT_D2Z *handle, int dir);
int _fftCUcomputeC2R(caFFT_C2R *handle, int dir);
int _fftCUcomputeZ2D(caFFT_Z2D *handle, int dir);
/**< compute wrappers */
}

#endif  // _CARMA_FFT_H_
