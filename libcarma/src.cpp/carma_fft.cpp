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

//! \file      carma_fft.cpp
//! \ingroup   libcarma
//! \class     carma_fft
//! \brief     this class provides the fft features to carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_obj.h>

/** These templates are used to select the proper cufft type
 from the T_in and T_out types. */

template <class T_in, class T_out>
cufftType carma_select_plan()
/**< Generic template for cufft type selection */
{
  return EXIT_FAILURE;
}

template <>
cufftType carma_select_plan<cuFloatComplex, cuFloatComplex>() {
  return CUFFT_C2C;
}

template <>
cufftType carma_select_plan<cufftReal, cuFloatComplex>() {
  return CUFFT_R2C;
}

template <>
cufftType carma_select_plan<cuFloatComplex, cufftReal>() {
  return CUFFT_C2R;
}

template <>
cufftType carma_select_plan<cuDoubleComplex, cuDoubleComplex>() {
  return CUFFT_Z2Z;
}

template <>
cufftType carma_select_plan<cufftDoubleReal, cuDoubleComplex>() {
  return CUFFT_D2Z;
}

template <>
cufftType carma_select_plan<cuDoubleComplex, cufftDoubleReal>() {
  return CUFFT_Z2D;
}

/** These templates are used to select the proper cufft executable
 from the T_in and T_out types. */

template <class T_in, class T_out>
cufftResult fft_compute(cufftHandle plan, T_in *d_input, T_out *d_output,
                        int dir)
/**< Generic template for cufft executable selection */
{
  return CUFFT_INVALID_VALUE;
}

template <>
cufftResult fft_compute<cuFloatComplex, cuFloatComplex>(
    cufftHandle plan, cuFloatComplex *d_input, cuFloatComplex *d_output,
    int dir) {
  return cufftExecC2C(plan, d_input, d_output, dir);
}

template <>
cufftResult fft_compute<cufftReal, cuFloatComplex>(cufftHandle plan,
                                                   cufftReal *d_input,
                                                   cuFloatComplex *d_output,
                                                   int dir) {
  return cufftExecR2C(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuFloatComplex, cufftReal>(cufftHandle plan,
                                                   cuFloatComplex *d_input,
                                                   cufftReal *d_output,
                                                   int dir) {
  return cufftExecC2R(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuDoubleComplex, cuDoubleComplex>(
    cufftHandle plan, cuDoubleComplex *d_input, cuDoubleComplex *d_output,
    int dir) {
  return cufftExecZ2Z(plan, d_input, d_output, dir);
}

template <>
cufftResult fft_compute<cufftDoubleReal, cuDoubleComplex>(
    cufftHandle plan, cufftDoubleReal *d_input, cuDoubleComplex *d_output,
    int dir) {
  return cufftExecD2Z(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuDoubleComplex, cufftDoubleReal>(
    cufftHandle plan, cuDoubleComplex *d_input, cufftDoubleReal *d_output,
    int dir) {
  return cufftExecZ2D(plan, d_input, d_output);
}

/** This is the carma_fft definition. */

template <class T_in, class T_out>
void carma_initfft(const long *dims_data, cufftHandle *plan, cufftType tPlan) {
  /** \brief carma_fft creator.
   * \param dims_data : the array size
   * \param size_data : =1 : 2D array, >1 : 3D array
   * \param inplace   : flag to select inplace transform
   */

  tPlan = carma_select_plan<T_in, T_out>();
  if (tPlan == -1) {
    DEBUG_TRACE("Wrong data type");
    throw "Wrong data type\n";
  }
  if (dims_data[0] == 1) /* Create a 1D FFT plan. */
    carmafftSafeCall(cufftPlan1d(plan, dims_data[1], tPlan, 1));
  else if (dims_data[0] == 2)
    /* Create a 2D FFT plan. */
    carmafftSafeCall(cufftPlan2d(plan, dims_data[1], dims_data[2], tPlan));
  else
  /* Create a 3D FFT plan. */ {
    int mdims[2];
    mdims[0] = (int)dims_data[1];
    mdims[1] = (int)dims_data[2];
    carmafftSafeCall(
        /*cufftPlan3d(plan,dims_data[1],dims_data[2],dims_data[3], tPlan));*/
        // cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C
        // ,(int)dims_data[3]));
        cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0, tPlan,
                      (int)dims_data[3]));
  }
}
template void carma_initfft<cuFloatComplex, cufftReal>(const long *dims_data,
                                                       cufftHandle *plan,
                                                       cufftType tPlan);
template void carma_initfft<cufftReal, cuFloatComplex>(const long *dims_data,
                                                       cufftHandle *plan,
                                                       cufftType tPlan);
template void carma_initfft<cuFloatComplex, cuFloatComplex>(
    const long *dims_data, cufftHandle *plan, cufftType tPlan);
template void carma_initfft<cuDoubleComplex, cufftDoubleReal>(
    const long *dims_data, cufftHandle *plan, cufftType tPlan);
template void carma_initfft<cufftDoubleReal, cuDoubleComplex>(
    const long *dims_data, cufftHandle *plan, cufftType tPlan);
template void carma_initfft<cuDoubleComplex, cuDoubleComplex>(
    const long *dims_data, cufftHandle *plan, cufftType tPlan);

template <class T_in, class T_out>
int carma_fft(T_in *input, T_out *output, int dir, cufftHandle plan) {
  /** \brief FFT computation.
   * \param input  : input array on the device
   * \param output : output array on the device
   * \param dir    : FFT direction
   *
   * this method computes the FFT of the input array in the specified direction
   * and stores the result in the output array
   */

  // CUFFT_FORWARD = -1 and CUFFT_INVERSE = 1 (cf cufft.h)
  // carmafftSafeCall( fft_compute(plan, (T_in*)input, (T_out*)output, dir *
  // CUFFT_FORWARD));
  carmafftSafeCall(
      fft_compute(plan, (T_in *)input, (T_out *)output, dir * CUFFT_FORWARD));
  return EXIT_SUCCESS;
}
template int carma_fft<cuFloatComplex, cufftReal>(cuFloatComplex *input,
                                                  cufftReal *output, int dir,
                                                  cufftHandle plan);
template int carma_fft<cufftReal, cuFloatComplex>(cufftReal *input,
                                                  cuFloatComplex *output,
                                                  int dir, cufftHandle plan);
template int carma_fft<cuFloatComplex, cuFloatComplex>(cuFloatComplex *input,
                                                       cuFloatComplex *output,
                                                       int dir,
                                                       cufftHandle plan);
template int carma_fft<cufftDoubleReal, cufftDoubleReal>(
    cufftDoubleReal *input, cufftDoubleReal *output, int dir, cufftHandle plan);
template int carma_fft<cufftDoubleReal, cuDoubleComplex>(
    cufftDoubleReal *input, cuDoubleComplex *output, int dir, cufftHandle plan);
template int carma_fft<cuDoubleComplex, cufftDoubleReal>(
    cuDoubleComplex *input, cufftDoubleReal *output, int dir, cufftHandle plan);
template int carma_fft<cuDoubleComplex, cuDoubleComplex>(
    cuDoubleComplex *input, cuDoubleComplex *output, int dir, cufftHandle plan);
