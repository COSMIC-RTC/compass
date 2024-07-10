// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_fft.cpp
//! \ingroup   libcarma
//! \class     CarmaFFT
//! \brief     this class provides the fft features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_obj.hpp>

/** These templates are used to select the proper cufft type
 from the T_in and T_out types. */

template <class T_in, class T_out>
cufftType carma_select_plan();
// /**< Generic template for cufft type selection */
// {
//   return EXIT_FAILURE;
// }

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
                        int32_t dir)
/**< Generic template for cufft executable selection */
{
  return CUFFT_INVALID_VALUE;
}

template <>
cufftResult fft_compute<cuFloatComplex, cuFloatComplex>(
    cufftHandle plan, cuFloatComplex *d_input, cuFloatComplex *d_output,
    int32_t dir) {
  return cufftExecC2C(plan, d_input, d_output, dir);
}

template <>
cufftResult fft_compute<cufftReal, cuFloatComplex>(cufftHandle plan,
                                                   cufftReal *d_input,
                                                   cuFloatComplex *d_output,
                                                   int32_t dir) {
  return cufftExecR2C(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuFloatComplex, cufftReal>(cufftHandle plan,
                                                   cuFloatComplex *d_input,
                                                   cufftReal *d_output,
                                                   int32_t dir) {
  return cufftExecC2R(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuDoubleComplex, cuDoubleComplex>(
    cufftHandle plan, cuDoubleComplex *d_input, cuDoubleComplex *d_output,
    int32_t dir) {
  return cufftExecZ2Z(plan, d_input, d_output, dir);
}

template <>
cufftResult fft_compute<cufftDoubleReal, cuDoubleComplex>(
    cufftHandle plan, cufftDoubleReal *d_input, cuDoubleComplex *d_output,
    int32_t dir) {
  return cufftExecD2Z(plan, d_input, d_output);
}

template <>
cufftResult fft_compute<cuDoubleComplex, cufftDoubleReal>(
    cufftHandle plan, cuDoubleComplex *d_input, cufftDoubleReal *d_output,
    int32_t dir) {
  return cufftExecZ2D(plan, d_input, d_output);
}

/** This is the CarmaFFT definition. */

template <class T_in, class T_out>
void carma_initfft(const int64_t *dims_data, cufftHandle *plan, cufftType type_plan) {
  /** \brief CarmaFFT creator.
   * \param dims_data : the array size
   * \param size_data : =1 : 2D array, >1 : 3D array
   * \param inplace   : flag to select inplace transform
   */

  type_plan = carma_select_plan<T_in, T_out>();
  if (dims_data[0] == 1) /* Create a 1D FFT plan. */
    carmafft_safe_call(cufftPlan1d(plan, dims_data[1], type_plan, 1));
  else if (dims_data[0] == 2)
    /* Create a 2D FFT plan. */
    carmafft_safe_call(cufftPlan2d(plan, dims_data[1], dims_data[2], type_plan));
  else
  /* Create a 3D FFT plan. */ {
    int32_t mdims[2];
    mdims[0] = (int32_t)dims_data[1];
    mdims[1] = (int32_t)dims_data[2];
    carmafft_safe_call(
        /*cufftPlan3d(plan,dims_data[1],dims_data[2],dims_data[3], type_plan));*/
        // cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C
        // ,(int32_t)dims_data[3]));
        cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0, type_plan,
                      (int32_t)dims_data[3]));
  }
}
template void carma_initfft<cuFloatComplex, cufftReal>(const int64_t *dims_data,
                                                       cufftHandle *plan,
                                                       cufftType type_plan);
template void carma_initfft<cufftReal, cuFloatComplex>(const int64_t *dims_data,
                                                       cufftHandle *plan,
                                                       cufftType type_plan);
template void carma_initfft<cuFloatComplex, cuFloatComplex>(
    const int64_t *dims_data, cufftHandle *plan, cufftType type_plan);
template void carma_initfft<cuDoubleComplex, cufftDoubleReal>(
    const int64_t *dims_data, cufftHandle *plan, cufftType type_plan);
template void carma_initfft<cufftDoubleReal, cuDoubleComplex>(
    const int64_t *dims_data, cufftHandle *plan, cufftType type_plan);
template void carma_initfft<cuDoubleComplex, cuDoubleComplex>(
    const int64_t *dims_data, cufftHandle *plan, cufftType type_plan);

template <class T_in, class T_out>
int32_t CarmaFFT(T_in *input, T_out *output, int32_t dir, cufftHandle plan) {
  /** \brief FFT computation.
   * \param input  : input array on the device
   * \param output : output array on the device
   * \param dir    : FFT direction
   *
   * this method computes the FFT of the input array in the specified direction
   * and stores the result in the output array
   */

  // CUFFT_FORWARD = -1 and CUFFT_INVERSE = 1 (cf cufft.h)
  // carmafft_safe_call( fft_compute(plan, (T_in*)input, (T_out*)output, dir *
  // CUFFT_FORWARD));
  carmafft_safe_call(
      fft_compute(plan, (T_in *)input, (T_out *)output, dir * CUFFT_FORWARD));
  return EXIT_SUCCESS;
}
template int32_t CarmaFFT<cuFloatComplex, cufftReal>(cuFloatComplex *input,
                                                  cufftReal *output, int32_t dir,
                                                  cufftHandle plan);
template int32_t CarmaFFT<cufftReal, cuFloatComplex>(cufftReal *input,
                                                  cuFloatComplex *output,
                                                  int32_t dir, cufftHandle plan);
template int32_t CarmaFFT<cuFloatComplex, cuFloatComplex>(cuFloatComplex *input,
                                                       cuFloatComplex *output,
                                                       int32_t dir,
                                                       cufftHandle plan);
template int32_t CarmaFFT<cufftDoubleReal, cufftDoubleReal>(
    cufftDoubleReal *input, cufftDoubleReal *output, int32_t dir, cufftHandle plan);
template int32_t CarmaFFT<cufftDoubleReal, cuDoubleComplex>(
    cufftDoubleReal *input, cuDoubleComplex *output, int32_t dir, cufftHandle plan);
template int32_t CarmaFFT<cuDoubleComplex, cufftDoubleReal>(
    cuDoubleComplex *input, cufftDoubleReal *output, int32_t dir, cufftHandle plan);
template int32_t CarmaFFT<cuDoubleComplex, cuDoubleComplex>(
    cuDoubleComplex *input, cuDoubleComplex *output, int32_t dir, cufftHandle plan);
