#include <carma_obj.h>
#include <carma_utils.h>

/** These templates are used to select the proper cufft type
from the T_in and T_out types. */

template<class T_in, class T_out> cufftType carma_select_plan()
/**< Generic template for cufft type selection */
{
  return EXIT_FAILURE;
}

template<> cufftType carma_select_plan<cuFloatComplex, cuFloatComplex>()
{
  return CUFFT_C2C;
}

template<> cufftType carma_select_plan<cufftReal, cuFloatComplex>()
{
  return CUFFT_R2C;
}

template<> cufftType carma_select_plan<cuFloatComplex, cufftReal>()
{
  return CUFFT_C2R;
}

template<> cufftType carma_select_plan<cuDoubleComplex, cuDoubleComplex>()
{
  return CUFFT_Z2Z;
}

template<> cufftType carma_select_plan<cufftDoubleReal, cuDoubleComplex>()
{
  return CUFFT_D2Z;
}

template<> cufftType carma_select_plan<cuDoubleComplex, cufftDoubleReal>()
{
  return CUFFT_Z2D;
}

/** These templates are used to select the proper cufft executable
from the T_in and T_out types. */

template<class T_in, class T_out> cufftResult fft_compute(cufftHandle plan, T_in *d_input,T_out *d_output, int dir)
/**< Generic template for cufft executable selection */
{
  return CUFFT_INVALID_VALUE;
}

template<> cufftResult 
fft_compute<cuFloatComplex, cuFloatComplex>(cufftHandle plan, cuFloatComplex *d_input,cuFloatComplex *d_output, int dir)
{
  return cufftExecC2C(plan, d_input, d_output, dir);
}

template<> cufftResult 
fft_compute<cufftReal, cuFloatComplex>(cufftHandle plan, cufftReal *d_input,cuFloatComplex *d_output, int dir)
{
  return cufftExecR2C(plan, d_input, d_output);
}

template<> cufftResult 
fft_compute<cuFloatComplex, cufftReal>(cufftHandle plan, cuFloatComplex *d_input,cufftReal *d_output, int dir)
{
  return cufftExecC2R(plan, d_input, d_output);
}

template<> cufftResult 
fft_compute<cuDoubleComplex, cuDoubleComplex>(cufftHandle plan, cuDoubleComplex *d_input,cuDoubleComplex *d_output, int dir)
{
  return cufftExecZ2Z(plan, d_input, d_output, dir);
}

template<> cufftResult 
fft_compute<cufftDoubleReal, cuDoubleComplex>(cufftHandle plan, cufftDoubleReal *d_input,cuDoubleComplex *d_output, int dir)
{
  return cufftExecD2Z(plan, d_input, d_output);
}

template<> cufftResult 
fft_compute<cuDoubleComplex, cufftDoubleReal>(cufftHandle plan, cuDoubleComplex *d_input, cufftDoubleReal *d_output, int dir)
{
  return cufftExecZ2D(plan, d_input, d_output);
}

/** This is the carma_fft definition. */

template<class T_in, class T_out> void carma_initfft(long *dims_data,cufftHandle *plan,cufftType tPlan){
  /** \brief carma_fft creator.
   * \param dims_data : the array size
   * \param size_data : =1 : 2D array, >1 : 3D array
   * \param inplace   : flag to select inplace transform
   */

  tPlan = carma_select_plan<T_in, T_out>();
  if(tPlan == -1) {
    fprintf(stderr, "Wrong data type\n");
    throw "Wrong data type\n";
  }
  
  if(dims_data[0]==2)
    /* Create a 2D FFT plan. */ 
    cufftSafeCall( cufftPlan2d(plan, dims_data[1], dims_data[2], tPlan));
  else
    /* Create a 3D FFT plan. */ {
    int mdims[2];
    mdims[0] = (int)dims_data[1];
    mdims[1] = (int)dims_data[2];
    cufftSafeCall( cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,(int)dims_data[3]));
  }
  
}
template void carma_initfft<float2, float>(long *dims_data,cufftHandle *plan,cufftType tPlan);
template void carma_initfft<float, float2>(long *dims_data,cufftHandle *plan,cufftType tPlan);
template void carma_initfft<float2, float2>(long *dims_data,cufftHandle *plan,cufftType tPlan);
template void carma_initfft<double2, double>(long *dims_data,cufftHandle *plan,cufftType tPlan);
template void carma_initfft<double, double2>(long *dims_data,cufftHandle *plan,cufftType tPlan);
template void carma_initfft<double2, double2>(long *dims_data,cufftHandle *plan,cufftType tPlan);


template<class T_in, class T_out> int carma_fft(T_in *input, T_out *output, int dir,cufftHandle plan){
  /** \brief FFT computation.
   * \param input  : input array on the device
   * \param output : output array on the device
   * \param dir    : FFT direction
   *
   * this method computes the FFT of the input array in the specified direction
   * and stores the result in the output array
   */

  //CUFFT_FORWARD = -1 and CUFFT_INVERSE = 1 (cf cufft.h)
  //cufftSafeCall( fft_compute(plan, (T_in*)input, (T_out*)output, dir * CUFFT_FORWARD));
  cufftSafeCall(fft_compute(plan, (T_in*)input, (T_out*)output, dir * CUFFT_FORWARD));
  return EXIT_SUCCESS;
}
template int carma_fft<float2, float>(float2 *input, float *output, int dir,cufftHandle plan);
template int carma_fft<float, float2>(float *input, float2 *output, int dir,cufftHandle plan);
template int carma_fft<float2, float2>(float2 *input, float2 *output, int dir,cufftHandle plan);
template int carma_fft<double, double>(double *input, double *output, int dir,cufftHandle plan);
template int carma_fft<double, double2>(double *input, double2 *output, int dir,cufftHandle plan);
template int carma_fft<double2, double>(double2 *input, double *output, int dir,cufftHandle plan);
template int carma_fft<double2, double2>(double2 *input, double2 *output, int dir,cufftHandle plan);

/*
template int carma_fft<cuDoubleComplex, double>(cuDoubleComplex *input, double *output, int dir,cufftHandle plan);
template int carma_fft<cuDoubleComplex,cuDoubleComplex >(cuDoubleComplex *input, cuDoubleComplex *output, int dir,cufftHandle plan);
template int carma_fft<double,cuDoubleComplex >(double *input, cuDoubleComplex *output, int dir,cufftHandle plan);
template int carma_fft<cuFloatComplex, float>(cuFloatComplex *input, float *output, int dir,cufftHandle plan);
template int carma_fft<cuFloatComplex,cuFloatComplex >(cuFloatComplex *input, cuFloatComplex *output, int dir,cufftHandle plan);
template int carma_fft<float,cuFloatComplex >(float *input, cuFloatComplex *output, int dir,cufftHandle plan);

 */
