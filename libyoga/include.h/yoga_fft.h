/**
 * \class yoga_fft
 *
 * \ingroup libyoga
 *
 * \brief this class provides wrappers to the cufft library
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */
#ifndef _YOGA_FFT_H_
#define _YOGA_FFT_H_

#include <cufft.h>
#include <yoga_obj.h>

using namespace std;

template<class T_in,class T_out>
class yoga_fft {

protected:
  yoga_obj<T_in> *d_input;///< Input data
  yoga_obj<T_out> *d_output;///< Output data
  cufftHandle plan;///< FFT plan
  cufftType tPlan;///< FFT plan type
  int inplace; ///< flag to select inplace transform or not (1 or 0)

 public:
  yoga_fft(long *dims_data, int inplace);
  ~yoga_fft();

  int host2device(T_in *data);
  int device2host(T_out *data);
  /**< Memory transfers both ways */

  int compute(int dir);
  /**< Compute on the class members */
  int compute(T_in *input, T_out *output, int dir);
  /**< Compute on any array */
};

typedef yoga_fft<cuFloatComplex,cuFloatComplex> yFFT_C2C;
typedef yoga_fft<cufftReal,cuFloatComplex> yFFT_R2C;
typedef yoga_fft<cuFloatComplex,cufftReal> yFFT_C2R;

typedef yoga_fft<cuDoubleComplex,cuDoubleComplex> yFFT_Z2Z;
typedef yoga_fft<cufftDoubleReal,cuDoubleComplex> yFFT_D2Z;
typedef yoga_fft<cuDoubleComplex,cufftDoubleReal> yFFT_Z2D;

extern "C" {

/** This is a collection of wrappers for yorick. */
  int _fftCUInitC2C(yFFT_C2C **handle, long *dims_data, int inplace);
  int _fftCUInitZ2Z(yFFT_Z2Z **handle, long *dims_data, int inplace);
  int _fftCUInitR2C(yFFT_R2C **handle, long *dims_data);
  int _fftCUInitD2Z(yFFT_D2Z **handle, long *dims_data);
  int _fftCUInitC2R(yFFT_C2R **handle, long *dims_data);
  int _fftCUInitZ2D(yFFT_Z2D **handle, long *dims_data);
  /**< Init wrappers */
  
  int _fftCUFreeC2C(yFFT_C2C **handle);
  int _fftCUFreeZ2Z(yFFT_Z2Z **handle);
  int _fftCUFreeR2C(yFFT_R2C **handle);
  int _fftCUFreeD2Z(yFFT_D2Z **handle);
  int _fftCUFreeC2R(yFFT_C2R **handle);
  int _fftCUFreeZ2D(yFFT_Z2D **handle);
  /**< Free wrappers */

  int _fftCUhost2deviceC2C(yFFT_C2C *handle, cuFloatComplex *data);
  int _fftCUhost2deviceZ2Z(yFFT_Z2Z *handle, cuDoubleComplex *data);
  int _fftCUhost2deviceR2C(yFFT_R2C *handle, cufftReal *data);
  int _fftCUhost2deviceD2Z(yFFT_D2Z *handle, cufftDoubleReal *data);
  int _fftCUhost2deviceC2R(yFFT_C2R *handle, cuFloatComplex *data);
  int _fftCUhost2deviceZ2D(yFFT_Z2D *handle, cuDoubleComplex *data);
  /**< host2device wrappers */

  int _fftCUdevice2hostC2C(yFFT_C2C *handle, cuFloatComplex *data);
  int _fftCUdevice2hostZ2Z(yFFT_Z2Z *handle, cuDoubleComplex *data);
  int _fftCUdevice2hostR2C(yFFT_R2C *handle, cuFloatComplex *data);
  int _fftCUdevice2hostD2Z(yFFT_D2Z *handle, cuDoubleComplex *data);
  int _fftCUdevice2hostC2R(yFFT_C2R *handle, cufftReal *data);
  int _fftCUdevice2hostZ2D(yFFT_Z2D *handle, cufftDoubleReal *data);
  /**< device2host wrappers */

  int _fftCUcomputeC2C(yFFT_C2C *handle, int dir);
  int _fftCUcomputeZ2Z(yFFT_Z2Z *handle, int dir);
  int _fftCUcomputeR2C(yFFT_R2C *handle, int dir);
  int _fftCUcomputeD2Z(yFFT_D2Z *handle, int dir);
  int _fftCUcomputeC2R(yFFT_C2R *handle, int dir);
  int _fftCUcomputeZ2D(yFFT_Z2D *handle, int dir);
  /**< compute wrappers */
};

#endif // _YOGA_FFT_H_

