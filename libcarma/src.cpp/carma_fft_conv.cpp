#include <carma_obj.h>
#include <string>
#include "convolutionFFT2D_common.h"

/*
 ____                      _           _____ _____ _____
 / ___|___  _ ____   _____ | |_   _____|  ___|  ___|_   _|
 | |   / _ \| '_ \ \ / / _ \| \ \ / / _ \ |_  | |_    | |
 | |__| (_) | | | \ V / (_) | |\ V /  __/  _| |  _|   | |
 \____\___/|_| |_|\_/ \___/|_| \_/ \___|_|   |_|     |_|

 */

int32_t carma_initfftconv(CarmaObjS *data_in, CarmaObjS *kernel_in, CarmaObjS *padded_data,
                      CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX) {
  float *odata = padded_data->get_o_data();
  cufftHandle *plan = padded_data->get_plan();  ///<  FFT plan

  if (odata == NULL)
    carma_safe_call(cudaMalloc(reinterpret_cast<void **>(&odata),
                             sizeof(float) * padded_data->get_nb_elements()));
  carma_safe_call(cudaMemset(padded_data->get_o_data(), 0,
                           padded_data->get_nb_elements() * sizeof(float)));
  if (data_in->get_dims(0) == 3) {
    int32_t mdims[2];
    mdims[0] = padded_data->get_dims(1);
    mdims[1] = padded_data->get_dims(2);
    carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_R2C, padded_data->get_dims(3)));
  } else {
    carmafft_safe_call(cufftPlan2d(plan, padded_data->get_dims(1),
                                 padded_data->get_dims(2), CUFFT_R2C));
  }
  // plan fwd
  cuFloatComplex *odata2 = padded_spectrum->get_o_data();
  if (odata2 == NULL)
    carma_safe_call(
        cudaMalloc(reinterpret_cast<void **>(&odata2),
                   sizeof(cuFloatComplex) * padded_spectrum->get_nb_elements()));
  carma_safe_call(cudaMemset(padded_spectrum->get_o_data(), 0,
                           padded_spectrum->get_nb_elements() * sizeof(float)));

  plan = padded_spectrum->get_plan();
  if (data_in->get_dims(0) == 3) {
    int32_t mdims[2];
    mdims[0] = padded_data->get_dims(1);
    mdims[1] = padded_data->get_dims(2);
    carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_C2R, padded_spectrum->get_dims(3)));
  } else {
    carmafft_safe_call(cufftPlan2d(plan, padded_data->get_dims(1),
                                 padded_data->get_dims(2), CUFFT_C2R));
  }
  // plan inv

  if (kernel_in->get_dims(0) == 3) {
    pad_kernel_3d(padded_data->get_o_data(), kernel_in->get_data(),
                padded_data->get_dims(1), padded_data->get_dims(2),
                kernel_in->get_dims(1), kernel_in->get_dims(2), kernelY, kernelX,
                kernel_in->get_dims(3));
  } else {
    pad_kernel(padded_data->get_o_data(), kernel_in->get_data(),
              padded_data->get_dims(1), padded_data->get_dims(2),
              kernel_in->get_dims(1), kernel_in->get_dims(2), kernelY, kernelX);
  }

  if (data_in->get_dims(0) == 3) {
    pad_data_clamp_to_border_3d(padded_data->get_data(), data_in->get_data(),
                           padded_data->get_dims(1), padded_data->get_dims(2),
                           data_in->get_dims(1), data_in->get_dims(2),
                           kernel_in->get_dims(1), kernel_in->get_dims(2),
                           kernelY, kernelX, data_in->get_dims(3));
  } else {
    pad_data_clamp_to_border(
        padded_data->get_data(), data_in->get_data(), padded_data->get_dims(1),
        padded_data->get_dims(2), data_in->get_dims(1), data_in->get_dims(2),
        kernel_in->get_dims(1), kernel_in->get_dims(2), kernelY, kernelX);
  }

  return EXIT_SUCCESS;
}

int32_t carma_fftconv(CarmaObjS *data_out, CarmaObjS *padded_data,
                  CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX) {
  CarmaFFT(padded_data->get_o_data(), padded_spectrum->get_o_data(), 1,
            *padded_data->get_plan());
  // keyword dir in fft is not relevant in this case
  CarmaFFT(padded_data->get_data(), padded_spectrum->get_data(), 1,
            *padded_data->get_plan());

  int32_t nim;
  nim = data_out->get_dims(0) == 3 ? data_out->get_dims(3) : 1;

  modulate_and_normalize((fComplex *)(padded_spectrum->get_data()),
                       (fComplex *)(padded_spectrum->get_o_data()),
                       padded_data->get_dims(1), padded_data->get_dims(2), 1,
                       nim);

  CarmaFFT(padded_spectrum->get_data(), padded_data->get_data(), 1,
            *padded_spectrum->get_plan());

  int32_t N = padded_data->get_dims(2) * padded_data->get_dims(1);
  int32_t n = data_out->get_dims(1) * data_out->get_dims(2);

  fftconv_unpad(data_out->get_data(), padded_data->get_data(),
                padded_data->get_dims(2), data_out->get_dims(1),
                data_out->get_dims(2), N, n, nim);

  return EXIT_SUCCESS;
}
