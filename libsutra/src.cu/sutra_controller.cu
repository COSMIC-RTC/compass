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

//! \file      sutra_controller.cu
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <math_constants.h>
#include <sutra_controller.hpp>
/*
  _  __                    _
 | |/ /___ _ __ _ __   ___| |___
 | ' // _ \ '__| '_ \ / _ \ / __|
 | . \  __/ |  | | | |  __/ \__ \
 |_|\_\___|_|  |_| |_|\___|_|___/

 */

__global__ void shift_krnl(float *data, int32_t offset, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = data[tid + offset * N];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void TT_filt_krnl(float *mat, int32_t n, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    tid % (n + 1) ? mat[tid] *= -1.0f : mat[tid] = (1.0f - mat[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fill_filtmat_krnl(float *filter, int32_t nactu, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    filter[tid] =
        tid % (nactu + 1) ? (float)-1. / nactu : (float)(1. - 1. / nactu);
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fill_cmat_krnl(float *cmat, float *wtt, float *Mtt, int64_t nact,
                               int64_t nslope, int64_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t i, j;
  while (tid < N) {
    i = tid / nact;
    j = tid - i * nact;
    // if(j < nact-2) cmat[tid] = 1.0f;
    cmat[tid] =
        (j < nact - 2) ? wtt[j + i * (nact - 2)] : Mtt[j - (nact - 2) + i * 2];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void do_statcov_krnl(float *statcov, float *xpos, float *ypos,
                                float norm, int64_t dim, int64_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t i, j;
  while (tid < N) {
    i = tid / dim;
    j = tid - i * dim;
    statcov[i * dim + j] =
        6.88 *
        powf(sqrtf((xpos[i] - xpos[j]) * (xpos[i] - xpos[j]) +
                   (ypos[i] - ypos[j]) * (ypos[i] - ypos[j])),
             5. / 3.) *
        norm;
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
__global__ void pupphase_krnl(T *o_data, float *i_data, int32_t *indx_pup, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    o_data[tid] = (T)i_data[indx_pup[tid]];
    tid += blockDim.x * gridDim.x;
  }
}

__device__ cuComplex exp_complex(cuComplex z) {
  cuComplex res;
  float t = expf(z.x);
  sincosf(z.y, &res.y, &res.x);
  res.x *= t;
  res.y *= t;

  return res;
}

__global__ void compute_Hcor_krnl(float *o_data, int32_t nrow, int32_t ncol, float Fs,
                                  float Te, float gmin, float gmax,
                                  float delay) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t j = tid / nrow;
  int32_t i = tid - j * nrow;

  float step = (Fs / 2.0f - (Fs / (2 * ncol))) / (ncol - 1);
  float f = Fs / (2 * ncol) + j * step;
  float G = gmin + i * (gmax - gmin) / (nrow - 1);
  cuFloatComplex pTe = make_cuFloatComplex(0.0f, 2 * CUDART_PI_F * f * Te);
  cuFloatComplex moins_pTe =
      make_cuFloatComplex(0.0f, -2 * CUDART_PI_F * f * Te);
  cuFloatComplex pTe2 = cuCmulf(pTe, pTe);
  cuFloatComplex UnMoinsepTe = make_cuFloatComplex(
      1.0f - exp_complex(moins_pTe).x, -exp_complex(moins_pTe).y);
  cuFloatComplex pdelay =
      make_cuFloatComplex(0.0f, -2 * CUDART_PI_F * f * Te * delay);

  cuFloatComplex res = cuCdivf(UnMoinsepTe, pTe2);
  cuFloatComplex Hbo = cuCmulf(res, exp_complex(pdelay));
  Hbo.x = 1 + G * Hbo.x;
  Hbo.y *= G;
  float mod = cuCabsf(Hbo);

  o_data[tid] = 1.0f / (mod * mod);
}

__global__ void absnormfft_krnl(cuFloatComplex *idata, float *odata, int32_t N,
                                float norm) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex cache;
  while (tid < N) {
    cache = idata[tid + 1];  // Reject first element (0 frequency)
    odata[tid] = (cache.x * cache.x + cache.y * cache.y) * norm;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void adjust_csrcol_krnl(int32_t *colind, int32_t *nnz, int32_t Nphi,
                                   int32_t nnztot) {
  int32_t tid = nnz[0] + threadIdx.x + blockIdx.x * blockDim.x;
  int32_t i = 1;
  int32_t N = nnz[0] + nnz[i];

  if (tid < nnztot) {
    while (tid > N) {
      i++;
      N += nnz[i];
    }
    __syncthreads();

    colind[tid] += i * Nphi;
  }
}

__global__ void adjust_csrrow_krnl(int32_t *rowind, int32_t *nact, int32_t *nnz,
                                   int32_t nact_tot) {
  int32_t tid = nact[0] + threadIdx.x + blockIdx.x * blockDim.x;
  int32_t i = 1;
  int32_t N = nact[0] + nact[1];

  if (tid < nact_tot) {
    while (tid > N) {
      i++;
      N += nact[i];
    }
    __syncthreads();

    rowind[tid] += nnz[i - 1];
  }
}

/*
  _                           _
 | |    __ _ _   _ _ __   ___| |__   ___ _ __ ___
 | |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
 | |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
 |_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/

 */

int32_t shift_buf(float *d_data, int32_t offset, int32_t N, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  shift_krnl<<<grid, threads>>>(d_data, offset, N);

  carma_check_msg("shift_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int32_t fill_filtmat(float *filter, int32_t nactu, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_filtmat_krnl<<<grid, threads>>>(filter, nactu, N);
  carma_check_msg("fill_filtmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
int32_t TT_filt(float *mat, int32_t n, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  int32_t N = n * n;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  TT_filt_krnl<<<grid, threads>>>(mat, n, N);
  carma_check_msg("TT_filt_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int32_t fill_cmat(float *cmat, float *wtt, float *Mtt, int32_t nactu, int32_t nslopes,
              CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  int32_t N = nactu * nslopes;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_cmat_krnl<<<grid, threads>>>(cmat, wtt, Mtt, nactu, nslopes, N);
  carma_check_msg("fill_cmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int32_t do_statmat(float *statcov, int64_t dim, float *xpos, float *ypos, float norm,
               CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  int32_t N = (dim * dim);
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  do_statcov_krnl<<<grid, threads>>>(statcov, xpos, ypos, norm, dim, N);
  carma_check_msg("do_statcov_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template <class T>
int32_t get_pupphase(T *o_data, float *i_data, int32_t *indx_pup, int32_t Nphi,
                 CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Nphi, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  pupphase_krnl<<<grid, threads>>>(o_data, i_data, indx_pup, Nphi);
  carma_check_msg("pupphase_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t get_pupphase<float>(float *o_data, float *i_data, int32_t *indx_pup,
                                 int32_t Nphi, CarmaDevice *device);
template int32_t get_pupphase<double>(double *o_data, float *i_data, int32_t *indx_pup,
                                  int32_t Nphi, CarmaDevice *device);

int32_t compute_Hcor_gpu(float *o_data, int32_t nrow, int32_t ncol, float Fs, float gmin,
                     float gmax, float delay, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nrow * ncol, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  compute_Hcor_krnl<<<grid, threads>>>(o_data, nrow, ncol, Fs, 1.0f / Fs, gmin,
                                       gmax, delay);
  carma_check_msg("compute_Hcor_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int32_t absnormfft(cuFloatComplex *idata, float *odata, int32_t N, float norm,
               CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  absnormfft_krnl<<<grid, threads>>>(idata, odata, N, norm);
  carma_check_msg("absnormfft_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int32_t adjust_csr_index(int32_t *rowind, int32_t *NNZ, int32_t *nact, int32_t nact_tot,
                     int32_t row_off, CarmaDevice *device) {
  int32_t N = nact_tot - row_off;
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid2(nb_blocks), threads2(nb_threads);

  adjust_csrrow_krnl<<<grid2, threads2>>>(rowind, nact, NNZ, nact_tot);

  return EXIT_SUCCESS;
}

template <typename T>
__global__ void convertVoltage_krnl(T *d_idata, uint16_t *d_odata, int32_t N,
                                    float volt_min, float volt_max, uint16_t val_max) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    d_odata[tid] =
        uint16_t((float(d_idata[tid]) - volt_min) / (volt_max - volt_min) * float(val_max));
    tid += blockDim.x * gridDim.x;
  }
}

template <typename Tin, typename Tout, std::enable_if_t<!std::is_same<Tin, Tout>::value, bool>>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int32_t N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  convertVoltage_krnl<<<grid, threads, 0, stream>>>(d_idata, d_odata, N, volt_min, volt_max,
                                         val_max);
  carma_check_msg("convertVoltage_krnl<<<>>> execution failed\n");
}

template void
    convert_to_voltage<float, uint16_t>(float *d_idata, uint16_t *d_odata, int32_t N,
                                      float volt_min, float volt_max, uint16_t val_max,
                                      CarmaDevice *device, cudaStream_t stream);

template <typename T>
__global__ void padCmat_krnl(T *idata, int32_t m, int32_t n, T *odata, int32_t m2, int32_t n2) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < m * n) {
    int32_t j = tid / m;
    int32_t i = tid - j * m;
    int32_t tid2 = i + j * m2;
    odata[tid2] = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pad_cmat(T *idata, int32_t m, int32_t n, T *odata, int32_t m2, int32_t n2,
              CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, m * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  padCmat_krnl<<<grid, threads>>>(idata, m, n, odata, m2, n2);
  carma_check_msg("padCmat_krnl<<<>>> execution failed\n");
}
