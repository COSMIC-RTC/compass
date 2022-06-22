// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller.cu
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <math_constants.h>
#include <sutra_controller.h>
/*
  _  __                    _
 | |/ /___ _ __ _ __   ___| |___
 | ' // _ \ '__| '_ \ / _ \ / __|
 | . \  __/ |  | | | |  __/ \__ \
 |_|\_\___|_|  |_| |_|\___|_|___/

 */

__global__ void shift_krnl(float *data, int offset, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = data[tid + offset * N];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void TT_filt_krnl(float *mat, int n, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    tid % (n + 1) ? mat[tid] *= -1.0f : mat[tid] = (1.0f - mat[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fill_filtmat_krnl(float *filter, int nactu, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    filter[tid] =
        tid % (nactu + 1) ? (float)-1. / nactu : (float)(1. - 1. / nactu);
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fill_cmat_krnl(float *cmat, float *wtt, float *Mtt, long nact,
                               long nslope, long N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j;
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
                                float norm, long dim, long N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j;
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
__global__ void pupphase_krnl(T *o_data, float *i_data, int *indx_pup, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
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

__global__ void compute_Hcor_krnl(float *o_data, int nrow, int ncol, float Fs,
                                  float Te, float gmin, float gmax,
                                  float delay) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int j = tid / nrow;
  int i = tid - j * nrow;

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

__global__ void absnormfft_krnl(cuFloatComplex *idata, float *odata, int N,
                                float norm) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex cache;
  while (tid < N) {
    cache = idata[tid + 1];  // Reject first element (0 frequency)
    odata[tid] = (cache.x * cache.x + cache.y * cache.y) * norm;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void adjust_csrcol_krnl(int *colind, int *nnz, int Nphi,
                                   int nnztot) {
  int tid = nnz[0] + threadIdx.x + blockIdx.x * blockDim.x;
  int i = 1;
  int N = nnz[0] + nnz[i];

  if (tid < nnztot) {
    while (tid > N) {
      i++;
      N += nnz[i];
    }
    __syncthreads();

    colind[tid] += i * Nphi;
  }
}

__global__ void adjust_csrrow_krnl(int *rowind, int *nact, int *nnz,
                                   int nact_tot) {
  int tid = nact[0] + threadIdx.x + blockIdx.x * blockDim.x;
  int i = 1;
  int N = nact[0] + nact[1];

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

int shift_buf(float *d_data, int offset, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  shift_krnl<<<grid, threads>>>(d_data, offset, N);

  carma_check_msg("shift_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int fill_filtmat(float *filter, int nactu, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_filtmat_krnl<<<grid, threads>>>(filter, nactu, N);
  carma_check_msg("fill_filtmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
int TT_filt(float *mat, int n, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  int N = n * n;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  TT_filt_krnl<<<grid, threads>>>(mat, n, N);
  carma_check_msg("TT_filt_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes,
              CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  int N = nactu * nslopes;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_cmat_krnl<<<grid, threads>>>(cmat, wtt, Mtt, nactu, nslopes, N);
  carma_check_msg("fill_cmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int do_statmat(float *statcov, long dim, float *xpos, float *ypos, float norm,
               CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  int N = (dim * dim);
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  do_statcov_krnl<<<grid, threads>>>(statcov, xpos, ypos, norm, dim, N);
  carma_check_msg("do_statcov_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template <class T>
int get_pupphase(T *o_data, float *i_data, int *indx_pup, int Nphi,
                 CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Nphi, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  pupphase_krnl<<<grid, threads>>>(o_data, i_data, indx_pup, Nphi);
  carma_check_msg("pupphase_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int get_pupphase<float>(float *o_data, float *i_data, int *indx_pup,
                                 int Nphi, CarmaDevice *device);
template int get_pupphase<double>(double *o_data, float *i_data, int *indx_pup,
                                  int Nphi, CarmaDevice *device);

int compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin,
                     float gmax, float delay, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nrow * ncol, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  compute_Hcor_krnl<<<grid, threads>>>(o_data, nrow, ncol, Fs, 1.0f / Fs, gmin,
                                       gmax, delay);
  carma_check_msg("compute_Hcor_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int absnormfft(cuFloatComplex *idata, float *odata, int N, float norm,
               CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  absnormfft_krnl<<<grid, threads>>>(idata, odata, N, norm);
  carma_check_msg("absnormfft_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int adjust_csr_index(int *rowind, int *NNZ, int *nact, int nact_tot,
                     int row_off, CarmaDevice *device) {
  int N = nact_tot - row_off;
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid2(nb_blocks), threads2(nb_threads);

  adjust_csrrow_krnl<<<grid2, threads2>>>(rowind, nact, NNZ, nact_tot);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
__global__ void convertVoltage_krnl(half *d_idata, float *d_odata, int N,
                                    float volt_min, float volt_max, uint16_t val_max) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    d_odata[tid] = __half2float(d_idata[tid]);
    tid += blockDim.x * gridDim.x;
  }
}
#endif

template <typename T>
__global__ void convertVoltage_krnl(T *d_idata, uint16_t *d_odata, int N,
                                    float volt_min, float volt_max, uint16_t val_max) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    d_odata[tid] =
        uint16_t((float(d_idata[tid]) - volt_min) / (volt_max - volt_min) * float(val_max));
    tid += blockDim.x * gridDim.x;
  }
}

template <typename Tin, typename Tout, std::enable_if_t<!std::is_same<Tin, Tout>::value, bool>>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  convertVoltage_krnl<<<grid, threads, 0, stream>>>(d_idata, d_odata, N, volt_min, volt_max,
                                         val_max);
  carma_check_msg("convertVoltage_krnl<<<>>> execution failed\n");
}

template void
    convert_to_voltage<float, uint16_t>(float *d_idata, uint16_t *d_odata, int N,
                                      float volt_min, float volt_max, uint16_t val_max,
                                      CarmaDevice *device, cudaStream_t stream);

#ifdef CAN_DO_HALF
template void
convert_to_voltage<half, float, std::enable_if<!std::is_same<half, float>::value>>(half *d_idata, float *d_odata, int N, float volt_min,
                              float volt_max, uint16_t val_max,
                              CarmaDevice *device, cudaStream_t stream);
template void
    convert_to_voltage<half, uint16_t, std::enable_if<!std::is_same<half, uint16_t>::value>>(half *d_idata, uint16_t *d_odata, int N,
                                     float volt_min, float volt_max, uint16_t val_max,
                                     CarmaDevice *device, cudaStream_t stream);
#endif

template <typename T>
__global__ void padCmat_krnl(T *idata, int m, int n, T *odata, int m2, int n2) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < m * n) {
    int j = tid / m;
    int i = tid - j * m;
    int tid2 = i + j * m2;
    odata[tid2] = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pad_cmat(T *idata, int m, int n, T *odata, int m2, int n2,
              CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, m * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  padCmat_krnl<<<grid, threads>>>(idata, m, n, odata, m2, n2);
  carma_check_msg("padCmat_krnl<<<>>> execution failed\n");
}

#ifdef CAN_DO_HALF
template void pad_cmat<half>(half *idata, int m, int n, half *odata, int m2,
                             int n2, CarmaDevice *device);
#endif
