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

//! \file      sutra_wfs.cu
//! \ingroup   libsutra
//! \class     SutraWfs
//! \brief     this class provides the wfs features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_utils.hpp>
#include <sutra_wfs.hpp>
#include "carma_utils.cuh"

__global__ void camplipup_krnl(cuFloatComplex *amplipup, float *phase,
                               float *offset, float *mask, float scale,
                               int32_t *istart, int32_t *jstart, int32_t *ivalid,
                               int32_t *jvalid, int32_t nphase, int32_t nphase2, int32_t npup,
                               int32_t Nfft, int32_t N, int32_t offset_phase) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x + offset_phase;

  while (tid < N) {
    int32_t nim = tid / nphase2;
    int32_t idim = tid - nim * nphase2;

    int32_t idimx = idim % nphase;  // nphase : size of the phase support in subaps
    int32_t idimy = idim / nphase;

    int32_t idphase = idimx + idimy * npup + istart[nim] + jstart[nim] * npup;

    // npup : size of the input phase screen

    int32_t idx = idimx + idimy * Nfft + nim * Nfft * Nfft;

    amplipup[idx].x =
        (cosf(phase[idphase] * scale - offset[idim])) * mask[idphase];
    amplipup[idx].y =
        (sinf(phase[idphase] * scale - offset[idim])) * mask[idphase];
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
                  float *mask, float scale, int32_t *istart, int32_t *jstart,
                  int32_t *ivalid, int32_t *jvalid, int32_t nphase, int32_t npup, int32_t Nfft,
                  int32_t Ntot, CarmaDevice *device, int32_t offset_phase = 0) {
  // here amplipup is a cube of data of size nfft x nfft x nsubap
  // phase is an array of size pupdiam x pupdiam
  // offset is an array of size pdiam x pdiam
  // mask is an array of size pupdiam x pupdiam
  // number of thread required : pdiam x pdiam x nsubap

  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  /*
     int32_t nb_threads = 0,nb_blocks = 0;
     get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);

     dim3 grid(nb_blocks), threads(nb_threads);
   */

  int32_t nphase2 = nphase * nphase;

  camplipup_krnl<<<grid, threads>>>(
      amplipup, phase, offset, mask, scale, istart, jstart, ivalid, jvalid,
      nphase, nphase2, npup, Nfft, Ntot, 0);  // offset_phase);
  carma_check_msg("camplipup_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void camplipup_krnl(cuFloatComplex *amplipup, cuFloatComplex *phase,
                               float *offset,
                               int32_t *istart, int32_t *jstart, int32_t *ivalid,
                               int32_t *jvalid, int32_t nphase, int32_t nphase2,
                               int32_t N, int32_t Nfft, int32_t Ntot) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int32_t nim = tid / nphase2;
    int32_t idim = tid - nim * nphase2;

    int32_t idimx = idim % nphase;  // nphase : size of the phase support in subaps
    int32_t idimy = idim / nphase;

    int32_t idphase = idimx + idimy * N + istart[nim] + jstart[nim] * N;

    // npup : size of the input phase screen

    int32_t idx = idimx + idimy * Nfft + nim * Nfft * Nfft;
    // Compute exp(1j*(phi - offset))
    float cos_phi = phase[idphase].x;
    float sin_phi = phase[idphase].y;
    float cos_off = cosf(offset[idim]);
    float sin_off = sinf(offset[idim]);

    // cos(phi - off) & normalize FFT
    amplipup[idx].x = (cos_phi * cos_off + sin_phi * sin_off) / N / N;
    // sin(phi - off) & normalize FFT
    amplipup[idx].y = (sin_phi * cos_off - cos_phi * sin_off) / N / N;

    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillcamplipup(cuFloatComplex *amplipup, cuFloatComplex *phase, float *offset,
                  int32_t *istart, int32_t *jstart,
                  int32_t *ivalid, int32_t *jvalid, int32_t nphase, int32_t N, int32_t Nfft,
                  int32_t Ntot, CarmaDevice *device) {

  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  /*
     int32_t nb_threads = 0,nb_blocks = 0;
     get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);

     dim3 grid(nb_blocks), threads(nb_threads);
   */

  int32_t nphase2 = nphase * nphase;

  camplipup_krnl<<<grid, threads>>>(
      amplipup, phase, offset, istart, jstart, ivalid, jvalid,
      nphase, nphase2, N, Nfft, Ntot);
  carma_check_msg("camplipup_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fsamplipup_krnl(cuFloatComplex *d_odata, float *idata,
                                float *mask, float scale,
                                int32_t Nfft, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N*N) {

    int32_t idimx = tid % N;  // nphase : size of the phase support in subaps
    int32_t idimy = tid / N;

    int32_t odx = idimx + idimy * Nfft;

    d_odata[odx].x =
        cosf(idata[tid] * scale) * mask[tid];
    d_odata[odx].y =
        sinf(idata[tid] * scale) * mask[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillfsamplipup(cuFloatComplex *d_odata, float *idata,
                    float *mask, float scale,
                    int32_t Nfft, int32_t N, CarmaDevice *device) {

  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, N * N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fsamplipup_krnl<<<grid, threads>>>(d_odata, idata, mask, scale, Nfft, N);
  carma_check_msg("fsamplipup_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void bimg_krnl(float *bimage, float *bcube, int32_t npix, int32_t npix2,
                          int32_t nsub, int32_t *ivalid, int32_t *jvalid, float alpha,
                          int32_t N) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    int32_t nim = tid / npix2;
    int32_t tidim = tid - nim * npix2;
    int32_t xim = tidim % npix;
    int32_t yim = tidim / npix;

    int32_t idbin = xim + yim * nsub + ivalid[nim] + jvalid[nim] * nsub;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];

    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillbinimg(float *bimage, float *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
               int32_t *ivalid, int32_t *jvalid, bool add, CarmaDevice *device) {
  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  float alpha;

  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  bimg_krnl<<<grid, threads>>>(bimage, bcube, npix, Npix, Nsub, ivalid, jvalid,
                               alpha, N);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

//
// __global__ void pyradd_krnl(pyradd(float *idata, float *odata, int32_t nelem) {
//  /*
//   indx is an array nrebin^2 * npix^2
//   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
//   of the subap Npix = npix x npix
//   */
//  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
//
//  while (tid < nelem) {
//	atomicAdd(odata+tid, idata[tid]);
//    tid += blockDim.x * gridDim.x;
//  }
// }
//
// int32_t pyradd(float *idata, float *odata, int32_t nelem, CarmaDevice *device) {
//  int32_t nb_threads = 0, nb_blocks = 0;
//  get_num_blocks_and_threads(device, nelem, nb_blocks, nb_threads);
//
//  dim3 grid(nb_blocks), threads(nb_threads);
//
//  pyradd_krnl<<<grid, threads>>>(idata, odata, nelem);
//
//  carma_check_msg("pyradd_krnl<<<>>> execution failed\n");
//
//  return EXIT_SUCCESS;
// }

__global__ void bimg_krnl_async(float *bimage, float *bcube, int32_t npix,
                                int32_t npix2, int32_t nsub, int32_t *ivalid, int32_t *jvalid,
                                float alpha, int32_t N, int32_t idstart) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  tid += idstart;

  while (tid < N) {
    int32_t nim = tid / npix2;
    int32_t tidim = tid - nim * npix2;
    int32_t xim = tidim % npix;
    int32_t yim = tidim / npix;

    int32_t idbin = xim + yim * nsub + ivalid[nim] + jvalid[nim] * nsub;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillbinimg_async(CarmaHostObj<float> *image_telemetry, float *bimage,
                     float *bcube, int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid,
                     int32_t *jvalid, int32_t nim, bool add, CarmaDevice *device) {
  float *hdata = image_telemetry->get_data();
  int32_t nstreams = image_telemetry->get_nb_streams();

  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
  dim3 grid(N / (nstreams * threads.x));
  float alpha;

  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int32_t i = 0; i < nstreams; i++) {
    bimg_krnl_async<<<grid, threads, 0, image_telemetry->get_cuda_stream(i)>>>(
        bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, alpha, N,
        i * N / nstreams);

    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
    cudaMemcpyAsync(&(hdata[i * nim / nstreams]), &(bimage[i * nim / nstreams]),
                    sizeof(float) * nim / nstreams, cudaMemcpyDeviceToHost,
                    image_telemetry->get_cuda_stream(i));
  }

  // cudaStreamSynchronize(image_telemetry->get_cuda_stream(nstreams-1));
  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int32_t fillbinimg_async(CarmaStreams *streams, CarmaObj<float> *bimage,
                     CarmaObj<float> *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
                     int32_t *ivalid, int32_t *jvalid, bool add, CarmaDevice *device) {
  float *g_image = bimage->get_data();
  float *g_cube = bcube->get_data();
  int32_t nstreams = streams->get_nb_streams();

  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
  dim3 grid(N / (nstreams * threads.x));

  float alpha;

  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int32_t i = 0; i < nstreams; i++) {
    bimg_krnl_async<<<grid, threads, 0, streams->get_stream(i)>>>(
        g_image, g_cube, npix, Npix, Nsub, ivalid, jvalid, alpha, N,
        i * N / nstreams);

    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
  }

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillbincube_krnl(float *bcube, cuFloatComplex *hrimage,
                                 int32_t *indxpix, int32_t Nfft, int32_t Npix, int32_t Nrebin,
                                 int32_t N) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t npix, nsubap, nrebin;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  // int32_t nim = tid / Npix;
  // int32_t tidim = tid - nim * Npix;
  // int32_t xim = tidim % 20;
  // int32_t yim = tidim / 20;
  while (tid < N) {
    nsubap = tid / Npix;
    npix = tid % Npix;

    // if (xim>=6 && xim<14 && yim>=6 && yim<14){
    for (int32_t i = 0; i < Nrebin; i++) {
      nrebin = indxpix[i + npix * Nrebin];
      bcube[tid] += hrimage[nrebin + Nfft * nsubap].x;
    }

    // }
    // else bcube[tid] = 0.;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillbincube(float *bcube, cuFloatComplex *hrimage, int32_t *indxpix, int32_t Nfft,
                int32_t Npix, int32_t Nrebin, int32_t Nsub, CarmaDevice *device) {
  int32_t N = Npix * Nsub;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  fillbincube_krnl<<<grid, threads>>>(bcube, hrimage, indxpix, Nfft, Npix,
                                      Nrebin, N);
  carma_check_msg("fillbincube_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillbincube_krnl_async(float *bcube, cuFloatComplex *hrimage,
                                       int32_t *indxpix, int32_t Nfft, int32_t Npix,
                                       int32_t Nrebin, int32_t N, int32_t idstart) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t npix, nsubap, nrebin;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  tid += idstart;

  while (tid < N) {
    nsubap = tid / Npix;
    npix = tid % Npix;

    for (int32_t i = 0; i < Nrebin; i++) {
      nrebin = indxpix[i + npix * Nrebin];
      bcube[tid] += hrimage[nrebin + Nfft * nsubap].x;
    }
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fillbincube_async(CarmaStreams *streams, float *bcube,
                      cuFloatComplex *hrimage, int32_t *indxpix, int32_t Nfft, int32_t Npix,
                      int32_t Nrebin, int32_t Nsub, CarmaDevice *device) {
  int32_t N = Npix * Nsub;
  int32_t nb_threads = 0, nb_blocks = 0;
  int32_t nstreams = streams->get_nb_streams();

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int32_t i = 0; i < nstreams; i++) {
    fillbincube_krnl_async<<<grid, threads, 0, streams->get_stream(i)>>>(
        bcube, hrimage, indxpix, Nfft, Npix, Nrebin, N, i * N / nstreams);
    carma_check_msg("fillbincubeasync_kernel<<<>>> execution failed\n");
  }
  return EXIT_SUCCESS;
}

__global__ void indexfill_krnl(cuFloatComplex *odata, cuFloatComplex *idata,
                               int32_t *indx, int32_t ntot, int32_t Ntot, int32_t N) {
  int32_t nim, idim;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / ntot;
    idim = tid - (nim * ntot);
    odata[indx[idim] + (nim * Ntot)].x = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t indexfill(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t *indx,
              int32_t nx, int32_t Nx, int32_t N, CarmaDevice *device) {
  int32_t ntot = nx * nx;
  int32_t Ntot = Nx * Nx;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  indexfill_krnl<<<grid, threads>>>(d_odata, d_idata, indx, ntot, Ntot, N);

  carma_check_msg("indexfill_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void conv_krnl(cuFloatComplex *odata, cuFloatComplex *idata, int32_t N,
                          int32_t n) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;
  int32_t nim;
  int32_t tidim;

  // int32_t nim = tid / n;
  // int32_t tidim = tid - nim * n;

  while (tid < N) {
    nim = tid / n;
    tidim = tid - nim * n;
    tmp.x = idata[tidim].x * odata[tid].x - idata[tidim].y * odata[tid].y;
    tmp.y = idata[tidim].y * odata[tid].x + idata[tidim].x * odata[tid].y;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t convolve_cube(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
                  int32_t n, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  conv_krnl<<<grid, threads>>>(d_odata, d_idata, N, n);

  carma_check_msg("conv_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template <class T>
__device__ void red_krnl(T *sdata, int32_t size, int32_t n) {
  if (!((size & (size - 1)) == 0)) {
    uint32_t s;

    if (size % 2 != 0)
      s = size / 2 + 1;
    else
      s = size / 2;
    uint32_t s_old = size;

    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        sdata[n] += sdata[n + s];
      }
      __syncthreads();
      s_old = s;
      s /= 2;

      if ((2 * s < s_old) && (s != 0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for (uint32_t s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        sdata[n] += sdata[n + s];
      }
      __syncthreads();
    }
  }
}

template <class T>
__global__ void reduce2(T *g_idata, T *g_odata, uint32_t n,
                        uint32_t nelem_thread = 1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;

  // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;

  for (int32_t cc = 0; cc < nelem_thread; cc++) {
    int32_t idim =
        tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;

    if (idim < n)
      sdata[tid] += g_idata[idim];
    else
      sdata[tid] += 0;
  }

  // sdata[tid] = (i < n) ? g_idata[i] : 0;

  __syncthreads();
  red_krnl(sdata, blockDim.x, tid);
  __syncthreads();

  //  if (tid == 0)
  //    printf("blockIdx.x %d \n", blockIdx.x);
  //  __syncthreads();
  //
  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[tid];
}

template <class T>
__global__ void reduce2(T *g_idata, T *g_odata, T thresh, uint32_t n,
                        uint32_t nelem_thread = 1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;

  // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;

  for (int32_t cc = 0; cc < nelem_thread; cc++) {
    int32_t idim =
        tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;

    if (idim < n) {
      if (g_idata[idim] > thresh)
        sdata[tid] += g_idata[idim];
      else
        sdata[tid] += 0;
    } else
      sdata[tid] += 0;
  }

  /*
     if (i < n)
     sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] : 0;
   */
  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class T>
__global__ void reduce2_new(T *g_idata, T *g_odata, T thresh, uint32_t n,
                            uint32_t nelem_thread = 1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;

  // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;

  for (int32_t cc = 0; cc < nelem_thread; cc++) {
    int32_t idim =
        tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;

    if (idim < n) {
      if (g_idata[idim] <= thresh) g_idata[idim] = 0;
      sdata[tid] += g_idata[idim];
    } else
      sdata[tid] += 0;
  }

  /*
     if (i < n)
     sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] : 0;
   */
  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class T>
__global__ void reduce2(T *g_idata, T *g_odata, T *weights, uint32_t n,
                        uint32_t nelem_thread = 1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;

  // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;

  for (int32_t cc = 0; cc < nelem_thread; cc++) {
    int32_t idim =
        tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;

    if (idim < n)
      sdata[tid] += g_idata[idim] * weights[idim];
    else
      sdata[tid] += 0;
  }

  //   sdata[tid] = (i < n) ? g_idata[i] * weights[i] : 0;

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[tid];
}

template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  CarmaDevice *device) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;

  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize = threads * sizeof(T);
  reduce2<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
                                              (uint32_t)size, nelem_thread);

  carma_check_msg("reduce2_kernel<<<>>> execution failed\n");
}

template void subap_reduce<float>(int32_t size, int32_t threads, int32_t blocks,
                                  float *d_idata, float *d_odata,
                                  CarmaDevice *device);
template void subap_reduce<double>(int32_t size, int32_t threads, int32_t blocks,
                                   double *d_idata, double *d_odata,
                                   CarmaDevice *device);

template <class T>
__global__ void reduce2_async(T *g_idata, T *g_odata, uint32_t n,
                              int32_t stream_offset) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  i += stream_offset * blockDim.x;

  // tid+=idstream*nb_blocks/nstreams;

  sdata[tid] = (i < n) ? g_idata[i] : 0;

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x + stream_offset] = sdata[0];
}

template <class T>
void subap_reduce_async(int32_t threads, int32_t blocks, CarmaStreams *streams,
                        T *d_idata, T *d_odata) {
  int32_t nstreams = streams->get_nb_streams();
  dim3 dimBlock(threads);
  dim3 dimGrid(blocks / nstreams);

  int32_t nbelem = threads * blocks;

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize = threads * sizeof(T);

  for (int32_t i = 0; i < nstreams; i++) {
    reduce2_async<T><<<dimGrid, dimBlock, smemSize, (*streams)[i]>>>(
        d_idata, d_odata, nbelem, i * blocks / nstreams);
  }

  carma_check_msg("reduce_kernel<<<>>> execution failed\n");
}

template void subap_reduce_async<float>(int32_t threads, int32_t blocks,
                                        CarmaStreams *streams, float *d_idata,
                                        float *d_odata);
template void subap_reduce_async<double>(int32_t threads, int32_t blocks,
                                         CarmaStreams *streams,
                                         double *d_idata, double *d_odata);

template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  T thresh, CarmaDevice *device) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;

  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize = threads * sizeof(T);
  reduce2<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, thresh, size,
                                              nelem_thread);
  carma_check_msg("reduce_kernel<<<>>> execution failed\n");
}

template void subap_reduce<float>(int32_t size, int32_t threads, int32_t blocks,
                                  float *d_idata, float *d_odata, float thresh,
                                  CarmaDevice *device);

template void subap_reduce<double>(int32_t size, int32_t threads, int32_t blocks,
                                   double *d_idata, double *d_odata,
                                   double thresh, CarmaDevice *device);

template <class T>
void subap_reduce_new(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                      T thresh, CarmaDevice *device) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;

  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize = threads * sizeof(T);
  reduce2_new<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, thresh,
                                                  size, nelem_thread);
  carma_check_msg("reduce_kernel<<<>>> execution failed\n");
}

template void subap_reduce_new<float>(int32_t size, int32_t threads, int32_t blocks,
                                      float *d_idata, float *d_odata,
                                      float thresh, CarmaDevice *device);

template void subap_reduce_new<double>(int32_t size, int32_t threads, int32_t blocks,
                                       double *d_idata, double *d_odata,
                                       double thresh, CarmaDevice *device);

template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  T *weights, CarmaDevice *device) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;

  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize = threads * sizeof(T);
  reduce2<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, weights, size,
                                              nelem_thread);
  carma_check_msg("reduce_kernel<<<>>> execution failed\n");
}

template void subap_reduce<float>(int32_t size, int32_t threads, int32_t blocks,
                                  float *d_idata, float *d_odata,
                                  float *weights, CarmaDevice *device);

template void subap_reduce<double>(int32_t size, int32_t threads, int32_t blocks,
                                   double *d_idata, double *d_odata,
                                   double *weights, CarmaDevice *device);

template <class T>
__global__ void reduce_phasex(T *g_idata, T *g_odata, int32_t *indx, uint32_t n,
                              T alpha) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * n * n;

  sdata[tid] = g_idata[indx[i + tid * n + n - 1]] - g_idata[indx[i + tid * n]];

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0] / n * alpha;
  }
}

template <class T>
__global__ void reduce_phasey(T *g_idata, T *g_odata, int32_t *indx, uint32_t n,
                              T alpha)

// FIXME
// full parallelization would require to use 4 x threads
// instead of doing 4 operations per threads
{
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * n * n + threadIdx.x;

  sdata[tid] = g_idata[indx[i + (n - 1) * n]] - g_idata[indx[i]];

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0] / n * alpha;
  }
}

template <class T>
__global__ void derive_phasex(T *g_idata, T *g_odata, int32_t *indx, T *mask,
                              T alpha, uint32_t n, uint32_t N,
                              float *fluxPerSub) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    if (tid % n == 0) {  // start of a column
      sdata[tid] = (g_idata[indx[i + 1]] - g_idata[indx[i]]) * mask[indx[i]];
    } else {
      if ((tid + 1) % n == 0) {  // end of a column
        sdata[tid] = (g_idata[indx[i]] - g_idata[indx[i - 1]]) * mask[indx[i]];
      } else
        sdata[tid] =
            (g_idata[indx[i + 1]] - g_idata[indx[i - 1]]) * mask[indx[i]];
    }
  }

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0] / n * alpha / fluxPerSub[blockIdx.x] / 2;
}

template <class T>
__global__ void derive_phasey(T *g_idata, T *g_odata, int32_t *indx, T *mask,
                              T alpha, uint32_t n, uint32_t N,
                              float *fluxPerSub) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    if (tid < n) {  // start of a column
      sdata[tid] = (g_idata[indx[i + n]] - g_idata[indx[i]]) * mask[indx[i]];
    } else {
      if (tid >= n * (n - 1)) {  // end of a column
        sdata[tid] = (g_idata[indx[i]] - g_idata[indx[i - n]]) * mask[indx[i]];
      } else
        sdata[tid] =
            (g_idata[indx[i + n]] - g_idata[indx[i - n]]) * mask[indx[i]];
    }
  }

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0] / n / 2 * alpha / fluxPerSub[blockIdx.x];
}

template <class T>
void phase_reduce(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t *indx,
                  T alpha) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  reduce_phasex<T>
      <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, indx, threads, alpha);

  carma_check_msg("reduce_phasex_kernel<<<>>> execution failed\n");

  reduce_phasey<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
                                                    indx, threads, alpha);

  carma_check_msg("reduce_phasey_kernel<<<>>> execution failed\n");
}

template void phase_reduce<float>(int32_t threads, int32_t blocks, float *d_idata,
                                  float *d_odata, int32_t *indx, float alpha);

template void phase_reduce<double>(int32_t threads, int32_t blocks, double *d_idata,
                                   double *d_odata, int32_t *indx, double alpha);

template <class T>
void phase_derive(int32_t size, int32_t threads, int32_t blocks, int32_t n, T *d_idata,
                  T *d_odata, int32_t *indx, T *mask, T alpha, float *fluxPerSub) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  derive_phasex<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, d_odata, indx, mask, alpha, n, size, fluxPerSub);

  carma_check_msg("phase_derivex_kernel<<<>>> execution failed\n");

  derive_phasey<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, &(d_odata[blocks]), indx, mask, alpha, n, size, fluxPerSub);

  carma_check_msg("phase_derivey_kernel<<<>>> execution failed\n");
}

template void phase_derive<float>(int32_t size, int32_t threads, int32_t blocks, int32_t n,
                                  float *d_idata, float *d_odata, int32_t *indx,
                                  float *mask, float alpha, float *fluxPerSub);

template void phase_derive<double>(int32_t size, int32_t threads, int32_t blocks, int32_t n,
                                   double *d_idata, double *d_odata, int32_t *indx,
                                   double *mask, double alpha,
                                   float *fluxPerSub);

template <class T>
__global__ void project_gather(T *g_idata, T *g_odata, int32_t *indx,
                               uint32_t N) {

  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    g_odata[i] = g_idata[indx[i]];
  }

}

template <class T>
void phase_project(int32_t nphase, int32_t nvalid, T *d_idata, T *d_odata, int32_t *indx,
                   T *d_ttprojmat, T *d_ttprojvec, CarmaDevice *device) {

  int32_t nb_blocks, nb_threads;
  int32_t size = nphase * nphase * nvalid;
  get_num_blocks_and_threads(device, size, nb_blocks, nb_threads);

  dim3 dimGridGather(nb_blocks, 1, 1);
  dim3 dimBlockGather(nb_threads, 1, 1);

  project_gather<T><<<dimGridGather,dimBlockGather>>>(
      d_idata, d_ttprojvec, indx, size);

  carma_check_msg("project_gather_kernel<<<>>> execution failed\n");

  // x-slope
  carma_gemm_strided_batched(device->get_cublas_handle(), 'n', 'n', 1, 1,
                            nphase * nphase, T(1.0f),
                            d_ttprojmat, 1, nphase * nphase,
                            d_ttprojvec, nphase * nphase, nphase * nphase,
                            T(0.0f), d_odata, 1, 1, nvalid);

  carma_check_msg("cublas_batched_gemm<<<>>> execution failed\n");

  // y-slope
  carma_gemm_strided_batched(device->get_cublas_handle(), 'n', 'n', 1, 1,
                    nphase * nphase, T(1.0f),
                    &(d_ttprojmat[nvalid*nphase*nphase]), 1, nphase * nphase,
                    d_ttprojvec, nphase * nphase, nphase * nphase,
                    T(0.0f), &(d_odata[nvalid]), 1, 1, nvalid);

  carma_check_msg("cublas_batched_gemm<<<>>> execution failed\n");

}

template void phase_project<float>(int32_t nphase, int32_t nvalid, float *d_idata,
                                  float *d_odata, int32_t *indx, float *d_ttprojmat,
                                  float *d_ttprojvec, CarmaDevice *device);

template void phase_project<double>(int32_t nphase, int32_t nvalid, double *d_idata,
                                double *d_odata, int32_t *indx, double *d_ttprojmat,
                                double *d_ttprojvec, CarmaDevice *device);

template <class Tout, class Tin>
__global__ void pyrgetpup_krnl(Tout *g_odata, Tin *g_idata, Tout *offsets,
                               Tin *pup, float lambda, uint32_t n) {
  // roll( pup * exp(i*phase) ) * offsets

  Tout *sdata = SharedMemory<Tout>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  int32_t x, y, xx, yy, i2;

  if (i < n * n / 2) {
    x = i % n;
    y = i / n;

    xx = (x + n / 2) % n;
    yy = (y + n / 2) % n;
    i2 = xx + yy * n;

    float mic2rad = (2 * CARMA_PI / lambda);

    sdata[2 * tid].x = cosf(g_idata[i] * mic2rad) * pup[i];
    sdata[2 * tid].y = sinf(g_idata[i] * mic2rad) * pup[i];

    sdata[2 * tid + 1].x = cosf(g_idata[i2] * mic2rad) * pup[i2];
    sdata[2 * tid + 1].y = sinf(g_idata[i2] * mic2rad) * pup[i2];
  }
  __syncthreads();

  if (i < n * n / 2) {
    g_odata[i].x = sdata[2 * tid + 1].x * offsets[i].x -
                   sdata[2 * tid + 1].y * offsets[i].y;
    g_odata[i].y = sdata[2 * tid + 1].x * offsets[i].y +
                   sdata[2 * tid + 1].y * offsets[i].x;

    g_odata[i2].x =
        sdata[2 * tid].x * offsets[i2].x - sdata[2 * tid].y * offsets[i2].y;
    g_odata[i2].y =
        sdata[2 * tid].x * offsets[i2].y + sdata[2 * tid].y * offsets[i2].x;
  }
}

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup,
                int32_t np, float lambda, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, np * np / 2, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  int32_t smemSize = 2 * nb_threads * sizeof(Tout);
  pyrgetpup_krnl<Tout, Tin><<<grid, threads, smemSize>>>(
      d_odata, d_idata, d_offsets, d_pup, lambda, np);

  carma_check_msg("pyrgetpup_kernel<<<>>> execution failed\n");
}

template void pyr_getpup<cuFloatComplex, float>(
    cuFloatComplex *d_odata, float *d_idata, cuFloatComplex *d_offsets,
    float *d_pup, int32_t np, float lambda, CarmaDevice *device);
template void pyr_getpup<cuDoubleComplex, double>(
    cuDoubleComplex *d_odata, double *d_idata, cuDoubleComplex *d_offsets,
    double *d_pup, int32_t np, float lambda, CarmaDevice *device);

__global__ void copyImginBinimg_krnl(float *binimg, int32_t *validsubsx,
                                     int32_t *validsubsy, int32_t Nb, float *img,
                                     int32_t *validx, int32_t *validy, int32_t Nim,
                                     int32_t Npix) {
  int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < Npix) {
    int32_t bind = validsubsx[tid] + validsubsy[tid] * Nb;
    int32_t iind = validx[tid] + validy[tid] * Nim;
    binimg[bind] = img[iind];

    tid += blockDim.x * gridDim.x;
  }
}

void copy_imgin_binimg(float *binimg, int32_t *validsubsx, int32_t *validsubsy, int32_t Nb,
                     float *img, int32_t *validx, int32_t *validy, int32_t Nim, int32_t Npix,
                     CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, Npix, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  copyImginBinimg_krnl<<<grid, threads>>>(binimg, validsubsx, validsubsy, Nb,
                                          img, validx, validy, Nim, Npix);
  carma_check_msg("copyImginBinimg_krnl<<<>>> execution failed\n");
}

/*
   ////////////////////////////////////////////////////////////////////////:
   New version of pyrgetpup for hr version of pyramid model
   ////////////////////////////////////////////////////////////////////////
 */

template <class Tout, class Tin>
__global__ void pyrgetpup_krnl(Tout *g_odata, Tin *g_idata, Tin *pup,
                               float lambda, float cx, float cy, uint32_t n,
                               uint32_t N) {
  // load shared mem
  // const uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    const int32_t x = i % n;
    const int32_t y = i / n;

    const float xx = (x - n / 2.0f);
    const float yy = (y - n / 2.0f);

    const float mic2rad = (2 * CARMA_PI / lambda);

    const float phi_modu = cx * xx + cy * yy;

    const int32_t offs = (N - n) / 2.0f;

    const int32_t i2 = x + offs + (y + offs) * N;

    g_odata[i2].x = cosf(g_idata[i] * mic2rad + phi_modu) * pup[i];
    g_odata[i2].y = sinf(g_idata[i] * mic2rad + phi_modu) * pup[i];

    i += blockDim.x * gridDim.x;
  }
}

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tin *d_pup, int32_t np, int32_t N,
                float lambda, float cx, float cy, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, np * np, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyrgetpup_krnl<Tout, Tin>
      <<<grid, threads>>>(d_odata, d_idata, d_pup, lambda, cx, cy, np, N);

  carma_check_msg("pyrgetpup_kernel<<<>>> execution failed\n");
}

template void pyr_getpup<cuFloatComplex, float>(cuFloatComplex *d_odata,
                                                float *d_idata, float *d_pup,
                                                int32_t np, int32_t N, float lambda,
                                                float cx, float cy,
                                                CarmaDevice *device);
template void pyr_getpup<cuDoubleComplex, double>(
    cuDoubleComplex *d_odata, double *d_idata, double *d_pup, int32_t np, int32_t N,
    float lambda, float cx, float cy, CarmaDevice *device);

template <class T>
__global__ void rollmod_krnl(T *g_odata, T *g_idata, T *g_mask, float cx,
                             float cy, uint32_t n, uint32_t ns,
                             uint32_t nim) {
  // roll( pup * exp(i*phase) ) * offsets

  // load shared mem
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  int32_t xx, yy, i2;

  if (i < ns * ns * nim) {
    int32_t n_im = i / (ns * ns);
    int32_t tidim = i - n_im * ns * ns;

    xx = tidim % ns;
    yy = tidim / ns;

    int32_t xx2 = 0;
    int32_t yy2 = 0;

    // cx *= 0;
    // cy *= 0;

    // here we re-order the quadrant to simulate a roll :
    // upper right is 1, upper left is 2, bottom right is 3 and bottom left is 4
    if (n_im == 0) {
      // first image : upper right
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy + (n - ns) - cy);
    }

    if (n_im == 1) {
      // second image : upper left
      xx2 = (xx - cx);
      yy2 = (yy + (n - ns) - cy);
    }

    if (n_im == 2) {
      // third image : lower right
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy - cy);
    }

    if (n_im == 3) {
      // fourth image : lower left
      xx2 = (xx - cx);
      yy2 = (yy - cy);
    }

    if (xx2 < 0)
      xx2 = xx2 + n;
    else
      xx2 = xx2 % n;

    if (yy2 < 0)
      yy2 = yy2 + n;
    else
      yy2 = yy2 % n;

    i2 = xx2 + yy2 * n;

    if (i2 < n * n) {
      T tmp1, tmp2;
      tmp1 = g_idata[i2];
      tmp2 = g_mask[tidim];
      g_odata[i].x = tmp1.x * tmp2.x - tmp1.y * tmp2.y;
      g_odata[i].y = tmp1.x * tmp2.y + tmp1.y * tmp2.x;
    }
  }
}

template <class T>
void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int32_t np,
                 int32_t ns, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns * 4, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  rollmod_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, d_mask, cx, cy, np, ns, 4);

  carma_check_msg("rollmod_kernel<<<>>> execution failed\n");
}

template void pyr_rollmod<cuFloatComplex>(cuFloatComplex *d_odata,
                                          cuFloatComplex *d_idata,
                                          cuFloatComplex *d_mask, float cx,
                                          float cy, int32_t np, int32_t ns,
                                          CarmaDevice *device);
template void pyr_rollmod<cuDoubleComplex>(cuDoubleComplex *d_odata,
                                           cuDoubleComplex *d_idata,
                                           cuDoubleComplex *d_mask, float cx,
                                           float cy, int32_t np, int32_t ns,
                                           CarmaDevice *device);

// ////////////////////////////////////////////////////////////
// ADDING PYR_ROLLMOD MODIFIED FOR ROOF-PRISM: ROOF_ROOLMOD //
// ////////////////////////////////////////////////////////////

template <class T>
__global__ void roof_rollmod_krnl(T *g_odata, T *g_idata, T *g_mask, float cx,
                                  float cy, uint32_t n, uint32_t ns,
                                  uint32_t nim) {
  // roll( pup * exp(i*phase) ) * offsets

  // load shared mem
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  int32_t xx, yy, i2;

  if (i < ns * ns * nim) {
    int32_t n_im = i / (ns * ns);
    int32_t tidim = i - n_im * ns * ns;

    xx = tidim % ns;
    yy = tidim / ns;

    int32_t xx2 = 0;
    int32_t yy2 = 0;

    // cx *= 0;
    // cy *= 0;

    // here we re-order the half-planes to simulate a roll :
    // top is 1 and bottom is 2 ; left is 3 and right is 4.
    if (n_im == 0) {
      // first image : top
      xx2 = (xx + (n - ns / 2));
      yy2 = (yy + (n - ns) - cy);
    }

    if (n_im == 1) {
      // second image : bottom
      xx2 = (xx + (n - ns / 2));
      yy2 = (yy - cy);
    }

    if (n_im == 2) {
      // third image : left
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy + (n - ns / 2));
    }

    if (n_im == 3) {
      // fourth image : right
      xx2 = (xx - cx);
      yy2 = (yy + (n - ns / 2));
    }

    if (xx2 < 0)
      xx2 = xx2 + n;
    else
      xx2 = xx2 % n;

    if (yy2 < 0)
      yy2 = yy2 + n;
    else
      yy2 = yy2 % n;

    i2 = xx2 + yy2 * n;

    if (i2 < n * n) {
      T tmp1, tmp2;
      tmp1 = g_idata[i2];
      tmp2 = g_mask[tidim];
      g_odata[i].x = tmp1.x * tmp2.x - tmp1.y * tmp2.y;
      g_odata[i].y = tmp1.x * tmp2.y + tmp1.y * tmp2.x;
    }
  }
}

template <class T>
void roof_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int32_t np,
                  int32_t ns, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns * 4, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  roof_rollmod_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, d_mask, cx, cy, np, ns, 4);

  carma_check_msg("roof_rollmod_kernel<<<>>> execution failed\n");
}

template void roof_rollmod<cuFloatComplex>(cuFloatComplex *d_odata,
                                           cuFloatComplex *d_idata,
                                           cuFloatComplex *d_mask, float cx,
                                           float cy, int32_t np, int32_t ns,
                                           CarmaDevice *device);
template void roof_rollmod<cuDoubleComplex>(cuDoubleComplex *d_odata,
                                            cuDoubleComplex *d_idata,
                                            cuDoubleComplex *d_mask, float cx,
                                            float cy, int32_t np, int32_t ns,
                                            CarmaDevice *device);

// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////

template <class T>
__global__ void fillbinpyr_krnl(T *g_odata, T *g_idata, uint32_t nrebin,
                                uint32_t n, uint32_t ns,
                                uint32_t nim) {
  // bin2d(hrimg)

  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  int32_t x, y;

  if (i < ns * ns * nim) {
    int32_t n_im = i / (ns * ns);
    int32_t tidim = i - n_im * (ns * ns);
    x = tidim % ns;
    y = tidim / ns;
    int32_t xim = x * nrebin;
    int32_t yim = y * nrebin;

    sdata[tid] = 0;

    for (int32_t cpty = 0; cpty < nrebin; cpty++) {
      for (int32_t cptx = 0; cptx < nrebin; cptx++) {
        int32_t tidim2 = (xim + cptx) + (yim + cpty) * n + n_im * n * n;
        sdata[tid] += g_idata[tidim2];
      }
    }

    //    sdata[tid] /= (nrebin * nrebin);
  }
  __syncthreads();

  if (i < ns * ns * nim) {
    g_odata[i] = sdata[tid];
  }
}

template <class T>
void pyr_fillbin(T *d_odata, T *d_idata, int32_t nrebin, int32_t np, int32_t ns, int32_t nim,
                 CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns * nim, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  int32_t smemSize = nb_threads * sizeof(T);
  fillbinpyr_krnl<T>
      <<<grid, threads, smemSize>>>(d_odata, d_idata, nrebin, np, ns, nim);

  carma_check_msg("pyrgetpup_kernel<<<>>> execution failed\n");
}

template void pyr_fillbin<float>(float *d_odata, float *d_idata, int32_t nrebin,
                                 int32_t np, int32_t ns, int32_t nim, CarmaDevice *device);
template void pyr_fillbin<double>(double *d_odata, double *d_idata, int32_t nrebin,
                                  int32_t np, int32_t ns, int32_t nim,
                                  CarmaDevice *device);

template <class T>
__global__ void pyr_bimg_krnl(T *bimage, const T *bcube, const int32_t nxsub,
                              const float alpha, const int32_t N) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int32_t nxsub2 = nxsub * nxsub;
  const int32_t nximg = (2 * nxsub + 3);

  while (tid < N) {
    const int32_t nim = tid / nxsub2;
    const int32_t tidim = tid - nim * nxsub2;
    const int32_t yim = tidim / nxsub;
    const int32_t xim = tidim - yim * nxsub;

    int32_t offset;

    switch (nim) {
      case 0:
        offset = 1;
        break;

      case 1:
        offset = 2 + nxsub;
        break;

      case 2:
        offset = 1 + (nximg * (1 + nxsub));
        break;

      default:
        offset = 1 + (nximg + 1) * (1 + nxsub);
        break;
    }
    const int32_t idbin = offset + xim + (yim + 1) * nximg;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t pyr_fillbinimg(T *bimage, const T *bcube, const int32_t nxsub, const bool add,
                   CarmaDevice *device) {
  const int32_t N = 4 * nxsub * nxsub;
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  T alpha;

  if (add)
    alpha = 1.0;
  else
    alpha = 0.0;

  pyr_bimg_krnl<<<grid, threads>>>(bimage, bcube, nxsub, alpha, N);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t pyr_fillbinimg<float>(float *bimage, const float *bcube,
                                   const int32_t nxsub, const bool add,
                                   CarmaDevice *device);
template int32_t pyr_fillbinimg<double>(double *bimage, const double *bcube,
                                    const int32_t nxsub, const bool add,
                                    CarmaDevice *device);

template <class T>
__global__ void pyr_bimg_krnl(T *oimage, const T *image, const int32_t n,
                              const int32_t rebin, const float alpha, const int32_t N) {
  /*
     indx is an array nrebin^2 * npix^2
     it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels
     of the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < n * n) {
    const int32_t x = tid % n;
    const int32_t y = tid / n;
    const int32_t xx = x * rebin;
    const int32_t yy = y * rebin;
    oimage[tid] *= alpha;

    for (int32_t ii = 0; ii < rebin; ii++) {
      for (int32_t jj = 0; jj < rebin; jj++) {
        const int32_t xim = xx + ii;
        const int32_t yim = yy + jj;
        const int32_t iim = xim + yim * N;
        oimage[tid] = oimage[tid] + image[iim];
      }
    }
    oimage[tid] /= (rebin * rebin);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t pyr_fillbinimg(T *oimage, const T *image, const int32_t n, const int32_t N,
                   const int32_t rebin, const bool add, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  T alpha;

  if (add)
    alpha = 1.0;
  else
    alpha = 0.0;

  pyr_bimg_krnl<<<grid, threads>>>(oimage, image, n, rebin, alpha, N);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t pyr_fillbinimg<float>(float *oimage, const float *image,
                                   const int32_t n, const int32_t N, const int32_t rebin,
                                   const bool add, CarmaDevice *device);
template int32_t pyr_fillbinimg<double>(double *oimage, const double *image,
                                    const int32_t n, const int32_t N, const int32_t rebin,
                                    const bool add, CarmaDevice *device);

// ////////////////////////////////////////////////////////////
// ADDING PYR_FILLBIN MODIFIED FOR ROOF-PRISM: ROOF_FILLBIN //
// ////////////////////////////////////////////////////////////

template <class T>
__global__ void fillbinroof_krnl(T *g_odata, T *g_idata, uint32_t nrebin,
                                 uint32_t n, uint32_t ns,
                                 uint32_t nim) {
  // bin2d(hrimg)

  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  int32_t x, y;

  if (i < ns * ns * nim) {
    int32_t n_im = i / (ns * ns);
    int32_t tidim = i - n_im * ns * ns;

    // / Reprendre ici...

    x = tidim % ns;
    y = tidim / ns;
    int32_t xim = x * nrebin;
    int32_t yim = y * nrebin;

    sdata[tid] = 0;

    for (int32_t cpty = 0; cpty < nrebin; cpty++) {
      for (int32_t cptx = 0; cptx < nrebin; cptx++) {
        int32_t tidim2 = (xim + cptx) + (yim + cpty) * n + n_im * n * n;
        sdata[tid] += g_idata[tidim2];
      }
    }
    sdata[tid] /= (nrebin * nrebin);
  }
  __syncthreads();

  if (i < ns * ns * nim) {
    g_odata[i] = sdata[tid];
  }
}

template <class T>
void roof_fillbin(T *d_odata, T *d_idata, int32_t nrebin, int32_t np, int32_t ns, int32_t nim,
                  CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns * nim, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  int32_t smemSize = nb_threads * sizeof(T);
  fillbinroof_krnl<T>
      <<<grid, threads, smemSize>>>(d_odata, d_idata, nrebin, np, ns, nim);

  carma_check_msg("pyrgetpup_kernel<<<>>> execution failed\n");
}

template void roof_fillbin<float>(float *d_odata, float *d_idata, int32_t nrebin,
                                  int32_t np, int32_t ns, int32_t nim,
                                  CarmaDevice *device);
template void roof_fillbin<double>(double *d_odata, double *d_idata, int32_t nrebin,
                                   int32_t np, int32_t ns, int32_t nim,
                                   CarmaDevice *device);

// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////

template <class Tout, class Tin>
__global__ void abspyr_krnl(Tout *g_odata, Tin *g_idata, uint32_t ns,
                            uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  int32_t x, y;

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  uint32_t i2 = x + y * ns;

  if (i < ns * ns) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_odata[i + cpt * ns * ns] =
          sqrtf(g_idata[i2 + cpt * ns * ns].x * g_idata[i2 + cpt * ns * ns].x +
                g_idata[i2 + cpt * ns * ns].y * g_idata[i2 + cpt * ns * ns].y);
      g_odata[i2 + cpt * ns * ns] =
          sqrtf(g_idata[i + cpt * ns * ns].x * g_idata[i + cpt * ns * ns].x +
                g_idata[i + cpt * ns * ns].y * g_idata[i + cpt * ns * ns].y);
    }
  }
}

template <class Tout, class Tin>
void pyr_abs(Tout *d_odata, Tin *d_idata, int32_t ns, int32_t nim,
             CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns / 2, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  abspyr_krnl<Tout, Tin><<<grid, threads>>>(d_odata, d_idata, ns, nim);

  carma_check_msg("abspyr_kernel<<<>>> execution failed\n");
}

template void pyr_abs<float, cuFloatComplex>(float *d_odata,
                                             cuFloatComplex *d_idata, int32_t ns,
                                             int32_t nim, CarmaDevice *device);
template void pyr_abs<double, cuDoubleComplex>(double *d_odata,
                                               cuDoubleComplex *d_idata, int32_t ns,
                                               int32_t nim, CarmaDevice *device);

template <class Tin, class Tout>
__global__ void abs2pyr_krnl(Tout *g_odata, Tin *g_idata, Tout fact,
                             uint32_t ns, uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  int32_t x, y;

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  uint32_t i2 = x + y * ns;

  // if (i2 < ns/2) i2 = i2 + ns/2 + ns/2*ns;

  Tin tmp1, tmp2;

  if (i < ns * ns / 2) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      tmp1 = g_idata[i + cpt * ns * ns];
      tmp2 = g_idata[i2 + cpt * ns * ns];

      g_odata[i + cpt * ns * ns] += (tmp2.x * tmp2.x + tmp2.y * tmp2.y) * fact;
      g_odata[i2 + cpt * ns * ns] += (tmp1.x * tmp1.x + tmp1.y * tmp1.y) * fact;
    }
  }
}

template <class Tin, class Tout>
void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int32_t ns, int32_t nim,
              CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns / 2, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  abs2pyr_krnl<Tin, Tout><<<grid, threads>>>(d_odata, d_idata, fact, ns, nim);

  carma_check_msg("abs2pyr_kernel<<<>>> execution failed\n");
}

template void pyr_abs2<cuFloatComplex, float>(float *d_odata,
                                              cuFloatComplex *d_idata,
                                              float fact, int32_t ns, int32_t nim,
                                              CarmaDevice *device);
template void pyr_abs2<cuDoubleComplex, double>(double *d_odata,
                                                cuDoubleComplex *d_idata,
                                                double fact, int32_t ns, int32_t nim,
                                                CarmaDevice *device);

// //////////////////////////////////////////////////////
// ADDING PYR_ABS2 MODIFIED FOR ROOF-PRISM: ROOF_ABS2 //
// //////////////////////////////////////////////////////

template <class Tin, class Tout>
__global__ void abs2roof_krnl(Tout *g_odata, Tin *g_idata, Tout fact,
                              uint32_t ns, uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  int32_t x, y;

  // The following four lines shifts the index to induce a roll
  // along both axes.

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  uint32_t i2 = x + y * ns;

  // if (i2 < ns/2) i2 = i2 + ns/2 + ns/2*ns;

  Tin tmp1, tmp2;

  if (i < ns * ns / 2) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      tmp1 = g_idata[i + cpt * ns * ns];
      tmp2 = g_idata[i2 + cpt * ns * ns];

      g_odata[i + cpt * ns * ns] += (tmp2.x * tmp2.x + tmp2.y * tmp2.y) * fact;
      g_odata[i2 + cpt * ns * ns] += (tmp1.x * tmp1.x + tmp1.y * tmp1.y) * fact;
    }
  }
}

template <class Tin, class Tout>
void roof_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int32_t ns, int32_t nim,
               CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, ns * ns / 2, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  abs2roof_krnl<Tin, Tout><<<grid, threads>>>(d_odata, d_idata, fact, ns, nim);

  carma_check_msg("abs2pyr_kernel<<<>>> execution failed\n");
}

template void roof_abs2<cuFloatComplex, float>(float *d_odata,
                                               cuFloatComplex *d_idata,
                                               float fact, int32_t ns, int32_t nim,
                                               CarmaDevice *device);
template void roof_abs2<cuDoubleComplex, double>(double *d_odata,
                                                 cuDoubleComplex *d_idata,
                                                 double fact, int32_t ns, int32_t nim,
                                                 CarmaDevice *device);

// //////////////////////////////////////////////////////
// //////////////////////////////////////////////////////
// //////////////////////////////////////////////////////

template <class Tout, class Tin>
__global__ void submask_krnl(Tout *g_odata, Tin *g_mask, uint32_t n) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    g_odata[i].x = g_odata[i].x * g_mask[i];
    g_odata[i].y = g_odata[i].y * g_mask[i];
    i += blockDim.x * gridDim.x;
  }
}

template <class Tout, class Tin>
void apply_submask(Tout *d_odata, Tin *d_mask, int32_t n, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  submask_krnl<Tout, Tin><<<grid, threads>>>(d_odata, d_mask, n);

  carma_check_msg("submask_kernel<<<>>> execution failed\n");
}

template void apply_submask<cuFloatComplex, float>(cuFloatComplex *d_odata,
                                                 float *d_mask, int32_t n,
                                                 CarmaDevice *device);
template void apply_submask<cuDoubleComplex, double>(cuDoubleComplex *d_odata,
                                                   double *d_idata, int32_t n,
                                                   CarmaDevice *device);

template <class T>
__global__ void submaskpyr_krnl(T *g_odata, T *g_mask, uint32_t n) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    T tmp;
    tmp.x = g_odata[i].x * g_mask[i].x - g_odata[i].y * g_mask[i].y;
    tmp.y = g_odata[i].y * g_mask[i].x + g_odata[i].x * g_mask[i].y;

    g_odata[i].x = tmp.x;
    g_odata[i].y = tmp.y;

    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void pyr_submaskpyr(T *d_odata, T *d_mask, int32_t n, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  submaskpyr_krnl<T><<<grid, threads>>>(d_odata, d_mask, n);

  carma_check_msg("submask_kernel<<<>>> execution failed\n");
}

template void pyr_submaskpyr<cuFloatComplex>(cuFloatComplex *d_odata,
                                             cuFloatComplex *d_mask, int32_t n,
                                             CarmaDevice *device);
template void pyr_submaskpyr<cuDoubleComplex>(cuDoubleComplex *d_odata,
                                              cuDoubleComplex *d_idata, int32_t n,
                                              CarmaDevice *device);

template <class T>
__global__ void submaskpyr_krnl(T *g_odata, float *g_mask, uint32_t n) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    T tmp;
    tmp.x = g_odata[i].x * g_mask[i] - g_odata[i].y * g_mask[i];
    tmp.y = g_odata[i].y * g_mask[i] + g_odata[i].x * g_mask[i];

    g_odata[i].x = tmp.x;
    g_odata[i].y = tmp.y;

    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void pyr_submaskpyr(T *d_odata, float *d_mask, int32_t n, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  submaskpyr_krnl<T><<<grid, threads>>>(d_odata, d_mask, n);

  carma_check_msg("submask_kernel<<<>>> execution failed\n");
}

template void pyr_submaskpyr<cuFloatComplex>(cuFloatComplex *d_odata,
                                             float *d_mask, int32_t n,
                                             CarmaDevice *device);
template void pyr_submaskpyr<cuDoubleComplex>(cuDoubleComplex *d_odata,
                                              float *d_idata, int32_t n,
                                              CarmaDevice *device);

template <class Tout, class Tin>
__global__ void submask3d_krnl(Tout *g_odata, Tin *g_mask, uint32_t n,
                               uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n * n) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_odata[i + cpt * n * n].x = g_odata[i + cpt * n * n].x * g_mask[i];
      g_odata[i + cpt * n * n].y = g_odata[i + cpt * n * n].y * g_mask[i];
    }
  }
}

template <class Tout, class Tin>
void pyr_submask3d(Tout *d_odata, Tin *d_mask, int32_t n, int32_t nim,
                   CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  submask3d_krnl<Tout, Tin><<<grid, threads>>>(d_odata, d_mask, n, nim);

  carma_check_msg("submask3d_kernel<<<>>> execution failed\n");
}

template void pyr_submask3d<cuFloatComplex, float>(cuFloatComplex *d_odata,
                                                   float *d_mask, int32_t n,
                                                   int32_t nim,
                                                   CarmaDevice *device);
template void pyr_submask3d<cuDoubleComplex, double>(cuDoubleComplex *d_odata,
                                                     double *d_idata, int32_t n,
                                                     int32_t nim,
                                                     CarmaDevice *device);

template <class T>
__global__ void intensities_krnl(T *g_odata, T *g_idata, int32_t *subindx,
                                 int32_t *subindy, uint32_t ns,
                                 uint32_t nvalid, uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int32_t i2;
    g_odata[i] = 0;

    for (int32_t cpt = 0; cpt < nim; cpt++) {
      i2 = subindx[i + cpt * nvalid] + subindy[i + cpt * nvalid] * ns;
      g_odata[i] += g_idata[i2];
    }
  }
}

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy, int32_t ns,
                     int32_t nvalid, int32_t nim, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  intensities_krnl<T>
      <<<grid, threads, 0, stream>>>(d_odata, d_idata, subindx, subindy, ns, nvalid, nim);

  carma_check_msg("intensities_kernel<<<>>> execution failed\n");
}

template void pyr_intensities<float>(float *d_odata, float *d_idata,
                                     int32_t *subindx, int32_t *subindy, int32_t ns,
                                     int32_t nvalid, int32_t nim, CarmaDevice *device, cudaStream_t stream);
template void pyr_intensities<double>(double *d_odata, double *d_idata,
                                      int32_t *subindx, int32_t *subindy, int32_t ns,
                                      int32_t nvalid, int32_t nim,
                                      CarmaDevice *device, cudaStream_t stream);

template <class T>
__global__ void intensities_krnl(T *g_odata, T *g_idata, int32_t *subindx,
                                 int32_t *subindy, uint32_t ns,
                                 uint32_t nvalid) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int32_t iq1 = subindx[i] + subindy[i] * ns;
    int32_t iq2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    int32_t iq3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    int32_t iq4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;
    g_odata[i] = g_idata[iq1] + g_idata[iq2] + g_idata[iq3] + g_idata[iq4];
  }
}

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy, int32_t ns,
                     int32_t nvalid, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  intensities_krnl<T>
      <<<grid, threads, 0, stream>>>(d_odata, d_idata, subindx, subindy, ns, nvalid);

  carma_check_msg("intensities_kernel<<<>>> execution failed\n");
}

template void pyr_intensities<float>(float *d_odata, float *d_idata,
                                     int32_t *subindx, int32_t *subindy, int32_t ns,
                                     int32_t nvalid, CarmaDevice *device, cudaStream_t stream);
template void pyr_intensities<double>(double *d_odata, double *d_idata,
                                      int32_t *subindx, int32_t *subindy, int32_t ns,
                                      int32_t nvalid, CarmaDevice *device, cudaStream_t stream);

// //////////////////////////////////////////////////////////
// ADDING PYR_SUBSUM MODIFIED FOR HR pyramid              //
// //////////////////////////////////////////////////////////

template <class T>
__global__ void intensities2_krnl(T *g_odata, T *g_idata, int32_t *subindx,
                                  int32_t *subindy, uint32_t ns,
                                  uint32_t nvalid) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int32_t iq1 = subindx[i] + subindy[i] * ns;
    int32_t iq2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    int32_t iq3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    int32_t iq4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;
    g_odata[i] = g_idata[iq1] + g_idata[iq2] + g_idata[iq3] + g_idata[iq4];
  }
}

template <class T>
void pyr_intensities2(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy,
                      int32_t ns, int32_t nvalid, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  intensities2_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns, nvalid);

  carma_check_msg("intensities_kernel<<<>>> execution failed\n");
}

template void pyr_intensities2<float>(float *d_odata, float *d_idata,
                                      int32_t *subindx, int32_t *subindy, int32_t ns,
                                      int32_t nvalid, CarmaDevice *device);
template void pyr_intensities2<double>(double *d_odata, double *d_idata,
                                       int32_t *subindx, int32_t *subindy, int32_t ns,
                                       int32_t nvalid, CarmaDevice *device);

// //////////////////////////////////////////////////////////
// ADDING PYR_SUBSUM MODIFIED FOR ROOF-PRISM: ROOF_SUBSUM //
// //////////////////////////////////////////////////////////

template <class T>
__global__ void roof_intensities_krnl(T *g_odata, T *g_idata, int32_t *subindx,
                                      int32_t *subindy, uint32_t ns,
                                      uint32_t nvalid, uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int32_t i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = 0;

    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_odata[i] += g_idata[i2 + cpt * ns * ns];
    }
  }
}

template <class T>
void roof_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy,
                      int32_t ns, int32_t nvalid, int32_t nim, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  roof_intensities_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns, nvalid, nim);

  carma_check_msg("intensities_kernel<<<>>> execution failed\n");
}

template void roof_intensities<float>(float *d_odata, float *d_idata,
                                      int32_t *subindx, int32_t *subindy, int32_t ns,
                                      int32_t nvalid, int32_t nim,
                                      CarmaDevice *device);
template void roof_intensities<double>(double *d_odata, double *d_idata,
                                       int32_t *subindx, int32_t *subindy, int32_t ns,
                                       int32_t nvalid, int32_t nim,
                                       CarmaDevice *device);

// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////

template <class T>
__global__ void pyrfact_krnl(T *g_data, T fact, uint32_t n,
                             uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n] = g_data[i + cpt * n * n] * fact;
    }
    i += blockDim.x * gridDim.x;
  }
}

__global__ void pyrfact_krnl(cuFloatComplex *g_data, float fact, uint32_t n,
                             uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n].x = g_data[i + cpt * n * n].x * fact;
      g_data[i + cpt * n * n].y = g_data[i + cpt * n * n].y * fact;
    }
    i += blockDim.x * gridDim.x;
  }
}

__global__ void pyrfact_krnl(float *g_data, float fact1, float *fact2,
                             uint32_t n, uint32_t nim) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int32_t cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n] = g_data[i + cpt * n * n] * fact1 / fact2[0];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void pyr_fact(T *d_data, T fact, int32_t n, int32_t nim, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyrfact_krnl<T><<<grid, threads>>>(d_data, fact, n, nim);

  carma_check_msg("pyrfact_kernel<<<>>> execution failed\n");
}

template void pyr_fact<float>(float *d_data, float fact, int32_t ns, int32_t nim,
                              CarmaDevice *device);
template void pyr_fact<double>(double *d_odata, double fact, int32_t ns, int32_t nim,
                               CarmaDevice *device);

void pyr_fact(cuFloatComplex *d_data, float fact, int32_t n, int32_t nim,
              CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyrfact_krnl<<<grid, threads>>>(d_data, fact, n, nim);

  carma_check_msg("pyrfact_kernel<<<>>> execution failed\n");
}

void pyr_fact(float *d_data, float fact1, float *fact2, int32_t n, int32_t nim,
              CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, n * n, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyrfact_krnl<<<grid, threads>>>(d_data, fact1, fact2, n, nim);

  carma_check_msg("pyrfact_kernel<<<>>> execution failed\n");
}

template <class T>
__global__ void digitalize_krnl(T *camimg, float *binimg, float *dark,
                                float *flat, int32_t max_flux_per_pix, int32_t max_pix_value,
                                int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  float tmp;
  while (tid < N) {
    tmp = (binimg[tid] * flat[tid] + dark[tid]) / (float)max_flux_per_pix;
    tmp = (tmp < 1) ? tmp : 1.f;
    camimg[tid] = (T)(tmp * max_pix_value);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t digitalize(T *camimg, float *binimg, float *dark, float *flat,
               int32_t max_flux_per_pix, int32_t max_pix_value, int32_t N,
               CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  digitalize_krnl<<<grid, threads>>>(camimg, binimg, dark, flat, max_flux_per_pix,
                                     max_pix_value, N);

  carma_check_msg("digitalize_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t digitalize<uint16_t>(uint16_t *camimg, float *binimg, float *dark,
                                  float *flat, int32_t max_flux_per_pix,
                                  int32_t max_pix_value, int32_t N, CarmaDevice *device);
/*
   __global__ void fillcamplipup_krnl(cuFloatComplex *amplipup, float
   *phase,float *offset, float *mask, int32_t *indx, int32_t Nfft, int32_t Npup, int32_t npup,
   int32_t N)
   {
   int32_t nim,idim,idimx,idimy,idx;
   int32_t tid   = threadIdx.x + blockIdx.x * blockDim.x;

   while (tid < N) {
   nim   = tid / Npup;
   idim  = tid - nim * Npup;
   idimx = idim % npup;
   idimy = idim / npup;
   idx   = idimx + idimy * Nfft + nim * Nfft * Nfft;
   amplipup[idx].x = (cosf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
   amplipup[idx].y = (sinf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
   tid  += blockDim.x * gridDim.x;
   }
   }

   int32_t fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
   float *mask, int32_t *indx, int32_t Nfft, int32_t Npup, int32_t Nsub, int32_t npup, CarmaDevice
   *device)
   // here amplipup is a cube of data of size nfft x nfft x nsubap
   // phase is an array of size pupdiam x pupdiam
   // offset is an array of size pdiam x pdiam
   // mask is an array of size pupdiam x pupdiam
   // indx is an array of size pdiam x pdiam x nsubap
   // number of thread required : pdiam x pdiam x nsubap
   // Npup = pdiam x pdiam
   {

   struct cudaDeviceProp deviceProperties;
   cudaGetDeviceProperties(&deviceProperties, device);

   int32_t nb_blocks,nb_threads;
   get_num_blocks_and_threads(device, Npup * Nsub, nb_blocks, nb_threads);
   dim3 grid(nb_blocks), threads(nb_threads);

   fillcamplipup_krnl<<<grid,
   threads>>>(amplipup,phase,offset,mask,indx,Nfft,Npup,npup,N);
   carma_check_msg("fillcamplipup_kernel<<<>>> execution failed\n");

   return EXIT_SUCCESS;
   }




 */
