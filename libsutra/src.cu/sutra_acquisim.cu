// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_acquisim.cu
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_acquisim.hpp>
#include <sutra_utils.hpp>

template <class T>
__global__ void bcube_krnl_2D(T *bimage, T *bcube, int32_t *num_ssp) {
  const uint32_t x =
      threadIdx.x;  // on lance npix thread en x avec npix nb pix par sous-pup
  const uint32_t y = threadIdx.y;  // idem en y

  const uint32_t deltaX =
      blockIdx.x;  // on lance nssp blocks en x, i.e. 80 dans le cas d'un 80x80
  const uint32_t deltaY = blockIdx.y;  // idem en y

  const uint32_t npix = blockDim.x;
  const uint32_t nssp = gridDim.x;

  const T input =
      bimage[(x + npix * deltaX) + nssp * npix * (y + npix * deltaY)];

  const uint32_t nim = deltaX + nssp * deltaY;

  if (num_ssp[nim])
    bcube[(num_ssp[nim] - 1) * npix * npix + x + y * npix] = input;
  // (par exemple en yorick : num_ssp = (*(y_wfs(1)._isvalid))(*)(cum)(2:) *
  // ((y_wfs(1)._isvalid))(*) )
}

template <class T>
int32_t fillbincube_2D(T *bimage, T *bcube, int32_t npix, int32_t nxsub, int32_t *num_ssp) {
  dim3 grid(nxsub, nxsub), threads(npix, npix);

  bcube_krnl_2D<T><<<grid, threads>>>(bimage, bcube, num_ssp);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t fillbincube_2D<float>(float *bimage, float *bcube, int32_t npix,
                                   int32_t nsub, int32_t *num_ssp);
template int32_t fillbincube_2D<double>(double *bimage, double *bcube, int32_t npix,
                                    int32_t nsub, int32_t *num_ssp);

template <class T>
__global__ void bcube_krnl(T *bimage, T *bcube, int32_t npix, int32_t npix2, int32_t nsub,
                           int32_t *ivalid, int32_t *jvalid, int32_t N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of
   the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    int32_t nim = tid / npix2;
    int32_t tidim = tid - nim * npix2;
    int32_t xim = tidim % npix;
    int32_t yim = tidim / npix;

    // int32_t idbin = xim + yim * nsub + ivalid[nim] * npix
    //           + jvalid[nim] * npix * nsub;
    int32_t idbin = xim + ivalid[nim] + (yim + jvalid[nim]) * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillbincube(T *bimage, T *bcube, int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid,
                int32_t *jvalid, CarmaDevice *device) {
  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  bcube_krnl<T>
      <<<grid, threads>>>(bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, N);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t fillbincube<float>(float *bimage, float *bcube, int32_t npix, int32_t nsub,
                                int32_t Nsub, int32_t *ivalid, int32_t *jvalid,
                                CarmaDevice *device);
template int32_t fillbincube<double>(double *bimage, double *bcube, int32_t npix,
                                 int32_t nsub, int32_t Nsub, int32_t *ivalid, int32_t *jvalid,
                                 CarmaDevice *device);

template <class T>
__global__ void bcube_krnl_async(T *bimage, T *bcube, int32_t npix, int32_t npix2,
                                 int32_t nsub, int32_t *ivalid, int32_t *jvalid, int32_t N,
                                 int32_t idstart) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of
   the subap Npix = npix x npix
   */
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += idstart;

  while (tid < N) {
    int32_t nim = tid / npix2;
    int32_t tidim = tid - nim * npix2;
    int32_t xim = tidim % npix;
    int32_t yim = tidim / npix;

    int32_t idbin =
        xim + yim * nsub + ivalid[nim] * npix + jvalid[nim] * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillbincube_async(CarmaHostObj<T> *image_telemetry, T *bimage, T *bcube,
                      int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid, int32_t *jvalid,
                      int32_t nim, CarmaDevice *device) {
  int32_t nstreams = image_telemetry->get_nb_streams();

  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
  dim3 grid(N / (nstreams * threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int32_t i = 0; i < nstreams; i++) {
    cudaMemcpyAsync(&(bimage[i * nim / nstreams]),
                    image_telemetry->get_data_at(i * nim / nstreams),
                    sizeof(float) * nim / nstreams, cudaMemcpyHostToDevice,
                    image_telemetry->get_cuda_stream(i));
    bcube_krnl_async<T>
        <<<grid, threads, 0, image_telemetry->get_cuda_stream(i)>>>(
            bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, N,
            i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
  }
  // cudaStreamSynchronize(image_telemetry->get_cuda_stream(nstreams-1));
  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t fillbincube_async<float>(CarmaHostObj<float> *image_telemetry,
                                      float *bimage, float *bcube, int32_t npix,
                                      int32_t nsub, int32_t Nsub, int32_t *ivalid,
                                      int32_t *jvalid, int32_t nim,
                                      CarmaDevice *device);
template int32_t fillbincube_async<double>(CarmaHostObj<double> *image_telemetry,
                                       double *bimage, double *bcube, int32_t npix,
                                       int32_t nsub, int32_t Nsub, int32_t *ivalid,
                                       int32_t *jvalid, int32_t nim,
                                       CarmaDevice *device);

template <class T>
int32_t fillbincube_async(CarmaStreams *streams, CarmaObj<T> *bimage,
                      CarmaObj<T> *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
                      int32_t *ivalid, int32_t *jvalid, CarmaDevice *device) {
  int32_t nstreams = streams->get_nb_streams();

  int32_t Npix = npix * npix;
  int32_t N = Npix * nsub;
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
  dim3 grid(N / (nstreams * threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int32_t i = 0; i < nstreams; i++) {
    bcube_krnl_async<T><<<grid, threads, 0, streams->get_stream(i)>>>(
        *bimage, *bcube, npix, Npix, Nsub, ivalid, jvalid, N, i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
  }

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t fillbincube_async<float>(CarmaStreams *streams,
                                      CarmaObj<float> *bimage,
                                      CarmaObj<float> *bcube, int32_t npix,
                                      int32_t nsub, int32_t Nsub, int32_t *ivalid,
                                      int32_t *jvalid, CarmaDevice *device);
template int32_t fillbincube_async<double>(CarmaStreams *streams,
                                       CarmaObj<double> *bimage,
                                       CarmaObj<double> *bcube, int32_t npix,
                                       int32_t nsub, int32_t Nsub, int32_t *ivalid,
                                       int32_t *jvalid, CarmaDevice *device);
