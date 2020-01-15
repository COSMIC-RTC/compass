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

//! \file      sutra_acquisim.cu
//! \ingroup   libsutra
//! \class     sutra_acquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_acquisim.h>
#include <sutra_utils.h>

template <class T>
__global__ void bcube_krnl_2D(T *bimage, T *bcube, int *num_ssp) {
  const unsigned int x =
      threadIdx.x;  // on lance npix thread en x avec npix nb pix par sous-pup
  const unsigned int y = threadIdx.y;  // idem en y

  const unsigned int deltaX =
      blockIdx.x;  // on lance nssp blocks en x, i.e. 80 dans le cas d'un 80x80
  const unsigned int deltaY = blockIdx.y;  // idem en y

  const unsigned int npix = blockDim.x;
  const unsigned int nssp = gridDim.x;

  const T input =
      bimage[(x + npix * deltaX) + nssp * npix * (y + npix * deltaY)];

  const unsigned int nim = deltaX + nssp * deltaY;

  if (num_ssp[nim])
    bcube[(num_ssp[nim] - 1) * npix * npix + x + y * npix] = input;
  // (par exemple en yorick : num_ssp = (*(y_wfs(1)._isvalid))(*)(cum)(2:) *
  // ((y_wfs(1)._isvalid))(*) )
}

template <class T>
int fillbincube_2D(T *bimage, T *bcube, int npix, int nxsub, int *num_ssp) {
  dim3 grid(nxsub, nxsub), threads(npix, npix);

  bcube_krnl_2D<T><<<grid, threads>>>(bimage, bcube, num_ssp);

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube_2D<float>(float *bimage, float *bcube, int npix,
                                   int nsub, int *num_ssp);
template int fillbincube_2D<double>(double *bimage, double *bcube, int npix,
                                    int nsub, int *num_ssp);

template <class T>
__global__ void bcube_krnl(T *bimage, T *bcube, int npix, int npix2, int nsub,
                           int *ivalid, int *jvalid, int N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of
   the subap Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    // int idbin = xim + yim * nsub + ivalid[nim] * npix
    //           + jvalid[nim] * npix * nsub;
    int idbin = xim + ivalid[nim] + (yim + jvalid[nim]) * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fillbincube(T *bimage, T *bcube, int npix, int nsub, int Nsub, int *ivalid,
                int *jvalid, carma_device *device) {
  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  bcube_krnl<T>
      <<<grid, threads>>>(bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, N);

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube<float>(float *bimage, float *bcube, int npix, int nsub,
                                int Nsub, int *ivalid, int *jvalid,
                                carma_device *device);
template int fillbincube<double>(double *bimage, double *bcube, int npix,
                                 int nsub, int Nsub, int *ivalid, int *jvalid,
                                 carma_device *device);

template <class T>
__global__ void bcube_krnl_async(T *bimage, T *bcube, int npix, int npix2,
                                 int nsub, int *ivalid, int *jvalid, int N,
                                 int idstart) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of
   the subap Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += idstart;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    int idbin =
        xim + yim * nsub + ivalid[nim] * npix + jvalid[nim] * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fillbincube_async(carma_host_obj<T> *image_telemetry, T *bimage, T *bcube,
                      int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
                      int nim, carma_device *device) {
  int nstreams = image_telemetry->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N / (nstreams * threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int i = 0; i < nstreams; i++) {
    cudaMemcpyAsync(&(bimage[i * nim / nstreams]),
                    image_telemetry->getDataAt(i * nim / nstreams),
                    sizeof(float) * nim / nstreams, cudaMemcpyHostToDevice,
                    image_telemetry->get_cudaStream_t(i));
    bcube_krnl_async<T>
        <<<grid, threads, 0, image_telemetry->get_cudaStream_t(i)>>>(
            bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, N,
            i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
  }
  // cudaStreamSynchronize(image_telemetry->get_cudaStream_t(nstreams-1));
  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube_async<float>(carma_host_obj<float> *image_telemetry,
                                      float *bimage, float *bcube, int npix,
                                      int nsub, int Nsub, int *ivalid,
                                      int *jvalid, int nim,
                                      carma_device *device);
template int fillbincube_async<double>(carma_host_obj<double> *image_telemetry,
                                       double *bimage, double *bcube, int npix,
                                       int nsub, int Nsub, int *ivalid,
                                       int *jvalid, int nim,
                                       carma_device *device);

template <class T>
int fillbincube_async(carma_streams *streams, carma_obj<T> *bimage,
                      carma_obj<T> *bcube, int npix, int nsub, int Nsub,
                      int *ivalid, int *jvalid, carma_device *device) {
  int nstreams = streams->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N / (nstreams * threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int i = 0; i < nstreams; i++) {
    bcube_krnl_async<T><<<grid, threads, 0, streams->get_stream(i)>>>(
        *bimage, *bcube, npix, Npix, Nsub, ivalid, jvalid, N, i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
    // will only
    //   commence executing when all previous CUDA calls in stream x have
    //   completed
  }

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube_async<float>(carma_streams *streams,
                                      carma_obj<float> *bimage,
                                      carma_obj<float> *bcube, int npix,
                                      int nsub, int Nsub, int *ivalid,
                                      int *jvalid, carma_device *device);
template int fillbincube_async<double>(carma_streams *streams,
                                       carma_obj<double> *bimage,
                                       carma_obj<double> *bcube, int npix,
                                       int nsub, int Nsub, int *ivalid,
                                       int *jvalid, carma_device *device);
