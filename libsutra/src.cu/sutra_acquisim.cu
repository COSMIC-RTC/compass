// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_acquisim.cu
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
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

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

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
                int *jvalid, CarmaDevice *device) {
  int Npix = npix * npix;
  int N = Npix * nsub;
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  bcube_krnl<T>
      <<<grid, threads>>>(bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, N);

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube<float>(float *bimage, float *bcube, int npix, int nsub,
                                int Nsub, int *ivalid, int *jvalid,
                                CarmaDevice *device);
template int fillbincube<double>(double *bimage, double *bcube, int npix,
                                 int nsub, int Nsub, int *ivalid, int *jvalid,
                                 CarmaDevice *device);

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
int fillbincube_async(CarmaHostObj<T> *image_telemetry, T *bimage, T *bcube,
                      int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
                      int nim, CarmaDevice *device) {
  int nstreams = image_telemetry->get_nb_streams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
  dim3 grid(N / (nstreams * threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion
  // of data
  for (int i = 0; i < nstreams; i++) {
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
template int fillbincube_async<float>(CarmaHostObj<float> *image_telemetry,
                                      float *bimage, float *bcube, int npix,
                                      int nsub, int Nsub, int *ivalid,
                                      int *jvalid, int nim,
                                      CarmaDevice *device);
template int fillbincube_async<double>(CarmaHostObj<double> *image_telemetry,
                                       double *bimage, double *bcube, int npix,
                                       int nsub, int Nsub, int *ivalid,
                                       int *jvalid, int nim,
                                       CarmaDevice *device);

template <class T>
int fillbincube_async(CarmaStreams *streams, CarmaObj<T> *bimage,
                      CarmaObj<T> *bcube, int npix, int nsub, int Nsub,
                      int *ivalid, int *jvalid, CarmaDevice *device) {
  int nstreams = streams->get_nb_streams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  // here nstreams should be : final image size / npix
  dim3 threads(nb_threads);
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

  carma_check_msg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int fillbincube_async<float>(CarmaStreams *streams,
                                      CarmaObj<float> *bimage,
                                      CarmaObj<float> *bcube, int npix,
                                      int nsub, int Nsub, int *ivalid,
                                      int *jvalid, CarmaDevice *device);
template int fillbincube_async<double>(CarmaStreams *streams,
                                       CarmaObj<double> *bimage,
                                       CarmaObj<double> *bcube, int npix,
                                       int nsub, int Nsub, int *ivalid,
                                       int *jvalid, CarmaDevice *device);
