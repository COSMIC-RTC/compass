#include <sutra_lgs.h>

texture<float, 3, cudaReadModeElementType> tex3;  // 3D texture

__device__ float
interp_transf(float *idata, float tidim, int npix, float pixsize, float norm,
    float shiftx, float deltah, int hmax) {
  int xref;
  float weightx;

  // determine the interpolated position
  if (npix % 2 == 0)
    tidim -= (npix / 2.0f + 0.5); // centered pix pos
  else
    tidim -= (npix / 2.0f); // centered pix pos

  tidim *= pixsize;     // centered pix size
  tidim /= norm;        // in meter 
  tidim += shiftx;      // in meter on profile centered on cog
  tidim /= deltah;      // in pixels on profile array

  if (tidim < 0) {
    xref = 0;
    weightx = 0.;
  } else {
    if (tidim < hmax - 2) {
      xref = (int) tidim;
      weightx = tidim - xref;
    } else {
      xref = hmax - 2;
      weightx = 1.;
    }
  }

  return (idata[xref] * (1 - weightx) + idata[xref + 1] * weightx);
}

__global__ void
iprof_krnl(cuFloatComplex *profout, float *profin, float *profinc, int npix,
    float *doffaxis, float hg, float pixsize, float h0, float deltah, int hmax,
    int Ntot) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int nim = tid / npix;
    float tidim = tid - nim * npix;
    float norm;

    norm = 206265.0f * doffaxis[nim] / (hg * hg);

    if (pixsize > deltah * norm) {
      profout[tid].x = interp_transf(profinc, tidim + 1, npix, pixsize, norm,
          hg - h0, deltah, hmax)
          - interp_transf(profinc, tidim, npix, pixsize, norm, hg - h0, deltah,
              hmax);
      profout[tid].y = 0.0f;
    } else {
      profout[tid].x = interp_transf(profin, tidim, npix, pixsize, norm,
          hg - h0, deltah, hmax);
      profout[tid].y = 0.0f;
    }

    tid += blockDim.x * gridDim.x;
  }
}

int
interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum, int npix,
    float *doffaxis, float hg, float pixsize, float h0, float deltah, int hmax,
    int Ntot, int device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, Ntot, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  iprof_krnl<<<grid, threads>>>(profout, prof1d, profcum, npix, doffaxis, hg,
      pixsize, h0, deltah, hmax, Ntot);
  cutilCheckMsg("iprof_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void
tftbeam_krnl(cuFloatComplex *profout, cuFloatComplex *fbeam, int N, int Ntot) {
  int idim;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    idim = tid % N;
    __shared__ cuFloatComplex tmp;
    tmp = profout[tid];
    profout[tid].x = tmp.x * fbeam[idim].x - tmp.y * fbeam[idim].y;
    profout[tid].y = tmp.y * fbeam[idim].x + tmp.x * fbeam[idim].y;
    tid += blockDim.x * gridDim.x;
  }
}

int
times_ftbeam(cuFloatComplex *profout, cuFloatComplex *ftbeam, int N, int Ntot,
    int device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, Ntot, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  tftbeam_krnl<<<grid, threads>>>(profout, ftbeam, N, Ntot);
  cutilCheckMsg("tftbeam_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void
rollbeamexp_krnl(float *imout, cuFloatComplex *iprof, float *beam, int N,
    int Ntot) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int nim = tid / (N * N);
    int tidim = tid % (N * N);
    int tidprof = tidim % N;
    int tidbeam = tidim / N;
    int x = (tidprof + (N / 2)) % N + nim * N; // for roll
    imout[tid] = abs(iprof[x].x) * beam[tidbeam];
    tid += blockDim.x * gridDim.x;
  }
}

int
rollbeamexp(float *imout, cuFloatComplex *iprof, float *beam, int N, int Ntot,
    int device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, Ntot, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  rollbeamexp_krnl<<<grid, threads>>>(imout, iprof, beam, N, Ntot);
  cutilCheckMsg("beamexp_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;

}

__global__ void
rotate_krnl(cuFloatComplex *odata, float *idata, int width, int height,
    float *theta, float center, int N, int Ntot)
// rotate a cube manually
    {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int nim = tid / N;
    int tidim = tid - nim * N;

    int y = tidim / width;
    int x = tidim - y * width;

    float ucent = width / 2 - center;
    float vcent = height / 2 - center;

    float u = x - ucent;
    float v = y - vcent;
    // transform coordinates
    float tu = u * cos(theta[nim]) + v * sin(theta[nim]) + ucent;
    float tv = -u * sin(theta[nim]) + v * cos(theta[nim]) + ucent;

    int uref, vref;
    float weightu, weightv;

    if (tu < 0) {
      uref = 0;
      weightu = 0.;
    } else {
      if (tu < width - 2) {
        uref = (int) tu;
        weightu = tu - uref;
      } else {
        uref = width - 2;
        weightu = 1.;
      }
    }

    if (tv < 0) {
      vref = 0;
      weightv = 0.;
    } else {
      if (tv < height - 2) {
        vref = (int) tv;
        weightv = tv - vref;
      } else {
        vref = height - 2;
        weightv = 1.;
      }
    }

    int ind1, ind2, ind3, ind4;
    ind1 = ind2 = ind3 = ind4 = nim * N;

    ind1 += uref + vref * width;
    ind2 += uref + 1 + vref * width;
    ind3 += uref + (vref + 1) * width;
    ind4 += uref + 1 + (vref + 1) * width;

    odata[tid].x = (idata[ind1] * (1 - weightu) + idata[ind2] * weightu)
        * (1 - weightv)
        + (idata[ind3] * (1 - weightu) + idata[ind4] * weightu) * weightv;
    odata[tid].y = 0.0f;

    tid += blockDim.x * gridDim.x;
  }
}

int
lgs_rotate(cuFloatComplex *odata, float *idata, int width, int height,
    float *theta, float center, int Ntot, int device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, Ntot, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  int N = width * height;

  rotate_krnl<<<grid, threads>>>(odata, idata, width, height, theta, center, N,
      Ntot);
  cutilCheckMsg("rotate_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void
rotate3d_krnl(cuFloatComplex *g_odata, int width, int height, int N,
    float *theta, float center, int Ntot)
// rotate a cube using 3d texture fetch
    {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < Ntot) {
    int nim = tid / N;
    int tidim = tid - nim * N;

    int y = tidim / width;
    int x = tidim - y * width;

    float ucent = width / 2 - center;
    float vcent = height / 2 - center;

    float u = x - ucent;
    float v = y - vcent;
    // transform coordinates
    float tu = u * cos(theta[nim]) + v * sin(theta[nim]) + ucent;
    float tv = -u * sin(theta[nim]) + v * cos(theta[nim]) + ucent;

    g_odata[tid].x = tex3D(tex3, tu + 0.5, tv + 0.5, nim + 0.5);
    g_odata[tid].y = 0.0f;
  }
}

int
rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
    cudaArray *d_array, cudaChannelFormatDesc channelDesc, int width,
    int height, float *theta, float center, int Ntot, int device) {
  tex3.normalized = false;
  tex3.filterMode = cudaFilterModeLinear;      // linear interpolation
  tex3.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
  tex3.addressMode[1] = cudaAddressModeClamp;
  tex3.addressMode[2] = cudaAddressModeClamp;

  // copy the data into the array
  cutilSafeCall(cudaMemcpy3D(&copyParams));
  // bind array to 3D texture
  cutilSafeCall(cudaBindTextureToArray(tex3, d_array, channelDesc));

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);

  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount * 8;
  int nThreads = (Ntot + nBlocks - 1) / nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (Ntot + nThreads - 1) / nThreads;
  }

  int N = width * height;

  dim3 grid(nBlocks), threads(nThreads);

  //cout << nBlocks << " " << nx / nBlocks << " " << ny / nBlocks << " " <<  nz / nBlocks <<endl;
  rotate3d_krnl<<<grid, threads, 0>>>(d_odata, width, height, N, theta, center,
      Ntot);

  cutilCheckMsg("rotate3d_krnl <<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
