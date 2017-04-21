#include <sutra_centroider.h>
#include "carma_utils.cuh"

/**********************************
  _  __                    _      *
 | |/ /___ _ __ _ __   ___| |___  *
 | ' // _ \ '__| '_ \ / _ \ / __| *
 | . \  __/ |  | | | |  __/ \__ \ *
 |_|\_\___|_|  |_| |_|\___|_|___/ *
                                  *
 **********************************/


template<class T>
__device__ void scanmax_krnl(T *sdata, int *values, int size, int n) {
  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          values[n] = n + s;
          sdata[n] = sdata[n + s];
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          values[n] = n + s;
          sdata[n] = sdata[n + s];
        }
      }
      __syncthreads();
    }
  }
}

template<class T>
__device__ inline void sortmax_krnl(T *sdata, unsigned int *values, int size,
    int n) {

  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
    }
  }
}

template<class T>
__device__ inline void sortmax_krnl(T *sdata, int *values, int size, int n) {

  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
    }
  }
}
template<class T>
__global__ void sortmax(T *g_idata, T *g_odata, unsigned int *values,
    int nmax, int Npix, int size, int nelem_thread) {
  extern __shared__ uint svalues[];
  T *sdata = (T*) &svalues[Npix];
/*
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = tid;
  sdata[tid] = g_idata[i];
*/
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int y = (tid / n) + 1;
  int idim;
  int sdim;

  for (int cc = 0; cc < nelem_thread; cc++) {
	 idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
	 sdim = tid * nelem_thread + cc;
	 if (idim < size){
	   sdata[sdim] = g_idata[idim];
	   svalues[sdim] = sdim;
	 }
  }

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {
	  for (int cpt = 0; cpt < nelem_thread ; cpt++){
		  sdim = tid * nelem_thread + cpt;
		  if (sdim >= cc)
			  sortmax_krnl(&(sdata[cc]), &(svalues[cc]), Npix - cc, sdim - cc);
		  __syncthreads();
	  }


  }
  for (int cpt = 0; cpt < nelem_thread ; cpt++){
	  sdim = tid * nelem_thread + cpt;
	  if (sdim < nmax) {
		g_odata[nmax * blockIdx.x + sdim] = sdata[sdim]; // - sdata[nmax - 1];
		values[nmax * blockIdx.x + sdim] = svalues[sdim];
	  }
  }
  /*
   __syncthreads();
   if ((blockIdx.x == 0) && (tid < nmax))
   printf("tid %d sdata %f \n",tid,g_odata[tid]);
   */
}

template<class T>
__device__ inline void sortmaxi_krnl(T *sdata, unsigned int *values, int size,
    int n) {

  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
    }
  }
}

template<class T>
__global__ void sortmaxi(T *g_idata, int *values, int nmax, int offx, int offy,
    int npix, int Npix, int Npix2) {
  extern __shared__ uint svalues[];
  T *sdata = (T*) &svalues[blockDim.x];

  // load shared mem
  unsigned int tid = threadIdx.x;

  svalues[tid] = tid;
  int nlig = tid / npix;
  int ncol = tid - nlig * npix;
  int idx = offx + ncol + (nlig + offy) * Npix + blockIdx.x * Npix2;
  sdata[tid] = g_idata[idx];

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {

    if (tid >= cc)
      sortmaxi_krnl(&(sdata[cc]), &(svalues[cc]), blockDim.x - cc, tid - cc);

    __syncthreads();

  }

  if (tid < nmax) {
    nlig = svalues[tid] / npix;
    ncol = svalues[tid] - nlig * npix;
    idx = offx + ncol + (nlig + offy) * Npix;
    values[nmax * blockIdx.x + tid] = idx;
  }
  __syncthreads();

}

template<class T>
__global__ void centroid_max(T *g_idata, T *g_odata, int n, int nmax, int nsub,
    T scale, T offset) {
  extern __shared__ unsigned int svalues[];
  T *sdata = (T*) &svalues[blockDim.x];
  T subsum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = tid;
  sdata[tid] = g_idata[i];

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {

    if (tid >= cc)
      sortmax_krnl(&(sdata[cc]), &(svalues[cc]), blockDim.x - cc, tid - cc);

    __syncthreads();

  }

  // at this point the nmax first elements of sdata are the nmax brightest
  // pixels and the nmax first elements of svalues are their positions

  // first copy the brightest values out for reduction
  //__syncthreads();
  // if(tid == 0)
  const T minim = sdata[nmax - 1];
  __syncthreads();
  if (tid < nmax - 1)
    sdata[tid] -= minim;       // BIG BUG if nmax > 32
  else
    sdata[tid] = 0.0f;

  __syncthreads();
  if ((tid >= nmax) && (tid < 2 * nmax)) {
    sdata[tid] = sdata[tid - nmax];
  }
  /*
   if(tid == 0){
   g_odata[blockIdx.x] = sdata[tid];
   g_odata[blockIdx.x + nsub] = minim;
   }
   */
  __syncthreads();

  //if (tid < nmax)
  reduce_krnl(sdata, nmax, tid); // BIG BUG le retour

  //__syncthreads();
  // get the sum per subap
  if (tid == 0)
    subsum = (abs(sdata[tid]) > 1.e-6 ? sdata[tid] : 0.0f);

  __syncthreads();

  // if (tid == 0)
//	  g_odata[blockIdx.x] = subsum;

  // put back the brightest pixels values
  if ((tid >= nmax) && (tid < 2 * nmax))
    sdata[tid - nmax] = sdata[tid];

  __syncthreads();

  // compute the centroid on the first part of the array
  if (tid < nmax)
    sdata[tid] *= ((svalues[tid] % n) + 1);
  // x centroid

  reduce_krnl(sdata, nmax, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x] = (
        subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);
  __syncthreads();
  // put back the brightest pixels values
  if ((tid >= nmax) && (tid < 2 * nmax))
    sdata[tid - nmax] = sdata[tid];

  __syncthreads();

  // compute the centroid on the first part of the array
  if (tid < nmax)
    sdata[tid] *= (svalues[tid] / n + 1);
  // y centroid

  reduce_krnl(sdata, nmax, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x + nsub] = (
        subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);

}

template<class T>
__global__ void centroid_max2(T *g_idata, T *g_odata, T *g_minim, int n,
    int nmax, int nsub, T scale, T offset) {
  extern __shared__ uint svalues[];
  T *sdata = (T*) &svalues[blockDim.x];
//  T subsum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = tid;
  sdata[tid] = g_idata[i];

  //__syncthreads();
  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {

    if (tid >= cc)
      sortmax_krnl(&(sdata[cc]), &(svalues[cc]), blockDim.x - cc, tid - cc);

    //__syncthreads();
    __syncthreads();
  }

  // at this point the nmax first elements of sdata are the nmax brightest
  // pixels and the nmax first elements of svalues are their positions

  // first copy the brightest values out for reduction
  //__syncthreads();
  __syncthreads();
  if (tid == 0) {
    g_minim[blockIdx.x] = sdata[nmax - 1];
  }
  //__syncthreads();
  __syncthreads();
  if (tid < nmax) {
    sdata[tid] -= g_minim[blockIdx.x];
  } else {
    sdata[tid] = 0.0f;
  }

  __syncthreads();
  if ((tid >= nmax) && (tid < 2 * nmax)) {
    sdata[tid] = sdata[tid - nmax];
  }

  //__syncthreads();
  __syncthreads();

  if (tid < nmax)
    g_minim[tid] = sdata[tid];
  __syncthreads();

  // if (tid == 0){
  reduce_krnl(g_minim, nmax, tid);
  /*
   for (int cc =1 ; cc<nmax ; cc++){
   g_minim[cc] = sdata[cc];
   sdata[tid] += g_minim[cc];
   __threadfence();
   }
   }
   */

  //__syncthreads();
  __syncthreads();
  // get the sum per subap
//  if (tid == 0)
//    subsum = (abs(sdata[tid]) > 1.e-6 ? sdata[tid] : 0.0f);

  //__syncthreads();
  __syncthreads();
  if (tid == 0) {
    g_odata[blockIdx.x] = g_minim[tid];
    g_odata[blockIdx.x + nsub] = 0.0f;
  }
  //__syncthreads();
  __syncthreads();
  /*
   // if (tid == 0)
   //	  g_odata[blockIdx.x] = subsum;

   // put back the brightest pixels values
   if ((tid >= nmax) && (tid < 2 * nmax))
   sdata[tid - nmax] = sdata[tid];

   __syncthreads();

   // compute the centroid on the first part of the array
   if (tid < nmax)
   sdata[tid] *= ((svalues[tid] % n) + 1);
   // x centroid

   reduce_krnl(sdata, nmax, tid);
   //__syncthreads();
   if (tid == 0)
   g_odata[blockIdx.x] =(subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);
   __syncthreads();
   // put back the brightest pixels values
   if ((tid >= nmax) && (tid < 2 * nmax))
   sdata[tid - nmax] = sdata[tid];

   __syncthreads();

   // compute the centroid on the first part of the array
   if (tid < nmax)
   sdata[tid] *= (svalues[tid] / n + 1);
   // y centroid

   reduce_krnl(sdata, nmax, tid);
   //__syncthreads();
   if (tid == 0)
   g_odata[blockIdx.x + nsub] = (subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);
   */
}
template<class T>
__global__ void centroid_bpix(int nsub, int n, T *g_idata, unsigned int *values,
    T *g_odata, T scale, T offset) {
  extern __shared__ uint svalues[];
  T *sdata = (T*) &svalues[blockDim.x];
  T subsum;
  // T minimum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = values[i];
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  __syncthreads();
  // get the sum per subap
  if (tid == 0)
    subsum = (abs(sdata[tid]) > 1.e-6 ? sdata[tid] : 0.0f);

  __syncthreads();

// Reload sdata
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= ((svalues[tid] % n) + 1);
  // x centroid
  __syncthreads();
  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x] = (
        subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);
  __syncthreads();
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= (svalues[tid] / n + 1);
  // y centroid
  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x + nsub] = (
        subsum != 0.0f ? ((sdata[tid] / subsum) - offset) * scale : 0.0f);

}
template<class T>
__global__ void centroidx_2D(T *g_idata, T *g_odata, T *alpha, unsigned int n,
    unsigned int N, T scale, T offset, unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  const unsigned int tid = threadIdx.x + blockIdx.y * blockDim.x;
  const unsigned int i = blockIdx.x * blockDim.x * blockDim.y + tid;
  //unsigned int x = (tid % n) + 1;
  sdata[tid] = (i < N) ? g_idata[i] * (threadIdx.x+1) : 0;

  __syncthreads();
  reduce_krnl(sdata, blockDim.x*blockDim.y, tid);
  __syncthreads();

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}
template<class T>
__global__ void centroidy_2D(T *g_idata, T *g_odata, T *alpha, unsigned int n,
    unsigned int N, T scale, T offset, unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  const unsigned int tid = threadIdx.x + blockIdx.y * blockDim.x;
  const unsigned int i = blockIdx.x * blockDim.x * blockDim.y + tid;
  //unsigned int x = (tid % n) + 1;
  sdata[tid] = (i < N) ? g_idata[i] * (threadIdx.y+1) : 0;

  __syncthreads();
  reduce_krnl(sdata, blockDim.x*blockDim.y, tid);
  __syncthreads();

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}
template
__global__ void centroidx_2D<float>(float *g_idata, float *g_odata, float *alpha, unsigned int n,
    unsigned int N, float scale, float offset, unsigned int nelem_thread);
template
__global__ void centroidy_2D<float>(float *g_idata, float *g_odata, float *alpha, unsigned int n,
    unsigned int N, float scale, float offset, unsigned int nelem_thread);

template<class T>
__global__ void centroidx(T *g_idata, T *g_odata, T *alpha, unsigned int n,
    unsigned int N, T scale, T offset, unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int x = (tid % n) + 1;
  unsigned int x;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * x;
    else
      sdata[tid] += 0;
  }

  //sdata[tid] = (i < N) ? g_idata[i] * x : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void centroidy(T *g_idata, T *g_odata, T *alpha, unsigned int n,
    unsigned int N, T scale, T offset, unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int y = (tid / n) + 1;
  unsigned int y;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    y = ((tid * nelem_thread + cc) / n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * y;
    else
      sdata[tid] += 0;
  }

  //sdata[tid] = (i < N) ? g_idata[i] * y : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void centroidx_async(T *g_idata, T *g_odata, T *alpha,
    unsigned int n, unsigned int N, T scale, T offset, int stream_offset) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int x = (tid % n) + 1;
  i += stream_offset * blockDim.x;

  sdata[tid] = (i < N) ? g_idata[i] * x : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x + stream_offset] = ((sdata[0] * 1.0
        / (alpha[blockIdx.x + stream_offset] + 1.e-6)) - offset) * scale;
}

template<class T>
__global__ void centroidy_async(T *g_idata, T *g_odata, T *alpha,
    unsigned int n, unsigned int N, T scale, T offset, int stream_offset) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = (tid / n) + 1;
  i += stream_offset * blockDim.x;

  sdata[tid] = (i < N) ? g_idata[i] * y : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x + stream_offset] = ((sdata[0] * 1.0
        / (alpha[blockIdx.x + stream_offset] + 1.e-6)) - offset) * scale;
}

template<class T>
void get_centroids_async(int threads, int blocks, int n, carma_streams *streams,
    T *d_idata, T *d_odata, T *alpha, T scale, T offset) {
  int nstreams = streams->get_nbStreams();
  int nbelem = threads * blocks;

  dim3 dimBlock(threads);
  dim3 dimGrid(blocks / nstreams);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  for (int i = 0; i < nstreams; i++) {
    centroidx_async<T> <<<dimGrid, dimBlock, smemSize, streams->get_stream(i)>>>(
        d_idata, d_odata, alpha, n, nbelem, scale, offset,
        i * blocks / nstreams);

    carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

    centroidy_async<T> <<<dimGrid, dimBlock, smemSize, streams->get_stream(i)>>>(
        d_idata, &(d_odata[blocks]), alpha, n, nbelem, scale, offset,
        i * blocks / nstreams);

    carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
  }
}

template void
get_centroids_async<float>(int threads, int blocks, int n,
    carma_streams *streams, float *d_idata, float *d_odata, float *alpha,
    float scale, float offset);

template<class T>
__global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T thresh,
    unsigned int n, unsigned int N, T scale, T offset,
    unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int x = (tid % n) + 1;
  unsigned int x;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += (g_idata[idim] > thresh) ? g_idata[idim] * x : 0;
    else
      sdata[tid] += 0;
  }

  // if (i < N)
  //   sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * x : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T thresh,
    unsigned int n, unsigned int N, T scale, T offset,
    unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int y = (tid / n) + 1;
  unsigned int y;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    y = ((tid * nelem_thread + cc) / n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += (g_idata[idim] > thresh) ? g_idata[idim] * y : 0;
    else
      sdata[tid] += 0;
  }

  // if (i < N)
  //   sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * y : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T *weights,
    unsigned int n, unsigned int N, T scale, T offset,
    unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int x = (tid % n) + 1;
  unsigned int x;
  int idim;
  sdata[tid] = 0.0f;
  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * x * weights[idim];
    else
      sdata[tid] += 0.0f;
  }
  __syncthreads();

  // sdata[tid] = (i < N) ? g_idata[i] * x * weights[i] : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);
  //if(tid == 0)
  //	printf("blockIdx %d sdata %f \n",blockIdx.x,sdata[tid]);
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T *weights,
    unsigned int n, unsigned int N, T scale, T offset,
    unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int y = (tid / n) + 1;
  unsigned int y;
  int idim;
  sdata[tid] = 0.0f;
  for (int cc = 0; cc < nelem_thread; cc++) {
    y = ((tid * nelem_thread + cc) / n) + 1;
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * y * weights[idim];
    else
      sdata[tid] += 0.0f;
  }

  //sdata[tid] = (i < N) ? g_idata[i] * y * weights[i] : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6))
        - offset) * scale;
}

template<class T>
__global__ void interp_parab(T *g_idata, T *g_centroids, int *g_values,
    T *g_matinterp, int sizex, int sizey, int nvalid, int Npix, int Npix2,
    T scale, T offset) {
  extern __shared__ T pidata[];
  T *scoeff = (T*) &pidata[blockDim.x];
  T *m_interp = (T*) &scoeff[6];

  int offy = g_values[blockIdx.x] / Npix;
  int offx = g_values[blockIdx.x] - offy * Npix;
  offx -= sizex / 2;
  offy -= sizey / 2;
  int nlig = threadIdx.x / sizex;
  int ncol = threadIdx.x - nlig * sizex;
  int idx = offx + ncol + (nlig + offy) * Npix + blockIdx.x * Npix2;

  // load shared mem
  pidata[threadIdx.x] = g_idata[idx];

  __syncthreads();

  for (int cc = 0; cc < 6; cc++)
    m_interp[cc * sizex * sizey + threadIdx.x] = g_matinterp[cc * sizex * sizey
        + threadIdx.x] * pidata[threadIdx.x];

  __syncthreads();

  // do reduction for each 6 coeffs
  if (threadIdx.x < 6) {
    scoeff[threadIdx.x] = 0.0f;
    for (int cc = 0; cc < sizex * sizey; cc++)
      scoeff[threadIdx.x] += m_interp[cc + threadIdx.x * sizex * sizey];
  }

  __syncthreads();

  // now retreive x0 and y0 from scoeff
  if (threadIdx.x < 2) {
    T denom = scoeff[2] * scoeff[2] - 4.0f * scoeff[1] * scoeff[0];
    if (denom == 0) {
      if (threadIdx.x == 0)
        g_centroids[blockIdx.x] = 0.0f;
      if (threadIdx.x == 1)
        g_centroids[blockIdx.x + nvalid] = 0.0f;
    } else {
      if (threadIdx.x == 0) {
        g_centroids[blockIdx.x] = (2.0f * scoeff[1] * scoeff[3]
            - scoeff[4] * scoeff[2]) / denom;
        int xm = (2 * offx + sizex);
        xm = ((xm & 1) == 0) ? xm / 2 : xm / 2 + 1; //(xm&1)==xm%2
        g_centroids[blockIdx.x] += (xm + 0.5 - (Npix + 1) / 4);
        g_centroids[blockIdx.x] = (g_centroids[blockIdx.x] - offset) * scale;
      }
      if (threadIdx.x == 1) {
        g_centroids[blockIdx.x + nvalid] = (2.0f * scoeff[0] * scoeff[4]
            - scoeff[3] * scoeff[2]) / denom;
        int ym = (2 * offy + sizey);
        ym = ((ym & 1) == 0) ? ym / 2 : ym / 2 + 1; //(ym&1)==ym%2
        g_centroids[blockIdx.x + nvalid] += (ym + 0.5 - (Npix + 1) / 4);
        g_centroids[blockIdx.x + nvalid] = (g_centroids[blockIdx.x + nvalid]
            - offset) * scale;
      }
    }
  }

  __syncthreads();
  /*
   if (threadIdx.x == 0) g_centroids[blockIdx.x] = (g_centroids[blockIdx.x] - offset) * scale;
   if (threadIdx.x == 1) g_centroids[blockIdx.x+nvalid] = (g_centroids[blockIdx.x+nvalid] - offset) * scale;
   */
}

__global__ void fillweights_krnl(float *d_out, float *weights, int Npix,
    int N) {
  int nim, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix;
    idx = tid - nim * Npix;
    d_out[tid] = weights[idx];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fillcorrcube_krnl(cuFloatComplex *d_out, float *d_in,
    int npix_in, int Npix_in, int npix_out, int Npix_out, int N) {
  int nim, npix, nlig, ncol, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix_in;
    npix = tid - nim * Npix_in;
    nlig = npix / npix_in;
    ncol = npix - nlig * npix_in;
    idx = nlig * npix_out + ncol + nim * Npix_out;
    d_out[idx].x = d_in[tid];
    d_out[idx].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fillcorrim_krnl(cuFloatComplex *d_out, float *d_in, int npix_in,
    int Npix_in, int npix_out, int Npix_out, int N) {
  int nim, npix, nlig, ncol, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix_in;
    npix = tid - nim * Npix_in;
    nlig = npix / npix_in;
    ncol = npix - nlig * npix_in;
    idx = nlig * npix_out + ncol + nim * Npix_out;
    d_out[idx].x = d_in[npix];
    d_out[idx].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void fillval_krnl(cuFloatComplex *d_out, float val, int npix_in,
    int npix_out, int N) {
  int nim, npix, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / npix_in;
    npix = tid % npix_in;
    idx = nim * npix_out + npix;
    d_out[idx].x = val;
    d_out[idx].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void corr_krnl(cuFloatComplex *odata, cuFloatComplex *idata, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ cuFloatComplex tmp;

  while (tid < N) {
    tmp.x = idata[tid].x * odata[tid].x + idata[tid].y * odata[tid].y;
    tmp.y = -1.0f * idata[tid].y * odata[tid].x + idata[tid].x * odata[tid].y;
    odata[tid].x = tmp.x;
    odata[tid].y = tmp.y;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void roll2real_krnl(float *odata, cuFloatComplex *idata, int n,
    int Npix, int N) {
  // here we need to roll and keep only the (2:,2:,) elements
  // startegy is to go through all elements of input array
  // if final index > 0 (in x & y) keep it in output with (idx-1,idy-1)
  // n is always odd because it is 2 x npix
  int nim, idpix, idx, idy, idt;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix;
    idpix = tid - nim * Npix;
    idy = idpix / n;
    idx = idpix - n * idy;

    idx = (idx + n / 2) % n;
    idy = (idy + n / 2) % n;

    if ((idx > 0) && (idy > 0)) {
      idt = idx - 1 + (idy - 1) * (n - 1) + nim * (n - 1) * (n - 1);
      odata[idt] = idata[tid].x;
    }
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void corrnorm_krnl(float *odata, float *idata, int Npix, int N) {

  int nim, idpix;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ float tmp;

  while (tid < N) {
    nim = tid / Npix;
    idpix = tid - nim * Npix;
    tmp = odata[tid] / idata[idpix];
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void convert_krnl(float *odata, float *idata, float offset,
    float scale, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (idata[tid] - offset) * scale;
    tid += blockDim.x * gridDim.x;
  }
}

/*
 _                           _
 | |    __ _ _   _ _ __   ___| |__   ___ _ __ ___
 | |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
 | |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
|_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/

 */
template<class T>
void subap_bpcentro(int threads, int blocks, int npix, T *d_idata,
    unsigned int *values, T *d_odata, T scale, T offset) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * (sizeof(T) + sizeof(uint));

  centroid_bpix<T> <<<dimGrid, dimBlock, smemSize>>>(blocks, npix, d_idata,
      values, d_odata, scale, offset);

  carmaCheckMsg("centroid_bpix<<<>>> execution failed\n");
}
template void
subap_bpcentro<float>(int threads, int blocks, int npix, float *d_idata,
    unsigned int *values, float *d_odata, float scale, float offset);
template void
subap_bpcentro<double>(int threads, int blocks, int npix, double *d_idata,
    unsigned int *values, double *d_odata, double scale, double offset);

template<class T>
void subap_sortmax(int threads, int blocks, T *d_idata, T *d_odata,
    unsigned int *values, int nmax, carma_device *device) {

	int maxThreads = device->get_properties().maxThreadsPerBlock;
	unsigned int nelem_thread = 1;
	while((threads/nelem_thread > maxThreads) || (threads % nelem_thread != 0)){
		nelem_thread++;
	}

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  size_t smemSize = threads * (sizeof(T) + sizeof(uint));
  sortmax<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, values, nmax, threads, threads*blocks, nelem_thread);

  carmaCheckMsg("sortmax_kernel<<<>>> execution failed\n");
}
template void
subap_sortmax<float>(int threads, int blocks, float *d_idata, float *d_odata,
    unsigned int *values, int nmax, carma_device *device);
template void
subap_sortmax<double>(int threads, int blocks, double *d_idata, double *d_odata,
    unsigned int *values, int nmax, carma_device *device);

template<class T>
void subap_centromax(int threads, int blocks, T *d_idata, T *d_odata, int npix,
    int nmax, T scale, T offset) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  /*  int smemSize =
   (threads <= 32) ?
   2 * threads * (sizeof(T) + sizeof(uint)) :
   threads * (sizeof(T) + sizeof(uint));*/
  int smemSize = 2 * threads * (sizeof(T) + sizeof(uint));

  centroid_max<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, npix,
      nmax, blocks, scale, offset);

  carmaCheckMsg("centroid_kernel<<<>>> execution failed\n");
}
template void
subap_centromax<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int npix, int nmax, float scale, float offset);
template void
subap_centromax<double>(int threads, int blocks, double *d_idata,
    double *d_odata, int npix, int nmax, double scale, double offset);

template<class T>
void subap_centromax2(int threads, int blocks, T *d_idata, T *d_odata,
    T *d_minim, int npix, int nmax, T scale, T offset) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  /*  int smemSize =
   (threads <= 32) ?
   2 * threads * (sizeof(T) + sizeof(uint)) :
   threads * (sizeof(T) + sizeof(uint));*/
  int smemSize = 2 * threads * (sizeof(T) + sizeof(uint));

  centroid_max2<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, d_minim,
      npix, nmax, blocks, scale, offset);

  carmaCheckMsg("centroid_kernel<<<>>> execution failed\n");
}
template void
subap_centromax2<float>(int threads, int blocks, float *d_idata, float *d_odata,
    float *d_minim, int npix, int nmax, float scale, float offset);
template void
subap_centromax2<double>(int threads, int blocks, double *d_idata,
    double *d_odata, double *d_minim, int npix, int nmax, double scale,
    double offset);

template<class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
    T *d_odata, T *alpha, T scale, T offset, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) || (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  centroidx<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, alpha, n,
      size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

  centroidy<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
      alpha, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void
get_centroids<float>(int size, int threads, int blocks, int n, float *d_idata,
    float *d_odata, float *alpha, float scale, float offset,
    carma_device *device);

template void
get_centroids<double>(int size, int threads, int blocks, int n, double *d_idata,
    double *d_odata, double *alpha, double scale, double offset,
    carma_device *device);

template<class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
    T *d_odata, T *alpha, T thresh, T scale, T offset, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) || (threads % nelem_thread != 0)) {
    nelem_thread++;
  }
  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  centroidx<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, alpha,
      thresh, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

  centroidy<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
      alpha, thresh, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void
get_centroids<float>(int size, int threads, int blocks, int n, float *d_idata,
    float *d_odata, float *alpha, float thresh, float scale, float offset,
    carma_device *device);

template void
get_centroids<double>(int size, int threads, int blocks, int n, double *d_idata,
    double *d_odata, double *alpha, double thresh, double scale, double offset,
    carma_device *device);

template<class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
    T *d_odata, T *alpha, T *weights, T scale, T offset, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) || (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  centroidx<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, alpha,
      weights, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

  centroidy<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
      alpha, weights, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void
get_centroids<float>(int size, int threads, int blocks, int n, float *d_idata,
    float *d_odata, float *alpha, float *weights, float scale, float offset,
    carma_device *device);

template void
get_centroids<double>(int size, int threads, int blocks, int n, double *d_idata,
    double *d_odata, double *alpha, double *weights, double scale,
    double offset, carma_device *device);

int fillweights(float *d_out, float *d_in, int npix, int N,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fillweights_krnl<<<grid, threads>>>(d_out, d_in, npix * npix, N);
  carmaCheckMsg("<<<fillweights_krnl>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fillcorr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out,
    int N, int nvalid, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  if (nvalid == 1) {
    //cout << "3d" << endl;
    fillcorrcube_krnl<<<grid, threads>>>(d_out, d_in, npix_in,
        npix_in * npix_in, npix_out, npix_out * npix_out, N);
  } else {
    fillcorrim_krnl<<<grid, threads>>>(d_out, d_in, npix_in, npix_in * npix_in,
        npix_out, npix_out * npix_out, N);
    //cout << "2d" << endl;
  }

  carmaCheckMsg("fillcorr_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int correl(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  corr_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  carmaCheckMsg("corr_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int roll2real(float *d_odata, cuFloatComplex *d_idata, int n, int Npix, int N,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  roll2real_krnl<<<grid, threads>>>(d_odata, d_idata, n, Npix, N);

  carmaCheckMsg("roll2real_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int corr_norm(float *d_odata, float *d_idata, int Npix, int N,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  corrnorm_krnl<<<grid, threads>>>(d_odata, d_idata, Npix, N);

  carmaCheckMsg("corrnorm_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out,
    int N, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fillval_krnl<<<grid, threads>>>(d_out, val, npix_in, npix_out, N);
  carmaCheckMsg("fillcorr_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template<class T>
void subap_sortmaxi(int threads, int blocks, T *d_idata, int *values, int nmax,
    int offx, int offy, int npix, int Npix)
// here idata is a [Npix,Npix,nvalid] array
// we want to get the [nvalid] max into subregions of [npix,npix] starting at [xoff,yoff]
// number of threads is npix * npix and number of blocks : nvalid
    {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ?
          2 * threads * (sizeof(T) + sizeof(uint)) :
          threads * (sizeof(T) + sizeof(uint));
  sortmaxi<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, values, nmax, offx,
      offy, npix, Npix, Npix * Npix);

  carmaCheckMsg("sortmaxi_kernel<<<>>> execution failed\n");
}
template void
subap_sortmaxi<float>(int threads, int blocks, float *d_idata, int *values,
    int nmax, int offx, int offy, int npix, int Npix);

template void
subap_sortmaxi<double>(int threads, int blocks, double *d_idata, int *values,
    int nmax, int offx, int offy, int npix, int Npix);

/*
 algorithm for parabolic interpolation
 we do the interpolation on nx x ny points
 we use a (nx * ny, 6) interp matrix
 we thus need nx * ny arrays in shared mem this is the number of threads
 we have nvalid blocks
 */

template<class T>
void subap_pinterp(int threads, int blocks, T *d_idata, int *values,
    T *d_centroids, T *d_matinterp, int sizex, int sizey, int nvalid, int Npix,
    T scale, T offset)
// here idata is a [Npix,Npix,nvalid] array
// we want to get the [nvalid] (x0,y0) into subregions of [sizex,sizey] around gvalue
// number of threads is sizex * sizey  and number of blocks : nvalid
    {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  int shsize = (threads + threads * 6 + 6);
  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds

  int smemSize = (shsize <= 32) ? 2 * shsize * sizeof(T) : shsize * sizeof(T);
  interp_parab<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_centroids,
      values, d_matinterp, sizex, sizey, nvalid, Npix, Npix * Npix, scale,
      offset);

  carmaCheckMsg("sortmaxi_kernel<<<>>> execution failed\n");
}

template void
subap_pinterp<float>(int threads, int blocks, float *d_idata, int *values,
    float *d_centroids, float *d_matinterp, int sizex, int sizey, int nvalid,
    int Npix, float scale, float offset);
/*
 template void subap_pinterp<double>(int threads, int blocks, double *d_idata,  int *values, double *d_centroids,
 double *d_matinterp, int sizex,int sizey, int nvalid, int Npix);
 */

int convert_centro(float *d_odata, float *d_idata, float offset, float scale,
    int N, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  convert_krnl<<<grid, threads>>>(d_odata, d_idata, offset, scale, N);

  carmaCheckMsg("convert_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template<class T>
__global__ void pyrslopes_krnl(T *g_odata, T *g_idata, int *subindx,
    int *subindy, T *subsum, unsigned int ns, unsigned int nvalid,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int i2 = subindx[i] + subindy[i] * ns;
    /*
     g_odata[i] = ((g_idata[i2 + ns * ns] + g_idata[i2 + 3 * ns * ns])
     - (g_idata[i2] + g_idata[i2 + 2 * ns * ns])) / subsum[0];
     g_odata[i + nvalid] = ((g_idata[i2 + 2 * ns * ns]
     + g_idata[i2 + 3 * ns * ns]) - (g_idata[i2] + g_idata[i2 + ns * ns]))
     / subsum[0];
     */
    g_odata[i] = ((g_idata[i2 + ns * ns] + g_idata[i2 + 3 * ns * ns])
        - (g_idata[i2] + g_idata[i2 + 2 * ns * ns])) / subsum[i];
    g_odata[i + nvalid] = ((g_idata[i2 + 2 * ns * ns]
        + g_idata[i2 + 3 * ns * ns]) - (g_idata[i2] + g_idata[i2 + ns * ns]))
        / subsum[i];

  }
}

template<class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
    int ns, int nvalid, int nim, carma_device *device) {
  //cout << "hello cu" << endl;

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrslopes_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
      subsum, ns, nvalid, nim);

  carmaCheckMsg("pyrslopes_kernel<<<>>> execution failed\n");
}

template void
pyr_slopes<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    float *subsum, int ns, int nvalid, int nim, carma_device *device);
template void
pyr_slopes<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
    double *subsum, int ns, int nvalid, int nim, carma_device *device);

template<class T, T fct_sin(T)>
__global__ void pyr2slopes_krnl(T *g_odata, T *g_idata, int *subindx,
    int *subindy, T *subsum, unsigned int ns, unsigned int nvalid, T scale, T valid_thresh) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  T tmp;
  const T cmin(-1);
  const T cmax(1);
  while (i < nvalid) {
    const int iq1 = subindx[i] + subindy[i] * ns;
    const int iq2 = subindx[i+nvalid] + subindy[i+nvalid] * ns;
    const int iq3 = subindx[i+2*nvalid] + subindy[i+2*nvalid] * ns;
    const int iq4 = subindx[i+3*nvalid] + subindy[i+3*nvalid] * ns;
    /*
     * g_odata[i] = ((g_idata[i2 + ns * ns] + g_idata[i2 + 3 * ns * ns])
     * - (g_idata[i2] + g_idata[i2 + 2 * ns * ns])) / subsum[0];
     * g_odata[i + nvalid] = ((g_idata[i2 + 2 * ns * ns]
     * + g_idata[i2 + 3 * ns * ns]) - (g_idata[i2] + g_idata[i2 + ns * ns]))
     * / subsum[0];
     */
    if (subsum[i]<valid_thresh) { // flux too low -> set slopes to 9
        g_odata[i] = 0;
        g_odata[i + nvalid] = 0;
    } else {
    	tmp = ((g_idata[iq1] + g_idata[iq4]) - (g_idata[iq2] + g_idata[iq3])) / subsum[i];
        tmp = carma_clip(tmp, cmin, cmax);            // clip unexpected values
    	g_odata[i] = scale*fct_sin(tmp/2.);           // fct_sin calculates the sine of the input argument × π .
        tmp = ((g_idata[iq1] + g_idata[iq3]) - (g_idata[iq2] + g_idata[iq4])) / subsum[i];
        tmp = carma_clip(tmp, cmin, cmax);            // clip unexpected values
        g_odata[i + nvalid] = scale*fct_sin(tmp/2.);  // fct_sin calculates the sine of the input argument × π .
    }
    i += blockDim.x * gridDim.x;
  }
}

template<class T, T fct_sin(T)>
void pyr2_slopes_full(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
    int ns, int nvalid, T scale, T valid_thresh, carma_device *device) {
  //cout << "hello cu" << endl;

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyr2slopes_krnl<T, fct_sin> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
      subsum, ns, nvalid,scale, valid_thresh);

  carmaCheckMsg("pyrslopes_kernel<<<>>> execution failed\n");
}

template<>
void pyr2_slopes<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
                        float *subsum, int ns, int nvalid, float scale,
                        float valid_thresh, carma_device *device){
    pyr2_slopes_full<float, sinpif>(d_odata, d_idata, subindx, subindy, subsum, ns, nvalid, scale, valid_thresh, device);
}
template <>
void pyr2_slopes<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
                         double *subsum, int ns, int nvalid, double scale, double valid_thresh, carma_device *device){
    pyr2_slopes_full<double, sinpi>(d_odata, d_idata, subindx, subindy, subsum, ns, nvalid, scale, valid_thresh, device);
}

////////////////////////////////////////////////////////////
// ADDING PYR_SLOPES MODIFIED FOR ROOF-PRISM: ROOF_SLOPES //
////////////////////////////////////////////////////////////

template<class T>
__global__ void roofslopes_krnl(T *g_odata, T *g_idata, int *subindx,
    int *subindy, T *subsum, unsigned int ns, unsigned int nvalid,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = (g_idata[i2 + ns * ns] - g_idata[i2]) / subsum[i];
    g_odata[i + nvalid] =
        (g_idata[i2 + 3 * ns * ns] - g_idata[i2 + 2 * ns * ns]) / subsum[i];
  }
}

template<class T>
void roof_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
    int ns, int nvalid, int nim, carma_device *device) {
  //cout << "hello cu" << endl;

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  roofslopes_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
      subsum, ns, nvalid, nim);

  carmaCheckMsg("roofslopes_kernel<<<>>> execution failed\n");
}

template void
roof_slopes<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    float *subsum, int ns, int nvalid, int nim, carma_device *device);
template void
roof_slopes<double>(double *d_odata, double *d_idata, int *subindx,
    int *subindy, double *subsum, int ns, int nvalid, int nim,
    carma_device *device);

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
