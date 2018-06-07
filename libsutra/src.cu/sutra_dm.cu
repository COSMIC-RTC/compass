#include <sutra_dm.h>
#include "carma_utils.cuh"

#include <cuda.h>

#include <stdio.h>

texture<float> t_cmdVector;

texture<float> t_influ;
texture<struct tuple_t<float>> t_inData;

texture<int> t_istart;
texture<int> t_npoints;

texture<int> t_influpos;

void bindTexture(float *cmdVector, float *influ, tuple_t<float> *inData,
                 int *istart, int *npoints, int *influpos, int n_npoints,
                 int n_influpos, int n_vector, int n_influ) {
  cudaBindTexture(NULL, t_cmdVector, cmdVector, n_vector);
#ifdef CHEAT_CODE
  cudaBindTexture(NULL, t_influ, influ, n_influ);
#else
  cudaBindTexture(NULL, t_influ, influ, n_influpos);
#endif

  int length = 0;

#if (COMPN == 2)
  length = n_influpos;
#elif (COMPN == 3)
  length = n_npoints * MAXSPOT;
#endif

  cudaBindTexture(NULL, t_inData, inData, length);

  cudaBindTexture(NULL, t_istart, istart, n_npoints + 1);
  cudaBindTexture(NULL, t_npoints, npoints, n_npoints);

  cudaBindTexture(NULL, t_influpos, influpos, n_influpos);
}

void unbindTexture() {
  cudaUnbindTexture(t_influ);
  cudaUnbindTexture(t_inData);

  cudaUnbindTexture(t_istart);
  cudaUnbindTexture(t_npoints);

  cudaUnbindTexture(t_influpos);
}

template <class T>
__device__ void reduce(T *vector, int size, int id) {
  int lVect = size, shift = size;

  while (shift > 1) {
    lVect >>= 1;  // lVect = lVect / 2;
    shift = CEIL(shift, 2);

    if (id < lVect) vector[id] += vector[id + shift];

    __syncthreads();
  }
}

#ifdef TEXTURE

#if (COMPN == 1)

__device__ inline struct doubleint getInterval(const int pos) {
  return {tex1Dfetch(t_istart, pos), tex1Dfetch(t_npoints, pos)};
}

template <class T>
__device__ inline T getData(const int pos, const T *cmdVector) {
  return cmdVector[tex1Dfetch(t_influpos, pos)] * tex1Dfetch(t_influ, pos);
}

#else

__device__ inline struct doubleint getInterval(const int pos) {
  const int start = tex1Dfetch(t_istart, pos);
  return {start, tex1Dfetch(t_istart, pos + 1) - start};
}

template <class T>
__device__ inline T getData(const int pos) {
  const struct tuple_t<T> data = tex1Dfetch(t_influ3, pos);
  return cmdVector[data.pos] * data.data;
}

#endif

#else

__device__ inline struct doubleint getInterval(const int pos,
                                               const int *iStart_t,
                                               const int *nbInflu_t) {
  return {iStart_t[pos], nbInflu_t[pos]};
}

__device__ inline struct doubleint getInterval(const int pos,
                                               const int *iStart_t) {
  const int start = iStart_t[pos];
  return {start, iStart_t[pos + 1] - start};
}

template <class T>
__device__ inline T getData(const int pos, const T *cmdVector,
                            const T *influData, const int *iPos) {
  return cmdVector[iPos[pos]] * influData[pos];
}

template <class T>
__device__ inline T getData(const int pos, const T *cmdVector,
                            const struct tuple_t<T> *inData) {
  const struct tuple_t<T> data = inData[pos];
  return cmdVector[data.pos] * data.data;
}

#endif

#ifdef TEXTURE
template <class T>
__global__ void compShape1(T *outData, const T *cmdVector, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id);

    T sum = 0;

    for (int pos = interval.start; pos < (interval.start + interval.nbInflu);
         ++pos)
      sum += getData(pos, cmdVector);

    outData[id] = sum;
  }
}

template <class T>
__global__ void compShape2(T *outData, const T *cmdVector, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id);

    T sum = 0;

    for (int pos = interval.start; pos < interval.start + interval.nbInflu;
         ++pos)
      sum += getData(pos, cmdVector);

    outData[id] = sum;
  }
}

template <class T>
__global__ void compShape3(T *outData, const T *cmdVector, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x, iStart = id * MAXSPOT;

  if (id < N) {
    T sum = 0;

    for (int pos = iStart; pos < iStart + MAXSPOT; ++pos) {
      sum += getData(pos, cmdVector);
    }

    outData[id] = sum;
  }
}

template <class T>
__global__ void sharedCompShape1(T *outData, const T *cmdVector, const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    const struct doubleint interval = getInterval(bId);

    if (tId < interval.nbInflu)
      tempVector[tId] = getData(tId + interval.start, cmdVector);

    __syncthreads();

    reduce(tempVector, interval.nbInflu, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}

template <class T>
__global__ void sharedCompShape2(T *outData, const T *cmdVector, const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    const struct doubleint interval = getInterval(bId);

    if (tId < interval.nbInflu)
      tempVector[tId] = getData(tId + interval.start, cmdVector);

    __syncthreads();

    reduce(tempVector, interval.nbInflu, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}

template <class T>
__global__ void sharedCompShape3(T *outData, const T *cmdVector, const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    tempVector[tId] = getData(tId + bId * MAXSPOT, cmdVector);

    __syncthreads();

    reduce(tempVector, MAXSPOT, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}

#else

template <class T>
__global__ void compShape1(T *outData, const T *cmdVector, const T *influData,
                           const int *iStart_t, const int *nbInflu_t,
                           const int *iPos, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id, iStart_t, nbInflu_t);

    T sum = 0;

    for (int pos = interval.start; pos < (interval.start + interval.nbInflu);
         ++pos)
      sum += getData(pos, cmdVector, influData, iPos);

    outData[id] = sum;
  }
}

template <class T>
__global__ void compShape2(T *outData, const T *cmdVector,
                           const struct tuple_t<T> *inData, const int *iStart_t,
                           const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id, iStart_t);

    T sum = 0;

    for (int pos = interval.start; pos < interval.start + interval.nbInflu;
         ++pos)
      sum += getData(pos, cmdVector, inData);

    outData[id] = sum;
  }
}

template <class T>
__global__ void compShape3(T *outData, const T *cmdVector,
                           const struct tuple_t<T> *inData, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x, iStart = id * MAXSPOT;

  if (id < N) {
    T sum = 0;

    for (int pos = iStart; pos < iStart + MAXSPOT; ++pos) {
      sum += getData(pos, cmdVector, inData);
    }

    outData[id] = sum;
  }
}

template <class T>
__global__ void sharedCompShape1(T *outData, const T *cmdVector,
                                 const T *influData, const int *iStart_t,
                                 const int *nbInflu_t, const int *iPos,
                                 const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    const struct doubleint interval = getInterval(bId, iStart_t, nbInflu_t);

    if (tId < interval.nbInflu)
      tempVector[tId] =
          getData(tId + interval.start, cmdVector, influData, iPos);

    __syncthreads();

    reduce(tempVector, interval.nbInflu, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}

template <class T>
__global__ void sharedCompShape2(T *outData, const T *cmdVector,
                                 const struct tuple_t<T> *inData,
                                 const int *iStart_t, const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    const struct doubleint interval = getInterval(bId, iStart_t);

    if (tId < interval.nbInflu)
      tempVector[tId] = getData(tId + interval.start, cmdVector, inData);

    __syncthreads();

    reduce(tempVector, interval.nbInflu, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}

template <class T>
__global__ void sharedCompShape3(T *outData, const T *cmdVector,
                                 const struct tuple_t<T> *inData, const int N) {
  const int tId = threadIdx.x, bId = blockIdx.x + (blockIdx.y * gridDim.x);

  T *tempVector = SharedMemory<T>();  // For pixel's value and centroid stockage

  if (bId < N) {
    tempVector[tId] = getData(tId + bId * MAXSPOT, cmdVector, inData);

    __syncthreads();

    reduce(tempVector, MAXSPOT, tId);

    if (tId == 0) outData[bId] = tempVector[0];
  }
}
#endif

#ifdef TEXTURE
template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const struct tuple_t<T> *inData, const int *iStart_t,
                   const int *nbInflu_t, const int *iPos, const int roiLength,
                   const dim3 threads, const dim3 blocks, const int shared) {
#ifdef REDUCTION

#if (COMPN == 1)
  sharedCompShape1<T><<<blocks, threads, shared>>>(outData, roiLength);
#elif (COMPN == 2)
  sharedCompShape2<T><<<blocks, threads, shared>>>(outData, roiLength);
#elif (COMPN == 3)
  sharedCompShape3<T><<<blocks, threads, shared>>>(outData, roiLength);
#endif

#else

#if (COMPN == 1)
  compShape1<T><<<blocks, threads>>>(outData, cmdVector, roiLength);
#elif (COMPN == 2)
  compShape2<T><<<blocks, threads>>>(outData, cmdVector, roiLength);
#elif (COMPN == 3)
  compShape3<T><<<blocks, threads>>>(outData, cmdVector, roiLength);
#endif

#endif

  carmaCheckMsg("comp_dmshape2<<<>>> execution failed\n");
}
#else
template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const struct tuple_t<T> *inData, const int *iStart_t,
                   const int *nbInflu_t, const int *iPos, const int roiLength,
                   const dim3 threads, const dim3 blocks, const int shared) {

#ifdef REDUCTION

#if (COMPN == 1)
  sharedCompShape1<T><<<blocks, threads, shared>>>(
      outData, cmdVector, influData, iStart_t, nbInflu_t, iPos, roiLength);
#elif (COMPN == 2)
  sharedCompShape2<T><<<blocks, threads, shared>>>(outData, cmdVector, inData,
                                                   iStart_t, roiLength);
#elif (COMPN == 3)
  sharedCompShape3<T>
      <<<blocks, threads, shared>>>(outData, cmdVector, inData, roiLength);
#endif

#else

#if (COMPN == 1)
  compShape1<T><<<blocks, threads>>>(outData, cmdVector, influData, iStart_t,
                                     nbInflu_t, iPos, roiLength);
#elif (COMPN == 2)
  compShape2<T>
      <<<blocks, threads>>>(outData, cmdVector, inData, iStart_t, roiLength);
#elif (COMPN == 3)
  compShape3<T><<<blocks, threads>>>(outData, cmdVector, inData, roiLength);
#endif

#endif

  carmaCheckMsg("comp_dmshape2<<<>>> execution failed\n");
}
#endif

template void comp_dmshape2<float>(float *outData, const float *cmdVector,
                                   const float *influData,
                                   const struct tuple_t<float> *inData,
                                   const int *iStart_t, const int *nbInflu_t,
                                   const int *iPos, const int roiLength,
                                   const dim3 threads, const dim3 blocks,
                                   const int shared);

template void comp_dmshape2<double>(double *outData, const double *cmdVector,
                                    const double *influData,
                                    const struct tuple_t<double> *inData,
                                    const int *iStart_t, const int *nbInflu_t,
                                    const int *iPos, const int roiLength,
                                    const dim3 threads, const dim3 blocks,
                                    const int shared);

template <class T>
__global__ void dmshape_krnl(T *g_idata, T *g_odata, int *pos, int *istart,
                             int *npts, T *comm, unsigned int n, int N) {
  // T *sdata = SharedMemory<T>();
  // load shared mem
  // unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int local_istart = istart[i];
    int local_npts = npts[i];

    // sdata[tid] = 0;
    T value = 0;

    if (local_npts > 0) {
      for (int cc = 0; cc < local_npts; cc++) {
        int lpos = pos[local_istart + cc];
        int ninflu = lpos / n;
        // sdata[tid] += comm[ninflu] * g_idata[lpos];
        value += comm[ninflu] * g_idata[lpos];
      }
    }
    g_odata[i] = value;
    i += blockDim.x * gridDim.x;
  }
  /*
  __syncthreads();

  if (i < N) {
  // write result for this block to global mem
  g_odata[i] = sdata[tid];
  }
   */
}

template <class T>
void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos,
                  int *istart, int *npts, T *comm, unsigned int n, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  // int smemSize =
  //    (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  dmshape_krnl<T><<<dimGrid, dimBlock /*, smemSize*/>>>(
      d_idata, d_odata, pos, istart, npts, comm, n, N);

  carmaCheckMsg("dmshape_kernel<<<>>> execution failed\n");
}

template void comp_dmshape<float>(int threads, int blocks, float *d_idata,
                                  float *d_odata, int *pos, int *istart,
                                  int *npts, float *comm, unsigned int n,
                                  int N);

template void comp_dmshape<double>(int threads, int blocks, double *d_idata,
                                   double *d_odata, int *pos, int *istart,
                                   int *npts, double *comm, unsigned int n,
                                   int N);

template <class T>
__global__ void oneactu_krnl(T *g_idata, T *g_odata, int nactu, T ampli,
                             int *xoff, int *yoff, int dim_im, int dim_influ,
                             int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int iy = i / dim_im;
    int ix = i - iy * dim_im;
    int ixactu = ix - xoff[nactu];
    int iyactu = iy - yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_influ) && (iyactu > -1) &&
        (iyactu < dim_influ)) {
      int tid = ixactu + iyactu * dim_influ + nactu * dim_influ * dim_influ;
      g_odata[i] = ampli * g_idata[tid];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
                                  int *xoff, int *yoff, int dim_im,
                                  int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;
    int ixactu = ix + xoff[nactu];
    int iyactu = iy + yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_im) && (iyactu > -1) &&
        (iyactu < dim_im)) {
      int tid = ixactu + iyactu * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * N];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  // oneactu_krnl<T><<< dimGrid, dimBlock >>>(d_idata,d_odata, nactu, ampli,
  // xoff, yoff, dim_im, dim_influ, N);
  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              xoff, yoff, dim_im, dim_influ, N);

  carmaCheckMsg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int threads, int blocks, float *d_idata,
                             float *d_odata, int nactu, float ampli, int *xoff,
                             int *yoff, int dim_im, int dim_influ, int N);

template void oneactu<double>(int threads, int blocks, double *d_idata,
                              double *d_odata, int nactu, double ampli,
                              int *xoff, int *yoff, int dim_im, int dim_influ,
                              int N);

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
                                  int dim_im, int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;

    // write result for this block to global mem
    if ((ix > -1) && (ix < dim_im) && (iy > -1) && (iy < dim_im)) {
      int tid = ix + iy * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * dim_influ * dim_influ];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  // oneactu_krnl<T><<< dimGrid, dimBlock >>>(d_idata,d_odata, nactu, ampli,
  // xoff, yoff, dim_im, dim_influ, N);
  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              dim_im, dim_influ, N);

  carmaCheckMsg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int threads, int blocks, float *d_idata,
                             float *d_odata, int nactu, float ampli, int dim_im,
                             int dim_influ, int N);

template void oneactu<double>(int threads, int blocks, double *d_idata,
                              double *d_odata, int nactu, double ampli,
                              int dim_im, int dim_influ, int N);

template <class T>
__global__ void fulldmshape_krnl(T *g_idata, T *g_odata, int ninflu,
                                 int diminflu, T *comm, int N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // TODO : replace if by while smartly
  if (i < N) {
    sdata[tid] = 0;

    for (int cc = 0; cc < ninflu; cc++) {
      sdata[tid] += comm[cc] * g_idata[i + cc * diminflu];
    }
  }
  __syncthreads();

  if (i < N) {
    // write result for this block to global mem
    g_odata[i] = sdata[tid];
  }
}

template <class T>
void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata,
                      int ninflu, int diminflu, T *comm, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  fulldmshape_krnl<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, ninflu,
                                                       diminflu, comm, N);

  carmaCheckMsg("fulldmshape_kernel<<<>>> execution failed\n");
}

template void comp_fulldmshape<float>(int threads, int blocks, float *d_idata,
                                      float *d_odata, int ninflu, int diminflu,
                                      float *comm, int N);

template void comp_fulldmshape<double>(int threads, int blocks, double *d_idata,
                                       double *d_odata, int ninflu,
                                       int diminflu, double *comm, int N);

template <class T>
__global__ void getIF_krnl(T *IF, float *dmshape, int *indx_pup, long nb_pts,
                           long column, long nb_col) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[indx_pup[tid]];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
__global__ void getIFfull_krnl(T *IF, float *dmshape, long nb_pts, long column,
                               long nb_col) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[tid];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
int getIF(T *IF, float *dmshape, int *indx_pup, long nb_pts, int column,
          long nb_col, int puponly, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, nb_pts, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);
  if (puponly)
    getIF_krnl<T>
        <<<grid, threads>>>(IF, dmshape, indx_pup, nb_pts, column, nb_col);
  else
    getIFfull_krnl<T><<<grid, threads>>>(IF, dmshape, nb_pts, column, nb_col);
  carmaCheckMsg("getIF_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int getIF<float>(float *IF, float *dmshape, int *indx_pup, long nb_pts,
                          int column, long nb_col, int puponly,
                          carma_device *device);
template int getIF<double>(double *IF, float *dmshape, int *indx_pup,
                           long nb_pts, int column, long nb_col, int puponly,
                           carma_device *device);

__global__ void do_statmat_krnl(float *statcov, float *xpos, float *ypos,
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
int dm_dostatmat(float *statcov, long dim, float *xpos, float *ypos, float norm,
                 carma_device *device) {
  int nthreads = 0, nblocks = 0;
  int N = (dim * dim);
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);
  do_statmat_krnl<<<grid, threads>>>(statcov, xpos, ypos, norm, dim, N);
  carmaCheckMsg("do_statcov_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fill_filtermat_krnl(float *filter, int nactu, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    filter[tid] =
        tid % (nactu + 1) ? (float)-1. / nactu : (float)(1. - 1. / nactu);
    tid += blockDim.x * gridDim.x;
  }
}

int fill_filtermat(float *filter, int nactu, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  fill_filtermat_krnl<<<grid, threads>>>(filter, nactu, N);
  carmaCheckMsg("fill_filtmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void multi_krnl(float *i_data, float gain, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid] * gain;
    tid += blockDim.x * gridDim.x;
  }
}

int multi_vect(float *d_data, float gain, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  multi_krnl<<<grid, threads>>>(d_data, gain, N);

  carmaCheckMsg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void fillpos_krnl(int *pos, int *xoff, int *yoff, int size,
                             int nactu, int dim, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    int pix = tid / nactu;
    int actu = tid - pix * nactu;

    if (pix == 0)
      pos[tid] = (xoff[actu] + size / 2) + (yoff[actu] + size / 2) * dim;
    if (pix == 1)
      pos[tid] = (xoff[actu] - 1 + size / 2) + (yoff[actu] + size / 2) * dim;
    if (pix == 2)
      pos[tid] = (xoff[actu] + size / 2) + (yoff[actu] - 1 + size / 2) * dim;
    if (pix == 3)
      pos[tid] =
          (xoff[actu] - 1 + size / 2) + (yoff[actu] - 1 + size / 2) * dim;
    tid += blockDim.x * gridDim.x;
  }
}

int fillpos(int threads, int blocks, int *pos, int *xoff, int *yoff, int size,
            int nactu, int dim, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  fillpos_krnl<<<dimGrid, dimBlock>>>(pos, xoff, yoff, size, nactu, dim, N);

  carmaCheckMsg("fillpos_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillmapactu_krnl(float *mapactu, float *comvec, int *pos,
                                 int nactu, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    int pix = tid / nactu;
    int actu = tid - pix * nactu;

    mapactu[pos[tid]] = comvec[actu] / 4.0f;
    tid += blockDim.x * gridDim.x;
  }
}

int fill_mapactu(int threads, int blocks, float *mapactu, int *pos,
                 float *comvec, int nactu, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  fillmapactu_krnl<<<dimGrid, dimBlock>>>(mapactu, comvec, pos, nactu, N);

  carmaCheckMsg("fillmapactu_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
