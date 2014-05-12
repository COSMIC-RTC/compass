/*
 * prana_rtc.cpp
 *
 *  Created on: May 9, 2014
 *      Author: sevin
 */

#define LOOP_GAIN 0.5f     // -0.1f for all except 40x40 -0.088888889f
#define ALLOC_VERBOSE        0
#if ALLOC_VERBOSE==1
#define ALLOC_TRACE(fmt, args...) printf(fmt, ## args)
#elif ALLOC_VERBOSE==2
#define ALLOC_TRACE(fmt, args...) printf("[%s@%d]:" fmt, __FUNCTION__, __LINE__, ## args)
#else
#define ALLOC_TRACE(fmt, args...) /* */
#endif

#define STAMP(fmt, args...) fprintf(stderr, "[%s@%d]:" fmt, __FUNCTION__, __LINE__, ## args)

#define USE_ASYNC_VERSION   0  //0:sync version, 1:async version, 2:"1StreamPerOps" version
#define FORCE_NB_SHOT       100 //0 to compute all the input
#define NGPU                1
//#define CHECK_RESULTS       0
#define MEMCPY_SLP_ASYNC    0
//#define USE_SYNCHRO         0
#define CTX_DESTROY         1 //bug with the profiler for the 60x60 case
#define USE_MANY_LOOPS      0
#define JITTER_OUTPUT       0

#define CHECK_CU_CALL       1
#define CHECK_CUDA_CALL     1
#define CHECK_CUBLAS_CALL   1
#define CHECK_SVIPC_CALL    1

#if NGPU>1 && !USE_ASYNC_VERSION
#error multiGPU version needs to use code with async version
#endif

#if USE_ASYNC_VERSION==0
#warning this code will not use streams, NSTREAM define to 1
#define NSTREAM             1
#elif USE_ASYNC_VERSION==1
#define NSTREAM             4
#else //USE_ASYNC_VERSION==2
#define NSTREAM             8 // cpy/cube/reduce/centrox/centroy/cmdx/cmdy/axpy
#define NCHUNK              4
#endif //USE_ASYNC_VERSION
#if CHECK_CU_CALL
#define CUCHECK(fct, err)                                                   \
  if ((err = fct) != CUDA_SUCCESS) {                                        \
        printf("ERROR %s@%d code: %d, exiting\n", __FILE__, __LINE__, err); \
        exit(EXIT_FAILURE);                                                 \
      }
#else
#define CUCHECK(fct, err) fct
#endif

#if CHECK_CUDA_CALL
#define CUDACHECK(fct, err)                                                 \
  if ((err = fct) != cudaSuccess) {                                         \
        printf("ERROR %s@%d code: %d, exiting\n", __FILE__, __LINE__, err); \
        exit(EXIT_FAILURE);                                                 \
      }
#else
#define CUDACHECK(fct, err) fct
#endif

#if CHECK_CUBLAS_CALL
#define CUBLASCHECK(fct, err)                                               \
  if ((err = fct) != CUBLAS_STATUS_SUCCESS) {                               \
        printf("ERROR %s@%d code: %d, exiting\n", __FILE__, __LINE__, err); \
        exit(EXIT_FAILURE);                                                 \
      }
#else
#define CUBLASCHECK(fct, err) fct
#endif


#include <prana_rtc.h>
#include <carma_utils.h>
#include <sched.h>
#include <stdio.h>
#include <string.h>
#include <sched.h>

unsigned int nextPow2(unsigned int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}


prana_rtc::prana_rtc(int nGPUs_, CUcontext *ctx_,
    cublasHandle_t cublas_handle_, int nPixImgX_, int nPixImgY_, int nSlp_,
    int nCmd_, int nSppX_, int nSppY_, short *h_iSppValid_, float *h_cmat_) :
    nGPUs(nGPUs_), ctx(ctx_), cublas_handle(cublas_handle_), nPixImgX(
        nPixImgX_), nPixImgY(nPixImgY_), nSlp(nSlp_), nCmd(nCmd_), nSppX(
        nSppX_), nSppY(nSppY_) {

#if USE_ASYNC_VERSION
  const int nStreams = NSTREAM; // number of streams for CUDA calls
#endif
#if USE_ASYNC_VERSION==2
  const int nChunk = NCHUNK; // number of streams for CUDA calls
#endif

  for(int id_gpu=0; id_gpu<nGPUs; id_gpu++){
    ALLOC_TRACE("> Attach context of the device %d\n",id_gpu);
    cuCtxPushCurrent(ctx_[id_gpu]);
  }

  CUCHECK(cuCtxSetCurrent(ctx[0]), err);

  ALLOC_TRACE("> initCUDA loading module: rtc-centroider.cubin\n");
  CUCHECK(cuModuleLoad(&cuModule_centroider, "/home/sevin/rtc-main-4yorick/rtc-centroider.cubin"), err);

  ALLOC_TRACE("> initCUDA loading module: rtc-acquisim.cubin\n");
  CUCHECK(cuModuleLoad(&cuModule_acquisim, "/home/sevin/rtc-main-4yorick/rtc-acquisim.cubin"), err);

  ALLOC_TRACE("> initCUDA loading module: rtc-sync.cubin\n");
  CUCHECK(cuModuleLoad(&cuModule_sync, "/home/sevin/rtc-main-4yorick/rtc-sync.cubin"), err);

#if USE_ASYNC_VERSION
  ALLOC_TRACE("> initCUDA loading function: <reduce2_async> from rtc-centroider.cubin\n");
  CUCHECK(
      cuModuleGetFunction(&reduce2, cuModule_centroider,
          "reduce2_async"), err);
  ALLOC_TRACE("> initCUDA loading function: <centroidx_async> from rtc-centroider.cubin\n");
  CUCHECK(
      cuModuleGetFunction(&centroidx, cuModule_centroider,
          "centroidx_async"), err);
  ALLOC_TRACE("> initCUDA loading function: <centroidy_async> from rtc-centroider.cubin\n");
  CUCHECK(
      cuModuleGetFunction(&centroidy, cuModule_centroider,
          "centroidy_async"), err);
  ALLOC_TRACE("> initCUDA loading function: <bcube_krnl> from rtc-acquisim.cubin\n");
  CUCHECK(
      cuModuleGetFunction(&bcube_krnl, cuModule_acquisim,
          "bcube_krnl"), err);
#else //!USE_ASYNC_VERSION
  ALLOC_TRACE(
      "> initCUDA loading function: <reduce2> from rtc-centroider.cubin\n");
  CUCHECK(cuModuleGetFunction(&reduce2, cuModule_centroider, "reduce2"), err);
  ALLOC_TRACE(
      "> initCUDA loading function: <centroidx> from rtc-centroider.cubin\n");
  CUCHECK(cuModuleGetFunction(&centroidx, cuModule_centroider, "centroidx"),
      err);
  ALLOC_TRACE(
      "> initCUDA loading function: <centroidy> from rtc-centroider.cubin\n");
  CUCHECK(cuModuleGetFunction(&centroidy, cuModule_centroider, "centroidy"),
      err);
  ALLOC_TRACE(
      "> initCUDA loading function: <bcube_krnl> from rtc-acquisim.cubin\n");
  CUCHECK(cuModuleGetFunction(&bcube_krnl, cuModule_acquisim, "bcube_krnl"),
      err);
#endif //USE_ASYNC_VERSION
  int nSppValid = nSlp / 2; // number of subpupils valid

  int nPixSppX = nPixImgX / nSppX; // number of pixels in a Spp in the first axis
  int nPixSppY = nPixImgY / nSppY; // number of pixels in a Spp in the second axis

  // allocate host memory
  size_t nbytes;


  // allocate device memory
  d_image.resize(nGPUs);
  d_cmat.resize(nGPUs);
  d_slopes.resize(nGPUs);
  d_subsum.resize(nGPUs);
  d_com.resize(nGPUs);
//  CUCHECK(cuMemAlloc(d_image, nGPUs*sizeof(CUdeviceptr)), err); // + 0xFFFF);
//  CUCHECK(cuMemAlloc(d_cmat, nGPUs*sizeof(CUdeviceptr)), err); // + 0xFFFF);
//  CUCHECK(cuMemAlloc(d_slopes, nGPUs*sizeof(CUdeviceptr)), err); // + 0xFFFF);
//  CUCHECK(cuMemAlloc(d_subsum, nGPUs*sizeof(CUdeviceptr)), err); // + 0xFFFF);
//  CUCHECK(cuMemAlloc(d_com, nGPUs*sizeof(CUdeviceptr)), err); // + 0xFFFF);

  /* Allocate buffer with extra space for 64kb alignment */
  // TODO: Modifier tous les cuCtxSetCurrent
  cuCtxSetCurrent(ctx[0]);
  nbytes = nPixImgX * nPixImgY * sizeof(float);
  ALLOC_TRACE("> initCUDA cuMemAlloc: d_image_full, %zu bytes\n", nbytes);
  CUCHECK(cuMemAlloc(&d_image_full, nbytes), err); // + 0xFFFF);
  if (d_image_full == 0) {
    printf("cuMemAlloc d_image_full failed\n");
    exit(EXIT_FAILURE);
  }
  CUCHECK(cuMemsetD32(d_image_full, 0, nPixImgX * nPixImgY), err);

  nbytes = nPixSppX * nPixSppY * nSppValid * sizeof(float);
  ALLOC_TRACE("> initCUDA cuMemAlloc: d_image_cube, %zu bytes\n", nbytes);
  CUCHECK(cuMemAlloc(&d_image_cube, nbytes), err); // + 0xFFFF);
  if (d_image_cube == 0) {
    printf("cuMemAlloc d_image_cube failed\n");
    exit(EXIT_FAILURE);
  }
  CUCHECK(cuMemsetD32(d_image_cube, 0, nPixSppX * nPixSppY * nSppValid), err);

  nbytes = nSppValid * 2 * sizeof(short);
  ALLOC_TRACE("> initCUDA cuMemAlloc: d_iSppValid, %zu bytes\n", nbytes);
  CUCHECK(cuMemAlloc(&d_iSppValid, nbytes), err); // + 0xFFFF);
  if (d_iSppValid == 0) {
    printf("cuMemAlloc d_iSppValid failed\n");
    exit(EXIT_FAILURE);
  }
  CUCHECK(cuMemcpyHtoD(d_iSppValid, h_iSppValid_, nbytes), err);

//  STAMP("iSppValid\n");
//  for (int aaa = 0; aaa < nSlp; aaa++)
//    fprintf(stderr, "%d ", hAligned_iSppValid[aaa]);
//  fprintf(stderr, "\n");

  for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
    cuCtxSetCurrent(ctx[id_gpu]);
    if (id_gpu == 0) {
      ALLOC_TRACE(
          "> initCUDA cuMemAlloc: d_image[%d], %d bytes : use d_image_cube\n",
          id_gpu, 0);
      d_image[id_gpu] = d_image_cube;
    } else {
      ALLOC_TRACE("> initCUDA cuMemAlloc: d_image[%d], %zu bytes\n", id_gpu,
          nbytes);
      nbytes = nPixSppX * nPixSppY * nSppValid / nGPUs * sizeof(float);
      CUCHECK(cuMemAlloc(&d_image[id_gpu], nbytes), err); // + 0xFFFF);
      if (d_image[id_gpu] == 0) {
        printf("cuMemAlloc d_image[%d] failed\n", id_gpu);
        exit(EXIT_FAILURE);
      }
      CUCHECK(
          cuMemsetD32(d_image[id_gpu], 0,
              nPixSppX * nPixSppY * nSppValid / nGPUs), err);
    }

    nbytes = nSlp * nCmd / nGPUs * sizeof(float);
    ALLOC_TRACE("> initCUDA cuMemAlloc: d_cmat[%d], %zu bytes\n", id_gpu,
        nbytes);
    CUCHECK(cuMemAlloc(&d_cmat[id_gpu], nbytes), err); // + 0xFFFF);
    if (d_cmat[id_gpu] == 0) {
      printf("cuMemAlloc d_cmat[%d] failed\n", id_gpu);
      exit(EXIT_FAILURE);
    }

    //CUCHECK(cuMemsetD32(d_cmat[id_gpu], 0, nSlp * nCmd / nGPUs), err);
    CUCHECK(
        cuMemcpyHtoD(d_cmat[id_gpu], &h_cmat_[id_gpu * nSlp * nCmd / nGPUs], nbytes),
        err);

    nbytes = nSlp / nGPUs * sizeof(float);
    ALLOC_TRACE("> initCUDA cuMemAlloc: d_slopes[%d], %zu bytes\n", id_gpu,
        nbytes);
    CUCHECK(cuMemAlloc(&d_slopes[id_gpu], nbytes), err); // + 0xFFFF);
    if (d_slopes[id_gpu] == 0) {
      printf("cuMemAlloc d_slopes[%d] failed\n", id_gpu);
      exit(EXIT_FAILURE);
    }
    CUCHECK(cuMemsetD32(d_slopes[id_gpu], 0, nSlp / nGPUs), err);

    nbytes = nSppValid / nGPUs * sizeof(float);
    ALLOC_TRACE("> initCUDA cuMemAlloc: d_subsum[%d], %zu bytes\n", id_gpu,
        nbytes);
    CUCHECK(cuMemAlloc(&d_subsum[id_gpu], nbytes), err); // + 0xFFFF);
    if (d_subsum[id_gpu] == 0) {
      printf("cuMemAlloc d_subsum[%d] failed\n", id_gpu);
      exit(EXIT_FAILURE);
    }
    CUCHECK(cuMemsetD32(d_subsum[id_gpu], 0, nSppValid / nGPUs), err);

    nbytes = nCmd * sizeof(float);
    ALLOC_TRACE("> initCUDA cuMemAlloc: d_com[%d], %zu bytes\n", id_gpu, nbytes);
    CUCHECK(cuMemAlloc(&d_com[id_gpu], nbytes), err); // + 0xFFFF);
    if (d_com[id_gpu] == 0) {
      printf("cuMemAlloc d_com[%d] failed\n", id_gpu);
      exit(EXIT_FAILURE);
    }
    CUCHECK(cuMemsetD32(d_com[id_gpu], 0, nCmd), err);
  }

  CUCHECK(cuCtxSetCurrent(ctx[0]), err);

  nbytes = nCmd * sizeof(float);
  ALLOC_TRACE("> initCUDA cuMemAlloc: d_com_full, %zu bytes\n", nbytes);
  CUCHECK(cuMemAlloc(&d_com_full, nbytes), err); // + 0xFFFF);
  if (d_com_full == 0) {
    printf("cuMemAlloc d_com_full failed\n");
    exit(EXIT_FAILURE);
  }
  CUCHECK(cuMemsetD32(d_com_full, 0, nCmd), err);

  //  cudaEvent_t start, stop;
  //  cudaEventCreate(&start);
  //  cudaEventCreate(&stop);
  CUCHECK(cuEventCreate(&start, CU_EVENT_DEFAULT), err);
  CUCHECK(cuEventCreate(&stop, CU_EVENT_DEFAULT), err);

#if JITTER_OUTPUT
  CUevent start_loop, end_loop;
  CUCHECK(cuEventCreate(&start_loop, CU_EVENT_DEFAULT), err);
  CUCHECK(cuEventCreate(&end_loop, CU_EVENT_DEFAULT), err);
#endif //JITTER_OUTPUT
#if USE_ASYNC_VERSION==1
  CUstream streamsCompute[nGPUs][nStreams];
  CUstream streamsCopy[nGPUs];
  CUevent eventsCompute[nGPUs][nStreams];
  for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
    CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
    CUCHECK(cuStreamCreate(&streamsCopy[id_gpu], CU_STREAM_NON_BLOCKING),
        err);
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      CUCHECK(
          cuEventCreate(&eventsCompute[id_gpu][id_stream],
              CU_EVENT_DISABLE_TIMING), err);
      CUCHECK(
          cuStreamCreate(&streamsCompute[id_gpu][id_stream],
              CU_STREAM_NON_BLOCKING), err);
    }
  }
#elif USE_ASYNC_VERSION==2
  CUevent next_frame;
  CUCHECK(cuEventCreate(&next_frame, CU_EVENT_DEFAULT), err);

  CUstream streamsCompute[nGPUs][nStreams];
  CUstream streamsCopy[nGPUs];
  for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
    CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
    CUCHECK(cuStreamCreate(&streamsCopy[id_gpu], CU_STREAM_NON_BLOCKING), err);
    CUCHECK(cuStreamCreate(&streamsCompute[id_gpu][0], CU_STREAM_NON_BLOCKING), err);
    for (int id_stream = 1; id_stream < nStreams; id_stream++) {
      CUCHECK(cuStreamCreate(&streamsCompute[id_gpu][id_stream], CU_STREAM_NON_BLOCKING), err);
    }
  }
#endif //USE_ASYNC_VERSION

#if JITTER_OUTPUT
  float *laps=(float*)malloc((nb_frames_replayed)*sizeof(float));
#endif //JITTER_OUTPUT
  CUCHECK(cuCtxSetCurrent(ctx[0]), err);
  CUCHECK(cuEventRecord(start, 0), err);

  start_compute = carma_create_barrier(1);
  stop_compute = carma_create_barrier(1);

  running = false;
}

prana_rtc::~prana_rtc() {
  // TODO Auto-generated destructor stub
  CUCHECK(cuEventDestroy(start), err);
  CUCHECK(cuEventDestroy(stop), err);

  ALLOC_TRACE("> unloadCUDA cuMemFree: d_image_full\n");
  CUCHECK(cuMemFree(d_image_full), err);
  ALLOC_TRACE("> unloadCUDA cuMemFree: d_image_cube\n");
  CUCHECK(cuMemFree(d_image_cube), err);
  ALLOC_TRACE("> unloadCUDA cuMemFree: d_iSppValid\n");
  CUCHECK(cuMemFree(d_iSppValid), err);
  for (int gpu = 0; gpu < nGPUs; gpu++) {
    CUCHECK(cuCtxSetCurrent(ctx[gpu]), err);
    if (gpu > 0) {
      ALLOC_TRACE("> unloadCUDA cuMemFree: d_image[%d]\n", gpu);
      CUCHECK(cuMemFree(d_image[gpu]), err);
    }
    ALLOC_TRACE("> unloadCUDA cuMemFree: d_cmat[%d]\n", gpu);
    CUCHECK(cuMemFree(d_cmat[gpu]), err);
    ALLOC_TRACE("> unloadCUDA cuMemFree: d_slopes[%d]\n", gpu);
    CUCHECK(cuMemFree(d_slopes[gpu]), err);
    ALLOC_TRACE("> unloadCUDA cuMemFree: d_subsum[%d]\n", gpu);
    CUCHECK(cuMemFree(d_subsum[gpu]), err);
    ALLOC_TRACE("> unloadCUDA cuMemFree: d_com[%d]\n", gpu);
    CUCHECK(cuMemFree(d_com[gpu]), err);
  }
}

void *prana_rtc::run_helper(void * data){
  return ((prana_rtc *)data)->run();
}

void* prana_rtc::run() {
  CUCHECK(cuCtxSetCurrent(ctx[0]), err);
  size_t nbytes;

  int nSppValid = nSlp / 2; // number of subpupils valid

  int nPixSppX = nPixImgX / nSppX; // number of pixels in a Spp in the first axis
  int nPixSppY = nPixImgY / nSppY; // number of pixels in a Spp in the second axis

  float alpha = -1.0f;
  float beta = 1.0f;

  int nPixSpp = nPixSppX * nPixSppY;
  int nPixRowSpp = nPixSpp * nSppX;
  int nPixValid = nPixSpp * nSppValid;

  int threadsPerBlock = nPixSpp; //number of pixels per subap;
  int blocksPerGrid = nSppValid;
  int nbelem = nPixValid; // threadsPerBlock * blocksPerGrid

// int smemSize = 2 * threadsPerBlock * sizeof(float);
  int smemSize =
      (threadsPerBlock <= 32) ?
          2 * threadsPerBlock * sizeof(float) : threadsPerBlock * sizeof(float);
//  float one = 1., s2 = nPixSppX / 2 + 0.5;
  float one = 1. / 3.4475925989, s2 = nPixSppX / 2 + 0.5; //sqrtf(threadsPerBlock);

  short *d_iSppValidX = (short*) d_iSppValid;
  short *d_iSppValidY = &d_iSppValidX[nSppValid];

  running = true;
  while(running) {
    nbytes = nPixImgX * nPixImgY * sizeof(float);
    //TODO: WARNING IMAGE
    //CUCHECK(cuMemcpyHtoD(d_image_full, hAligned_image, nbytes), err);
    //DEBUG_TRACE("wait\n");
    carma_wait4barrier(&start_compute);
    if(!running) break;
    //DEBUG_TRACE("go\n");

#if JITTER_OUTPUT
    CUCHECK(cuEventRecord(start_loop, 0), err);
#endif //JITTER_OUTPUT
#if USE_ASYNC_VERSION==1 // classic async version 1 stream per part of the problem
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      int offset = id_stream * nPixImgX * nPixImgY / nStreams;
      nbytes = nPixImgX * nPixImgY / nStreams * sizeof(float);

      //TODO: WARNING! same as before
      CUCHECK(
          cuMemcpyHtoDAsync(d_image_full + id_stream * nbytes,
              &h_image[offset],
              nbytes, streamsCompute[0][id_stream]), err);
    }

    for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
      CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
      CUCHECK(
          cuMemsetD32Async(d_com[id_gpu], 0, nCmd,
              streamsCopy[id_gpu]), err);
    }

    //compute...
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      //    fillbinimg(this->d_binimg->getData(), this->d_bincube->getData(), this->npix, this->nvalid,this->npix*this->nxsub,
      //               this->d_validsubsx->getData(),this->d_validsubsy->getData(),(this->noise > 0),this->device);

      //    int fillbincube(float *bimage, float *bcube, int npix, int nsub, int Nsub, short *ivalid, short *jvalid, int nblocks, int nthreads);
      //    int Npix = npix * npix;
      //    int N = Npix * nsub;
      //
      //    dim3 grid(nblocks), threads(nthreads);
      //
      //    bcube_krnl<<<grid, threads>>>(bimage,bcube,npix,Npix,Nsub,ivalid,jvalid,N);
      int index_cube = id_stream * nPixValid / nStreams;
      nbytes = nPixValid * sizeof(float);

      void* args_bcube_krnl[] = {&d_image_full, &d_image_cube, &nPixSppX,
        &nPixSpp, &nPixImgX, &d_iSppValidX, &d_iSppValidY,
        &nPixValid, &index_cube};
      CUCHECK(
          cuLaunchKernel(bcube_krnl_async, blocksPerGrid / nStreams,
              1, 1, threadsPerBlock, 1, 1, 0,
              streamsCompute[0][id_stream], args_bcube_krnl, 0),
          err);

      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        //STAMP("--> stream:%d, gpu:%d\n", id_stream, id_gpu);
        // Invoke kernel
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);

        if (id_gpu > 0) {
          nbytes = nPixValid / nGPUs * sizeof(float);
          CUCHECK(
              cuMemcpyDtoDAsync(d_image[id_gpu],
                  d_image_cube + id_gpu * nbytes, nbytes,
                  streamsCompute[id_gpu][id_stream]), err);
        }

#if USE_MANY_LOOPS
      }
    }

    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //compute_slopes(dimGrid, dimBlock, smemSize, streamsCompute[g], i, d_image[id_gpu], d_subsum[g], threads,
        //               blocks, nstreams[id_gpu], d_slopes[g], nslp, nbelem);
        //void compute_slopes(dim3 dimGrid, dim3 dimBlock, int smemSize, cudaStream_t* streamsCompute, int i, float* d_image, float* d_subsum, int threads,
        //                    int blocks, int nstreams, float* d_slopes, int nslp, int nbelem) {
        //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_subsum, threads * blocks, i * blocks / nstreams);
        //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nstreams);
        //entroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nStreams);

        //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_subsum, threads * blocks, i * blocks / nStreams);
        int index_reduce = id_stream * nSppValid / nStreams / nGPUs;
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
        void* args_reduce2_async[] = {&d_image_ptr, &d_subsum_ptr,
          &nbelem, &index_reduce};
        CUCHECK(
            cuLaunchKernel(reduce2_async, blocksPerGrid / nStreams,
                1, 1, threadsPerBlock, 1, 1, smemSize,
                streamsCompute[id_gpu][id_stream],
                args_reduce2_async, 0), err);
#if USE_MANY_LOOPS
      }
    }
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
        //        sqrtf(threads), i * blocks / nStreams);
        //index_centroidx =  id_stream*nSppValid/nStreams/nGPUs;
        int index_centroidx = id_stream * nSppValid / nStreams / nGPUs;
        CUdeviceptr d_slopesx = d_slopes[id_gpu];
#if USE_MANY_LOOPS
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
#endif
        void* args_centroidx_async[] = {&d_image[id_gpu], &d_slopesx,
          &d_subsum_ptr, &nPixSppX, &nbelem, &one, &s2,
          &index_centroidx};
        CUCHECK(
            cuLaunchKernel(centroidx_async,
                blocksPerGrid / nStreams, 1, 1, threadsPerBlock,
                1, 1, smemSize,
                streamsCompute[id_gpu][id_stream],
                args_centroidx_async, 0), err);
#if USE_MANY_LOOPS
      }
    }
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //centroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
        //        sqrtf(threads), i * blocks / nStreams);
        //index_centroidy = nSppValid/nGPUs;//id_stream * blocksPerGrid / nStreams;
        int index_centroidy = id_stream * nSppValid / nStreams / nGPUs;
        CUdeviceptr d_slopesy = d_slopes[id_gpu]
        + nSppValid * sizeof(float);
#if USE_MANY_LOOPS
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
#endif
        void* args_centroidy_async[] = {&d_image[id_gpu], &d_slopesy,
          &d_subsum_ptr, &nPixSppX, &nbelem, &one, &s2,
          &index_centroidy};
        CUCHECK(
            cuLaunchKernel(centroidy_async,
                blocksPerGrid / nStreams, 1, 1, threadsPerBlock,
                1, 1, smemSize,
                streamsCompute[id_gpu][id_stream],
                args_centroidy_async, 0), err);
#if USE_MANY_LOOPS
      }
    }
#endif
    int sub_size = nSppValid / nGPUs / nStreams;
#if USE_MANY_LOOPS
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
        CUBLASCHECK(
            cublasSetStream(cublas_handle, streamsCompute[id_gpu][id_stream]),
            status_cublas);

        /////////////////////////////////
        // calcul de la contribution aux commandes des blocks de pentes
        //////////////////////////////////
        sub_size = nSppValid / nGPUs / nStreams;
        int istart = id_stream * sub_size;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_cmat = (float*) (d_cmat[id_gpu]);//+ nCmd * nbytes);
        float *sub_slp = (float*) (d_slopes[id_gpu]);// + nbytes);
        float *sub_cmd = (float*) d_com[id_gpu];
        alpha = -1.0f;
        beta = 1.0f;
        CUBLASCHECK(
            cublasSgemv(cublas_handle, CUBLAS_OP_N, nCmd, sub_size, &alpha,
                sub_cmat, nCmd, sub_slp, 1, &beta, sub_cmd, 1), status_cublas);

#if USE_MANY_LOOPS
      }
    }
    for (int id_stream = 0; id_stream < nStreams; id_stream++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        cuCtxSetCurrent(ctx[id_gpu]);
#endif
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
        CUBLASCHECK(
            cublasSetStream(cublas_handle, streamsCompute[id_gpu][id_stream]),
            status_cublas);

#if USE_MANY_LOOPS
        int istart = id_stream * sub_size + nSppValid;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_slp = (float*) (d_slopes[id_gpu] + nbytes);
        float *sub_cmd = (float*) d_com[id_gpu];
#else
        istart += nSppValid;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        sub_slp = (float*) (d_slopes[id_gpu] + nbytes);
        sub_cmd = (float*) d_com[id_gpu];
#endif
        CUBLASCHECK(
            cublasSgemv(cublas_handle, CUBLAS_OP_N, nCmd, sub_size, &alpha,
                sub_cmat, nCmd, sub_slp, 1, &beta, sub_cmd, 1), status_cublas);
      }
    }

    //cuMemsetD32Async done at the start of the loop
    //cuMemsetD32Async(h_com + idx_frame * nCmd * sizeof(float), 0, nCmd, streamsCompute[0][nStreams - 1]);
    //memset(&hAligned_com[idx_frame * nCmd], 0, nCmd*sizeof(float));

    for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
      CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
      CUBLASCHECK(
          cublasSetStream(cublas_handle,
              streamsCompute[id_gpu][nStreams - 1]),
          status_cublas);
      alpha = LOOP_GAIN;
      //Add the previous command vector the new computed value
      CUBLASCHECK(
          cublasSaxpy(cublas_handle, nCmd, &alpha,
              (float* ) d_com[id_gpu], 1,
              hAligned_com, 1), status_cublas);
    }

    /////////////////////////////////
    // copy back des commandes
    //////////////////////////////////
    //cuMemcpyDtoHAsync(hAligned_com, (float*)d_com, ncmd * sizeof(float), cudaMemcpyDeviceToHost, streamsCompute[nstreams-1]));

#elif USE_ASYNC_VERSION==2 // alternative async version 1 stream per operation
    CUCHECK(cuEventRecord(next_frame, streamsCompute[nGPUs-1][5]), err);
    CUCHECK(cuEventSynchronize(next_frame), err);
    //CUCHECK(cuStreamWaitEvent(streamsCompute[0][0], next_frame, 0), err);

    nbytes = nPixImgX * nPixImgY / nChunk * sizeof(float);
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      int offset = id_chunk * nPixImgX * nPixImgY / nChunk;
      //TODO: WARNING IMAGE
      CUCHECK(cuMemcpyHtoDAsync(d_image_full + id_chunk * nbytes, &h_image[offset], nbytes, streamsCompute[0][0]), err);
    }

    for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
      CUCHECK(cuMemsetD32Async(d_com[id_gpu], 0, nCmd, streamsCopy[id_gpu]), err);
    }

    //compute...
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        //fprintf(stderr, "--> stream:%d, gpu:%d\n", id_chunk,id_gpu);
        // Invoke kernel
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);

//    fillbinimg(this->d_binimg->getData(), this->d_bincube->getData(), this->npix, this->nvalid,this->npix*this->nxsub,
//               this->d_validsubsx->getData(),this->d_validsubsy->getData(),(this->noise > 0),this->device);

//    int fillbincube(float *bimage, float *bcube, int npix, int nsub, int Nsub, short *ivalid, short *jvalid, int nblocks, int nthreads);
//    int Npix = npix * npix;
//    int N = Npix * nsub;
//
//    dim3 grid(nblocks), threads(nthreads);
//
//    bcube_krnl<<<grid, threads>>>(bimage,bcube,npix,Npix,Nsub,ivalid,jvalid,N);
        int index_cube = id_chunk * nPixValid / nChunk;
        nbytes = nPixValid / nGPUs * sizeof(float);

        void* args_bcube_krnl[] = {&d_image_full, &d_image_cube, &nPixSppX, &nPixSpp, &nPixImgX, &d_iSppValidX, &d_iSppValidY, &nPixValid, &index_cube};
        CUCHECK(cuLaunchKernel(bcube_krnl_async, blocksPerGrid / nChunk, 1, 1, threadsPerBlock, 1, 1, 0, streamsCompute[id_gpu][1], args_bcube_krnl, 0), err);

#if USE_MANY_LOOPS
      }
    }

    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //compute_slopes(dimGrid, dimBlock, smemSize, streamsCompute[g], i, d_image[id_gpu], d_subsum[g], threads,
        //               blocks, nstreams[id_gpu], d_slopes[g], nslp, nbelem);
        //void compute_slopes(dim3 dimGrid, dim3 dimBlock, int smemSize, cudaStream_t* streamsCompute, int i, float* d_image, float* d_subsum, int threads,
        //                    int blocks, int nstreams, float* d_slopes, int nslp, int nbelem) {
        //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_subsum, threads * blocks, i * blocks / nstreams);
        //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nstreams);
        //entroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nStreams);

        //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_subsum, threads * blocks, i * blocks / nStreams);
        int index_reduce = id_chunk * nSppValid / nChunk / nGPUs;
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
        void* args_reduce2_async[] = {&d_image_ptr, &d_subsum_ptr, &nbelem, &index_reduce};
        CUCHECK(cuLaunchKernel(reduce2_async, blocksPerGrid / nChunk, 1, 1, threadsPerBlock, 1, 1, smemSize, streamsCompute[id_gpu][2], args_reduce2_async, 0), err);

#if USE_MANY_LOOPS
      }
    }
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
        //        sqrtf(threads), i * blocks / nStreams);
        //index_centroidx =  id_chunk*nSppValid/nStreams/nGPUs;
        int index_centroidx = id_chunk * nSppValid / nChunk / nGPUs;
        CUdeviceptr d_slopesx = d_slopes[id_gpu];
#if USE_MANY_LOOPS
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
#endif
        void* args_centroidx_async[] = {&d_image[id_gpu], &d_slopesx, &d_subsum_ptr, &nPixSppX, &nbelem, &one, &s2, &index_centroidx};
        CUCHECK(cuLaunchKernel(centroidx_async, blocksPerGrid / nChunk, 1, 1, threadsPerBlock, 1, 1, smemSize, streamsCompute[id_gpu][3], args_centroidx_async, 0), err);

#if USE_MANY_LOOPS
      }
    }
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        //centroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
        //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
        //        sqrtf(threads), i * blocks / nStreams);
        //index_centroidy = nSppValid/nGPUs;//id_chunk * blocksPerGrid / nStreams;
        int index_centroidy = id_chunk * nSppValid / nChunk / nGPUs;
        CUdeviceptr d_slopesy = d_slopes[id_gpu] + nSppValid * sizeof(float);
#if USE_MANY_LOOPS
        CUdeviceptr d_image_ptr = d_image[id_gpu];
        CUdeviceptr d_subsum_ptr = d_subsum[id_gpu];
#endif
        void* args_centroidy_async[] = {&d_image[id_gpu], &d_slopesy, &d_subsum_ptr, &nPixSppX, &nbelem, &one, &s2, &index_centroidy};
        CUCHECK(cuLaunchKernel(centroidy_async, blocksPerGrid / nChunk, 1, 1, threadsPerBlock, 1, 1, smemSize, streamsCompute[id_gpu][4], args_centroidy_async, 0), err);

#if USE_MANY_LOOPS
      }
    }
#endif
    int sub_size = nSppValid / nGPUs / nChunk;
#if USE_MANY_LOOPS
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        CUBLASCHECK(cublasSetStream(cublas_handle, streamsCompute[id_gpu][5]), status_cublas);

        /////////////////////////////////
        // calcul de la contribution aux commandes des blocks de pentes
        //////////////////////////////////
        sub_size = nSppValid / nGPUs / nChunk;
        int istart = id_chunk * sub_size;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_slp = (float*) (d_slopes[id_gpu] + nbytes);
        float *sub_cmd = (float*) d_com[id_gpu];
        alpha = -1.0f;
        beta = 1.0f;
        CUBLASCHECK(cublasSgemv(cublas_handle, CUBLAS_OP_N, nCmd, sub_size, &alpha, sub_cmat, nCmd, sub_slp, 1, &beta, sub_cmd, 1), status_cublas);

#if USE_MANY_LOOPS
      }
    }
    for (int id_chunk = 0; id_chunk < nChunk; id_chunk++) {
      for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
        CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
#endif
        CUBLASCHECK(cublasSetStream(cublas_handle, streamsCompute[id_gpu][5]), status_cublas);

#if USE_MANY_LOOPS
        int istart = id_chunk * sub_size + nSppValid;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        float *sub_slp = (float*) (d_slopes[id_gpu] + nbytes);
        float *sub_cmd = (float*) d_com[id_gpu];
#else
        istart += nSppValid;
        nbytes = istart * sizeof(float);
        //float *sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        sub_cmat = (float*) (d_cmat[id_gpu] + nCmd * nbytes);
        sub_slp = (float*) (d_slopes[id_gpu] + nbytes);
        sub_cmd = (float*) d_com[id_gpu];
#endif
        CUBLASCHECK(cublasSgemv(cublas_handle, CUBLAS_OP_N, nCmd, sub_size, &alpha, sub_cmat, nCmd, sub_slp, 1, &beta, sub_cmd, 1), status_cublas);
      }
    }

    //cuMemsetD32Async done at the start of the loop
    //cuMemsetD32Async(h_com + idx_frame * nCmd * sizeof(float), 0, nCmd, streamsCompute[id_gpu][6]);
    //memset(&hAligned_com[idx_frame * nCmd], 0, nCmd*sizeof(float));

    for (int id_gpu = 0; id_gpu < nGPUs; id_gpu++) {
      CUCHECK(cuCtxSetCurrent(ctx[id_gpu]), err);
      CUBLASCHECK(cublasSetStream(cublas_handle, streamsCompute[id_gpu][5]), status_cublas);
      alpha = LOOP_GAIN;
      //Add the previous command vector the new computed value
      CUBLASCHECK(cublasSaxpy(cublas_handle, nCmd, &alpha, (float* ) d_com[id_gpu], 1, hAligned_com, 1), status_cublas);
    }

    /////////////////////////////////
    // copy back des commandes
    //////////////////////////////////
    //cuMemcpyDtoHAsync(hAligned_com, (float*)d_com, ncmd * sizeof(float), cudaMemcpyDeviceToHost, streamsCompute[nstreams-1]));

#else //sync version

//    fillbinimg(this->d_binimg->getData(), this->d_bincube->getData(), this->npix, this->nvalid,this->npix*this->nxsub,
//               this->d_validsubsx->getData(),this->d_validsubsy->getData(),(this->noise > 0),this->device);

//    int fillbincube(float *bimage, float *bcube, int npix, int nsub, int Nsub, short *ivalid, short *jvalid, int nblocks, int nthreads);
//    int Npix = npix * npix;
//    int N = Npix * nsub;
//
//    dim3 grid(nblocks), threads(nthreads);
//
//    bcube_krnl<<<grid, threads>>>(bimage,bcube,npix,Npix,Nsub,ivalid,jvalid,N);

    void* args_bcube_krnl[] = { &d_image_full, &d_image_cube, &nPixSppX,
        &nPixSpp, &nPixImgX, &d_iSppValidX, &d_iSppValidY, &nPixValid };
    CUCHECK(
        cuLaunchKernel(bcube_krnl, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1,
            0, 0, args_bcube_krnl, 0), err);

//    cuMemcpyDtoH((void*) hAligned_image, d_image_full,
//        nPixImgX * nPixImgY * sizeof(float));
//    STAMP("loop d_image %d\ng_image=[ [ %f", id_frame, hAligned_image[0]);
//    for (int aaa = 1; aaa < nPixImgX * nPixImgY; aaa++) {
//
//      if (aaa % nPixImgX == 0)
//        fprintf(stderr, "],\n[%f", hAligned_image[aaa]);
//      else
//        fprintf(stderr, ", %f", hAligned_image[aaa]);
//    }
//    fprintf(stderr, "] ]\n");

//        cuMemcpyDtoH((void*) hAligned_image, d_image_cube,
//            nPixValid * sizeof(float));
//        STAMP("loop d_image_cube %d\ng_image_cube=[ %f", id_frame, hAligned_image[0]);
//        for (int aaa = 1; aaa < nPixValid; aaa++) {
//          fprintf(stderr, ", %f", hAligned_image[aaa]);
//          if (aaa % nPixSpp == 0)
//            fprintf(stderr, "\n");
//        }
//        fprintf(stderr, "] \n");

   nbytes = nPixValid / nGPUs * sizeof(float);
    for (int id_gpu = 1; id_gpu < nGPUs; id_gpu++) {
      //for (int i = 0; i < nStreams; i++) {
      //cuCtxSetCurrent(ctx[id_gpu]);
      CUCHECK(
          cuMemcpyDtoD(d_image[id_gpu], d_image_cube + id_gpu * nbytes, nbytes),
          err);
      //}
    }

    //compute...
    for (int g = 0; g < nGPUs; g++) {
      //fprintf(stderr, "--> stream:%d, gpu:%d\n", i,g);
      // Invoke kernel
      CUCHECK(cuCtxSetCurrent(ctx[g]), err);

      //compute_slopes(dimGrid, dimBlock, smemSize, streamsCompute[g], i, d_image[g], d_subsum[g], threads,
      //               blocks, nstreams[g], d_slopes[g], nslp, nbelem);
      //void compute_slopes(dim3 dimGrid, dim3 dimBlock, int smemSize, cudaStream_t* streamsCompute, int i, float* d_image, float* d_subsum, int threads,
      //                    int blocks, int nstreams, float* d_slopes, int nslp, int nbelem) {
      //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_subsum, threads * blocks, i * blocks / nstreams);
      //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nstreams);
      //entroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(d_image, d_slopes, d_subsum, nslp, nbelem, 1.0f, sqrtf(threads), i * blocks / nStreams);

      //reduce2_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
      //        d_image, d_subsum, threads * blocks, i * blocks / nStreams);

      void* args_reduce2[] = { &d_image[g], &d_subsum[g], &nbelem };
      CUCHECK(
          cuLaunchKernel(reduce2, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1,
              smemSize, 0, args_reduce2, 0), err);
//        void* args_reduce2[] = { &d_image[g], &d_subsum[g], &nbelem };
//        CUCHECK(cuLaunchKernel(reduce2, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1, smemSize, 0, args_reduce2, 0), err);

//      cuMemcpyDtoH((void*) hAligned_image, d_subsum[g],
//          nSppValid * sizeof(float));
//      STAMP("loop d_subsum %d\ng_subsum=[ %f", id_frame, hAligned_image[0]);
//      for (int aaa = 1; aaa < nSppValid; aaa++)
//        fprintf(stderr, ", %f", hAligned_image[aaa]);
//      fprintf(stderr, "]\n");

      //centroidx_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
      //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
      //        sqrtf(threads), i * blocks / nStreams);
      void* args_centroidx[] = { &d_image[g], &d_slopes[g], &d_subsum[g],
          &nPixSppX, &nbelem, &one, &s2 };
      CUCHECK(
          cuLaunchKernel(centroidx, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1,
              smemSize, 0, args_centroidx, 0), err);
//        void* args_centroidx[] = { &d_image[g], &d_slopes[g], &d_subsum[g], &nPixSppX, &nbelem, &one, &s2};
//        CUCHECK(cuLaunchKernel(centroidx,
//            blocksPerGrid, 1, 1, threadsPerBlock, 1, 1, smemSize, 0,
//            args_centroidx, 0), err);

      //centroidy_async<<<dimGrid, dimBlock, smemSize, streamsCompute[i]>>>(
      //        d_image, d_slopes, d_subsum, nSlp, nbelem, 1.0f,
      //        sqrtf(threads), i * blocks / nStreams);
      CUdeviceptr d_slopesy = d_slopes[g] + nSppValid * sizeof(float);
      void* args_centroidy[] = { &d_image[g], &d_slopesy, &d_subsum[g],
          &nPixSppX, &nbelem, &one, &s2 };
      CUCHECK(
          cuLaunchKernel(centroidy, blocksPerGrid, 1, 1, threadsPerBlock, 1, 1,
              smemSize, 0, args_centroidy, 0), err);
//        CUdeviceptr d_slopesy= d_slopes[g]+nSppValid*sizeof(float);
//        void* args_centroidy[] = { &d_image[g], &d_slopesy, &d_subsum[g], &nPixSppX, &nbelem, &one, &s2 };
//        CUCHECK(cuLaunchKernel(centroidy,
//            blocksPerGrid, 1, 1, threadsPerBlock, 1, 1, smemSize, 0,
//            args_centroidy, 0), err);

//      cuMemcpyDtoH((void*) hAligned_image, d_slopes[g],
//          nSlp * sizeof(float));
//      STAMP("loop slopes %d\ng_slp=[ %f", id_frame, hAligned_image[0]);
//      for (int aaa = 1; aaa < nSlp; aaa++)
//        fprintf(stderr, ", %f", hAligned_image[aaa]);
//      fprintf(stderr, "]\n");

      CUCHECK(cuMemsetD32(d_com[g], 0, nCmd), err);

      int sub_size = nSlp / nGPUs;
      nbytes = sub_size * sizeof(float);
      /////////////////////////////////
      // calcul de la contribution aux commandes des blocks de pentes
      //////////////////////////////////
      float *sub_cmat = (float*) d_cmat[g];
      float *sub_slp = (float*) d_slopes[g];
      float *sub_cmd = (float*) d_com[g];
      alpha = -1.0f;
      beta = 1.0f;
      CUBLASCHECK(
          cublasSgemv(cublas_handle, CUBLAS_OP_N, nCmd, sub_size, &alpha, sub_cmat, nCmd, sub_slp, 1, &beta, sub_cmd, 1),
          status_cublas);
    }

    alpha = LOOP_GAIN;
    for (int g = 0; g < nGPUs; g++) {
      CUCHECK(cuCtxSetCurrent(ctx[g]), err);
      CUBLASCHECK(
          cublasSaxpy(cublas_handle, nCmd, &alpha, (float* ) d_com[g], 1,
              (float*)d_com_full, 1), status_cublas);
    }

    /////////////////////////////////
    // copy back des commandes
    //////////////////////////////////
    //cuMemcpyDtoHAsync(hAligned_com, (float*)d_com, ncmd * sizeof(float), cudaMemcpyDeviceToHost, streamsCompute[nstreams-1]));
#endif //USE_ASYNC_VERSION
    //cerr << "done " << idx_frame << endl;
    //DEBUG_TRACE("inc\n");
    carma_increment_barrier(&stop_compute);
    //DEBUG_TRACE("done\n");
    /*
     switch(err) {
     case cudaErrorInvalidValue:
     cout << "cudaMemcpyFromSymbol : cudaErrorInvalidValue" << endl;
     break;
     case cudaErrorInvalidSymbol:
     cout << "cudaMemcpyFromSymbol : cudaErrorInvalidSymbol" << endl;
     break;
     case cudaErrorInvalidDevicePointer:
     cout << "cudaMemcpyFromSymbol : cudaErrorInvalidDevicePointer" << endl;
     break;
     case cudaErrorInvalidMemcpyDirection:
     cout << "cudaMemcpyFromSymbol : cudaErrorInvalidMemcpyDirection" << endl;
     break;
     }
     */

#if JITTER_OUTPUT
    CUCHECK(cuCtxSetCurrent(ctx[0]), err);
    CUCHECK(cuEventRecord(end_loop, 0), err);
    CUCHECK(cuEventSynchronize(end_loop), err);
    CUCHECK(cuEventElapsedTime(&laps[id_frame], start_loop, end_loop), err);
#endif //JITTER_OUTPUT
    //if(newFrame) break;
    //else usleep(100000);
    //cuCtxSynchronize();
  }
  return 0L;
}

void set_affinity(int cpu){
  struct sched_param param_sched;
  param_sched.sched_priority = sched_get_priority_max(SCHED_FIFO);

  if (sched_setscheduler(0, SCHED_FIFO, &param_sched) != 0) {
    printf("soucis de scheduling RT\n");
  }

  // set thread affinity to core 0
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(0, &cpu_set);
  if (sched_setaffinity(0, sizeof(cpu_set_t), &cpu_set) != 0) {
    printf("soucis de setaffinity\n");
  }
}
void prana_rtc::start_rtc() {
  //DEBUG_TRACE("here\n");
  rtc_loop = carma_start_thread(&prana_rtc::run_helper, (void*)this);
  //DEBUG_TRACE("here\n");
}

void prana_rtc::stop_rtc() {
  running = false;
  //DEBUG_TRACE("inc\n");
  carma_increment_barrier(&start_compute);
  //DEBUG_TRACE("end\n");
  carma_end_thread(rtc_loop);
  //DEBUG_TRACE("done\n");
}

void prana_rtc::set_matcom(CUdeviceptr d_matcom_, CUcontext ctx_) {
}

void prana_rtc::set_image(CUdeviceptr d_image_, CUcontext ctx_) {
  CUCHECK(
      cuMemcpyPeer(d_image_full, ctx[0], d_image_, ctx_, nPixImgX*nPixImgY*sizeof(float)),
      err);
  carma_increment_barrier(&start_compute);
}

void prana_rtc::get_commands(CUdeviceptr d_commands_, CUcontext ctx_) {
  carma_wait4barrier(&stop_compute);
  CUCHECK(
      cuMemcpyPeer(d_commands_, ctx_, d_com_full, ctx[0], nCmd*sizeof(float)),
      err);
}
