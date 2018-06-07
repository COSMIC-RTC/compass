#include <sutra_target.h>

// declare texture reference for 2D float texture
texture<float, 2, cudaReadModeElementType> tex2d;

// declare shared memory arrays
// extern __shared__ cuFloatComplex cachec[];
// extern __shared__ float cache[];

__global__ void texraytrace_krnl(float *g_odata, int nx, int ny, float xoff,
                                 float yoff) {
  // calculate normalized texture coordinates
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

  if ((x < nx) && (y < ny)) {
    float xref = x + xoff;
    float yref = y + yoff;
    // read from texture and write to global memory
    g_odata[y * nx + x] += tex2D(tex2d, xref + 0.5, yref + 0.5);
  }
}

int target_texraytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                       int Ny, float xoff, float yoff, int Ntot,
                       cudaChannelFormatDesc channelDesc,
                       carma_device *device) {
  tex2d.addressMode[0] = cudaAddressModeClamp;
  tex2d.addressMode[1] = cudaAddressModeClamp;
  tex2d.filterMode = cudaFilterModeLinear;
  tex2d.normalized = false;

  // size_t offset = 0;
  carmaSafeCall(cudaBindTexture2D(0, tex2d, d_idata, channelDesc, Nx, Ny,
                                  sizeof(float) * Nx));
  /*
   FIX ME !
   code sometimes code crashes here after several launches...
   this is a known alignment issue !
   see below ...
   to solve this d_idata should be allocated via cudaMallocPitch
   need to add this option in carma_obj

   "The problem here is that the memory bound to the 2D texture does not have
   the proper alignment restrictions. Both the base offset of the texture
   memory, and the pitch, have certain HW dependant alignment restrictions.
   However, currently in the CUDA API, we only expose the base offset
   restriction as a device property, and not the pitch restriction.

   The pitch restriction will be addressed in a future CUDA release. Meanwhile,
   it's recommended that apps use cudaMallocPitch() when allocating pitched
   memory, so that the driver takes care of satisfying all restrictions."
   */

  int nBlocks = device->get_properties().multiProcessorCount * 2;

  dim3 dimBlock(nBlocks, nBlocks, 1);
  dim3 dimGrid(nx / dimBlock.x + 1, ny / dimBlock.y + 1, 1);

  texraytrace_krnl<<<dimGrid, dimBlock, 0>>>(d_odata, nx, ny, xoff, yoff);

  carmaCheckMsg("texraytrace_krnl <<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__device__ void generic_raytrace(float *odata, float *idata, int nx, int ny,
                                 float xoff, float yoff, float G, float thetaML,
                                 float dx, float dy, int Nx, int blockSize,
                                 int istart) {
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  y += istart;

  int tido;

  int iref, jref;  // tidi;
  float x1 = ((x - nx / 2.0f) * G);
  float y1 = ((y - ny / 2.0f) * G);

  float x2 = (cosf(thetaML) * x1 - sinf(thetaML) * y1) + nx / 2.0f + xoff;
  float y2 = (sinf(thetaML) * x1 + cosf(thetaML) * y1) + ny / 2.0f + yoff;

  float xref = x2 + dx;
  float yref = y2 + dy;

  float xshift, yshift, wx1, wx2, wy1, wy2;

  iref = (int)xref;
  jref = (int)yref;

  /*
  // Utilisation de la memoire partagee (valable dans le cas ou (Nx >=
  nx+xoff+1) et (Ny=Nx >= ny+yoff+1) )

  tidi = iref + jref * Nx;

  if (tidi < Nx * Nx) {

    // copie des elements idata vers la memeoire partagee (variable cache) du
  bloc,
    // pour cache[0:(blockDim.x+1)-2 ; (blockDim.x+1):2*(blockDim.x+1)-2 ; ... ;
  ((blockDim.y+1)-1)*(blockDim.x+1):(blockDim.y+1)*(blockDim.x+1)-2]
    cache[threadIdx.x + threadIdx.y * (blockDim.x+1)] = idata[tidi];

    // si le thread est le dernier element du bloc suivant l'axe x, on copie des
  elements idata vers la memeoire partagee (variable cache) du bloc,
    // pour cache[(blockDim.x+1)-1 ; 2*(blockDim.x+1)-1 ; ... ;
  blockDim.y*(blockDim.x+1)-1] if(threadIdx.x == blockDim.x-1) cache[threadIdx.x
  + threadIdx.y * (blockDim.x+1) + 1] = idata[tidi+1]; if((threadIdx.y ==
  blockDim.y-1) )
    {
       // si le thread est le dernier element du bloc suivant l'axe y, on copie
  des elements idata vers la memeoire partagee (variable cache) du bloc,
       // pour
  cache[((blockDim.y+1)-1)*(blockDim.x+1):(blockDim.y+1)*(blockDim.x+1)-2]
       cache[threadIdx.x + threadIdx.y * (blockDim.x+1) + (blockDim.x + 1)] =
  idata[tidi+Nx];

       // si le thread est le dernier element du bloc suivant l'axe y, on copie
  de l'element idata vers la memeoire partagee (variable cache) du bloc,
       // pour cache[(blockDim.y+1)*(blockDim.x+1)-1]
       if(threadIdx.x == blockDim.x-1)
          cache[threadIdx.x + threadIdx.y * (blockDim.x+1) + (blockDim.x+1) +
  1]= idata[tidi+Nx+1];
    }
  }

  __syncthreads();*/

  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref - iref;
    yshift = yref - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    /*odata[tido] += (wx1 * wy1 * cache[threadIdx.x + threadIdx.y *
       (blockDim.x+1)]
        + wx2 * wy1 * cache[threadIdx.x + 1 + threadIdx.y * (blockDim.x+1)]
        + wx1 * wy2 * cache[threadIdx.x + (threadIdx.y + 1) * (blockDim.x+1)]
        + wx2 * wy2 * cache[threadIdx.x + 1 + (threadIdx.y + 1) *
       (blockDim.x+1)]);*/

    if ((iref + 1 < Nx) && (jref + 1 < Nx)) {
      odata[tido] += (wx1 * wy1 * idata[iref + jref * Nx] +
                      wx2 * wy1 * idata[iref + 1 + jref * Nx] +
                      wx1 * wy2 * idata[iref + (jref + 1) * Nx] +
                      wx2 * wy2 * idata[iref + 1 + (jref + 1) * Nx]);
    } else {
      odata[tido] += 0.0f;
    }
  }
}

__device__ void lgs_raytrace(float *odata, float *idata, int nx, int ny,
                             float xoff, float yoff, int Nx, int blockSize,
                             int istart, float delta) {
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  y += istart;

  int tido;

  int iref, jref;  //, tidi;
  float xref = x * delta + xoff;
  float yref = y * delta + yoff;

  float xshift, yshift, wx1, wx2, wy1, wy2;

  iref = (int)xref;
  jref = (int)yref;
  // tidi = iref + jref * Nx;

  // ATTENTION : following optimization had beed commented because of being not
  // so efficient...

  // if ((x < nx) && (y < ny)) {
  /*
    if (tidi < Nx * Nx)
      cache[threadIdx.x + threadIdx.y * (blockSize+1)] = idata[tidi];
    //}
    // Remplissage des bords du cache :
    if ((threadIdx.x == blockSize-1) || (threadIdx.y == blockSize-1)){
          if ((threadIdx.x == blockSize-1)){
                  int bx = threadIdx.x + 1;
                  int xx = bx + blockIdx.x * (blockDim.x);
                  float bxref = xx * delta + xoff;
                  int biref = (int) bxref;
                  int tidib = biref + jref * Nx;
                  if (tidib < Nx * Nx)
                          cache[bx + threadIdx.y * (blockSize+1)] =
    idata[tidib]; else cache[bx + threadIdx.y * (blockSize+1)] = 0.0f;
          }
          if ((threadIdx.y == blockSize-1)){
                  int by = threadIdx.y + 1;
                  int yy = by + blockIdx.y * (blockDim.y);
                  float byref = yy * delta + yoff;
                  int bjref = (int) byref;
                  int tidib = iref + bjref * Nx;
                  if (tidib < Nx * Nx)
                          cache[threadIdx.x + by * (blockSize+1)] =
    idata[tidib]; else cache[threadIdx.x + by * (blockSize+1)] = 0.0f;
          }
          if ((threadIdx.x == blockSize-1) && (threadIdx.y == blockSize-1)){
                  int bx = threadIdx.x + 1;
                  int by = threadIdx.y + 1;
                  int xx = bx + blockIdx.x * (blockDim.x);
                  int yy = by + blockIdx.y * (blockDim.y);
                  float bxref = xx * delta + xoff;
                  float byref = yy * delta + yoff;
                  int biref = (int) bxref;
                  int bjref = (int) byref;
                  int tidib = biref + bjref * Nx;
                  if (tidib < Nx * Nx)
                          cache[bx + by * (blockSize+1)] = idata[tidib];
                  else cache[bx + by * (blockSize+1)] = 0.0f;
          }
    }

    __syncthreads();
  */
  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref - iref;
    yshift = yref - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    odata[tido] += (wx1 * wy1 * idata[iref + jref * Nx] +
                    wx2 * wy1 * idata[iref + 1 + jref * Nx] +
                    wx1 * wy2 * idata[iref + (jref + 1) * Nx] +
                    wx2 * wy2 * idata[iref + 1 + (jref + 1) * Nx]);
    /*
        odata[tido] += (wx1 * wy1 * cache[threadIdx.x + threadIdx.y *
       (blockSize+1)]
                + wx2 * wy1 * cache[threadIdx.x + 1 + threadIdx.y *
       (blockSize+1)]
                + wx1 * wy2 * cache[threadIdx.x + (threadIdx.y + 1) *
       (blockSize+1)]
                + wx2 * wy2 * cache[threadIdx.x + 1 + (threadIdx.y + 1) *
       (blockSize+1)]);*/
  }
}

__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,
                              float xoff, float yoff, float G, float thetaML,
                              float dx, float dy, int Nx, int blockSize) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx,
                   blockSize, 0);
}

__global__ void raytrace_lgs_krnl(float *odata, float *idata, int nx, int ny,
                                  float xoff, float yoff, int Nx, int blocksize,
                                  float delta) {
  lgs_raytrace(odata, idata, nx, ny, xoff, yoff, Nx, blocksize, 0, delta);
}

__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,
                              float xoff, float yoff, int Nx, int blockSize,
                              int istart) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, 1.0f, 0.0f, 0.0f, 0.0f, Nx,
                   blockSize, istart);
}

int target_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                    float xoff, float yoff, float G, float thetaML, float dx,
                    float dy, int block_size) {
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = ny + block_size - ny % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  raytrace_krnl<<<blocks, threads, smemSize>>>(
      d_odata, d_idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx, block_size);

  carmaCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int target_lgs_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                        float xoff, float yoff, float delta, int block_size) {
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = ny + block_size - ny % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);
  raytrace_lgs_krnl<<<blocks, threads, smemSize>>>(
      d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, delta);

  carmaCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int target_raytrace_async(carma_streams *streams, float *d_odata,
                          float *d_idata, int nx, int ny, int Nx, float xoff,
                          float yoff, int block_size) {
  int nstreams = streams->get_nbStreams();

  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize, streams->get_stream(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size);

  carmaCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int target_raytrace_async(carma_host_obj<float> *phase_telemetry,
                          float *d_odata, float *d_idata, int nx, int ny,
                          int Nx, float xoff, float yoff, int block_size) {
  float *hdata = phase_telemetry->getData();
  int nstreams = phase_telemetry->get_nbStreams();

  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize,
                    phase_telemetry->get_cudaStream_t(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size);

  carmaCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
  // will only
  //   commence executing when all previous CUDA calls in stream x have
  //   completed
  int delta = block_size * nx;
  for (int i = 0; i < nstreams; i++) {
    int nbcopy = nx * ny - (i + 1) * block_size * nx;
    if (nbcopy > 0) {
      nbcopy = delta;
    } else
      nbcopy = delta + nbcopy;
    cudaMemcpyAsync(&(hdata[i * block_size * nx]),
                    &(d_odata[i * block_size * nx]), sizeof(float) * nbcopy,
                    cudaMemcpyDeviceToHost,
                    phase_telemetry->get_cudaStream_t(i));
  }

  return EXIT_SUCCESS;
}

__global__ void fillamplikrnl(cuFloatComplex *amplipup, float *phase,
                              float *mask, float scale, int puponly, int nx,
                              int Np, int Nx) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Np) {
    int nlinep = tid / nx;
    int ncol = tid - nlinep * nx;
    int nim = ncol + nlinep * Nx;

    if (puponly == 1) {
      amplipup[nim].x = mask[tid];
      amplipup[nim].y = 0.0f;
    } else {
      amplipup[nim].x = cosf(-scale * phase[tid]) * mask[tid];
      amplipup[nim].y = sinf(-scale * phase[tid]) * mask[tid];
    }
    tid += blockDim.x * gridDim.x;
  }
}

int fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask,
                  float scale, int puponly, int nx, int ny, int Nx,
                  carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nx * ny, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fillamplikrnl<<<grid, threads>>>(amplipup, phase, mask, scale, puponly, nx,
                                   nx * ny, Nx);
  carmaCheckMsg("fillamplikrnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

/*
 __global__ void fillampli_krnl(cuFloatComplex *odata, float *idata, float
 *mask, int nx, int ny, int Nx, int blockSize)
 {
 int x = threadIdx.x + blockIdx.x * blockDim.x;
 int y = threadIdx.y + blockIdx.y * blockDim.y;
 //int tid = x + y *blockDim.x * gridDim.x;
 int tid = x + y * nx;

 if ((x < nx) && (y < ny)) {
 (cachec[threadIdx.x + threadIdx.y * blockSize]).x =
 0.0f;//cosf(idata[tid]);//*mask[tid]; (cachec[threadIdx.x + threadIdx.y *
 blockSize]).y =  0.0f;//sinf(idata[tid]);//*mask[tid]; } else {
 (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  0.0f;
 (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
 }

 __syncthreads();

 tid = x + y * Nx;

 if ((x < nx) && (y < ny)) {
 odata[tid] = cachec[threadIdx.x + threadIdx.y * blockSize];
 }
 }



 int fillampli(cuFloatComplex *d_odata,float *d_idata, float *mask,int nx, int
 ny, int Nx, carma_device *device)
 {

 struct cudaDeviceProp deviceProperties = ;
 cudaGetDeviceProperties(&deviceProperties, device);
 // FIX ME !!!!!!!!
 int blockSize = 16;

 int nnx = nx + blockSize - nx%blockSize; // find next multiple of BLOCK_SZ
 int nny = ny + blockSize - ny%blockSize;
 dim3 blocks(nnx/blockSize,nny/blockSize), threads(blockSize,blockSize);

 int smemSize = (blockSize + 1) * (blockSize + 1) * sizeof(float);

 carmaCheckMsg("fillampli_kernel<<<>>> execution failed\n");
 fillampli_krnl<<<blocks, threads, smemSize>>>(d_odata, d_idata, mask, nx,ny,Nx,
 blockSize);

 return EXIT_SUCCESS;
 }

 __global__ void fillpupil_krnl(cuFloatComplex *odata, float *mask, int nx, int
 ny, int Nx, int blockSize)
 {
 int x = threadIdx.x + blockIdx.x * blockDim.x;
 int y = threadIdx.y + blockIdx.y * blockDim.y;
 //int tid = x + y *blockDim.x * gridDim.x;
 int tid = x + y * nx;

 if ((x < nx) && (y < ny)) {
 (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  mask[tid];
 (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
 } else {
 (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  0.0f;
 (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
 }

 __syncthreads();

 tid = x + y * Nx;

 if ((x < nx) && (y < ny)) {
 odata[tid] = cachec[threadIdx.x + threadIdx.y * blockSize];
 }
 }


 int fillpupil(cuFloatComplex *d_odata,float *mask,int nx, int ny, int Nx,
 carma_device *device)
 {
 struct cudaDeviceProp deviceProperties;
 cudaGetDeviceProperties(&deviceProperties, device);
 // FIX ME !!!!!!!!
 int blockSize = 16;

 int nnx = nx + blockSize - nx%blockSize; // find next multiple of BLOCK_SZ
 int nny = ny + blockSize - ny%blockSize;
 dim3 blocks(nnx/blockSize,nny/blockSize), threads(blockSize,blockSize);

 int smemSize = blockSize * blockSize * sizeof(float);

 fillpupil_krnl<<<blocks, threads,smemSize>>>(d_odata, mask,
 nx,ny,Nx,blockSize);

 carmaCheckMsg("fillpupil_kernel<<<>>> execution failed\n");
 return EXIT_SUCCESS;
 }


 __global__ void texraytrace_krnl2( float* g_odata, int nx,  int ny, float
 *xref, float *yref)
 {
 // calculate normalized texture coordinates
 unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
 unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;

 if ((x < nx) && (y < ny))
 // read from texture and write to global memory
 g_odata[y*nx + x] += tex2D(tex,xref[x]+0.5,yref[y]+0.5);
 }


 __global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,float
 xoff, float yoff, int Nx,int blockSize)
 {
 extern __shared__ float cache[];

 int x = threadIdx.x + blockIdx.x * blockDim.x;
 int y = threadIdx.y + blockIdx.y * blockDim.y;
 int tido;

 int iref,jref,tidi;
 int xref = x + xoff;
 int yref = y + xoff;

 float xshift,yshift,wx1,wx2,wy1,wy2;

 if ((x < nx) && (y < ny)) {
 iref = (int)xref;
 jref = (int)yref;
 tidi = iref + jref * Nx;
 cache[threadIdx.x + threadIdx.y * blockSize] =  idata[tidi];
 }

 if ((x == nx-1) && (y < ny-1)) {
 cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
 if (threadIdx.y == blockSize-1) {
 cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
 cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
 }
 }

 if ((x < nx-1) && (y == ny-1)) {
 cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
 if (threadIdx.x == blockSize-1) {
 cache[threadIdx.x + 1 + threadIdx.y * blockSize]   = idata[tidi+1];
 cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
 }
 }

 if ((x == nx-1) && (y == ny-1)) {
 cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
 cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
 cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
 }

 __syncthreads();

 if ((x < nx) && (y < ny)) {
 if ((threadIdx.x == blockSize-1) && (threadIdx.y == blockSize-1)) {
 cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
 cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
 cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
 } else {
 if (threadIdx.x == blockSize-1)
 cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
 if (threadIdx.y == blockSize-1)
 cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
 }
 }
 __syncthreads();

 if ((x < nx) && (y < ny)) {
 tido = x + y * nx;

 xshift = xref - iref;
 yshift = yref - jref;

 wx1 = (1.0f - xshift);
 wx2 = xshift;
 wy1 = (1.0f - yshift);
 wy2 = yshift;

 odata[tido] += (wx1 * wy1 * cache[threadIdx.x + threadIdx.y * blockSize] +
 wx2 * wy1 * cache[threadIdx.x + 1 + threadIdx.y * blockSize] +
 wx1 * wy2 * cache[threadIdx.x + (threadIdx.y+1) * blockSize] +
 wx2 * wy2 * cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize]);
 }
 }



 */
