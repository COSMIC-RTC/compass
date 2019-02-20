#include <sutra_centroider_pyr.h>
#include <carma_utils.cuh>

template <class T>
__global__ void pyrslopes_krnl(T *g_odata, T *g_idata, int *subindx,
                               int *subindy, float *intensities, unsigned int ns,
                               unsigned int nvalid, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int i1 = subindx[i] + subindy[i] * ns;
    int i2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    int i3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    int i4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;

    g_odata[i] = ((g_idata[i2] + g_idata[i4]) - (g_idata[i1] + g_idata[i3])) /
                 intensities[i];
    g_odata[i + nvalid] =
        ((g_idata[i3] + g_idata[i4]) - (g_idata[i1] + g_idata[i2])) /
        intensities[i];
  }
}

template <class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy,
                float *intensities, int ns, int nvalid, int nim,
                carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrslopes_krnl<T><<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
                                       intensities, ns, nvalid, nim);

  carmaCheckMsg("pyrslopes_kernel<<<>>> execution failed\n");
}

template void pyr_slopes<float>(float *d_odata, float *d_idata, int *subindx,
                                int *subindy, float *intensities, int ns,
                                int nvalid, int nim, carma_device *device);
template void pyr_slopes<double>(double *d_odata, double *d_idata, int *subindx,
                                 int *subindy, float *intensities, int ns,
                                 int nvalid, int nim, carma_device *device);

template <class T, T fct_sin(T)>
__global__ void pyr2slopes_krnl(T *g_odata, T *ref, T *g_idata, int *subindx,
                                int *subindy, float *intensities, unsigned int ns,
                                unsigned int nvalid, float scale,
                                T valid_thresh, int do_sin) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  T tmp;
  const T cmin(-1);
  const T cmax(1);
  while (i < nvalid) {
    const int iq1 = subindx[i] + subindy[i] * ns;
    const int iq2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    const int iq3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    const int iq4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;

    if (intensities[i] < valid_thresh) {  // flux too low -> set slopes to 9
      g_odata[i] = 0;
      g_odata[i + nvalid] = 0;
    } else {
      tmp = ((g_idata[iq1] + g_idata[iq4]) - (g_idata[iq2] + g_idata[iq3])) /
            intensities[i];
      tmp = carma_clip(tmp, cmin, cmax);  // clip unexpected values
      if (do_sin) {
        g_odata[i] =
            scale *
            fct_sin(
                tmp /
                2.);  // fct_sin calculates the sine of the input argument × π .
      } else {
        g_odata[i] = scale * tmp;
      }
      tmp = ((g_idata[iq1] + g_idata[iq3]) - (g_idata[iq2] + g_idata[iq4])) /
            intensities[i];
      tmp = carma_clip(tmp, cmin, cmax);  // clip unexpected values
      if (do_sin) {
        g_odata[i + nvalid] =
            scale *
            fct_sin(
                tmp /
                2.);  // fct_sin calculates the sine of the input argument × π .
      } else {
        g_odata[i + nvalid] = scale * tmp;
      }
    }
    g_odata[i] -= ref[i];
    g_odata[i + nvalid] -= ref[i];
    i += blockDim.x * gridDim.x;
  }
}

template <class T, T fct_sin(T)>
void pyr2_slopes_full(T *d_odata, T *ref, T *d_idata, int *subindx,
                      int *subindy, float *intensities, int ns, int nvalid,
                      float scale, T valid_thresh, int do_sin,
                      carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyr2slopes_krnl<T, fct_sin>
      <<<grid, threads>>>(d_odata, ref, d_idata, subindx, subindy, intensities,
                          ns, nvalid, scale, valid_thresh, do_sin);

  carmaCheckMsg("pyrslopes_kernel<<<>>> execution failed\n");
}

template <>
void pyr2_slopes<float>(float *d_odata, float *ref, float *d_idata,
                        int *subindx, int *subindy, float *intensities, int ns,
                        int nvalid, float scale, float valid_thresh, int do_sin,
                        carma_device *device) {
  pyr2_slopes_full<float, sinpif>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, device);
}
template <>
void pyr2_slopes<double>(double *d_odata, double *ref, double *d_idata,
                         int *subindx, int *subindy, float *intensities,
                         int ns, int nvalid, float scale, double valid_thresh,
                         int do_sin, carma_device *device) {
  pyr2_slopes_full<double, sinpi>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, device);
}
