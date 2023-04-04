// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pyr.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderPyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <carma_utils.cuh>
#include "sutra_centroider_utils.cuh"
#include <sutra_centroider_pyr.h>

template <class T>
__global__ void pyrslopes_krnl(T *g_odata, T *g_idata, int *subindx,
                               int *subindy, float *intensities,
                               sutra::SlopesIndex si,
                               unsigned int ns, unsigned int nvalid,
                               unsigned int nim) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nvalid) {
        int i1 = subindx[i] + subindy[i] * ns;
        int i2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
        int i3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
        int i4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;

        g_odata[si.x(i)] = ((g_idata[i2] + g_idata[i4]) - (g_idata[i1] + g_idata[i3])) / intensities[i];
        g_odata[si.y(i)] = ((g_idata[i3] + g_idata[i4]) - (g_idata[i1] + g_idata[i2])) / intensities[i];
    }
}

template <class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy,
                float *intensities, int ns, int nvalid, int nim,
                SlopeOrder slope_order, CarmaDevice *device) {
  // cout << "hello cu" << endl;

    int nb_blocks, nb_threads;
    get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
    dim3 grid(nb_blocks), threads(nb_threads);

    sutra::SlopesIndex si{nvalid, slope_order};


    pyrslopes_krnl<T><<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
                                         intensities, si, ns, nvalid, nim);

    carma_check_msg("pyrslopes_kernel<<<>>> execution failed\n");
}

template void pyr_slopes<float> (float *d_odata, float *d_idata, int *subindx,
                                 int *subindy, float *intensities, int ns,
                                 int nvalid, int nim, SlopeOrder slope_order,
                                 CarmaDevice *device);
template void pyr_slopes<double>(double *d_odata, double *d_idata, int *subindx,
                                 int *subindy, float *intensities, int ns,
                                 int nvalid, int nim, SlopeOrder slope_order,
                                 CarmaDevice *device);

template <class T, T fct_sin(T)>
__global__ void pyr2slopes_krnl(T *g_odata, T *ref, T *g_idata, int *subindx,
                                int *subindy, float *intensities,
                                sutra::SlopesIndex si,
                                unsigned int ns, unsigned int nvalid,
                                float scale, T valid_thresh, int do_sin) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  T tmp;
  const T cmin(-1);
  const T cmax(1);
  while (i < nvalid) {
    const int iq1 = subindx[i] + subindy[i] * ns;
    const int iq2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    const int iq3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    const int iq4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;

    if (intensities[i] < valid_thresh) { // flux too low -> set slopes to 9
        g_odata[si.x(i)] = 0;
        g_odata[si.y(i)] = 0;
    } else {
        tmp = ((g_idata[iq1] + g_idata[iq4]) - (g_idata[iq2] + g_idata[iq3])) /
              intensities[i];
        //tmp = carma_clip(tmp, cmin, cmax); // clip unexpected values
        if (do_sin) {
            g_odata[si.x(i)] = scale * fct_sin( tmp / 2. ); // fct_sin calculates the sine of the input argument × π .
        } else {
            g_odata[si.x(i)] = scale * tmp;
        }
        tmp = ((g_idata[iq1] + g_idata[iq3]) - (g_idata[iq2] + g_idata[iq4])) /
              intensities[i];
        //tmp = carma_clip(tmp, cmin, cmax); // clip unexpected values
        if (do_sin) {
            g_odata[si.y(i)] = scale * fct_sin( tmp / 2. ); // fct_sin calculates the sine of the input argument × π .
        } else {
            g_odata[si.y(i)] = scale * tmp;
        }
    }
    g_odata[si.x(i)] -= ref[si.x(i)];
    g_odata[si.y(i)] -= ref[si.y(i)];
    i += blockDim.x * gridDim.x;
  }
}

template <class T, T fct_sin(T)>
void pyr2_slopes_full(T *d_odata, T *ref, T *d_idata, int *subindx,
                      int *subindy, float *intensities, int ns, int nvalid,
                      float scale, T valid_thresh, int do_sin,
                      SlopeOrder slope_order, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  sutra::SlopesIndex si{nvalid, slope_order};

  pyr2slopes_krnl<T, fct_sin>
      <<<grid, threads>>>(d_odata, ref, d_idata, subindx, subindy, intensities,
                          si,
                          ns, nvalid, scale, valid_thresh, do_sin);

  carma_check_msg("pyrslopes_kernel<<<>>> execution failed\n");
}

template <>
void pyr2_slopes<float>(float *d_odata, float *ref, float *d_idata,
                        int *subindx, int *subindy, float *intensities, int ns,
                        int nvalid, float scale, float valid_thresh, int do_sin,
                        SlopeOrder slope_order, CarmaDevice *device) {
  pyr2_slopes_full<float, sinpif>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, slope_order, device);
}
template <>
void pyr2_slopes<double>(double *d_odata, double *ref, double *d_idata,
                         int *subindx, int *subindy, float *intensities, int ns,
                         int nvalid, float scale, double valid_thresh, int do_sin,
                         SlopeOrder slope_order, CarmaDevice *device) {
  pyr2_slopes_full<double, sinpi>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, slope_order, device);
}
