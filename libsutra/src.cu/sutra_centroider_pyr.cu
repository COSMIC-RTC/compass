// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_centroider_pyr.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderPyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <carma_utils.cuh>
#include "sutra_centroider_utils.cuh"
#include <sutra_centroider_pyr.hpp>

template <class T, T fct_sin(T)>
__global__ void pyr2slopes_krnl(T *g_odata, T *ref, T *g_idata, int32_t *subindx,
                                int32_t *subindy, float *intensities,
                                sutra::SlopesIndex si,
                                uint32_t ns, uint32_t nvalid,
                                float scale, T valid_thresh, int32_t do_sin) {

  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  T tmp;
  const T cmin(-1);
  const T cmax(1);
  while (i < nvalid) {
    const int32_t iq1 = subindx[i] + subindy[i] * ns;
    const int32_t iq2 = subindx[i + nvalid] + subindy[i + nvalid] * ns;
    const int32_t iq3 = subindx[i + 2 * nvalid] + subindy[i + 2 * nvalid] * ns;
    const int32_t iq4 = subindx[i + 3 * nvalid] + subindy[i + 3 * nvalid] * ns;

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
void pyr_slopes_full(T *d_odata, T *ref, T *d_idata, int32_t *subindx,
                      int32_t *subindy, float *intensities, int32_t ns, int32_t nvalid,
                      float scale, T valid_thresh, int32_t do_sin,
                      SlopeOrder slope_order, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  sutra::SlopesIndex si{nvalid, slope_order};

  pyr2slopes_krnl<T, fct_sin>
      <<<grid, threads, 0, stream>>>(d_odata, ref, d_idata, subindx, subindy, intensities,
                          si,
                          ns, nvalid, scale, valid_thresh, do_sin);

  carma_check_msg("pyrslopes_kernel<<<>>> execution failed\n");
}

template <>
void pyr_slopes<float>(float *d_odata, float *ref, float *d_idata,
                        int32_t *subindx, int32_t *subindy, float *intensities, int32_t ns,
                        int32_t nvalid, float scale, float valid_thresh, int32_t do_sin,
                        SlopeOrder slope_order, CarmaDevice *device, cudaStream_t stream) {
  pyr_slopes_full<float, sinpif>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, slope_order, device, stream);
}
template <>
void pyr_slopes<double>(double *d_odata, double *ref, double *d_idata,
                         int32_t *subindx, int32_t *subindy, float *intensities, int32_t ns,
                         int32_t nvalid, float scale, double valid_thresh, int32_t do_sin,
                         SlopeOrder slope_order, CarmaDevice *device, cudaStream_t stream) {
  pyr_slopes_full<double, sinpi>(d_odata, ref, d_idata, subindx, subindy,
                                  intensities, ns, nvalid, scale, valid_thresh,
                                  do_sin, slope_order, device, stream);
}
