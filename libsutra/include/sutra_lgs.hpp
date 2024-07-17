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

//! \file      sutra_lgs.hpp
//! \ingroup   libsutra
//! \class     SutraLGS
//! \brief     this class provides the lgs features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_LGS_H_
#define _SUTRA_LGS_H_

#include <map>
#include <carma.hpp>
#include <carma_obj.hpp>
#include <sutra_utils.hpp>

using std::map;
using std::pair;
using std::vector;

class SutraLGS {
 public:
  int32_t device;
  int64_t nvalid;
  int64_t npix;
  int64_t nmaxhr;
  float hg;
  float h0;
  float deltah;
  float pixsize;
  int64_t nprof;

  cufftHandle *ftlgskern_plan;
  CarmaObj<float> *d_doffaxis;
  CarmaObj<float> *d_azimuth;
  CarmaObj<float> *d_prof1d;
  CarmaObj<float> *d_profcum;
  CarmaObj<cuFloatComplex> *d_prof2d;
  CarmaObj<float> *d_beam;
  CarmaObj<cuFloatComplex> *d_ftbeam;
  CarmaObj<float> *d_lgskern;
  CarmaObj<cuFloatComplex> *d_ftlgskern;

  CarmaContext *current_context;
  /*
   cudaArray                *d_spotarray;
   cudaChannelFormatDesc    channel_desc;
   cudaMemcpy3DParms        copyParams;
   */

 public:
  SutraLGS(CarmaContext *context, CarmaObj<float> *d_lgskern,
            CarmaObj<cuFloatComplex> *d_ftlgskern,
            map<vector<int32_t>, cufftHandle *> ftlgskern_plans, int64_t nvalid,
            int64_t npix, int64_t nmaxhr);
  ~SutraLGS();

  int32_t lgs_init(int32_t nprof, float hg, float h0, float deltah, float pixsie,
               float *doffaxis, float *prof1d, float *profcum, float *beam,
               cuFloatComplex *ftbeam, float *azimuth);
  int32_t load_prof(float *prof1d, float *profcum, float hg, float h0,
                float deltah);
  int32_t lgs_update(CarmaDevice *device);
  int32_t lgs_makespot(CarmaDevice *device, int32_t nin);
  int32_t load_kernels(float *h_lgskern, CarmaDevice *device);
};

// General utilities
int32_t interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum,
                int32_t npix, float *doffaxis, float hg, float pixsize, float h0,
                float deltah, int32_t hmax, int32_t Ntot, CarmaDevice *device);
int32_t times_ftbeam(cuFloatComplex *profout, cuFloatComplex *fbeam, int32_t N,
                 int32_t Ntot, CarmaDevice *device);
int32_t roll_beam_exp(float *imout, cuFloatComplex *iprof, float *beam, int32_t N,
                int32_t Ntot, CarmaDevice *device);
int32_t lgs_rotate(cuFloatComplex *odata, float *idata, int32_t width, int32_t height,
               float *theta, float center, int32_t Ntot, CarmaDevice *device);
// int32_t rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
//              cudaArray *d_array, cudaChannelFormatDesc channel_desc, int32_t width,
//              int32_t height, float *theta, float center, int32_t Ntot,
//              CarmaDevice *device);

#endif  // _SUTRA_LGS_H_
