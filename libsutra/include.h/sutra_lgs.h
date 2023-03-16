// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_lgs.h
//! \ingroup   libsutra
//! \class     SutraLGS
//! \brief     this class provides the lgs features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _SUTRA_LGS_H_
#define _SUTRA_LGS_H_

#include <carma.h>
#include <carma_obj.h>
#include <sutra_utils.h>

using std::map;
using std::pair;
using std::vector;

class SutraLGS {
 public:
  int device;
  long nvalid;
  long npix;
  long nmaxhr;
  float hg;
  float h0;
  float deltah;
  float pixsize;
  long nprof;

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
            map<vector<int>, cufftHandle *> ftlgskern_plans, long nvalid,
            long npix, long nmaxhr);
  ~SutraLGS();

  int lgs_init(int nprof, float hg, float h0, float deltah, float pixsie,
               float *doffaxis, float *prof1d, float *profcum, float *beam,
               cuFloatComplex *ftbeam, float *azimuth);
  int load_prof(float *prof1d, float *profcum, float hg, float h0,
                float deltah);
  int lgs_update(CarmaDevice *device);
  int lgs_makespot(CarmaDevice *device, int nin);
  int load_kernels(float *h_lgskern, CarmaDevice *device);
};

// General utilities
int interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum,
                int npix, float *doffaxis, float hg, float pixsize, float h0,
                float deltah, int hmax, int Ntot, CarmaDevice *device);
int times_ftbeam(cuFloatComplex *profout, cuFloatComplex *fbeam, int N,
                 int Ntot, CarmaDevice *device);
int roll_beam_exp(float *imout, cuFloatComplex *iprof, float *beam, int N,
                int Ntot, CarmaDevice *device);
int lgs_rotate(cuFloatComplex *odata, float *idata, int width, int height,
               float *theta, float center, int Ntot, CarmaDevice *device);
// int rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
//              cudaArray *d_array, cudaChannelFormatDesc channel_desc, int width,
//              int height, float *theta, float center, int Ntot,
//              CarmaDevice *device);

#endif  // _SUTRA_LGS_H_
