// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_lgs.h
//! \ingroup   libsutra
//! \class     sutra_lgs
//! \brief     this class provides the lgs features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_LGS_H_
#define _SUTRA_LGS_H_

#include <carma.h>
#include <carma_obj.h>
#include <sutra_utils.h>

using std::map;
using std::pair;
using std::vector;

class sutra_lgs {
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
  carma_obj<float> *d_doffaxis;
  carma_obj<float> *d_azimuth;
  carma_obj<float> *d_prof1d;
  carma_obj<float> *d_profcum;
  carma_obj<cuFloatComplex> *d_prof2d;
  carma_obj<float> *d_beam;
  carma_obj<cuFloatComplex> *d_ftbeam;
  carma_obj<float> *d_lgskern;
  carma_obj<cuFloatComplex> *d_ftlgskern;

  carma_context *current_context;
  /*
   cudaArray                *d_spotarray;
   cudaChannelFormatDesc    channelDesc;
   cudaMemcpy3DParms        copyParams;
   */

 public:
  sutra_lgs(carma_context *context, carma_obj<float> *d_lgskern,
            carma_obj<cuFloatComplex> *d_ftlgskern,
            map<vector<int>, cufftHandle *> ftlgskern_plans, long nvalid,
            long npix, long nmaxhr);
  ~sutra_lgs();

  int lgs_init(int nprof, float hg, float h0, float deltah, float pixsie,
               float *doffaxis, float *prof1d, float *profcum, float *beam,
               cuFloatComplex *ftbeam, float *azimuth);
  int load_prof(float *prof1d, float *profcum, float hg, float h0,
                float deltah);
  int lgs_update(carma_device *device);
  int lgs_makespot(carma_device *device, int nin);
  int load_kernels(float *h_lgskern, carma_device *device);
};

// General utilities
int interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum,
                int npix, float *doffaxis, float hg, float pixsize, float h0,
                float deltah, int hmax, int Ntot, carma_device *device);
int times_ftbeam(cuFloatComplex *profout, cuFloatComplex *fbeam, int N,
                 int Ntot, carma_device *device);
int rollbeamexp(float *imout, cuFloatComplex *iprof, float *beam, int N,
                int Ntot, carma_device *device);
int lgs_rotate(cuFloatComplex *odata, float *idata, int width, int height,
               float *theta, float center, int Ntot, carma_device *device);
int rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
             cudaArray *d_array, cudaChannelFormatDesc channelDesc, int width,
             int height, float *theta, float center, int Ntot,
             carma_device *device);

#endif  // _SUTRA_LGS_H_
