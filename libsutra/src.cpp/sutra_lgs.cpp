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

//! \file      sutra_lgs.cpp
//! \ingroup   libsutra
//! \class     SutraLGS
//! \brief     this class provides the lgs features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_lgs.hpp>
#include <sutra_utils.hpp>
#include <sutra_wfs.hpp>

SutraLGS::SutraLGS(CarmaContext *context, CarmaObj<float> *d_lgskern,
                     CarmaObj<cuFloatComplex> *d_ftlgskern,
                     map<vector<int32_t>, cufftHandle *> ftlgskern_plans,
                     int64_t nvalid, int64_t npix, int64_t nmaxhr) {
  this->current_context = context;
  this->device = current_context->get_active_device();

  this->nvalid = nvalid;
  this->npix = npix;
  this->nprof = 0;
  this->nmaxhr = nmaxhr;
  this->d_doffaxis = this->d_prof1d = this->d_profcum = 0;
  this->hg = this->pixsize = this->h0 = this->deltah = 0;

  int64_t *dims_data1 = new int64_t[2];
  dims_data1[0] = 1;
  int64_t *dims_data2 = new int64_t[3];
  dims_data2[0] = 2;
  int64_t *dims_data3 = new int64_t[4];
  dims_data3[0] = 3;

  dims_data1[1] = npix;
  this->d_beam = new CarmaObj<float>(context, dims_data1);
  this->d_ftbeam = new CarmaObj<cuFloatComplex>(context, dims_data1);

  dims_data1[1] = nvalid;
  this->d_azimuth = new CarmaObj<float>(context, dims_data1);

  dims_data2[1] = npix;
  dims_data2[2] = nvalid;
  this->d_prof2d = new CarmaObj<cuFloatComplex>(context, dims_data2);

  int32_t mdims[2];
  mdims[0] = (int32_t)dims_data2[1];

  cufftHandle *plan = this->d_prof2d->get_plan();  ///< FFT plan
  carmafft_safe_call(cufftPlanMany(plan, 1, mdims, NULL, 1, 0, NULL, 1, 0,
                                 CUFFT_C2C, (int32_t)dims_data2[2]));

  dims_data3[1] = npix;
  dims_data3[2] = npix;
  dims_data3[3] = nmaxhr;
  // this->d_lgskern = new CarmaObj<float>(context, dims_data3);
  // this->d_ftlgskern = new CarmaObj<cuFloatComplex>(context, dims_data3);
  this->d_lgskern = d_lgskern;
  this->d_ftlgskern = d_ftlgskern;

  mdims[0] = (int32_t)dims_data3[1];
  mdims[1] = (int32_t)dims_data3[2];

  vector<int32_t> vdims(dims_data3 + 1, dims_data3 + 4);

  if (ftlgskern_plans.find(vdims) == ftlgskern_plans.end()) {
    // DEBUG_TRACE("Creating FFT plan : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);print_mem_info();
    cufftHandle *plan = (cufftHandle *)malloc(
        sizeof(cufftHandle));  // = this->d_camplipup->get_plan(); ///< FFT plan
    carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_C2C, (int32_t)dims_data3[3]));

    ftlgskern_plans.insert(pair<vector<int32_t>, cufftHandle *>(vdims, plan));

    this->ftlgskern_plan = plan;
    // DEBUG_TRACE("FFT plan created");print_mem_info();
  } else {
    this->ftlgskern_plan = ftlgskern_plans.at(vdims);
    // DEBUG_TRACE("FFT plan already exists : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);
  }

  // plan = this->d_ftlgskern->get_plan();
  // carmafft_safe_call(
  //     cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
  //     (int32_t)dims_data3[3]));

  /*
   cudaExtent volumeSize = make_cudaExtent(npix,npix,nvalid);

   this->channel_desc = cudaCreateChannelDesc(32, 0, 0, 0,
   cudaChannelFormatKindFloat);
   carma_safe_call( cudaMalloc3DArray(&(this->d_spotarray), &(this->channel_desc),
   volumeSize) );

   // prepare 3d cppy
   this->copyParams.srcPtr =
   make_cudaPitchedPtr((void*)(this->d_lgskern->get_data()),
   volumeSize.width*sizeof(float),
   volumeSize.width, volumeSize.height);
   copyParams.dstArray = this->d_spotarray;
   copyParams.extent   = volumeSize;
   copyParams.kind     = cudaMemcpyDeviceToDevice;
   */
  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;
}

SutraLGS::~SutraLGS() {
  current_context->set_active_device(device, 1);
  delete this->d_doffaxis;
  delete this->d_azimuth;
  delete this->d_prof1d;
  delete this->d_profcum;
  delete this->d_prof2d;
  delete this->d_beam;
  delete this->d_ftbeam;
  // delete this->d_lgskern;
  // delete this->d_ftlgskern;

  // delete this->current_context;

  // carma_safe_call(cudaFreeArray(this->d_spotarray));
}

int32_t SutraLGS::lgs_init(int32_t nprof, float hg, float h0, float deltah,
                        float pixsize, float *doffaxis, float *prof1d,
                        float *profcum, float *beam, cuFloatComplex *ftbeam,
                        float *azimuth) {
  current_context->set_active_device(device, 1);
  this->nprof = nprof;
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;
  this->pixsize = pixsize;

  int64_t *dims_data1 = new int64_t[2];
  dims_data1[0] = 1;
  dims_data1[1] = this->nvalid;
  this->d_doffaxis = new CarmaObj<float>(current_context, dims_data1);
  this->d_doffaxis->host2device(doffaxis);

  dims_data1[1] = nprof + 1;
  this->d_prof1d = new CarmaObj<float>(current_context, dims_data1);
  this->d_profcum = new CarmaObj<float>(current_context, dims_data1);

  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->d_beam->host2device(beam);
  this->d_ftbeam->host2device(ftbeam);
  this->d_azimuth->host2device(azimuth);

  delete[] dims_data1;
  this->lgs_update(current_context->get_device(this->device));
  this->lgs_makespot(current_context->get_device(this->device), 0);

  return EXIT_SUCCESS;
}

int32_t SutraLGS::load_prof(float *prof1d, float *profcum, float hg, float h0,
                         float deltah) {
  current_context->set_active_device(device, 1);
  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;

  return EXIT_SUCCESS;
}

int32_t SutraLGS::lgs_update(CarmaDevice *device) {
  interp_prof(this->d_prof2d->get_data(), this->d_prof1d->get_data(),
              this->d_profcum->get_data(), this->npix,
              this->d_doffaxis->get_data(), this->hg, this->pixsize, this->h0,
              this->deltah, this->nprof, this->d_prof2d->get_nb_elements(), device);

  // convolution by beam
  // do fft on prof2d
  CarmaFFT(this->d_prof2d->get_data(), this->d_prof2d->get_data(), 1,
            *this->d_prof2d->get_plan());

  // mult by beamft
  times_ftbeam(this->d_prof2d->get_data(), this->d_ftbeam->get_data(), this->npix,
               this->d_prof2d->get_nb_elements(), device);

  // fft back
  CarmaFFT(this->d_prof2d->get_data(), this->d_prof2d->get_data(), -1,
            *this->d_prof2d->get_plan());

  return EXIT_SUCCESS;
}

int32_t SutraLGS::lgs_makespot(CarmaDevice *device, int32_t nin) {
  carma_safe_call(cudaMemset(this->d_lgskern->get_data(), 0,
                           sizeof(float) * this->d_lgskern->get_nb_elements()));
  // build final image
  // get abs of real and roll
  cuFloatComplex *data = this->d_prof2d->get_data();
  roll_beam_exp(
      this->d_lgskern->get_data(), &(data[nin]), this->d_beam->get_data(),
      this->npix,
      this->npix * this->npix * this->nmaxhr /*this->d_lgskern->get_nb_elements()*/,
      device);

  // rotate image and fill kernels ft
  float *data2 = this->d_azimuth->get_data();
  lgs_rotate(this->d_ftlgskern->get_data(), this->d_lgskern->get_data(),
             this->npix, this->npix, &(data2[nin / this->npix]), 0.0f,
             this->npix * this->npix * this->nmaxhr, device);

  // same with textures
  /*
   rotate3d(this->d_ftlgskern->get_data(),this->copyParams,this->d_spotarray,this->channel_desc,
   this->npix,this->npix,this->d_azimuth->get_data(),0.0f,this->npix*this->npix*this->nvalid,device);
   */

  // prepare for wfs code
  CarmaFFT(this->d_ftlgskern->get_data(), this->d_ftlgskern->get_data(), 1,
            *this->ftlgskern_plan);

  return EXIT_SUCCESS;
}

int32_t SutraLGS::load_kernels(float *h_lgskern, CarmaDevice *device) {
  this->d_lgskern->host2device(h_lgskern);
  cfillrealp(
      this->d_ftlgskern->get_data(), this->d_lgskern->get_data(),
      /*this->d_ftlgskern->get_nb_elements()*/ this->npix * this->npix * this->nmaxhr,
      device);
  CarmaFFT(this->d_ftlgskern->get_data(), this->d_ftlgskern->get_data(), 1,
            *this->ftlgskern_plan);

  return EXIT_SUCCESS;
}
