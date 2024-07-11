// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_dm.hpp
//! \ingroup   libsutra
//! \class     SutraDm
//! \brief     this class provides the dm features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_DM_H_
#define _SUTRA_DM_H_

#include <carma_utils.hpp>
#include <sutra_kl.hpp>
#include <sutra_phase.hpp>
#include <sutra_utils.hpp>
#include <map>

#include <cuda.h>

#define CHEAT_CODE
// #define COMPN 1  // 1, 2 or 3
// //#define REDUCTION
// //#define TEXTURE
#define BLOCKSIZE 512
#define CEIL(a, b) ((a) + (b)-1) / (b)
#define MAXSPOT 16
#define PIXELPERBLOCK 10

using std::string;
using std::vector;

//  ██████╗ ███╗   ███╗
//  ██╔══██╗████╗ ████║
//  ██║  ██║██╔████╔██║
//  ██║  ██║██║╚██╔╝██║
//  ██████╔╝██║ ╚═╝ ██║
//  ╚═════╝ ╚═╝     ╚═╝
//

class SutraDm {
 public:
  int32_t device;
  string type;
  float altitude;
  int64_t nactus;
  int64_t influsize;
  int64_t dim;
  float push4imat;
  float volt_min;
  float volt_max;
  float dx;
  float dy;
  float thetaML;
  float G;

  uint16_t val_max;

  SutraPhase *d_shape;

  CarmaObj<float> *d_com;

  CarmaObj<float> *d_influ;  // if relevant

  CarmaObj<int32_t> *d_istart;
  CarmaObj<int32_t> *d_npoints;

  CarmaObj<int32_t> *d_influpos;

  // pzt
  CarmaObj<int32_t> *d_xoff;
  CarmaObj<int32_t> *d_yoff;
  CarmaObj<float> *d_KLbasis;
  // CarmaSparseObj<float> *d_IFsparse;
  // CarmaObj<float> *d_comdouble;
  // CarmaObj<double> *d_shapedouble;

  SutraKL *d_kl;

  CarmaContext *current_context;
  cublasHandle_t cublas_handle() { return current_context->get_cublas_handle(); }
  cusparseHandle_t cusparse_handle() {
    return current_context->get_cusparse_handle();
  }

 public:
  SutraDm(CarmaContext *context, const char *type, float altitude, int64_t dim,
           int64_t nactus, int64_t influsize, int64_t ninflupos, int64_t n_npoints,
           float push4imat, int64_t nord, float dx, float dy, float thetaML, float G, int32_t device);
  ~SutraDm();

  int32_t nact();
  int32_t pzt_loadarrays(float *influ, int32_t *influpos, int32_t *npoints, int32_t *istart,
                     int32_t *xoff, int32_t *yoff);
  int32_t kl_loadarrays(float *rabas, float *azbas, int32_t *ord, float *cr, float *cp);
  int32_t tt_loadarrays(float *influ);
  int32_t reset_shape();
  int32_t comp_shape();

  int32_t comp_shape(uint16_t *comm);
  int32_t comp_shape(float *comm);
  int32_t comp_oneactu(int32_t nactu, float ampli);
  // Florian features
  int32_t kl_floloadarrays(float *covmat, float *filter, float *evals, float *bas);

  template <class T>
  int32_t get_IF(T *IF, int32_t *indx_pup, int64_t nb_pts, float ampli);
  template <class T>
  int32_t get_IF_sparse(CarmaSparseObj<T> *&d_IFsparse, int32_t *indx_pup,
                    int64_t nb_pts, float ampli, int32_t puponly);

  int32_t do_geomat(float *d_geocov, float *d_IF, int64_t n_pts);

  template <class T>
  int32_t do_geomat_from_sparse(T *d_geocov, CarmaSparseObj<T> *d_IFsparse);

  int32_t DDiago(CarmaObj<float> *d_statcov, CarmaObj<float> *d_geocov);
  int32_t compute_KLbasis(float *xpos, float *ypos, int32_t *indx, int64_t dim, float norm,
                      float ampli);
  int32_t piston_filt(CarmaObj<float> *d_statcov);
  int32_t set_registration(float dx, float dy, float thetaML, float G);
};

//  ██████╗ ███╗   ███╗███████╗
//  ██╔══██╗████╗ ████║██╔════╝
//  ██║  ██║██╔████╔██║███████╗
//  ██║  ██║██║╚██╔╝██║╚════██║
//  ██████╔╝██║ ╚═╝ ██║███████║
//  ╚═════╝ ╚═╝     ╚═╝╚══════╝
//

class SutraDms {
 public:
  vector<SutraDm *> d_dms;

 public:
  SutraDms();
  ~SutraDms();

  int32_t add_dm(CarmaContext *context, const char *type, float alt, int64_t dim,
             int64_t nactus, int64_t influsize, int64_t ninflupos, int64_t n_npoints,
             float push4imat, int64_t nord, float dx, float dy, float thetaML, float G,  int32_t device);
  int32_t add_dm(CarmaContext *context, const char *type, float alt, int64_t dim,
             int64_t nactus, int64_t influsize, int64_t ninflupos, int64_t n_npoints,
             float push4imat, int64_t nord, int32_t device);
  int32_t insert_dm(CarmaContext *context, const char *type, float alt, int64_t dim,
                int64_t nactus, int64_t influsize, int64_t ninflupos, int64_t n_npoints,
                float push4imat, int64_t nord, float dx, float dy, float thetaML, float G, int32_t device, int32_t idx);
  int32_t remove_dm(int32_t idx);

  int32_t ndm() { return d_dms.size(); };
  int32_t nact_total();
};

template <class T>
void comp_dmshape(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t *pos,
                  int32_t *istart, int32_t *npts, T *comm, uint32_t n, int32_t N);

template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const int32_t *iStart_t, const int32_t *iPos, const int32_t roiLength,
                   const dim3 threads, const dim3 blocks, const int32_t shared);

template <class T>
void oneactu(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t nactu,
             T ampli, int32_t *xoff, int32_t *yoff, int32_t dim_im, int32_t dim_influ, int32_t N);
template <class T>
void oneactu(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t nactu,
             T ampli, int32_t dim_im, int32_t dim_influ, int32_t N);
template <class T>
void comp_fulldmshape(int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                      int32_t nactus, int32_t diminflu, T *comm, int32_t N);

template <class T>
int32_t getIF(T *IF, float *dmshape, int32_t *indx_pup, int64_t nb_pts, int32_t column,
          int64_t nb_col, int32_t puponly, CarmaDevice *device);
int32_t dm_dostatmat(float *d_statcov, int64_t Nkl, float *d_xpos, float *d_ypos,
                 float norm, CarmaDevice *device);
int32_t fill_filtermat(float *filter, int32_t nactu, int32_t N, CarmaDevice *device);
int32_t find_nnz(float *d_data, int32_t N, CarmaDevice *device);
int32_t convertToCom(uint16_t *volts, float *com, int32_t N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device);

#endif  // _SUTRA_DM_H_
