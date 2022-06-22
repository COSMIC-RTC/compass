// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_dm.h
//! \ingroup   libsutra
//! \class     SutraDm
//! \brief     this class provides the dm features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#ifndef _SUTRA_DM_H_
#define _SUTRA_DM_H_

#include <carma_utils.h>
#include <sutra_kl.h>
#include <sutra_phase.h>
#include <sutra_utils.h>
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
  int device;
  string type;
  float altitude;
  long nactus;
  long influsize;
  long dim;
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

  CarmaObj<int> *d_istart;
  CarmaObj<int> *d_npoints;

  CarmaObj<int> *d_influpos;

  // pzt
  CarmaObj<int> *d_xoff;
  CarmaObj<int> *d_yoff;
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
  SutraDm(CarmaContext *context, const char *type, float altitude, long dim,
           long nactus, long influsize, long ninflupos, long n_npoints,
           float push4imat, long nord, float dx, float dy, float thetaML, float G, int device);
  ~SutraDm();

  int nact();
  int pzt_loadarrays(float *influ, int *influpos, int *npoints, int *istart,
                     int *xoff, int *yoff);
  int kl_loadarrays(float *rabas, float *azbas, int *ord, float *cr, float *cp);
  int tt_loadarrays(float *influ);
  int reset_shape();
  int comp_shape();

  int comp_shape(uint16_t *comm);
  int comp_shape(float *comm);
  int comp_oneactu(int nactu, float ampli);
  // Florian features
  int kl_floloadarrays(float *covmat, float *filter, float *evals, float *bas);

  template <class T>
  int get_IF(T *IF, int *indx_pup, long nb_pts, float ampli);
  template <class T>
  int get_IF_sparse(CarmaSparseObj<T> *&d_IFsparse, int *indx_pup,
                    long nb_pts, float ampli, int puponly);

  int do_geomat(float *d_geocov, float *d_IF, long n_pts);

  template <class T>
  int do_geomat_from_sparse(T *d_geocov, CarmaSparseObj<T> *d_IFsparse);

  int DDiago(CarmaObj<float> *d_statcov, CarmaObj<float> *d_geocov);
  int compute_KLbasis(float *xpos, float *ypos, int *indx, long dim, float norm,
                      float ampli);
  int piston_filt(CarmaObj<float> *d_statcov);
  int set_registration(float dx, float dy, float thetaML, float G);
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

  int add_dm(CarmaContext *context, const char *type, float alt, long dim,
             long nactus, long influsize, long ninflupos, long n_npoints,
             float push4imat, long nord, float dx, float dy, float thetaML, float G,  int device);
  int add_dm(CarmaContext *context, const char *type, float alt, long dim,
             long nactus, long influsize, long ninflupos, long n_npoints,
             float push4imat, long nord, int device);
  int insert_dm(CarmaContext *context, const char *type, float alt, long dim,
                long nactus, long influsize, long ninflupos, long n_npoints,
                float push4imat, long nord, float dx, float dy, float thetaML, float G, int device, int idx);
  int remove_dm(int idx);

  int ndm() { return d_dms.size(); };
  int nact_total();
};

template <class T>
void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos,
                  int *istart, int *npts, T *comm, unsigned int n, int N);

template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const int *iStart_t, const int *iPos, const int roiLength,
                   const dim3 threads, const dim3 blocks, const int shared);

template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N);
template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int dim_im, int dim_influ, int N);
template <class T>
void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata,
                      int nactus, int diminflu, T *comm, int N);

template <class T>
int getIF(T *IF, float *dmshape, int *indx_pup, long nb_pts, int column,
          long nb_col, int puponly, CarmaDevice *device);
int dm_dostatmat(float *d_statcov, long Nkl, float *d_xpos, float *d_ypos,
                 float norm, CarmaDevice *device);
int fill_filtermat(float *filter, int nactu, int N, CarmaDevice *device);
int find_nnz(float *d_data, int N, CarmaDevice *device);
int convertToCom(uint16_t *volts, float *com, int N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device);

#endif  // _SUTRA_DM_H_
