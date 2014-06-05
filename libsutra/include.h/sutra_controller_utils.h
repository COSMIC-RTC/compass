#ifndef SUTRA_CONTROLLER_UTILS_H_
#define SUTRA_CONTROLLER_UTILS_H_

#include <sutra_turbu.h>
#include <sutra_wfs.h>
#include <cstdio>

struct gtomo_struct {
  float DiamTel;
  float obs;
  float pasDPHI;
  long Nw;
  long Nx;
  long *Nssp;
  float *diamPup;
  float *XPup;
  float *YPup;
  float *thetaML;
  float *alphaX;
  float *alphaY;
  float *GsAlt;

  float lgs_cst;
  float spot_width;
  float lgs_depth;
  float lgs_alt;
  int    nlgs;

  long   *ioff_d;
  long   *Nssp_d;

  float *alphaX_d;
  float *alphaY_d;
  float *GsAlt_d;
  float *diamPup_d;
  float *thetaML_d;
  float *X_d;
  float *Y_d;
  float *XPup_d;
  float *YPup_d;

  long    max_Nl0;
  long   *indexL0_d;
  float *L0diff_d;
  float *h_d;
  float *cn2_d;
  float *tabDPHI_d;
  float *u_d;
  float *v_d;
  float *sspSizeL_d;

  cudaStream_t matcov_stream;
};

void process_err(cudaError_t e, const char* str);
void matts_gpu_gb(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_atmos *atmos, sutra_sensors *sensors, float diamTel, float cobs);
void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu);
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors, sutra_atmos *atmos, float *L0, float *cn2,float *alphaX, float *alphaY);
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors, float *alphaX, float *alphaY);
void matcov_gpu_3(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void matcov_gpu_4(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct gtomo_struct *tomo_gpu, sutra_atmos *atmos, sutra_sensors *sensors, float *alphaX, float *alphaY);
void generateXY(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors);
void tab_dphi_gpu_gb(float *tab_dphi, struct gtomo_struct *tomo_gpu, long Ndphi, float *L0diff_d, int Nl0, float convert);
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, long Nlayer, long Nw, long Nsubap);

#endif // SUTRA_CONTROLLER_UTILS_H_
