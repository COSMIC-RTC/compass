#ifndef SUTRA_CONTROLLER_UTILS_H_
#define SUTRA_CONTROLLER_UTILS_H_

#include <sutra_turbu.h>
#include <sutra_wfs.h>
#include <sutra_dm.h>
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
  int   nlgs;

  long  *ioff_d;
  long  *Nssp_d;

  float *alphaX_d;
  float *alphaY_d;
  float *GsAlt_d;
  float *diamPup_d;
  float *thetaML_d;
  float *X_d;
  float *Y_d;
  float *XPup_d;
  float *YPup_d;

  long   max_Nl0;
  long  *indexL0_d;
  float *L0diff_d;
  float *h_d;
  float *cn2_d;
  float *tabDPHI_d;
  float *u_d;
  float *v_d;
  float *sspSizeL_d;

  cudaStream_t matcov_stream;
};

struct cphim_struct{
	  float DiamTel;
	  long Ndphi;
	  float pasDPHI;
	  long int_npts;
	  long Nw;
	  long Nx;
	  long Ndm;
	  long Nactu;
	  long *Nssp;
	  float *diamPup;
	  float *XPup;
	  float *YPup;
	  float *thetaML;
	  float *alphaX;
	  float *alphaY;
	  float *GsAlt;
	  float x0;
	  float y0;
	  float dx;

	  float lgs_cst;
	  float spot_width;
	  float lgs_depth;
	  float lgs_alt;
	  int   nlgs;

	  long  *ioff_d;
	  long  *Nssp_d;

	  float *alphaX_d;
	  float *alphaY_d;
	  float *GsAlt_d;
	  float *diamPup_d;
	  float *thetaML_d;
	  float *X_d;
	  float *Y_d;
	  float *XPup_d;
	  float *YPup_d;

	  long   max_Nl0;
	  long  *indexL0_d;
	  float *L0diff_d;
	  float *h_d;
	  float *cn2_d;
	  float *tabDPHI_d;
	  float *tab_int_x;
	  float *tab_int_y;
	  float *xact_d;
	  float *yact_d;
	  float *u_d;
	  float *v_d;
	  float *sspSizeL_d;

	  cudaStream_t cphim_stream;
};

void process_err(cudaError_t e, const char* str);
void matts_gpu_gb(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_atmos *atmos, sutra_sensors *sensors, float diamTel, float cobs);
void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu);
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors, sutra_atmos *atmos, float *L0, float *cn2,float *alphaX, float *alphaY);
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors, float *alphaX, float *alphaY);
void update_cphim_atm(struct cphim_struct *cphim_struct, sutra_sensors *sensors, sutra_atmos *atmos, float *L0, float *cn2, float *alphaX, float *alphaY);
void update_cphim_sys(struct cphim_struct *cphim_struct, sutra_sensors *sensors, float *alphaX, float *alphaY, float *xactu, float *yactu, float *X, float *Y);
void matcov_gpu_3(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void matcov_gpu_4(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct gtomo_struct *tomo_gpu, sutra_atmos *atmos, sutra_sensors *sensors, float *alphaX, float *alphaY);
void CPHIM(float* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct cphim_struct *cphim_struct, sutra_atmos *atmos, sutra_sensors *sensors, float* alphaX, float *alphaY, float k2, carma_device *device);
void generateXY(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors);
void tab_dphi_gpu_gb(float *tab_dphi, struct gtomo_struct *tomo_gpu, long Ndphi, float *L0diff_d, int Nl0, float convert);
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, long Nlayer, long Nw, long Nsubap);
void sub_pos_cphim(struct cphim_struct *cphim_struct, long Nlayer, long Nw, long Nsubap);
void tab_u831J0(struct cphim_struct *cphim_gpu, float tmin, float tmax, carma_device *device);
void cuda_zcen(float *idata, float *odata, int N, carma_device *device);
void cumsum(float *odata, float *idata, int N);
void init_cphim_struct(struct cphim_struct *cphim_struct, sutra_atmos *atmos, sutra_sensors *sensors, sutra_dms *dms, float diamTel);
void free_cphim_struct(struct cphim_struct *cphim_struct);

#endif // SUTRA_CONTROLLER_UTILS_H_
