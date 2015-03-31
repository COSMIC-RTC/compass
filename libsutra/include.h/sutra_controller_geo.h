#ifndef SUTRA_CONTROLLER_GEO_H_
#define SUTRA_CONTROLLER_GEO_H_

#include <sutra_controller.h>

class sutra_controller_geo: public sutra_controller {
public:
  float gain;
  long Nphi;

  carma_obj<float> *d_gain;
  carma_obj<float> *d_proj;
  carma_obj<double> *d_phi;
  carma_obj<int> *d_indx_pup;
  carma_sparse_obj<double> *d_IFsparse;
  carma_obj<float> *d_geocov;
  carma_obj<double> *d_compdouble;
  carma_obj<float> *d_compfloat;
  //carma_obj<float> *d_cenbuff; // centroids circular buffer

public:
  sutra_controller_geo(carma_context *context, long nactu, long Nphi,
      float delay, sutra_dms *dms, char **type, float *alt, int ndm);
  sutra_controller_geo(const sutra_controller_geo& controller);
  ~sutra_controller_geo();

  string
  get_type();

  cusparseHandle_t cusparse_handle() {
		return current_context->get_cusparseHandle();
	}

  int
  set_gain(float gain);
  int
  load_mgain(float *mgain);
  int
  comp_dphi(sutra_source *target);
  int
  comp_com();
  int
  init_proj(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup);
  int
  init_proj_sparse(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup);


};



#endif /* SUTRA_CONTROLLER_GEO_H_ */
