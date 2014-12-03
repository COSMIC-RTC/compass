#ifndef SUTRA_CONTROLLER_GEO_H_
#define SUTRA_CONTROLLER_GEO_H_

#include <sutra_controller.h>

class sutra_controller_geo: public sutra_controller {
public:
  int delay;
  float gain;
  long Nphi;

  carma_obj<float> *d_gain;
  carma_obj<float> *d_proj;
  carma_obj<float> *d_phi;
  carma_obj<int> *d_indx_pup;
  //carma_obj<float> *d_cenbuff; // centroids circular buffer

public:
  sutra_controller_geo(carma_context *context, long nactu, long Nphi,
      long delay);
  sutra_controller_geo(const sutra_controller_geo& controller);
  ~sutra_controller_geo();

  string
  get_type();

  int
  set_gain(float gain);
  int
  set_delay(int delay);
  int
  load_mgain(float *mgain);
  int
  comp_dphi(sutra_source *target);
  int
  comp_com();
  int
  init_proj(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup);


};



#endif /* SUTRA_CONTROLLER_GEO_H_ */
