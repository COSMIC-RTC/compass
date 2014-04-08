#ifndef _sutra_controller_cured_H_
#define _sutra_controller_cured_H_

#include <sutra_controller.h>

class sutra_controller_cured: public sutra_controller {
public:
  float gain;

  // data for CuReD */
  carma_host_obj<float> *h_centroids;
  carma_host_obj<float> *h_err;

  // structures needed to run CuReD */
  sysCure* h_syscure;
  parCure* h_parcure;

public:
  sutra_controller_cured(carma_context *context, long nvalid, long nactu,
      long delay);
  sutra_controller_cured(const sutra_controller_cured& controller);
  ~sutra_controller_cured();

  string
  get_type();
  int
  set_gain(float gain);

  int
  comp_com();

  int
  init_cured(int nxsubs, int *isvalid);
};

#endif // _sutra_controller_cured_H_
