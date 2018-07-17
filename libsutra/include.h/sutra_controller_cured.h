#ifndef _sutra_controller_cured_H_
#define _sutra_controller_cured_H_

#include <sutra_controller.h>

class sutra_controller_cured : public sutra_controller {
 public:
  float gain;
  int ndivs;     // number of subdivision levels for cured
  bool tt_flag;  // flag for separate tt

  // data for CuReD */
  carma_host_obj<float> *h_centroids;
  carma_host_obj<float> *h_err;
  carma_obj<float> *d_err;      // current error
  carma_obj<float> *d_cenbuff;  // centroids circular buffer

  // data for CuReD */
  carma_obj<float> *d_imat;

  // structures needed to run CuReD */
  // sysCure* h_syscure;
  void *h_syscure;
  // parCure* h_parcure;
  void *h_parcure;

 public:
  sutra_controller_cured(carma_context *context, long nvalid, long nactu,
                         float delay, sutra_dms *dms, int *idx_dms, int ndm);
  sutra_controller_cured(const sutra_controller_cured &controller);
  ~sutra_controller_cured();

  string get_type() { return "cured"; }
  int set_gain(float gain);

  int comp_com();

  int init_cured(int nxsubs, int *isvalid, int ndivs, int tt);
  int frame_delay();
  int set_delay(float delay);
};

#endif  // _sutra_controller_cured_H_
