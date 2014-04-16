#ifndef _sutra_controller_cured_H_
#define _sutra_controller_cured_H_

#include <sutra_controller.h>

class sutra_controller_cured: public sutra_controller {
public:
  float gain;
  int   ndivs; //number of subdivision levels for cured

  // data for CuReD */
  carma_host_obj<float> *h_centroids;
  carma_host_obj<float> *h_err;
  carma_obj<float> *d_err; // current error

  // data for CuReD */
  carma_obj<float> *d_imat;

  // structures needed to run CuReD */
  //sysCure* h_syscure;
  void* h_syscure;
  //parCure* h_parcure;
  void* h_parcure;

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
    init_cured(int nxsubs, int *isvalid, int ndivs);
};

#endif // _sutra_controller_cured_H_
