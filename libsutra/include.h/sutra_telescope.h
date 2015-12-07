#ifndef _SUTRA_TELESCOPE_H_
#define _SUTRA_TELESCOPE_H_

#include <map>
#include <string>
#include <vector>
#include <carma.h>
#include <sutra_phase.h>

using namespace std;
using std::string;

class sutra_telescope {
public:

  carma_context *current_context;
  int device; // device #
    
  long pup_size; // size of pupil
  long num_eleme_pup; // number of points in the pupil
  
  carma_obj<float> *d_pupil; // the pupil mask
  carma_obj<float> *d_phase_ab_M1; // the phase aberration for M1
  
  long pup_size_m; // size of pupil
    
  carma_obj<float> *d_pupil_m; // the pupil mask
  carma_obj<float> *d_phase_ab_M1_m; // the phase aberration for M1
    
public:
  sutra_telescope(carma_context *context, long pup_size, long num_eleme_pup, float *pupil, float *phase_ab_M1, long pup_size_m, float *pupil_m, float *phase_ab_m1_m);
  ~sutra_telescope();
  
    
};

#endif // _SUTRA_TELESCOPE_H_