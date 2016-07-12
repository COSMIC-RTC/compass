#ifndef _SUTRA_PHASE_H_
#define _SUTRA_PHASE_H_

// this is the generic class for a phase
// contains a yoga obj for the phase screen itself
// the screen size (assumed square)
// the following is only initialized on demand :
// an array of zernike coeffs on which the phase can be decomposed
// an array of zernike polynomials (carma_object)
// a matrix to decompose the phase on zernike coeffs

#include <carma.h>
#include <carma_obj.h>

class sutra_phase {
public:

  carma_context *current_context;
  int device;

  carma_obj<float> *d_screen;
  long screen_size;
  float *zerCoeff;
  carma_obj<float> *zernikes;
  carma_obj<float> *mat;

public:
  sutra_phase(carma_context *current_context, long size);
  sutra_phase(const sutra_phase& phase);
  ~sutra_phase();

};

#endif // _SUTRA_PHASE_H_
