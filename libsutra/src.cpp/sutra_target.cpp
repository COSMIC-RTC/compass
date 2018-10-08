#include <sutra_target.h>

sutra_target::sutra_target(carma_context *context, sutra_telescope *yTelescope,
                           int ntargets, float *xpos, float *ypos,
                           float *lambda, float *mag, float zerop, long *sizes,
                           int Npts, int device) {
  this->ntargets = ntargets;

  for (int i = 0; i < ntargets; i++) {
    d_targets.push_back(new sutra_source(context, xpos[i], ypos[i], lambda[i],
                                         mag[i], zerop, sizes[i], "target",
                                         yTelescope->d_pupil, Npts, device));
  }
}

sutra_target::~sutra_target() {
  //  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
  while ((this->d_targets).size() > 0) {
    delete this->d_targets.back();
    d_targets.pop_back();
  }
}
