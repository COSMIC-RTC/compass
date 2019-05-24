/**
 * \file sutra_target.h
 *
 * \class sutra_target
 *
 * \ingroup libsutra
 *
 * \brief this class provides the target features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_TARGET_H_
#define _SUTRA_TARGET_H_

#include <sutra_source.h>
#include <map>
#include <string>
#include <vector>

using std::vector;

typedef std::pair<std::string, int> type_screen;

class sutra_target {
 public:
  int ntargets;
  vector<sutra_source *> d_targets;

 public:
  sutra_target(carma_context *context, sutra_telescope *d_tel, int ntargets,
               float *xpos, float *ypos, float *lambda, float *mag, float zerop,
               long *sizes, int Npts, int device);
  ~sutra_target();

  int get_phase(int ntarget, float *dest);
};

#endif  // _SUTRA_TARGET_H_
