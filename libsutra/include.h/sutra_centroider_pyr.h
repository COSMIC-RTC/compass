// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pyr.h
//! \ingroup   libsutra
//! \class     sutra_centroider_pyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_CENTROIDER_PYR_H_
#define _SUTRA_CENTROIDER_PYR_H_

#include <sutra_centroider.h>

struct Method_CoG {
  bool isLocal = false;
  bool isSinus = true;

  Method_CoG(bool isLocal_ = false, bool isSinus_ = true)
      : isLocal(isLocal_), isSinus(isSinus_) {}

  /** Method_CoG(int method)
   * where method is
   *        0: sinus global
   *        1: nosinus global
   *        2: sinus local)
   *        3: nosinus local
   **/
  Method_CoG(uint8_t method) : isLocal(method > 1), isSinus(!(method % 2)) {}

  static const char *str(const struct Method_CoG &method) {
    if (method.isSinus && !method.isLocal) return "sinus global";     // 0
    if (!method.isSinus && !method.isLocal) return "nosinus global";  // 1
    if (method.isSinus && method.isLocal) return "sinus local";       // 2
    if (!method.isSinus && method.isLocal) return "nosinus local";    // 3
    throw "method unknown";
  };
};

template <class Tin, class T>
class sutra_centroider_pyr : public sutra_centroider<Tin, T> {
 private:
  string pyr_type;

 public:
  sutra_centroider_pyr(carma_context *context, sutra_wfs *wfs, long nvalid,
                       float offset, float scale, bool filter_TT, int device);
  sutra_centroider_pyr(const sutra_centroider_pyr &centroider);
  ~sutra_centroider_pyr();

  string get_type();
  int set_valid_thresh(T valid_thresh);
  T get_valid_thresh();

  int set_method(Method_CoG method);
  Method_CoG get_method();
  string get_method_str();

  int get_pyr(float *cube, float *intensities, T *centroids, int *subindx,
              int *subindy, int nvalid, int ns, int nim);
  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();

 private:
  T valid_thresh;
  Method_CoG method;
};

template <class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy,
                float *intensities, int ns, int nvalid, int nim,
                carma_device *device);

template <class T>
void pyr2_slopes(T *d_odata, T *ref, T *d_idata, int *subindx, int *subindy,
                 float *intensities, int ns, int nvalid, float scale,
                 T valid_thresh, int do_sin, carma_device *device);
#endif  // _SUTRA_CENTROIDER_PYR_H_
