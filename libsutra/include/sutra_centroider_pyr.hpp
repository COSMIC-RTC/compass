// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pyr.hpp
//! \ingroup   libsutra
//! \class     SutraCentroiderPyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CENTROIDER_PYR_H_
#define _SUTRA_CENTROIDER_PYR_H_

#include <sutra_centroider.hpp>

struct Method_CoG {
  bool is_local = false;
  bool is_sinus = true;

  Method_CoG(bool isLocal_ = false, bool isSinus_ = true)
      : is_local(isLocal_), is_sinus(isSinus_) {}

  /** Method_CoG(int32_t method)
   * where method is
   *        0: sinus global
   *        1: nosinus global
   *        2: sinus local)
   *        3: nosinus local
   **/
  Method_CoG(uint8_t method) : is_local(method > 1), is_sinus(!(method % 2)) {}

  static const char *str(const struct Method_CoG &method) {
    if (method.is_sinus && !method.is_local) return "sinus global";     // 0
    if (!method.is_sinus && !method.is_local) return "nosinus global";  // 1
    if (method.is_sinus && method.is_local) return "sinus local";       // 2
    if (!method.is_sinus && method.is_local) return "nosinus local";    // 3
    throw "method unknown";
  };
};

template <class Tin, class T>
class SutraCentroiderPyr : public SutraCentroider<Tin, T> {
 private:
  string pyr_type;

 public:
  SutraCentroiderPyr(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                       float offset, float scale, bool filter_TT, int32_t device);
  SutraCentroiderPyr(const SutraCentroiderPyr &centroider);
  ~SutraCentroiderPyr();

  string get_type();
  int32_t set_valid_thresh(T valid_thresh);
  T get_valid_thresh();

  int32_t set_method(Method_CoG method);
  Method_CoG get_method();
  string get_method_str();

  int32_t get_pyr(float *cube, float *intensities, T *centroids, int32_t *subindx,
              int32_t *subindy, int32_t nvalid, int32_t ns, int32_t nim, cudaStream_t stream=0);
  int32_t get_cog(float *cube, float *intensities, T *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, T *slopes, bool noise);
  int32_t get_cog();

 private:
  T valid_thresh;
  Method_CoG method;
};

template <class T>
void pyr_slopes(T *d_odata, T *ref, T *d_idata, int32_t *subindx, int32_t *subindy,
                 float *intensities, int32_t ns, int32_t nvalid, float scale,
                 T valid_thresh, int32_t do_sin, SlopeOrder slope_order,
                 CarmaDevice *device, cudaStream_t stream=0);
#endif  // _SUTRA_CENTROIDER_PYR_H_
