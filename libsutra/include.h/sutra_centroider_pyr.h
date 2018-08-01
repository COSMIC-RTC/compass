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
   *        0: nosinus global
   *        1: sinus global
   *        2: nosinus local
   *        3: sinus local)
   **/
  Method_CoG(uint8_t method) : isLocal(method > 1), isSinus(!(method % 2)) {}

  static const char *str(const struct Method_CoG &method) {
    if (method.isSinus && method.isLocal) return "sinus local";
    if (!method.isSinus && method.isLocal) return "nosinus local";
    if (method.isSinus && !method.isLocal) return "sinus global";
    if (!method.isSinus && !method.isLocal) return "nosinus global";
    throw "method unknown";
  };
};

class sutra_centroider_pyr : public sutra_centroider {
 private:
  string pyr_type;

 public:
  sutra_centroider_pyr(carma_context *context, sutra_wfs *wfs, long nvalid,
                       float offset, float scale, int device);
  sutra_centroider_pyr(const sutra_centroider_pyr &centroider);
  ~sutra_centroider_pyr();

  string get_type();
  int set_valid_thresh(float valid_thresh);
  float get_valid_thresh();

  int set_method(Method_CoG method);
  Method_CoG get_method();
  string get_method_str();

  int get_pyr(float *cube, float *subsum, float *centroids, int *subindx,
              int *subindy, int nvalid, int ns, int nim);
  int get_cog(carma_streams *streams, float *cube, float *subsum,
              float *centroids, int nvalid, int npix, int ntot);
  int get_cog(float *subsum, float *slopes, bool noise);
  int get_cog();

 private:
  float valid_thresh;
  Method_CoG method;
};

template <class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                int ns, int nvalid, int nim, carma_device *device);

template <class T>
void pyr2_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                 int ns, int nvalid, T scale, T valid_thresh, int do_sin,
                 carma_device *device);

#endif  // _SUTRA_CENTROIDER_PYR_H_
