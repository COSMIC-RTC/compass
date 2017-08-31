#ifndef _SUTRA_CENTROIDER_PYR_H_
#define _SUTRA_CENTROIDER_PYR_H_

#include <sutra_centroider.h>

struct Method_CoG {
  enum Flags : uint8_t {Sinus=0x01, Local=0x02, Other=0x04};

  static const char* str(uint8_t method) {
    if (method>=Other) return "method unknown";
    if (  method&Sinus  &&   method&Local)  return "sinus local";
    if (~(method&Sinus) &&   method&Local)  return "nosinus local";
    if (  method&Sinus  && ~(method&Local)) return "sinus global";
    if (~(method&Sinus) && ~(method&Local)) return "nosinus global";
    throw "method unknown";
  };
};

class sutra_centroider_pyr: public sutra_centroider {
 private:
  string pyr_type;

 public:
  sutra_centroider_pyr(carma_context *context, sutra_sensors *sensors, int nwfs, long nvalid,
                       float offset, float scale, int device);
  sutra_centroider_pyr(const sutra_centroider_pyr& centroider);
  ~sutra_centroider_pyr();

  string get_type();
  int set_valid_thresh(float valid_thresh);
  float get_valid_thresh();

  int set_method(uint8_t method);
  uint8_t get_method();
  string get_method_str();

  int get_pyr(float *cube, float *subsum, float *centroids, int *subindx,
              int *subindy, int nvalid, int ns, int nim);
  int get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
              int nvalid, int npix, int ntot);
  int get_cog(float *subsum, float *slopes, bool noise);
  int get_cog();

 private:
  float valid_thresh;
  uint8_t method;
};

#endif // _SUTRA_CENTROIDER_PYR_H_
