#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
#include <string>

class sutra_centroider {
 public:
  int device;
  sutra_wfs *wfs;
  int nvalid;
  int nslopes;

  float offset;
  float scale;

  carma_context *current_context;

  carma_obj<float> *d_bincube;
  carma_obj<float> *d_subsum;
  carma_obj<float> *d_img;
  carma_obj<float> *d_img_raw;
  carma_obj<float> *d_dark;
  carma_obj<float> *d_flat;
  carma_obj<int> *d_validx;
  carma_obj<int> *d_validy;

 protected:
  sutra_centroider(carma_context *context, sutra_wfs *wfs, long nvalid,
                   float offset, float scale, int device);

 public:
  virtual ~sutra_centroider();
  int set_scale(float scale);
  int set_dark(float *dark, int n);
  int set_flat(float *flat, int n);
  int calibrate_img(bool save_raw = false);
  int load_validpos(int *ivalid, int *jvalid, int N);
  int fill_bincube(int npix);
  int load_img(float *img, int n);
  int load_pyrimg(float *img, int n);
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }
  int load_img_gpu(float *img, int n);

  virtual string get_type() = 0;

  virtual int get_cog(carma_streams *streams, float *cube, float *subsum,
                      float *centroids, int nvalid, int npix, int ntot) = 0;
  virtual int get_cog(float *subsum, float *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};

template <class T>
int convert_centro(T *d_odata, T *d_idata, T offset, T scale, int N,
                   carma_device *device);

#endif  // _SUTRA_CENTROIDER_H_
