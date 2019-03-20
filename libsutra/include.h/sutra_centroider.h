#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
#include <string>

template <class Tin, class Tout>
class sutra_centroider {
 public:
  int device;
  sutra_wfs *wfs;
  int nvalid;
  int nslopes;
  int npix;
  int nxsub;

  float offset;
  float scale;

  carma_context *current_context;

  carma_obj<Tout> *d_bincube;
  carma_obj<float> *d_intensities;
  carma_obj<Tout> *d_centroids_ref;  // ref centroids
  carma_obj<float> *d_img;
  carma_obj<Tin> *d_img_raw;
  carma_obj<float> *d_dark;
  carma_obj<float> *d_flat;
  carma_obj<int> *d_validx;
  carma_obj<int> *d_validy;
  carma_obj<int> *d_validMask;

 protected:
  sutra_centroider(carma_context *context, sutra_wfs *wfs, long nvalid,
                   float offset, float scale, int device);

 public:
  virtual ~sutra_centroider();
  int set_scale(float scale);
  int set_dark(float *dark, int n);
  int set_flat(float *flat, int n);
  int set_centroids_ref(float *centroids_ref);
  int calibrate_img();
  int load_validpos(int *ivalid, int *jvalid, int N);
  int set_npix(int npix);
  int set_nxsub(int nxsub);
  int load_img(Tin *img, int n);
  int load_img(Tin *img, int n, int location);
  int load_img(carma_obj<Tin> *img);
  int get_validMask();
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }

  virtual string get_type() = 0;

  virtual int get_cog(float *img, float *intensities, Tout *centroids,
                      int nvalid, int npix, int ntot) = 0;
  virtual int get_cog(float *intensities, Tout *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};
template <class Tin>
int calibration(Tin *img_raw, float *img_cal, float *dark, float *flat, int N,
                carma_device *device);

template <typename T>
int convert_centro(T *d_odata, T *d_idata, float offset, float scale, int N,
                   carma_device *device);
int fill_validMask(int size, int npix, int blocks, int *d_validMask,
                   int *validx, int *validy, carma_device *device);

#endif  // _SUTRA_CENTROIDER_H_
