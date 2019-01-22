#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
#include <string>

template <class T>
class sutra_centroider {
 public:
  int device;
  sutra_wfs *wfs;
  int nvalid;
  int nslopes;
  int npix;
  int nxsub;

  T offset;
  T scale;

  carma_context *current_context;

  carma_obj<T> *d_bincube;
  carma_obj<T> *d_intensities;
  carma_obj<T> *d_centroids_ref;  // ref centroids
  carma_obj<T> *d_img;
  carma_obj<T> *d_img_raw;
  carma_obj<T> *d_dark;
  carma_obj<T> *d_flat;
  carma_obj<int> *d_validx;
  carma_obj<int> *d_validy;

 protected:
  sutra_centroider(carma_context *context, sutra_wfs *wfs, long nvalid,
                   T offset, T scale, int device);

 public:
  virtual ~sutra_centroider();
  int set_scale(T scale);
  int set_dark(T *dark, int n);
  int set_flat(T *flat, int n);
  int set_centroids_ref(T *centroids_ref);
  int calibrate_img(bool save_raw = false);
  int load_validpos(int *ivalid, int *jvalid, int N);
  int set_npix(int npix);
  int set_nxsub(int nxsub);
  int load_img(T *img, int n);
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }

  virtual string get_type() = 0;

  virtual int get_cog(T *cube, T *intensities, T *centroids, int nvalid,
                      int npix, int ntot) = 0;
  virtual int get_cog(T *intensities, T *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};

template <class T>
int convert_centro(T *d_odata, T *d_idata, T offset, T scale, int N,
                   carma_device *device);

#endif  // _SUTRA_CENTROIDER_H_
