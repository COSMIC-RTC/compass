#include <sutra_centroider_cog.h>
#include <string>

sutra_centroider_cog::sutra_centroider_cog(carma_context *context,
                                           sutra_wfs *wfs, long nvalid,
                                           float offset, float scale,
                                           int device)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  this->nslopes = 2 * nvalid;
}

sutra_centroider_cog::~sutra_centroider_cog() {}

string sutra_centroider_cog::get_type() { return "cog"; }

int sutra_centroider_cog::get_cog(carma_streams *streams, float *cube,
                                  float *ref, float *centroids, int nvalid,
                                  int npix, int ntot) {
  current_context->set_activeDevice(device, 1);
  // simple cog
  int nstreams;
  if (streams != nullptr)
    nstreams = streams->get_nbStreams();
  else
    nstreams = 0;
  // fprintf(stderr, "\n[%s@%d]: nstreams=%d\n", __FILE__, __LINE__, nstreams);
  if (nstreams > 1) {
    // fprintf(stderr, "\n[%s@%d]: i=%d istart=%d npix=%d nvalid=%d\n",
    // __FILE__, __LINE__, i, istart, npix, nvalid);
    subap_reduce_async(npix * npix, nvalid, streams, cube, ref);
    // fprintf(stderr, "\n[%s@%d] I'm here\n", __FILE__, __LINE__);
    get_centroids_async(npix * npix, nvalid, npix, streams, cube, centroids,
                        ref, this->scale, this->offset);
    // fprintf(stderr, "\n[%s@%d] I'm here\n", __FILE__, __LINE__);
    // streams->wait_all_streams();
  } else {
    // subap_reduce(ntot, (npix * npix), nvalid, cube, ref,
    //              current_context->get_device(device));
    get_centroids(ntot, (npix * npix), nvalid, npix, cube, centroids, ref,
                  this->d_validx->getData(), this->d_validy->getData(),
                  this->scale, this->offset,
                  current_context->get_device(device));
  }
  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(float *ref, float *slopes, bool noise) {
  if (this->wfs != nullptr) {
    if (noise || wfs->roket == false) {
      return this->get_cog(wfs->streams, *(wfs->d_binimg), ref, slopes,
                           wfs->nvalid_tot, wfs->npix,
                           wfs->d_binimg->getDims()[1]);
    } else {
      return this->get_cog(wfs->streams, *(wfs->d_bincube_notnoisy), ref,
                           slopes, wfs->nvalid_tot, wfs->npix,
                           wfs->d_binimg->getDims()[1]);
    }
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_cog::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
