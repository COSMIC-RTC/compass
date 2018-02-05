#include <sutra_centroider_pyr.h>
#include <iostream>
#include <string>

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context,
    sutra_sensors *sensors, int nwfs, long nvalid, float offset,
    float scale, int device) : sutra_centroider(context, sensors, nwfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);
  if(sensors != nullptr) {
    this->pyr_type = sensors->d_wfs[nwfs]->type;
  } else {
    this->pyr_type = "pyrhr";
  }
  this->valid_thresh = 1e-4;

  // centroider method by default sin_global
  this->method = Method_CoG(false, true);
}

sutra_centroider_pyr::~sutra_centroider_pyr() {
}

string sutra_centroider_pyr::get_type() {
  return this->pyr_type;
}

int sutra_centroider_pyr::set_valid_thresh(float valid_thresh) {
  this->valid_thresh = valid_thresh;
  return EXIT_SUCCESS;
}
float sutra_centroider_pyr::get_valid_thresh() {
  return this->valid_thresh;
}

int sutra_centroider_pyr::set_method(Method_CoG method) {
  this->method = method;
  return EXIT_SUCCESS;
}

Method_CoG sutra_centroider_pyr::get_method() {
  return this->method;
}

string sutra_centroider_pyr::get_method_str() {
  return Method_CoG::str(this->method);
}

int sutra_centroider_pyr::get_cog(carma_streams *streams, float *cube,
                                  float *subsum, float *centroids, int nvalid, int npix, int ntot) {
  // TODO(Implement sutra_centroider_pyr::get_cog)

  return get_pyr(cube, subsum, centroids, this->d_validx->getData(), this->d_validy->getData(), this->nvalid, this->d_bincube->getDims(1), 4);
}

int sutra_centroider_pyr::get_pyr(float *cube, float *subsum, float *centroids,
                                  int *subindx, int *subindy, int nvalid, int ns, int nim) {
  current_context->set_activeDevice(device, 1);
  if (this->pyr_type == "pyr" || this->pyr_type == "roof") {
    DEBUG_TRACE("Not working anymore");
    throw "Not working anymore";
    // pyr_subsum(subsum, cube, subindx, subindy, ns, nvalid, nim,
    //         this->current_context->get_device(device));
    //
    // pyr_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim,
    //         this->current_context->get_device(device));
  } else if (this->pyr_type == "pyrhr") {
    pyr_subsum(subsum, cube, subindx, subindy, ns, nvalid,
               this->current_context->get_device(device));

    if( !(this->method.isLocal) ) {
      // if we are using a global method
      // DEBUG_TRACE("Global : %s", Method_CoG::str(this->method));
      int blocks, threads;
      this->current_context->set_activeDevice(device,1);
      sumGetNumBlocksAndThreads(nvalid,
                                this->current_context->get_device(device), blocks, threads);
      long dims[2] = {1, nvalid};
      carma_obj<float> tmp_obj(current_context, dims);
      reduce(nvalid, threads, blocks, subsum, tmp_obj.getData());

      fillvalues(subsum, tmp_obj.getData(), nvalid,
                 this->current_context->get_device(device));
      // } else {
      //     DEBUG_TRACE("Local : %s", Method_CoG::str(this->method));
    }

    // if(this->method.isSinus){  // if we are using a global method
    //     DEBUG_TRACE("Sinus : %s", Method_CoG::str(this->method));
    // } else {
    //     DEBUG_TRACE("NoSinus : %s", Method_CoG::str(this->method));
    // }

    pyr2_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, this->scale,
                this->valid_thresh, this->method.isSinus,  // if we are using a sin method
                this->current_context->get_device(device));
  } else {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(float *subsum, float *slopes, bool noise) {
  if(this->wfs != nullptr) {
    if (this->pyr_type == "pyr" || this->pyr_type == "roof")
      return this->get_pyr(*(wfs->d_bincube), subsum, slopes,
                           *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
                           wfs->nfft / wfs->nrebin, 4);
    else if (this->pyr_type == "pyrhr") {
      if(noise || wfs->error_budget == false) {
        return this->get_pyr(*(wfs->d_binimg), subsum, slopes,
                             *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
                             wfs->nfft / wfs->nrebin, 4);
      } else
        return this->get_pyr(*(wfs->d_binimg_notnoisy), subsum, slopes,
                             *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
                             wfs->nfft / wfs->nrebin, 4);
    }
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_pyr::get_cog() {
  if(this->wfs != nullptr)
    return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes),true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;

}
