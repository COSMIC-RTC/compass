#include <sutra_wfs_pyr_roof.h>
#include <sutra_ao_utils.h>
#include <carma_utils.h>

sutra_wfs_pyr_roof::sutra_wfs_pyr_roof(carma_context *context, sutra_telescope *d_tel, sutra_sensors *sensors, long nxsub,
    long nvalid, long npix, long nphase, long nrebin, long nfft, long ntot,
    long npup, float pdiam, float nphotons, float nphot4imat, int lgs, int device) : sutra_wfs_pyr(context, d_tel, sensors, nxsub,
    nvalid, npix, nphase, nrebin, nfft, ntot, npup, pdiam, nphotons, nphot4imat, lgs, device){
  this->type = "roof";
}

sutra_wfs_pyr_roof::~sutra_wfs_pyr_roof(){

}

//////////////////////////////
// ROOF WAVEFRONT SENSOR //
//////////////////////////////

// It starts by looking for the type of sensor. By default it assumes
// a pyramid wfs. The pyramid can also be explicitely asked for, or
// a roof prism can be asked for as well.

int sutra_wfs_pyr_roof::comp_generic() {

  current_context->set_activeDevice(device,1);
  //___________________________________________________________________
  //  ROOF SENSOR
  //PYR_GETPUP: reads pupil & phase and computes rolled electric field
  pyr_getpup(this->d_camplipup->getData(),
      this->d_gs->d_phase->d_screen->getData(), this->d_phalfxy->getData(),
      this->d_pupil->getData(), this->ntot, this->d_gs->lambda,
      this->current_context->get_device(device));

  carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), -1,
      *this->d_camplipup->getPlan());

  pyr_submask(this->d_camplifoc->getData(), this->d_submask->getData(),
      this->ntot, this->current_context->get_device(device));

  carmaSafeCall(
      cudaMemset(this->d_hrimg->getData(), 0,
          sizeof(float) * this->d_hrimg->getNbElem()));

  //this->npup = 1;
  for (int cpt = 0; cpt < this->npup; cpt++) {
    // modulation loop
    // computes the high resolution images
    carmaSafeCall(
        cudaMemset(this->d_fttotim->getData(), 0,
            sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));

    roof_rollmod(this->d_fttotim->getData(), this->d_camplifoc->getData(),
        this->d_poffsets->getData(), (this->pyr_cx->getData())[cpt],
        (this->pyr_cy->getData())[cpt], this->ntot, this->nfft,
        this->current_context->get_device(device));
    //
    //    pyr_rollmod(this->d_fttotim->getData(),this->d_camplifoc->getData(), this->d_poffsets->getData(),0,
    //    0,this->ntot , this->nfft, this->current_context->get_device(device));
    //

    carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
        *this->d_fttotim->getPlan());

    float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
    //if (cpt == this->npup-1) fact = fact / this->npup;

    roof_abs2(this->d_hrimg->getData(), this->d_fttotim->getData(), fact,
        this->nfft, 4, this->current_context->get_device(device));
  }

  if (this->noise > 0) {
    this->d_bincube->prng('N', this->noise);
  } else
    carmaSafeCall(
        cudaMemset(this->d_bincube->getData(), 0,
            sizeof(float) * this->d_bincube->getNbElem()));

  roof_fillbin(this->d_bincube->getData(), this->d_hrimg->getData(),
      this->nrebin, this->nfft, this->nfft / this->nrebin, 4,
      this->current_context->get_device(device));

  pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
      this->d_validsubsx->getData(), this->d_validsubsy->getData(),
      this->nfft / this->nrebin, this->nvalid, 4,
      this->current_context->get_device(device));

  int blocks, threads;
  getNumBlocksAndThreads(this->current_context->get_device(device),
      this->nvalid, blocks, threads);
  reduce(this->nvalid, threads, blocks, this->d_subsum->getData(),
      this->d_subsum->getData());

  pyr_fact(this->d_bincube->getData(), this->nphot, this->d_subsum->getData(),
      this->nfft / this->nrebin, 4, this->current_context->get_device(device));

  pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
      this->d_validsubsx->getData(), this->d_validsubsy->getData(),
      this->nfft / this->nrebin, this->nvalid, 4,
      this->current_context->get_device(device));

  return EXIT_SUCCESS;

}

int sutra_wfs_pyr_roof::comp_image() {

  current_context->set_activeDevice(device,1);
  int result = comp_generic();
  return result;
}
