#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_ao_utils.h>
#include <carma_utils.h>

sutra_wfs_pyr_pyrhr::sutra_wfs_pyr_pyrhr(carma_context *context,
                                         sutra_telescope *d_tel,
                                         sutra_sensors *sensors, long nxsub,
                                         long nvalid, long npix, long nphase,
                                         long nrebin, long nfft, long ntot,
                                         long npup, float pdiam, float nphotons,
                                         float nphot4imat, int lgs, int device) :
    sutra_wfs_pyr(context, d_tel, sensors, nxsub, nvalid, npix, nphase, nrebin,
                  nfft, ntot, npup, pdiam, nphotons, nphot4imat, lgs, device,
                  static_cast<char*>("pyr")) {
  this->type = "pyr";

}

sutra_wfs_pyr_pyrhr::~sutra_wfs_pyr_pyrhr() {

}

int sutra_wfs_pyr_pyrhr::wfs_initarrays(cuFloatComplex *halfxy, int *cx,
                                        int *cy, float *sincar, int *validsubsx,
                                        int *validsubsy) {
  current_context->set_activeDevice(device,1);
  this->d_phalfxy->host2device(halfxy);
  this->pyr_cx->fill_from(cx);
  this->pyr_cy->fill_from(cy);
  this->d_sincar->host2device(sincar);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);

  return EXIT_SUCCESS;
}

//////////////////////////////
// PYRAMID WAVEFRONT SENSOR //
//////////////////////////////

// It starts by looking for the type of sensor. By default it assumes
// a pyramid wfs. The pyramid can also be explicitely asked for, or
// a roof prism can be asked for as well.

int sutra_wfs_pyr_pyrhr::comp_generic() {
  /*
   //___________________________________________________________________
   //  PYRAMID SENSOR MODEL

   This code generates pupil images as seen from behind a pyramid wavefront sensor
   algorithm:
   for (i=0;i<mod_pts;i++) {
   get phase and multiply by exp(i*modu)
   do fft
   apply field stop
   myltiply by exp(i*pyramid) where pyramid is the pyramid shape
   do fft-1
   do abs2
   add to previous modulation image
   }
   do fft
   multiply by sinc (pixes transfer function)
   take 1 pixels over nrebin pixels in the image
   normalize
   add noise
   */
  current_context->set_activeDevice(device,1);

  carmaSafeCall(
      cudaMemset(this->d_hrimg->getData(), 0,
                 sizeof(float) * this->d_hrimg->getNbElem()));
  //this->npup = 1;
  for (int cpt = 0; cpt < this->npup; cpt++) {
    //int cpt = 0;
    carmaSafeCall(
        cudaMemset(this->d_camplipup->getData(), 0,
                   2 * sizeof(float) * this->d_camplipup->getNbElem()));

    pyr_getpup(this->d_camplipup->getData(),
               this->d_gs->d_phase->d_screen->getData(),
               this->d_pupil->getData(), this->ntot, this->nfft,
               this->d_gs->lambda, (this->pyr_cx->getData())[cpt],
               (this->pyr_cy->getData())[cpt],
               this->current_context->get_device(device));

    carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), -1,
              *this->d_camplipup->getPlan());
    /*
     pyr_submask(this->d_camplifoc->getData(), this->d_submask->getData(),
     this->nfft, this->current_context->get_device(device));
     */
    pyr_submaskpyr(this->d_camplifoc->getData(), this->d_phalfxy->getData(),
                   this->nfft, this->current_context->get_device(device));

    carma_fft(this->d_camplifoc->getData(), this->d_fttotim->getData(), 1,
              *this->d_camplipup->getPlan());
    //float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
    float fact = 1.0f;

    abs2(this->d_hrimg->getData(), this->d_fttotim->getData(),
         this->nfft * this->nfft, fact,
         this->current_context->get_device(device));

  }

  cfillrealp(this->d_fttotim->getData(), this->d_hrimg->getData(),
             this->d_hrimg->getNbElem(),
             this->current_context->get_device(device));

  carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
            *this->d_camplipup->getPlan());

  pyr_submask(this->d_fttotim->getData(), this->d_sincar->getData(), this->nfft,
              this->current_context->get_device(device));

  carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
            *this->d_camplipup->getPlan());

  cgetrealp(this->d_hrimg->getData(), this->d_fttotim->getData(),
            this->d_hrimg->getNbElem(),
            this->current_context->get_device(device));

  pyr_fillbinimg(this->d_binimg->getData(), this->d_hrimg->getData(),
                 this->nfft / this->nrebin, this->nfft, this->nrebin, false,
                 this->current_context->get_device(device));

  pyr_subsum(this->d_subsum->getData(), this->d_binimg->getData(),
             this->d_validsubsx->getData(), this->d_validsubsy->getData(),
             this->nfft / this->nrebin, this->nvalid,
             this->current_context->get_device(device));

  int blocks, threads;
//  getNumBlocksAndThreads(current_context->get_device(device), this->nvalid,
//      blocks, threads);
  sumGetNumBlocksAndThreads(this->nvalid, device, blocks, threads);

  reduce(this->nvalid, threads, blocks, this->d_subsum->getData(),
         this->d_subsum->getData());

  pyr_fact(this->d_binimg->getData(), this->nphot, this->d_subsum->getData(),
           this->nfft / this->nrebin, 1,
           this->current_context->get_device(device));

  // add noise
  if (this->noise > -1) {
    //cout << "adding poisson noise" << endl;
    this->d_bincube->prng('P');
  }
  if (this->noise > 0) {
    //cout << "adding detector noise" << endl;
    this->d_bincube->prng('N', this->noise, 1.0f);
  }

  pyr_subsum(this->d_subsum->getData(), this->d_binimg->getData(),
             this->d_validsubsx->getData(), this->d_validsubsy->getData(),
             this->nfft / this->nrebin, this->nvalid,
             this->current_context->get_device(device));

  return EXIT_SUCCESS;

}

int sutra_wfs_pyr_pyrhr::comp_image() {

  current_context->set_activeDevice(device,1);
  int result = comp_generic();
  return result;
}
