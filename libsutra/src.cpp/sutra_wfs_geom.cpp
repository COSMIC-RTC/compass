#include <carma_utils.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>
#include <sutra_wfs_geom.h>

sutra_wfs_geom::sutra_wfs_geom(carma_context *context, sutra_telescope *d_tel,
                               long nxsub, long nvalid, long nphase, long npup,
                               float pdiam, int device)
    : sutra_wfs(context, d_tel, nullptr, "geo", nxsub, nvalid, 0, nphase, 0, 0,
                0, npup, pdiam, 0, 0, false, device) {
  context->set_activeDevice(device, 1);

  this->nstreams = 1;  // nvalid/10;
  this->streams = new carma_streams(nstreams);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new carma_obj<float>(context, dims_data1);

  dims_data2[1] = nphase;
  dims_data2[2] = nphase;
  this->d_offsets = new carma_obj<float>(context, dims_data2);

  dims_data1[1] = nvalid;
  this->d_subsum = new carma_obj<float>(context, dims_data1);

  this->d_fluxPerSub = new carma_obj<float>(context, dims_data1);
  this->d_validsubsx = new carma_obj<int>(context, dims_data1);
  this->d_validsubsy = new carma_obj<int>(context, dims_data1);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new carma_obj<int>(context, dims_data2);
  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;
}

sutra_wfs_geom::~sutra_wfs_geom() {
  current_context->set_activeDevice(device, 1);
  if (this->type != "sh" && this->d_camplipup != 0L) delete this->d_camplipup;
  if (this->type != "sh" && this->d_camplifoc != 0L) delete this->d_camplifoc;

  if (this->type != "sh" && this->d_fttotim != 0L) delete this->d_fttotim;

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_subsum != 0L) delete this->d_subsum;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;

  if (this->lgs) delete this->d_gs->d_lgs;

  delete this->d_gs;

  delete this->streams;

  // delete this->current_context;
}

int sutra_wfs_geom::wfs_initarrays(int *phasemap, float *offsets,
                                   float *fluxPerSub, int *validsubsx,
                                   int *validsubsy) {
  current_context->set_activeDevice(device, 1);
  this->d_phasemap->host2device(phasemap);
  this->d_offsets->host2device(offsets);
  this->d_fluxPerSub->host2device(fluxPerSub);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);

  return EXIT_SUCCESS;
}

int sutra_wfs_geom::slopes_geom(int type, float *slopes) {
  current_context->set_activeDevice(device, 1);
  /*
   normalization notes :
   σ² = 0.17 (λ/D)^2 (D/r_0)^(5/3) , σ² en radians d'angle
   σ = sqrt(0.17 (λ/D)^2 (D/r_0)^(5/3)) * 206265 , σ en secondes

   // computing subaperture phase difference at edges

   todo : integrale( x * phase ) / integrale (x^2);
   with x = span(-0.5,0.5,npixels)(,-:1:npixels) * subap_diam * 2 * pi / lambda
   / 0.206265
   */
  if (type == 0) {
    // this is to convert in arcsec
    //> 206265* 0.000001/ 2 / 3.14159265 = 0.0328281
    // it would have been the case if the phase was given in radiants
    // but it is given in microns so normalization factor is
    // just 206265* 0.000001 = 0.206265

    // float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd;
    phase_reduce(this->nphase, this->nvalid,
                 this->d_gs->d_phase->d_screen->getData(), slopes,
                 this->d_phasemap->getData(), alpha);
  }

  if (type == 1) {
    // float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd;
    phase_derive(this->nphase * this->nphase * this->nvalid,
                 this->nphase * this->nphase, this->nvalid, this->nphase,
                 this->d_gs->d_phase->d_screen->getData(), slopes,
                 this->d_phasemap->getData(), this->d_pupil->getData(), alpha,
                 this->d_fluxPerSub->getData());
  }

  return EXIT_SUCCESS;
}

int sutra_wfs_geom::slopes_geom(int type) {
  this->slopes_geom(type, this->d_slopes->getData());

  return EXIT_SUCCESS;
}
