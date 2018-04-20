#include <sutra_gamora.h>

sutra_gamora::sutra_gamora(carma_context *context, int device, char *type, int nactus,
                           int nmodes, int niter, float *IFvalue, int *IFrowind, int *IFcolind, int IFnz,
                           float *TT, float *pupil, int size, int Npts, float scale, float *Btt, float *covmodes) {

  this->current_context = context;

  const int ngpu = context->get_ndevice();
  DEBUG_TRACE("using GAMORA with %d GPUs", ngpu);
  if(ngpu == 1)
    this->device = device;
  else {
    int devices[ngpu];
    for (int i = 0; i < ngpu; i++) {
      devices[i] = i;
    }
    this->device = devices[0];
  }

  this->current_context->set_activeDevice(this->device,1);

  this->scale = scale; // 2*pi/lambda, lambda in microns
  this->Npts = Npts; // number of pupil points
  this->nactus = nactus;
  this->nmodes = nmodes;
  this->niter = niter;
  this->size = size; // pupdiam in pixels
  this->d_term1 = NULL;
  this->d_term2 = NULL;
  this->d_newmodek = NULL;
  this->d_otftel = NULL;
  this->d_otfVii = NULL;
  this->d_Btt = NULL;
  this->d_covmodes = NULL;
  this->d_mask = NULL;
  this->d_err = NULL;
  this->d_amplipup = NULL;
  this->d_wherephase = NULL;
  this->d_phase = NULL;
  this->d_psf = NULL;
  this->d_IF = NULL;
  this->d_TT = NULL;

  int mradix = 2;
  int fft_size = pow(mradix, (long) (logf(2 * size) / logf(mradix)) + 1);

  long dims_data2[3] = {2,niter,nactus};
  long dims_data1[2] = {1,Npts};

  if(strcmp(type,"roket") == 0) {
    // Command error
    this->d_err = new carma_obj<float>(this->current_context,dims_data2);

  }

  int *wherephase;
  wherephase = (int*)malloc(Npts*sizeof(int));
  int cpt = 0;
  // Phase point index in spupil
  for (int cc=0; cc < size*size; cc++) {
    if(pupil[cc]>0) {
      wherephase[cpt] = cc;
      cpt += 1;
    }
  }
  this->d_wherephase = new carma_obj<int>(this->current_context,dims_data1,wherephase);
  this->d_phase = new carma_obj<float>(this->current_context,dims_data1);

  //dims_data2[1] = size;
  //dims_data2[2] = size;

  //spupil
  //this->d_pupil = new carma_obj<float>(this->current_context,dims_data2,pupil);
  // FFT good size
  dims_data2[1] = fft_size;
  dims_data2[2] = fft_size;

  this->d_psf = new carma_obj<float>(this->current_context, dims_data2);
  this->d_amplipup = new carma_obj<cuFloatComplex>(this->current_context, dims_data2);

  cufftHandle *plan = this->d_amplipup->getPlan(); ///< FFT plan
  carmafftSafeCall(
    cufftPlan2d(plan, this->d_amplipup->getDims(1),
                this->d_amplipup->getDims(2), CUFFT_C2C));

  dims_data1[1] = IFnz;

  carma_obj<float> d_val(this->current_context, dims_data1, IFvalue);
  carma_obj<int> d_row(this->current_context, dims_data1, IFrowind);
  dims_data1[1] = this->nactus-2 + 1;
  carma_obj<int> d_col(this->current_context, dims_data1, IFcolind);
  dims_data2[1] = nactus-2;
  dims_data2[2] = Npts;
  this->d_IF = new carma_sparse_obj<float>(this->current_context, dims_data2,
      d_val.getData(), d_row.getData(), d_col.getData(), IFnz, false);
  dims_data2[1] = 2;
  dims_data2[2] = Npts;
  this->d_TT = new carma_obj<float>(this->current_context, dims_data2, TT);

  if(strcmp(type,"Vii") == 0) {
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size;
    this->d_term1 = new carma_obj<float>(this->current_context, dims_data2);
    this->d_term2 = new carma_obj<float>(this->current_context, dims_data2);
    this->d_Dphi = new carma_obj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_newmodek = new carma_obj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_pupfft = new carma_obj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_otftel = new carma_obj<float>(this->current_context, dims_data2);
    this->d_otfVii = new carma_obj<float>(this->current_context, dims_data2);
    this->d_mask = new carma_obj<float>(this->current_context, dims_data2);
    dims_data2[1] = nactus;
    dims_data2[2] = nmodes;
    this->d_Btt = new carma_obj<float>(this->current_context, dims_data2, Btt);
    dims_data2[1] = nmodes;
    this->d_covmodes = new carma_obj<float>(this->current_context, dims_data2, covmodes);
    dims_data1[1] = nmodes;
    this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);
    //FFT plans
    carmafftSafeCall(
      cufftPlan2d(this->d_pupfft->getPlan(), this->d_pupfft->getDims(1),
                  this->d_pupfft->getDims(2), CUFFT_C2C));
    carmafftSafeCall(
      cufftPlan2d(this->d_newmodek->getPlan(), this->d_newmodek->getDims(1),
                  this->d_newmodek->getDims(2), CUFFT_C2C));
    carmafftSafeCall(
      cufftPlan2d(this->d_Dphi->getPlan(), this->d_Dphi->getDims(1),
                  this->d_Dphi->getDims(2), CUFFT_C2C));

    this->d_amplipup_ngpu.push_back(this->d_amplipup);
    this->d_newmodek_ngpu.push_back(this->d_newmodek);
    this->d_Btt_ngpu.push_back(this->d_Btt);
    this->d_covmodes_ngpu.push_back(this->d_covmodes);
    this->d_term1_ngpu.push_back(this->d_term1);
    this->d_term2_ngpu.push_back(this->d_term2);
    this->d_IF_ngpu.push_back(this->d_IF);
    this->d_TT_ngpu.push_back(this->d_TT);
    this->d_phase_ngpu.push_back(this->d_phase);
    this->d_wherephase_ngpu.push_back(this->d_wherephase);
    this->d_pupfft_ngpu.push_back(this->d_pupfft);
    this->d_Dphi_ngpu.push_back(this->d_Dphi);

    for (int d=1; d < ngpu; d++) {
      current_context->set_activeDevice(d,1);
      dims_data2[1] = fft_size;
      dims_data2[2] = fft_size;
      d_amplipup_ngpu.push_back(new carma_obj<cuFloatComplex>(this->current_context, dims_data2));
      cufftHandle *plan = this->d_amplipup_ngpu[d]->getPlan(); ///< FFT plan
      carmafftSafeCall(
        cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));
      d_term1_ngpu.push_back(new carma_obj<float>(this->current_context, dims_data2));
      d_term2_ngpu.push_back(new carma_obj<float>(this->current_context, dims_data2));
      d_Dphi_ngpu.push_back(new carma_obj<cuFloatComplex>(this->current_context, dims_data2));
      d_newmodek_ngpu.push_back(new carma_obj<cuFloatComplex>(this->current_context, dims_data2));
      d_pupfft_ngpu.push_back(new carma_obj<cuFloatComplex>(this->current_context, dims_data2));
      carmafftSafeCall(
        cufftPlan2d(this->d_pupfft_ngpu[d]->getPlan(),
                    dims_data2[1], dims_data2[2], CUFFT_C2C));
      carmafftSafeCall(
        cufftPlan2d(this->d_newmodek_ngpu[d]->getPlan(),
                    dims_data2[1], dims_data2[2], CUFFT_C2C));
      carmafftSafeCall(
        cufftPlan2d(this->d_Dphi_ngpu[d]->getPlan(),
                    dims_data2[1], dims_data2[2], CUFFT_C2C));
      dims_data2[1] = nactus;
      dims_data2[2] = nmodes;
      d_Btt_ngpu.push_back(new carma_obj<float>(this->current_context, dims_data2, Btt));
      dims_data2[1] = nmodes;
      d_covmodes_ngpu.push_back(new carma_obj<float>(this->current_context, dims_data2, covmodes));
      dims_data2[1] = nactus-2;
      dims_data2[2] = Npts;

      dims_data1[1] = IFnz;
      carma_obj<float> d_val_tmp(this->current_context, dims_data1, IFvalue);
      carma_obj<int> d_row_tmp(this->current_context, dims_data1, IFrowind);
      dims_data1[1] = this->nactus-2 + 1;
      carma_obj<int> d_col_tmp(this->current_context, dims_data1, IFcolind);

      d_IF_ngpu.push_back(new carma_sparse_obj<float>(this->current_context, dims_data2,
                          d_val_tmp.getData(), d_row_tmp.getData(), d_col_tmp.getData(), IFnz, false));

      dims_data2[1] = 2;
      dims_data2[2] = Npts;
      d_TT_ngpu.push_back(new carma_obj<float>(this->current_context, dims_data2, TT));
      dims_data1[1] = Npts;
      d_wherephase_ngpu.push_back(new carma_obj<int>(this->current_context,dims_data1,wherephase));
      d_phase_ngpu.push_back(new carma_obj<float>(this->current_context,dims_data1));


    }


  }



}

sutra_gamora::~sutra_gamora() {
  this->current_context->set_activeDevice(this->device,1);
  if(this->d_err)
    delete this->d_err;

  if(this->d_amplipup) {
    for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
           this->d_amplipup_ngpu.begin(); this->d_amplipup_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
  }

  if(this->d_psf)
    delete this->d_psf;

  if(this->d_pupfft) {
    for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
           this->d_pupfft_ngpu.begin(); this->d_pupfft_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
  }

  if(this->d_phase) {
    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_phase_ngpu.begin(); this->d_phase_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
  }

  if(this->d_wherephase) {
    for (std::vector<carma_obj<int> *>::iterator it =
           this->d_wherephase_ngpu.begin(); this->d_wherephase_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
  }

  if(this->d_IF) {
    for (std::vector<carma_sparse_obj<float> *>::iterator it =
           this->d_IF_ngpu.begin(); this->d_IF_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
  }

  if(this->d_term1) {
    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_term1_ngpu.begin(); this->d_term1_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }

    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_term2_ngpu.begin(); this->d_term2_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
    for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
           this->d_newmodek_ngpu.begin(); this->d_newmodek_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
    delete this->d_otftel;
    delete this->d_otfVii;
    delete this->d_mask;
    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_Btt_ngpu.begin(); this->d_Btt_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_covmodes_ngpu.begin(); this->d_covmodes_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
    delete this->h_eigenvals;
    for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
           this->d_Dphi_ngpu.begin(); this->d_Dphi_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }
    for (std::vector<carma_obj<float> *>::iterator it =
           this->d_TT_ngpu.begin(); this->d_TT_ngpu.end() != it;
         ++it) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      delete *it;
    }

  }
}

int sutra_gamora::psf_rec_roket(float *err) {
  this->current_context->set_activeDevice(this->device,1);
  // Get the error command buffer
  this->d_err->host2device(err);
  // set psf to 0
  carmaSafeCall(
    cudaMemset(this->d_psf->getData(), 0,
               sizeof(float) * this->d_psf->getNbElem()));

  for(int cc = 0; cc<this->niter; cc++) {
    // set amplipup to 0
    carmaSafeCall(
      cudaMemset(this->d_amplipup->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup->getNbElem()));
    // Apply iter #cc on the DM to get residual phase
    carma_gemv(this->current_context->get_cusparseHandle(),'t',1.0f,this->d_IF,this->d_err->getDataAt(cc*this->nactus),0.0f,this->d_phase->getData());
    carma_gemv(this->current_context->get_cublasHandle(),'t',this->d_TT->getDims(1),this->d_TT->getDims(2),
               1.0f,this->d_TT->getData(),this->d_TT->getDims(1),
               this->d_err->getDataAt(cc*this->nactus+(this->nactus-2)),1,1.0f,this->d_phase->getData(),1);

    // complex amplitude in the pupil put in the fft support
    fill_amplipup(this->d_amplipup->getData(), this->d_phase->getData(), this->d_wherephase->getData(),
                  this->scale, this->Npts,  this->size, this->d_amplipup->getDims()[1], 0, this->current_context->get_device(this->device));
    // complex amplitude in the focal plane
    carma_fft(this->d_amplipup->getData(), this->d_amplipup->getData(), 1,
              *this->d_amplipup->getPlan());
    // take square modulus and add it a the psf LE
    cumulpsf(this->d_psf->getData(), this->d_amplipup->getData(),
             this->d_psf->getDims(1) * this->d_psf->getDims(2), current_context->get_device(device));

    printf("\rComputing and stacking %d PSFs : %d%%",this->niter,(cc*100/this->niter));
  }
  // rescale because of fft and number of iterations
  this->d_psf->scale(1.0f/(this->d_wherephase->getDims(1)*this->d_wherephase->getDims(1)),1);
  this->d_psf->scale(1.0f/this->niter,1);
  //DEBUG_TRACE("%d %d",this->Npts*this->Npts*this->niter,this->niter);
  return EXIT_SUCCESS;
}

int sutra_gamora::psf_rec_Vii() {
// Telescope OTF computation and mask
  // Get the pupil
  this->current_context->set_activeDevice(this->device,1);

  printf("Computing Telescope OTF and corresponding mask...\n");
  carmaSafeCall(
    cudaMemset(this->d_pupfft->getData(), 0,
               sizeof(cuFloatComplex) * this->d_pupfft->getNbElem()));
  carmaSafeCall(
    cudaMemset(this->d_otftel->getData(), 0,
               sizeof(float) * this->d_otftel->getNbElem()));

  fill_amplipup(this->d_pupfft->getData(), this->d_term1->getData(), this->d_wherephase->getData(),
                1.0f, this->Npts,  this->size, this->d_pupfft->getDims(1), 1, this->current_context->get_device(this->device));
  // compute fft(pupil)

  carma_fft(this->d_pupfft->getData(), this->d_pupfft->getData(), 1,
            *this->d_pupfft->getPlan());
  // compute pupfft * conjpupfft as abs(pupfft)**2 and store it in
  // the real part of d_amplipup
  abs2complex(this->d_amplipup->getData(),this->d_pupfft->getData(),
              this->d_pupfft->getNbElem(),this->current_context->get_device(this->device));
  // compute ifft(pupfft*conjpupfft)

  carma_fft(this->d_amplipup->getData(), this->d_amplipup->getData(), -1,
            *this->d_amplipup->getPlan());
  ifftscale(this->d_amplipup->getData(),1.0f/this->d_amplipup->getNbElem(),
            this->d_amplipup->getNbElem(),this->current_context->get_device(this->device));
  // Get telescope OTF as real part of amplipup
  real(this->d_otftel->getData(),this->d_amplipup->getData(),
       this->d_amplipup->getNbElem(),this->current_context->get_device(this->device));
  // Deduce the mask from telescope OTF
  fill_mask(this->d_mask->getData(),this->d_otftel->getData(),
            this->d_mask->getNbElem(), this->Npts, this->current_context->get_device(this->device));
  printf("Done\n");

// OTF|| computation with Vii algorithm
  carmaSafeCall(
    cudaMemset(this->d_Dphi->getData(), 0,
               sizeof(cuFloatComplex) * this->d_Dphi->getNbElem()));
  carmaSafeCall(
    cudaMemset(this->d_otfVii->getData(), 0,
               sizeof(float) * this->d_otfVii->getNbElem()));

  // Diagonalisation of covmodes
  carma_syevd<float,1>('V', this->d_covmodes, this->h_eigenvals);

  for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
         this->d_pupfft_ngpu.begin(); this->d_pupfft_ngpu.end() != it;
       ++it) {
    carma_obj<cuFloatComplex> *tmp_pupfft = this->d_pupfft;
    if (*it != tmp_pupfft) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      (*it)->copyFrom(tmp_pupfft->getData(), tmp_pupfft->getNbElem());
    }
  }
  for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
         this->d_Dphi_ngpu.begin(); this->d_Dphi_ngpu.end() != it;
       ++it) {
    if (*it != this->d_Dphi) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      carmaSafeCall(
        cudaMemset((*it)->getData(), 0,
                   sizeof(float) * (*it)->getNbElem()));
    }
  }
  for (std::vector<carma_obj<float> *>::iterator it =
         this->d_covmodes_ngpu.begin(); this->d_covmodes_ngpu.end() != it;
       ++it) {
    carma_obj<float> *tmp_cov = this->d_covmodes;
    if (*it != tmp_cov) {
      current_context->set_activeDevice((*it)->getDevice(),1);
      (*it)->copyFrom(tmp_cov->getData(), tmp_cov->getNbElem());
    }
  }


  // Loop on the modes to compute OTF from Vii on the fly
  std::cout << "Computing Vii: " << std::endl;
  auto progress = carma_utils::ProgressBar(this->nmodes);
  for (int k = 0; k < this->nmodes; k++) {
    compute_Dphi_on_mode_k(k);
    progress.update();
    //printf("\rComputing OTF with %d Vii : %d%%",this->nmodes,(k*100/this->nmodes));
  }
  progress.finish();
  this->current_context->set_activeDevice(this->device,1);

  carma_obj<cuFloatComplex> *tmp_vector = new carma_obj<cuFloatComplex>(current_context,
      this->d_Dphi->getDims());

  cuFloatComplex alpha;
  alpha.x = 1.0f;
  alpha.y = 0.0f;
  for (std::vector<carma_obj<cuFloatComplex> *>::iterator it =
         this->d_Dphi_ngpu.begin(); this->d_Dphi_ngpu.end() != it; ++it) {
    if (*it != d_Dphi) {
      if (current_context->canP2P(d_Dphi->getDevice(),
                                  (*it)->getDevice())) {
        d_Dphi->axpy(alpha, (*it), 1, 1);
      } else {
        tmp_vector->copyFrom((*it)->getData(), (*it)->getNbElem());
        d_Dphi->axpy(alpha, tmp_vector, 1, 1);
      }
    }
  }
  delete tmp_vector;

  // ifft(2*Dphi)
  carma_fft(this->d_Dphi->getData(), this->d_Dphi->getData(), -1,
            *this->d_Dphi->getPlan());
  ifftscale(this->d_Dphi->getData(),1.0f/this->d_Dphi->getNbElem(),
            this->d_Dphi->getNbElem(),this->current_context->get_device(this->device));

  // OTF = exp(-0.5*real(ifft(2*Dphi)) * mask * scale**2 / otftel) * mask
  computeOTFvii(this->d_otfVii->getData(),this->d_Dphi->getData(),
                this->d_otftel->getData(),this->d_mask->getData(),
                this->scale,this->d_otfVii->getNbElem(),
                this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}


void sutra_gamora::compute_Dphi_on_mode_k(int k) {
  this->current_context->set_activeDevice(this->device,1);
  int ngpu = d_pupfft_ngpu.size();
  if(ngpu < 2) {
    carmaSafeCall(
      cudaMemset(this->d_amplipup->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup->getNbElem()));
    carmaSafeCall(
      cudaMemset(this->d_newmodek->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_newmodek->getNbElem()));
    // Get the new mode shape from d_U, IF and Btt
    carma_gemv<float>(this->current_context->get_cublasHandle(), 'n', this->nactus, this->nmodes, 1.0f,
                      this->d_Btt->getData(), this->nactus,
                      this->d_covmodes->getDataAt(k*this->nmodes),1, 0.0f,
                      this->d_term1->getData(), 1);
    carma_gemv(this->current_context->get_cusparseHandle(),'t',1.0f,
               this->d_IF,this->d_term1->getData(),0.0f,
               this->d_phase->getData());

    carma_gemv(this->current_context->get_cublasHandle(),'t',this->d_TT->getDims(1),this->d_TT->getDims(2),
               1.0f,this->d_TT->getData(),this->d_TT->getDims(1),
               this->d_term1->getDataAt(this->nactus-2),1,1.0f,this->d_phase->getData(),1);
    fill_amplipup(this->d_newmodek->getData(), this->d_phase->getData(), this->d_wherephase->getData(),
                  1.0f, this->Npts,  this->size, this->d_newmodek->getDims(1), 2, this->current_context->get_device(this->device));
    // Compute term2 = abs(fft(newmode))**2

    carma_fft(this->d_newmodek->getData(), this->d_amplipup->getData(), 1,
              *this->d_amplipup->getPlan());
    modulus2(this->d_term2->getData(),this->d_amplipup->getData(),
             this->d_term2->getNbElem(),this->current_context->get_device(this->device));
    //Compute term1 = real(fft(newmodek**2)*conjpupfft)
    pow2(this->d_newmodek->getData(),this->d_newmodek->getData(),
         this->d_newmodek->getNbElem(),this->current_context->get_device(this->device));
    carma_fft(this->d_newmodek->getData(), this->d_newmodek->getData(), 1,
              *this->d_newmodek->getPlan());
    fill_term1(this->d_term1->getData(),this->d_newmodek->getData(), this->d_pupfft->getData(),
               this->d_term1->getNbElem(),this->current_context->get_device(this->device));
    // Dphi += (term1 - term2) * eigenvals[k]
    add2Dphi(this->d_Dphi->getData(),this->d_term1->getData(), this->d_term2->getData(),
             this->h_eigenvals->getData()[k],this->d_Dphi->getNbElem(),this->current_context->get_device(this->device));
  } else {
    int cur_device = k % ngpu;
    this->current_context->set_activeDevice(cur_device,1);
    carmaSafeCall(
      cudaMemset(this->d_amplipup_ngpu[cur_device]->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup_ngpu[cur_device]->getNbElem()));
    carmaSafeCall(
      cudaMemset(this->d_newmodek_ngpu[cur_device]->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_newmodek_ngpu[cur_device]->getNbElem()));
    // Get the new mode shape from d_U, IF and Btt
    carma_gemv<float>(this->current_context->get_cublasHandle(), 'n', this->nactus, this->nmodes, 1.0f,
                      this->d_Btt_ngpu[cur_device]->getData(), this->nactus,
                      this->d_covmodes_ngpu[cur_device]->getDataAt(k*this->nmodes),1, 0.0f,
                      this->d_term1_ngpu[cur_device]->getData(), 1);

    carma_gemv(this->current_context->get_cusparseHandle(),'t',1.0f,
               this->d_IF_ngpu[cur_device],this->d_term1_ngpu[cur_device]->getData(),0.0f,
               this->d_phase_ngpu[cur_device]->getData());

    carma_gemv(this->current_context->get_cublasHandle(),'t',this->d_TT_ngpu[cur_device]->getDims(1),this->d_TT_ngpu[cur_device]->getDims(2),
               1.0f,this->d_TT_ngpu[cur_device]->getData(),this->d_TT_ngpu[cur_device]->getDims(1),
               this->d_term1_ngpu[cur_device]->getDataAt(this->nactus-2),1,1.0f,this->d_phase_ngpu[cur_device]->getData(),1);

    fill_amplipup(this->d_newmodek_ngpu[cur_device]->getData(), this->d_phase_ngpu[cur_device]->getData(), this->d_wherephase_ngpu[cur_device]->getData(),
                  1.0f, this->Npts,  this->size, this->d_newmodek_ngpu[cur_device]->getDims(1), 2, this->current_context->get_device(cur_device));
    // Compute term2 = abs(fft(newmode))**2

    carma_fft(this->d_newmodek_ngpu[cur_device]->getData(), this->d_amplipup_ngpu[cur_device]->getData(), 1,
              *this->d_amplipup_ngpu[cur_device]->getPlan());

    modulus2(this->d_term2_ngpu[cur_device]->getData(),this->d_amplipup_ngpu[cur_device]->getData(),
             this->d_term2_ngpu[cur_device]->getNbElem(),this->current_context->get_device(cur_device));
    //Compute term1 = real(fft(newmodek**2)*conjpupfft)
    pow2(this->d_newmodek_ngpu[cur_device]->getData(),this->d_newmodek_ngpu[cur_device]->getData(),
         this->d_newmodek_ngpu[cur_device]->getNbElem(),this->current_context->get_device(cur_device));

    carma_fft(this->d_newmodek_ngpu[cur_device]->getData(), this->d_newmodek_ngpu[cur_device]->getData(), 1,
              *this->d_newmodek_ngpu[cur_device]->getPlan());

    fill_term1(this->d_term1_ngpu[cur_device]->getData(),this->d_newmodek_ngpu[cur_device]->getData(), this->d_pupfft_ngpu[cur_device]->getData(),
               this->d_term1_ngpu[cur_device]->getNbElem(),this->current_context->get_device(cur_device));
    // Dphi += (term1 - term2) * eigenvals[k]
    add2Dphi(this->d_Dphi_ngpu[cur_device]->getData(),this->d_term1_ngpu[cur_device]->getData(), this->d_term2_ngpu[cur_device]->getData(),
             this->h_eigenvals->getData()[k],this->d_Dphi_ngpu[cur_device]->getNbElem(),this->current_context->get_device(cur_device));

  }
}
