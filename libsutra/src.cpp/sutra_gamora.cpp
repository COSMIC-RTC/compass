#include <sutra_gamora.h>

sutra_gamora::sutra_gamora(carma_context *context, int device, char *type, int nactus,
            int nmodes, int niter, float *IFvalue, int *IFrowind, int *IFcolind, int IFnz,
            float *TT, float *pupil, int size, int Npts, float scale, float *Btt, float *covmodes){

    this->current_context = context;
    this->device = device;
    this->current_context->set_activeDevice(device,1);

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

    if(strcmp(type,"roket") == 0){
        // Command error
        this->d_err = new carma_obj<float>(this->current_context,dims_data2);

    }

    if(strcmp(type,"roket") == 0 || strcmp(type,"Vii") == 0){

        int wherephase[Npts];
        int cpt = 0;
        // Phase point index in spupil
        for (int cc=0; cc < size*size;cc++){
            if(pupil[cc]>0){
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
    }

    if(strcmp(type,"Vii") == 0){
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
    }

}

sutra_gamora::~sutra_gamora(){
    this->current_context->set_activeDevice(this->device,1);
    if(this->d_err)
        delete this->d_err;
    if(this->d_amplipup)
        delete this->d_amplipup;
    if(this->d_psf)
        delete this->d_psf;
    if(this->d_pupfft)
        delete this->d_pupfft;
    if(this->d_phase)
        delete this->d_phase;
    if(this->d_wherephase)
        delete this->d_wherephase;
    if(this->d_IF)
        delete this->d_IF;
    if(this->d_term1){
        delete this->d_term1;
        delete this->d_term2;
        delete this->d_newmodek;
        delete this->d_otftel;
        delete this->d_otfVii;
        delete this->d_mask;
        delete this->d_Btt;
        delete this->d_covmodes;
        delete this->h_eigenvals;
        delete this->d_Dphi;
    }
}

int sutra_gamora::psf_rec_roket(float *err){
    this->current_context->set_activeDevice(this->device,1);
    // Get the error command buffer
    this->d_err->host2device(err);
    // set psf to 0
    carmaSafeCall(
        cudaMemset(this->d_psf->getData(), 0,
            sizeof(float) * this->d_psf->getNbElem()));

    for(int cc = 0; cc<this->niter; cc++){
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
                        this->scale * (-1), this->Npts,  this->size, this->d_amplipup->getDims()[1], 0, this->current_context->get_device(this->device));
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

int sutra_gamora::psf_rec_Vii(){
// Telescope OTF computation and mask
    // Get the pupil
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
    // Loop on the modes to compute OTF from Vii on the fly
    for (int k=0 ; k < this->nmodes ; k++){
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

        printf("\rComputing OTF with %d Vii : %d%%",this->nmodes,(k*100/this->nmodes));

    }

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
