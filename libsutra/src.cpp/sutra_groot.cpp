#include <sutra_groot.h>

template<class T_data>
sutra_groot<T_data>::sutra_groot(carma_context *context, int device, int nactus,
            int nlayers, T_data gsangle, T_data *vdt,
            T_data *Htheta, T_data *L0, T_data *winddir, T_data *scale,
            T_data *pzt2tt, T_data *TTPfilter, T_data *Nact,
            T_data *xpos, T_data *ypos, T_data fc) {

    this->current_context = context;
    this->device = device;
    this->current_context->set_activeDevice(device,1);

    this->nactus = nactus;
    this->nlayers = nlayers;
    this->gsangle = gsangle;
    this->fc = fc;

    long dims_data1[2] = {1,nlayers};
    this->h_vdt = new carma_host_obj<T_data>(dims_data1, vdt, MA_PAGELOCK);
    this->h_Htheta = new carma_host_obj<T_data>(dims_data1, Htheta, MA_PAGELOCK);
    this->h_L0 = new carma_host_obj<T_data>(dims_data1, L0, MA_PAGELOCK);
    this->h_winddir = new carma_host_obj<T_data>(dims_data1, winddir, MA_PAGELOCK);
    this->h_scale = new carma_host_obj<T_data>(dims_data1, scale, MA_PAGELOCK);

    dims_data1[1] = nactus;
    this->d_xpos = new carma_obj<T_data>(context, dims_data1, xpos);
    this->d_ypos = new carma_obj<T_data>(context, dims_data1, ypos);

    long dims_data2[3] = {2, nactus, nactus};
    this->d_Nact = new carma_obj<T_data>(context, dims_data2, Nact);
    this->d_Cerr = new carma_obj<T_data>(context, dims_data2);
    this->d_TTPfilter = new carma_obj<T_data>(context, dims_data2, TTPfilter);
    dims_data2[1] = 2;
    dims_data2[2] = 2;
    this->d_TT = new carma_obj<T_data>(context, dims_data2);

    dims_data2[1] = 2;
    dims_data2[2] = nactus;
    this->d_pzt2tt = new carma_obj<T_data>(context, dims_data2, pzt2tt);

    int N = 10000;
    dims_data1[1] = N;

    this->d_tab_int_x = new carma_obj<T_data>(context, dims_data1);
    this->d_tab_int_y = new carma_obj<T_data>(context, dims_data1);

    tab_u831J0(this->d_tab_int_x->getData(), this->d_tab_int_y->getData(), N,
                this->current_context->get_device(device));

    printf("I am Groot\n");
}

template
sutra_groot<float>::sutra_groot(carma_context *context, int device, int nactus,
            int nlayers, float gsangle, float *vdt,
            float *Htheta, float *L0, float *winddir, float *scale,
            float *pzt2tt, float *TTPfilter, float *Nact,
            float *xpos, float *ypos, float fc);
template
sutra_groot<double>::sutra_groot(carma_context *context, int device, int nactus,
            int nlayers, double gsangle, double *vdt,
            double *Htheta, double *L0, double *winddir, double *scale,
            double *pzt2tt, double *TTPfilter, double *Nact,
            double *xpos, double *ypos, double fc);


template<class T_data>
sutra_groot<T_data>::~sutra_groot(){
    this->current_context->set_activeDevice(this->device,1);
    delete this->h_vdt;
    delete this->h_Htheta;
    delete this->h_L0;
    delete this->h_winddir;
    delete this->h_scale;
    delete this->d_xpos;
    delete this->d_ypos;
    delete this->d_Nact;
    delete this->d_TT;
    delete this->d_Cerr;
    delete this->d_pzt2tt;
    delete this->d_TTPfilter;
}

template
sutra_groot<float>::~sutra_groot();
template
sutra_groot<double>::~sutra_groot();

template<class T_data>
int sutra_groot<T_data>::compute_Cerr(){

    this->current_context->set_activeDevice(this->device,1);

    carmaSafeCall(
        cudaMemset(this->d_Cerr->getData(), 0,
            sizeof(T_data) * this->d_Cerr->getNbElem()));
    printf("Computing Cerr...\n");

    for(int l = 0; l < this->nlayers; l++){
        compute_Cerr_layer(this->d_Cerr->getData(), this->nactus, this->d_tab_int_x->getData(),
                            this->d_tab_int_y->getData(), this->d_xpos->getData(),
                            this->d_ypos->getData(), (*this->h_vdt)[l],
                            (*this->h_Htheta)[l], (*this->h_L0)[l], this->fc, (*this->h_winddir)[l],
                            this->gsangle, (*this->h_scale)[l], this->d_tab_int_y->getNbElem(),
                            this->current_context->get_device(device));
    }
    add_transpose(this->d_Cerr->getData(), this->nactus, this->current_context->get_device(device));
    printf("Done\n");

    printf("Applying coupling matrix...\n");
    // Coupling matrix filter
    //carma_potri(this->d_Nact);
    carma_obj<T_data> d_tmp(this->d_Cerr);
    carma_gemm(cublas_handle(), 'n', 'n', this->nactus,
                this->nactus, this->nactus, (T_data)1.,
                this->d_Nact->getData(), this->nactus,
                this->d_Cerr->getData(), this->nactus, (T_data)0.0,
                d_tmp.getData(), this->nactus);
    carma_gemm(cublas_handle(), 'n', 'n', this->nactus,
                this->nactus, this->nactus, (T_data)1.0,
                d_tmp.getData(), this->nactus,
                this->d_Nact->getData(), this->nactus, (T_data)0.0,
                this->d_Cerr->getData(), this->nactus);
    printf("Done\n");

    // Tip-tilt component
    printf("Computing TT component...\n");
    carma_obj<T_data> d_tmp2(this->d_pzt2tt);
    carma_gemm(cublas_handle(), 'n', 'n', 2,
                this->nactus, this->nactus, (T_data)1.0,
                this->d_pzt2tt->getData(), 2,
                this->d_Cerr->getData(), this->nactus, (T_data)0.0,
                d_tmp2.getData(), 2);
    carma_gemm(cublas_handle(), 'n', 't', 2,
                2, this->nactus, (T_data)1.0,
                d_tmp2.getData(), 2,
                this->d_pzt2tt->getData(), 2, (T_data)0.0,
                this->d_TT->getData(), 2);
    printf("Done\n");

    // Filtering TT + piston from Cerr
    printf("Filtering TT + piston from Cerr...\n");
    carma_gemm(cublas_handle(), 'n', 'n', this->nactus,
                this->nactus, this->nactus, (T_data)1.0,
                this->d_TTPfilter->getData(), this->nactus,
                this->d_Cerr->getData(), this->nactus, (T_data)0.0,
                d_tmp.getData(), this->nactus);
    carma_gemm(cublas_handle(), 'n', 't', this->nactus,
                this->nactus, this->nactus, (T_data)1.0,
                d_tmp.getData(), this->nactus,
                this->d_TTPfilter->getData(), this->nactus, (T_data)0.0,
                this->d_Cerr->getData(), this->nactus);
    printf("Done\n");

    return EXIT_SUCCESS;

}

template
int sutra_groot<float>::compute_Cerr();
template
int sutra_groot<double>::compute_Cerr();
