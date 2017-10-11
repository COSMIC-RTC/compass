#include <sutra_dm.h>
#include <stdio.h>

//#include <sutra_dm.cuh>


sutra_dm::sutra_dm(carma_context *context, const char* type, long dim,
                   long ninflu, long influsize, long ninflupos, long n_npoints,
                   float push4imat, long nord, int device) {

  this->d_influ = NULL;
  this->d_influ2 = NULL;
  this->d_influ3 = NULL;
  this->d_influpos = NULL;
  this->d_influpos2 = NULL;
  this->d_npoints = NULL;
  this->d_istart = NULL;
  this->d_xoff = NULL;
  this->d_yoff = NULL;
  this->d_coeffs = NULL;
  this->d_mask = NULL;
  this->d_zr = NULL;
  this->d_ztheta = NULL;
  this->d_KLbasis = NULL;
  this->d_pos = NULL;
  //this->d_IFsparse = NULL;
  //this->d_commdouble = NULL;
  //this->d_shapedouble = NULL;

  this->current_context = context;
  this->device = device;
  current_context->set_activeDevice(device,1);
  this->ninflu = ninflu;
  this->dim = dim;
  this->influsize = influsize;
  this->d_shape = new sutra_phase(context, dim);
  this->type = type;
  this->push4imat = push4imat;

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = ninflu;
  this->d_comm = new carma_obj<float>(context, dims_data1);

  if (strcmp(type, "kl") != 0) {
    long *dims_data3 = new long[4];
    dims_data3[0] = 3;
    dims_data3[1] = influsize;
    dims_data3[2] = influsize;
    dims_data3[3] = ninflu;
    this->d_influ = new carma_obj<float>(context, dims_data3);
    delete[] dims_data3;

    dims_data1[1] = ninflu;
    this->d_xoff = new carma_obj<int>(context, dims_data1);
    this->d_yoff = new carma_obj<int>(context, dims_data1);
  }

  if (strcmp(type, "pzt") == 0) {
    // Convolution version
    /*
     dims_data1[1] = 4 * ninflu;
     this->d_pos = new carma_obj<int>(context, dims_data1);
     long dims_data2[3] = {2,this->dim,this->dim};
     this->d_kernconv = new carma_obj<float>(context,dims_data2);
     this->d_mapactu = new carma_obj<float>(context,dims_data2);
     dims_data2[2] = this->dim/2 + 1;
     this->d_ftkernconv = new carma_obj<cuFloatComplex>(context,dims_data2);
     this->d_ftmapactu = new carma_obj<cuFloatComplex>(context,dims_data2);
     cufftHandle *plan = this->d_ftkernconv->getPlan();
     carmafftSafeCall(cufftPlan2d(plan, dims_data2[1], dims_data2[1], CUFFT_R2C));
     plan = this->d_ftmapactu->getPlan();
     carmafftSafeCall(cufftPlan2d(plan, dims_data2[1], dims_data2[1], CUFFT_C2R));
    	 */
    //this->d_IFsparse = 0L;
    //this->d_commdouble = new carma_obj<float>(context, dims_data1);
    dims_data1[1] = ninflupos;
    this->d_influpos = new carma_obj<int>(context, dims_data1);
    this->d_influpos2 = new carma_obj<int>(context, dims_data1);
    this->d_influ2 = new carma_obj<float>(context, dims_data1);

    int length = 0;
#if(COMPN == 2)
    length = ninflupos;
#elif(COMPN == 3)
    length = n_npoints * MAXSPOT;
#endif
    dims_data1[1] = length;
    this->d_influ3 = new carma_obj<struct tuple_t<float>>(context, dims_data1);
    //this->d_shapedouble = new carma_obj<double>(context,this->d_shape->d_screen->getDims());

    dims_data1[1] = n_npoints;
    this->d_npoints = new carma_obj<int>(context, dims_data1);
    dims_data1[1] = n_npoints + 1;
    this->d_istart = new carma_obj<int>(context, dims_data1);
  }
  if (strcmp(type, "kl") == 0) {
    this->d_kl = new sutra_kl(context, influsize, ninflupos, n_npoints, ninflu, nord,
                              device);
  }
  // Florian features
  /*
  if (strcmp(type, "flo_kl") == 0) {
    this->d_kl = new sutra_kl(context, influsize, ninflupos, n_npoints, ninflu,
                              device);
  }
*/
#ifdef TEXTURE
#ifdef CHEAT_CODE
  bindTexture(d_comm->getData(), d_influ->getData(), d_influ3->getData(), d_istart->getData(), d_npoints->getData(), d_influpos->getData());
#else
  bindTexture(d_comm->getData(), d_influ2->getData(), d_influ3->getData(), d_istart->getData(), d_npoints->getData(), d_influpos->getData());
#endif
#endif

}

sutra_dm::~sutra_dm() {
  //delete this->current_context;
  current_context->set_activeDevice(device,1);

  delete this->d_shape;
  delete this->d_comm;

  if (this->d_influ != NULL)
    delete this->d_influ;
  if (this->d_influ2 != NULL)
    delete this->d_influ2;
  if (this->d_influ3 != NULL)
    delete this->d_influ3;
  if (this->d_influpos != NULL)
    delete this->d_influpos;
  if (this->d_influpos != NULL)
    delete this->d_influpos2;
  if (this->d_npoints != NULL)
    delete this->d_npoints;
  if (this->d_istart != NULL)
    delete this->d_istart;
  if (this->d_xoff != NULL)
    delete this->d_xoff;
  if (this->d_yoff != NULL)
    delete this->d_yoff;
  if (this->d_coeffs != NULL)
    delete this->d_coeffs;
  if (this->d_mask != NULL)
    delete this->d_mask;
  if (this->d_zr != NULL)
    delete this->d_zr;
  if (this->d_ztheta != NULL)
    delete this->d_ztheta;
  if (this->d_KLbasis != NULL)
    delete this->d_KLbasis;
  if (this->d_pos != NULL)
    delete this->d_pos;
#ifdef TEXTURE
  unbindTexture();
#endif
  /*
  if (this->d_IFsparse != NULL)
  delete this->d_IFsparse;
  if (this->d_commdouble != NULL)
  delete this->d_commdouble;
  if (this->d_shapedouble != NULL)
  delete this->d_shapedouble;
   */
}

int sutra_dm::nact() {
  return this->ninflu;
}
/*
int sutra_dm::prepare_convolve() {

	current_context->set_activeDevice(device,1);
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(current_context->get_device(device),
			this->d_pos->getNbElem(), nblocks, nthreads);
	fillpos(nthreads, nblocks, this->d_pos->getData(), this->d_xoff->getData(),
			this->d_yoff->getData(), this->influsize, this->ninflu, this->dim,
			this->d_pos->getNbElem());
	carma_fft(this->d_kernconv->getData(), this->d_ftkernconv->getData(), 1,
			*this->d_ftkernconv->getPlan());

	return EXIT_SUCCESS;
}
*/
int sutra_dm::pzt_loadarrays(float *influ,float *influ2, struct tuple_t<float> *influ3,  int *influpos, int *influpos2, int *npoints,
                             int *istart, int *xoff, int *yoff) {
  current_context->set_activeDevice(device,1);
  this->d_influ->host2device(influ);
  this->d_influ2->host2device(influ2);
#if(COMPN != 1)
  this->d_influ3->host2device(influ3);
#endif
  this->d_xoff->host2device(xoff);
  this->d_yoff->host2device(yoff);
  this->d_istart->host2device(istart);
  this->d_influpos->host2device(influpos);
  this->d_influpos2->host2device(influpos2);
  this->d_npoints->host2device(npoints);
  //this->d_kernconv->host2device(kernconv);
  //this->prepare_convolve();

  /*
  int *osef;
  this->get_IF_sparse(this->d_IFsparse,osef,this->d_shape->d_screen->getNbElem(),1.0f,0);
   */

  return EXIT_SUCCESS;
}

int sutra_dm::kl_loadarrays(float *rabas, float *azbas, int *ord, float *cr,
                            float *cp) {
  current_context->set_activeDevice(device,1);
  this->d_kl->d_rabas->host2device(rabas);
  this->d_kl->d_azbas->host2device(azbas);
  this->d_kl->h_ord->fill_from(ord);
  this->d_kl->d_ord->host2device(ord);
  this->d_kl->d_cr->host2device(cr);
  this->d_kl->d_cp->host2device(cp);

  return EXIT_SUCCESS;
}

int sutra_dm::reset_shape() {
  current_context->set_activeDevice(device,1);

  carmaSafeCall(
    cudaMemset(this->d_shape->d_screen->getData(), 0,
               sizeof(float) * this->d_shape->d_screen->getNbElem()));

  return EXIT_SUCCESS;
}

int sutra_dm::comp_shape() {
  return this->comp_shape(this->d_comm->getData());
}

#ifdef CHEAT_CODE
int sutra_dm::comp_shape(float *comvec) {
  current_context->set_activeDevice(device,1);
  this->reset_shape();

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(current_context->get_device(device),
                         this->d_shape->d_screen->getNbElem(), nblocks,
                         nthreads);

  if (this->type == "pzt") {
    comp_dmshape(nthreads, nblocks, this->d_influ->getData(),
                 this->d_shape->d_screen->getData(),
                 this->d_influpos->getData(), this->d_istart->getData(),
                 this->d_npoints->getData(), comvec,
                 this->influsize * this->influsize,
                 this->d_shape->d_screen->getNbElem());
  }


  if (this->type == "tt")
    comp_fulldmshape(nthreads, nblocks, this->d_influ->getData(),
                     this->d_shape->d_screen->getData(), this->ninflu,
                     this->influsize * this->influsize, comvec,
                     this->d_shape->d_screen->getNbElem());

  if (this->type == "kl") {
    int xoff = (int) ((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)
                      / 2.0f);
    int yoff = xoff;
    this->d_kl->do_combi(comvec, this->d_shape->d_screen->getData(),
                         this->d_shape->d_screen->getDims()[1], xoff, yoff);
  }
  return EXIT_SUCCESS;
}

#else

int sutra_dm::comp_shape(float *comvec) {
  current_context->set_activeDevice(device,1);
  this->reset_shape();

#ifdef REDUCTION
  dim3 threads(MAXSPOT);
  dim3 blocks(this->d_shape->screen_size, this->d_shape->screen_size);
  int shared = MAXSPOT * sizeof(float);
#else
  dim3 threads(BLOCKSIZE);
  dim3 blocks(CEIL(this->d_shape->d_screen->getNbElem() << 2, threads.x));
  int shared = 0;
#endif

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(current_context->get_device(device),
                         this->d_shape->d_screen->getNbElem(), nblocks,
                         nthreads);

  if (this->type == "pzt") {

    comp_dmshape2<float>(this->d_shape->d_screen->getData(),
                         comvec, this->d_influ2->getData(),
#if(COMPN != 1)
                         this->d_influ3->getData(),
#else
                         NULL,
#endif
                         this->d_istart->getData(),this->d_npoints->getData(),this->d_influpos2->getData(),
                         this->d_shape->d_screen->getNbElem(),
                         threads, blocks, shared
                        );
  }


  if (this->type == "tt")
    comp_fulldmshape(nthreads, nblocks, this->d_influ->getData(),
                     this->d_shape->d_screen->getData(), this->ninflu,
                     this->influsize * this->influsize, comvec,
                     this->d_shape->d_screen->getNbElem());

  if (this->type == "kl") {
    int xoff = (int) ((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)
                      / 2.0f);
    int yoff = xoff;
    this->d_kl->do_combi(comvec, this->d_shape->d_screen->getData(),
                         this->d_shape->d_screen->getDims()[1], xoff, yoff);
  }
  return EXIT_SUCCESS;
}

#endif

int sutra_dm::comp_oneactu(int nactu, float ampli) {
  current_context->set_activeDevice(device,1);
  this->reset_shape();
  int nthreads = 0, nblocks = 0;
  //getNumBlocksAndThreads(this->device,this->dim * this->dim, nblocks, nthreads);
  getNumBlocksAndThreads(current_context->get_device(device),
                         this->influsize * this->influsize, nblocks, nthreads);
  if (this->type == "pzt")
    oneactu(nthreads, nblocks, this->d_influ->getData(),
            this->d_shape->d_screen->getData(), nactu, ampli,
            this->d_xoff->getData(), this->d_yoff->getData(), this->dim,
            this->influsize, this->influsize * this->influsize);
  if (this->type == "tt")
    oneactu(nthreads, nblocks, this->d_influ->getData(),
            this->d_shape->d_screen->getData(), nactu, ampli, this->dim,
            this->influsize, this->influsize * this->influsize);
  if (this->type == "kl") {
    int xoff = (int) ((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)
                      / 2.0f);
    int yoff = xoff;
    this->d_kl->do_compute(ampli, this->d_shape->d_screen->getData(), nactu,
                           this->d_shape->d_screen->getDims()[1], xoff, yoff);
  }

  return EXIT_SUCCESS;
}

template<class T>
int sutra_dm::get_IF(T *IF, int *indx_pup, long nb_pts, float ampli) {

  for (int i = 0; i < this->ninflu; i++) {
    this->comp_oneactu(i, ampli);
    getIF<T>(IF, this->d_shape->d_screen->getData(), indx_pup, nb_pts, i,
             this->ninflu, 1, current_context->get_device(device));
  }

  this->reset_shape();

  return EXIT_SUCCESS;
}
template int
sutra_dm::get_IF<float>(float *IF, int *indx_pup, long nb_pts, float ampli);
template int
sutra_dm::get_IF<double>(double *IF, int *indx_pup, long nb_pts, float ampli);

template<class T>
int sutra_dm::get_IF_sparse(carma_sparse_obj<T> *&d_IFsparse, int *indx_pup,
                            long nb_pts, float ampli, int puponly) {

  current_context->set_activeDevice(device,1);
  int nnz_tot = 0;
  float *values[this->ninflu];
  int *colind[this->ninflu];
  int NZ[this->ninflu];
  long dims_data2[3] = { 2, 1, nb_pts };
  carma_obj<T> d_IF(current_context, dims_data2);
  carma_sparse_obj<T> *d_IFsparse_vec;

  std::cout << "Computing IF sparse..." << std::endl;
  for (int i = 0; i < this->ninflu; i++) {
    //Compute and store IF for actu i in d_IF
    this->comp_oneactu(i, ampli);
    getIF<T>(d_IF.getData(), this->d_shape->d_screen->getData(), indx_pup,
             nb_pts, 0, this->ninflu, puponly,
             this->current_context->get_device(device));
    //CUsparse d_IF
    d_IFsparse_vec = new carma_sparse_obj<T>(&d_IF);
    // Retrieve nnz, values and colind from d_IFsparse_vec, stored on CPU
    //DEBUG_TRACE("nnz : %d \n",d_IFsparse_vec->getNzElem());

    NZ[i] = d_IFsparse_vec->getNzElem();
    values[i] = (float*) malloc(NZ[i] * sizeof(T));
    colind[i] = (int*) malloc(NZ[i] * sizeof(int));

    carmaSafeCall(
      cudaMemcpyAsync(values[i], d_IFsparse_vec->getData(), sizeof(T) * NZ[i],
                      cudaMemcpyDeviceToHost));
    carmaSafeCall(
      cudaMemcpyAsync(colind[i], d_IFsparse_vec->d_colind,
                      sizeof(int) * NZ[i], cudaMemcpyDeviceToHost));

    nnz_tot += NZ[i];

    delete d_IFsparse_vec;
  }
  //Reconstruction of d_data, d_colind, d_rowind for IFsparse
  long dims_data[2] = { 1, nnz_tot };
  carma_obj<T> d_val(current_context, dims_data);
  carma_obj<int> d_col(current_context, dims_data);
  dims_data[1] = this->ninflu + 1;
  carma_obj<int> d_row(current_context, dims_data);
  int cpt[this->ninflu + 1];
  cpt[0] = 0;

  for (int i = 0; i < this->ninflu; i++) {
    carmaSafeCall(
      cudaMemcpyAsync(d_val.getDataAt(cpt[i]), values[i], sizeof(T) * NZ[i],
                      cudaMemcpyHostToDevice));
    carmaSafeCall(
      cudaMemcpyAsync(d_col.getDataAt(cpt[i]), colind[i], sizeof(int) * NZ[i],
                      cudaMemcpyHostToDevice));
    cpt[i + 1] = cpt[i] + NZ[i];
  }
  carmaSafeCall(
    cudaMemcpyAsync(d_row.getData(), cpt, sizeof(int) * (this->ninflu + 1),
                    cudaMemcpyHostToDevice));
  dims_data2[1] = this->ninflu;

  d_IFsparse = new carma_sparse_obj<T>(current_context, dims_data2,
                                       d_val.getData(), d_col.getData(),
                                       d_row.getData(), nnz_tot, false);

  return EXIT_SUCCESS;
}
template int
sutra_dm::get_IF_sparse<float>(carma_sparse_obj<float> *&d_IFsparse,
                               int *indx_pup, long nb_pts, float ampli,
                               int puponly);
template int
sutra_dm::get_IF_sparse<double>(carma_sparse_obj<double> *&d_IFsparse,
                                int *indx_pup, long nb_pts, float ampli,
                                int puponly);

int sutra_dm::compute_KLbasis(float *xpos, float *ypos, int *indx, long dim,
                              float norm, float ampli) {

  current_context->set_activeDevice(device,1);
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = this->ninflu;
  dims_data[2] = this->ninflu;
  carma_obj<float> *d_statcov = new carma_obj<float>(current_context,
      dims_data);
  this->d_KLbasis = new carma_obj<float>(current_context, dims_data);
  long dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = dim;
  carma_obj<int> *d_indx = new carma_obj<int>(this->current_context,
      dims_data2);

  //Compute the statistic matrix from actuators positions & Kolmogorov statistic
  dims_data2[1] = this->ninflu;
  carma_obj<float> *d_xpos = new carma_obj<float>(current_context, dims_data2);
  carma_obj<float> *d_ypos = new carma_obj<float>(current_context, dims_data2);

  d_xpos->host2device(xpos);
  d_ypos->host2device(ypos);
  dm_dostatmat(d_statcov->getData(), this->ninflu, d_xpos->getData(),
               d_ypos->getData(), norm, current_context->get_device(device));

  delete d_xpos;
  delete d_ypos;

  // Compute and apply piston filter
  this->piston_filt(d_statcov);

  dims_data[1] = this->ninflu;
  dims_data[2] = this->ninflu;
  carma_obj<float> *d_geocov = new carma_obj<float>(current_context, dims_data);
  d_indx->host2device(indx);

  //Sparse version for geomat
  carma_sparse_obj<double> *d_IFsparsepup;
  carma_obj<double> *d_geodouble = new carma_obj<double>(current_context,
      d_geocov->getDims());
  this->get_IF_sparse<double>(d_IFsparsepup, d_indx->getData(), dim, ampli, 1);
  this->do_geomatFromSparse(d_geodouble->getData(), d_IFsparsepup);
  doubletofloat(d_geodouble->getData(), d_geocov->getData(),
                d_geocov->getNbElem(), current_context->get_device(device));

  delete d_geodouble;
  delete d_IFsparsepup;
  //Dense version for geomat
  /*
  dims_data[1] = dim;
  dims_data[2] = this->ninflu;

  carma_obj<float> *d_IF = new carma_obj<float>(this->current_context, dims_data);
  // Get influence functions of the DM
  this->get_IF(d_IF->getData(),d_indx->getData(),dim,ampli);
  // Compute geometric matrix (to be CUsparsed)
  this->do_geomat(d_geocov->getData(),d_IF->getData(),dim);
   */

  // Double diagonalisation to obtain KL basis on actuators
  this->DDiago(d_statcov, d_geocov);

  delete d_geocov;
  delete d_indx;

  return EXIT_SUCCESS;
}

template<class T>
int sutra_dm::do_geomatFromSparse(T *d_geocov,
                                  carma_sparse_obj<T> *d_IFsparse) {

  current_context->set_activeDevice(device,1);
  carma_sparse_obj<T> *d_tmp = new carma_sparse_obj<T>(this->current_context);

  carma_gemm<T>(cusparse_handle(), 'n', 't', d_IFsparse, d_IFsparse, d_tmp);
  carma_csr2dense<T>(d_tmp, d_geocov);

  delete d_tmp;
  return EXIT_SUCCESS;
}
template int
sutra_dm::do_geomatFromSparse<float>(float *d_geocov,
                                     carma_sparse_obj<float> *d_IFsparse);
template int
sutra_dm::do_geomatFromSparse<double>(double *d_geocov,
                                      carma_sparse_obj<double> *d_IFsparse);

int sutra_dm::do_geomat(float *d_geocov, float *d_IF, long n_pts) {
  current_context->set_activeDevice(device,1);
  carma_gemm(cublas_handle(), 't', 'n', this->ninflu, this->ninflu, n_pts, 1.0f,
             d_IF, n_pts, d_IF, n_pts, 0.0f, d_geocov, this->ninflu);

  return EXIT_SUCCESS;
}

int sutra_dm::piston_filt(carma_obj<float> *d_statcov) {
  current_context->set_activeDevice(device,1);
  long Nmod = d_statcov->getDims()[1];
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  carma_obj<float> *d_F = new carma_obj<float>(current_context, dims_data);
  carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);

  int N = d_statcov->getDims()[1] * d_statcov->getDims()[1];
  fill_filtermat(d_F->getData(), Nmod, N, current_context->get_device(device));

  carma_gemm(cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f, d_F->getData(),
             Nmod, d_statcov->getData(), Nmod, 0.0f, d_tmp->getData(), Nmod);
  carma_gemm(cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
             d_tmp->getData(), Nmod, d_F->getData(), Nmod, 0.0f,
             d_statcov->getData(), Nmod);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

int sutra_dm::set_comkl(float *comvec) {
  current_context->set_activeDevice(device,1);
  if (this->d_KLbasis == NULL) {
    std::cerr << "KLbasis has to be computed before calling this function" << std::endl;
    return EXIT_FAILURE;
  } else {
    long dims_data[2] = { 1, this->ninflu };
    carma_obj<float> d_comkl(current_context, dims_data, comvec);

    carma_gemv(cublas_handle(), 'n', ninflu, ninflu, 1.0f,
               this->d_KLbasis->getData(), this->d_KLbasis->getDims()[1],
               d_comkl.getData(), 1, 0.0f, this->d_comm->getData(), 1);

    return EXIT_SUCCESS;
  }
}

int sutra_dm::DDiago(carma_obj<float> *d_statcov, carma_obj<float> *d_geocov) {
  current_context->set_activeDevice(device,1);
  const long dims_data[3] = { 2, this->ninflu, this->ninflu };
  carma_obj<float> *d_M1 = new carma_obj<float>(current_context, dims_data);
  carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
  carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context, dims_data);

  const long dims_data2[2] = { 1, this->ninflu };
  carma_obj<float> *d_eigenvals = new carma_obj<float>(current_context,
      dims_data2);
  carma_obj<float> *d_eigenvals_sqrt = new carma_obj<float>(current_context,
      dims_data2);
  carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,
      dims_data2);
  carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2,
      MA_PAGELOCK);
  carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(
    dims_data2, MA_PAGELOCK);
  carma_host_obj<float> *h_eigenvals_sqrt = new carma_host_obj<float>(
    dims_data2, MA_PAGELOCK);

  // 1. SVdec(geocov,U) --> Ut * geocov * U = D������
  carma_syevd<float, 1>('V', d_geocov, h_eigenvals);

  d_eigenvals->host2device(*h_eigenvals);
  for (int i = 0; i < this->ninflu; i++) {
    h_eigenvals_sqrt->getData()[i] = sqrt(h_eigenvals->getData()[i]); // D = sqrt(D������)
    h_eigenvals_inv->getData()[i] = 1. / sqrt(h_eigenvals->getData()[i]); // D��������������� = 1/sqrt(D������)
  }
  d_eigenvals_sqrt->host2device(*h_eigenvals_sqrt);
  d_eigenvals_inv->host2device(*h_eigenvals_inv);

  // 2. M��������������� = sqrt(eigenvals) * Ut : here, we have transpose(M���������������)
  /*
  carma_dgmm<float>(cublas_handle(),CUBLAS_SIDE_RIGHT,this->ninflu,this->ninflu,
  d_geocov->getData(), this->ninflu, d_eigenvals_inv->getData(),1,
  d_M1->getData(), this->ninflu);*/

  carma_dgmm<float>(cublas_handle(), CUBLAS_SIDE_RIGHT, this->ninflu,
                    this->ninflu, d_geocov->getData(), this->ninflu,
                    d_eigenvals_sqrt->getData(), 1, d_M1->getData(),
                    this->ninflu);

  // 3. C' = M��������������� * statcov * M���������������t
  carma_gemm<float>(cublas_handle(), 't', 'n', ninflu, ninflu, ninflu, 1.0f,
                    d_M1->getData(), ninflu, d_statcov->getData(), ninflu, 0.0f,
                    d_tmp->getData(), ninflu);

  carma_gemm<float>(cublas_handle(), 'n', 'n', ninflu, ninflu, ninflu, 1.0f,
                    d_tmp->getData(), ninflu, d_M1->getData(), ninflu, 0.0f,
                    d_tmp2->getData(), ninflu);

  // 4. SVdec(C',A)
  carma_syevd<float, 1>('V', d_tmp2, h_eigenvals);

  // 5. M = U * D���������������
  carma_dgmm<float>(cublas_handle(), CUBLAS_SIDE_RIGHT, this->ninflu,
                    this->ninflu, d_geocov->getData(), this->ninflu,
                    d_eigenvals_inv->getData(), 1, d_tmp->getData(),
                    this->ninflu);

  // 6. B = M * A;
  carma_gemm<float>(cublas_handle(), 'n', 'n', ninflu, ninflu, ninflu, 1.0f,
                    d_tmp->getData(), ninflu, d_tmp2->getData(), ninflu, 0.0f,
                    d_KLbasis->getData(), ninflu);

  delete d_M1;
  delete d_tmp;
  delete d_tmp2;
  delete d_eigenvals;
  delete d_eigenvals_sqrt;
  delete d_eigenvals_inv;
  delete h_eigenvals;
  delete h_eigenvals_sqrt;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
sutra_dms::sutra_dms(int ndm) {
}

sutra_dms::~sutra_dms() {
  for (std::vector<sutra_dm *>::iterator it = this->d_dms.begin();
       this->d_dms.end() != it; it++) {
    delete *it;
  }
  this->d_dms.clear();
  this->d_type.clear();
}

int sutra_dms::add_dm(carma_context *context, const char* type, float alt,
                      long dim, long ninflu, long influsize, long ninflupos,
                      long n_npoints, float push4imat, long nord, int device) {
  if(get_inddm(type, alt)>0)
    return EXIT_FAILURE;
  d_dms.push_back(
    new sutra_dm(context, type, dim, ninflu, influsize, ninflupos, n_npoints,
                 push4imat, nord, device));
  d_type.push_back(std::make_pair(type, alt));

  return EXIT_SUCCESS;
}

int sutra_dms::get_inddm(string type, float alt) {
  type_screen dm_type = make_pair(type, alt);
  std::vector<type_screen>::iterator it;

  //  it = std::find(d_type.begin(), d_type.end(), dm_type);
  //  if (it != d_type.end())
  //    return it-d_type.begin();

  for (size_t idx = 0; idx < this->d_dms.size(); idx++) {
    if( this->d_type[idx].first==dm_type.first && this->d_type[idx].second==dm_type.second)
      return idx;
  }

  return -1;
}

int sutra_dms::remove_dm(string type, float alt) {
  int pos = get_inddm(type, alt);
  if(pos<0)
    return EXIT_FAILURE;
  sutra_dm *dm=d_dms.at(pos);
  d_dms.erase(d_dms.begin()+pos);
  d_type.erase(d_type.begin()+pos);
  delete dm;
  return EXIT_SUCCESS;
}

// Florian features
int sutra_dm::kl_floloadarrays(float *covmat, float *filter, float *evals,
                               float *bas) {
  current_context->set_activeDevice(device,1);
  this->d_kl->d_covmat->host2device(covmat);
  this->d_kl->d_filter->host2device(filter);
  this->d_kl->d_bas->host2device(bas);
  this->d_kl->d_evals->host2device(evals);

  return EXIT_SUCCESS;
}

int sutra_dms::nact_total() {
  int nact=0;
  for (size_t idx = 0; idx < this->d_dms.size(); idx++) {
    nact += d_dms[idx]->ninflu;
  }
  return nact;
}
