#include <sutra_dm.h>

sutra_dm::sutra_dm(carma_context *context, const char* type, long dim,
    long ninflu, long influsize, long ninflupos, long n_npoints,
    float push4imat, int device) {

  this->d_influ = NULL;
  this->d_influpos = NULL;
  this->d_npoints = NULL;
  this->d_istart = NULL;
  this->d_xoff = NULL;
  this->d_yoff = NULL;
  this->d_coeffs = NULL;
  this->d_mask = NULL;
  this->d_zr = NULL;
  this->d_ztheta = NULL;
  this->d_KLbasis = NULL;

  this->current_context = context;
  this->ninflu = ninflu;
  this->dim = dim;
  this->influsize = influsize;
  this->d_shape = new sutra_phase(context, dim);
  this->type = type;
  this->push4imat = push4imat;
  this->device = device;

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
    dims_data1[1] = ninflupos;
    this->d_influpos = new carma_obj<int>(context, dims_data1);

    dims_data1[1] = n_npoints;
    this->d_npoints = new carma_obj<int>(context, dims_data1);
    this->d_istart = new carma_obj<int>(context, dims_data1);
  }
  if (strcmp(type, "kl") == 0) {
    this->d_kl = new sutra_kl(context, influsize, ninflupos, n_npoints, ninflu,
        device);
  }
  // Florian features
  if (strcmp(type, "flo_kl") == 0) {
    this->d_kl = new sutra_kl(context, influsize, ninflupos, n_npoints, ninflu,
        device);
  }

}

sutra_dm::~sutra_dm() {
  //delete this->current_context;

  delete this->d_shape;
  delete this->d_comm;

  if (this->d_influ != NULL)
    delete this->d_influ;
  if (this->d_influpos != NULL)
    delete this->d_influpos;
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
}

int sutra_dm::pzt_loadarrays(float *influ, int *influpos, int *npoints,
    int *istart, int *xoff, int *yoff) {
  this->d_influ->host2device(influ);
  this->d_xoff->host2device(xoff);
  this->d_yoff->host2device(yoff);
  this->d_istart->host2device(istart);
  this->d_influpos->host2device(influpos);
  this->d_npoints->host2device(npoints);

  return EXIT_SUCCESS;
}

int sutra_dm::kl_loadarrays(float *rabas, float *azbas, int *ord, float *cr,
    float *cp) {
  this->d_kl->d_rabas->host2device(rabas);
  this->d_kl->d_azbas->host2device(azbas);
  this->d_kl->h_ord->fill_from(ord);
  this->d_kl->d_ord->host2device(ord);
  this->d_kl->d_cr->host2device(cr);
  this->d_kl->d_cp->host2device(cp);

  return EXIT_SUCCESS;
}

int sutra_dm::reset_shape() {
  current_context->set_activeDeviceForce(device);

  cutilSafeCall(
      cudaMemset(this->d_shape->d_screen->getData(), 0,
          sizeof(float) * this->d_shape->d_screen->getNbElem()));

  return EXIT_SUCCESS;
}

int sutra_dm::comp_shape(float *comvec) {
  current_context->set_activeDeviceForce(device);
  this->reset_shape();

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(current_context->get_device(device), this->d_shape->d_screen->getNbElem(),
      nblocks, nthreads);

  if (this->type == "pzt")
    comp_dmshape(nthreads, nblocks, this->d_influ->getData(),
        this->d_shape->d_screen->getData(), this->d_influpos->getData(),
        this->d_istart->getData(), this->d_npoints->getData(), comvec,
        this->influsize * this->influsize,
        this->d_shape->d_screen->getNbElem());

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

int sutra_dm::comp_shape() {
  return this->comp_shape(this->d_comm->getData());
}

int sutra_dm::comp_oneactu(int nactu, float ampli) {
  this->reset_shape();
  int nthreads = 0, nblocks = 0;
  //getNumBlocksAndThreads(this->device,this->dim * this->dim, nblocks, nthreads);
  getNumBlocksAndThreads(current_context->get_device(device), this->influsize * this->influsize,
      nblocks, nthreads);
  if (this->type == "pzt")
    oneactu(nthreads, nblocks, this->d_influ->getData(),
        this->d_shape->d_screen->getData(), nactu, ampli,
        this->d_xoff->getData(), this->d_yoff->getData(), this->dim,
        this->influsize, this->d_shape->d_screen->getNbElem());
  if (this->type == "tt")
    oneactu(nthreads, nblocks, this->d_influ->getData(),
        this->d_shape->d_screen->getData(), nactu, ampli, this->dim,
        this->influsize, this->d_shape->d_screen->getNbElem());
  if (this->type == "kl") {
    int xoff = (int) ((this->d_shape->d_screen->getDims()[1] - this->d_kl->dim)
        / 2.0f);
    int yoff = xoff;
    this->d_kl->do_compute(ampli, this->d_shape->d_screen->getData(), nactu,
        this->d_shape->d_screen->getDims()[1], xoff, yoff);
  }

  return EXIT_SUCCESS;
}

int
sutra_dm::get_IF(float *IF, int *indx_pup, long nb_pts, float ampli){

	for (int i=0 ; i<this->ninflu ; i++){
		this->comp_oneactu(i,ampli);
		getIF(IF,this->d_shape->d_screen->getData(),indx_pup,nb_pts,i,this->ninflu,current_context->get_device(device));
	}

	this->reset_shape();

	return EXIT_SUCCESS;
}

int
sutra_dm::get_IF_sparse(carma_sparse_obj<float> *&d_IFsparse, int *indx_pup, long nb_pts, float ampli){

	int nnz_tot = 0;
	float *values[this->ninflu];
	int *colind[this->ninflu];
	int NZ[this->ninflu];
	long dims_data2[3] = {2,1,nb_pts};
	carma_obj<float> d_IF(current_context,dims_data2);
	carma_sparse_obj<float> *d_IFsparse_vec;

	cout << "Computing IF sparse..." << endl;
	for(int i=0 ; i<this->ninflu ; i++){
		//Compute and store IF for actu i in d_IF
		this->comp_oneactu(i,ampli);
		getIF(d_IF.getData(),this->d_shape->d_screen->getData(),indx_pup,nb_pts,0,this->ninflu,this->current_context->get_device(device));
		//CUsparse d_IF
		d_IFsparse_vec = new carma_sparse_obj<float>(&d_IF);
		// Retrieve nnz, values and colind from d_IFsparse_vec, stored on CPU
		//DEBUG_TRACE("nnz : %d \n",d_IFsparse_vec->getNzElem());
		NZ[i] = d_IFsparse_vec->getNzElem();
		values[i] = (float*)malloc(NZ[i]*sizeof(float));
		colind[i] = (int*)malloc(NZ[i]*sizeof(int));

		cutilSafeCall(
		      cudaMemcpy(values[i], d_IFsparse_vec->getData(), sizeof(float) * NZ[i],
		          cudaMemcpyDeviceToHost));
		cutilSafeCall(
			  cudaMemcpy(colind[i], d_IFsparse_vec->d_colind, sizeof(int) * NZ[i],
				  cudaMemcpyDeviceToHost));

		nnz_tot += NZ[i];

		delete d_IFsparse_vec;
	}
	//Reconstruction of d_data, d_colind, d_rowind for IFsparse
	long dims_data[2] = {1,nnz_tot};
	carma_obj<float> d_val(current_context,dims_data);
	carma_obj<int> d_col(current_context,dims_data);
	dims_data[1] = this->ninflu + 1;
	carma_obj<int> d_row(current_context,dims_data);
	int cpt[this->ninflu + 1];
	cpt[0] = 0;

	for(int i=0 ; i<this->ninflu ; i++){
		cutilSafeCall(
			  cudaMemcpyAsync(d_val.getData(cpt[i]), values[i], sizeof(float) * NZ[i],
				  cudaMemcpyHostToDevice));
		cutilSafeCall(
			  cudaMemcpyAsync(d_col.getData(cpt[i]), colind[i], sizeof(int) * NZ[i],
				  cudaMemcpyHostToDevice));
		cpt[i+1] = cpt[i] + NZ[i];
	}
	cutilSafeCall(
		  cudaMemcpy(d_row.getData(), cpt, sizeof(int) * (this->ninflu + 1),
			  cudaMemcpyHostToDevice));
	dims_data2[1] = this->ninflu;
	d_IFsparse = new carma_sparse_obj<float>(current_context, dims_data2,
			d_val.getData(), d_col.getData(), d_row.getData(),nnz_tot,false);

	return EXIT_SUCCESS;
}

int
sutra_dm::compute_KLbasis(float *xpos, float *ypos, int *indx, long dim, float norm, float ampli){

	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = this->ninflu;
	dims_data[2] = this->ninflu;
	carma_obj<float> *d_statcov = new carma_obj<float>(current_context, dims_data);
	this->d_KLbasis = new carma_obj<float>(current_context, dims_data);
	long dims_data2[2];
	dims_data2[0] = 1;
	dims_data2[1] = dim;
	carma_obj<int> *d_indx = new carma_obj<int>(this->current_context, dims_data2);

	//Compute the statistic matrix from actuators positions & Kolmogorov statistic
	carma_obj<float> *d_xpos = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_ypos = new carma_obj<float>(current_context, dims_data);

	d_xpos->host2device(xpos);
	d_ypos->host2device(ypos);
	dm_dostatmat(d_statcov->getData(), this->ninflu,d_xpos->getData(), d_ypos->getData(), norm, current_context->get_device(device));
	// Compute and apply piston filter

	this->piston_filt(d_statcov);

	dims_data[1] = this->ninflu;
	dims_data[2] = this->ninflu;
	carma_obj<float> *d_geocov = new carma_obj<float>(current_context, dims_data);
	d_indx->host2device(indx);

	//Sparse version for geomat
	carma_sparse_obj<float> *d_IFsparse;
	this->get_IF_sparse(d_IFsparse,d_indx->getData(),dim,1.0f);
	this->do_geomatFromSparse(d_geocov->getData(),d_IFsparse,ampli);

	delete d_IFsparse;

	//Dense version for geomat
	/*
	dims_data[1] = dim;
	dims_data[2] = this->ninflu;
	carma_obj<float> *d_IF = new carma_obj<float>(this->current_context, dims_data);
	// Get influence functions of the DM
	this->get_IF(d_IF->getData(),d_indx->getData(),dim,1.0f);
	// Compute geometric matrix (to be CUsparsed)
	this->do_geomat(d_geocov->getData(),d_IF->getData(),dim,ampli);
*/
	// Double diagonalisation to obtain KL basis on actuators
	this->DDiago(d_statcov,d_geocov);

	delete d_geocov;
	delete d_indx;
	delete d_xpos;
	delete d_ypos;

	return EXIT_SUCCESS;
}

int
sutra_dm::do_geomatFromSparse(float *d_geocov, carma_sparse_obj<float> *d_IFsparse, float ampli){

	carma_sparse_obj<float> *d_tmp = new carma_sparse_obj<float>(this->current_context);

	carma_gemm<float>(cusparse_handle(),'n','t',d_IFsparse,d_IFsparse,d_tmp);
	carma_csr2dense<float>(d_tmp, d_geocov);
	multi_vect(d_geocov,ampli,this->ninflu*this->ninflu,this->current_context->get_device(device));

	delete d_tmp;
	return EXIT_SUCCESS;
}
int
sutra_dm::do_geomat(float *d_geocov, float *d_IF, long n_pts, float ampli){
	carma_gemm(cublas_handle(), 't', 'n', this->ninflu, this->ninflu, n_pts, 1.0f,
	      d_IF, n_pts, d_IF, n_pts, 0.0f,
	      d_geocov, this->ninflu);
	multi_vect(d_geocov,ampli,this->ninflu*this->ninflu,this->current_context->get_device(device));

	return EXIT_SUCCESS;
}

int
sutra_dm::piston_filt(carma_obj <float> *d_statcov){
	long Nmod = d_statcov->getDims()[1];
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = Nmod;
	dims_data[2] = Nmod;
	carma_obj<float> *d_F = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);

	int N = d_statcov->getDims()[1] * d_statcov->getDims()[1];
	fill_filtermat(d_F->getData(),Nmod, N, current_context->get_device(device));

	carma_gemm(cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
			  d_F->getData(), Nmod, d_statcov->getData(), Nmod, 0.0f,
			  d_tmp->getData(), Nmod);
	carma_gemm(cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
			      d_tmp->getData(), Nmod, d_F->getData(), Nmod, 0.0f,
			      d_statcov->getData(), Nmod);

	delete d_tmp;
	delete d_F;

	return EXIT_SUCCESS;
}

int
sutra_dm::set_comkl(float *comvec){
	if(this->d_KLbasis == NULL){
		cerr << "KLbasis has to be computed before calling this function" << endl;
		return EXIT_FAILURE;
	}
	else{
		long dims_data[2] = {1,this->ninflu};
		carma_obj<float> d_comkl(current_context,dims_data,comvec);

		carma_gemv(cublas_handle(),'n', ninflu, ninflu, 1.0f, this->d_KLbasis->getData(),this->d_KLbasis->getDims()[1],
					d_comkl.getData(),1, 0.0f, this->d_comm->getData(),1);

		return EXIT_SUCCESS;
	}
}

int
sutra_dm::DDiago(carma_obj<float> *d_statcov, carma_obj<float> *d_geocov){
	const long dims_data[3] = {2, this->ninflu, this->ninflu};
	carma_obj<float> *d_M1 = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context,dims_data);

	const long dims_data2[2]={1,this->ninflu};
	carma_obj<float> *d_eigenvals = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_sqrt = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,dims_data2);
	carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_sqrt = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

	// 1. SVdec(geocov,U) --> Ut * geocov * U = D²
	carma_syevd<float,1>('V', d_geocov, h_eigenvals);

	d_eigenvals->host2device(*h_eigenvals);
	for (int i=0 ; i<this->ninflu ; i++){
		h_eigenvals_sqrt->getData()[i] = sqrt(h_eigenvals->getData()[i]); // D = sqrt(D²)
		h_eigenvals_inv->getData()[i] = 1./sqrt(h_eigenvals->getData()[i]);// D⁻¹ = 1/sqrt(D²)
	}
	d_eigenvals_sqrt->host2device(*h_eigenvals_sqrt);
	d_eigenvals_inv->host2device(*h_eigenvals_inv);

	// 2. M⁻¹ = sqrt(eigenvals) * Ut : here, we have transpose(M⁻¹)
	/*
	carma_dgmm<float>(cublas_handle(),CUBLAS_SIDE_RIGHT,this->ninflu,this->ninflu,
					d_geocov->getData(), this->ninflu, d_eigenvals_inv->getData(),1,
					d_M1->getData(), this->ninflu);*/

	carma_dgmm<float>(cublas_handle(),CUBLAS_SIDE_RIGHT,this->ninflu,this->ninflu,
					d_geocov->getData(), this->ninflu, d_eigenvals_sqrt->getData(),1,
					d_M1->getData(), this->ninflu);

	// 3. C' = M⁻¹ * statcov * M⁻¹t
	carma_gemm<float>(cublas_handle(), 't', 'n', ninflu, ninflu, ninflu, 1.0f,
		      d_M1->getData(), ninflu, d_statcov->getData(), ninflu, 0.0f,
		      d_tmp->getData(), ninflu);

	carma_gemm<float>(cublas_handle(), 'n', 'n', ninflu, ninflu, ninflu, 1.0f,
			      d_tmp->getData(), ninflu, d_M1->getData(), ninflu, 0.0f,
			      d_tmp2->getData(), ninflu);

	// 4. SVdec(C',A)
	carma_syevd<float,1>('V', d_tmp2, h_eigenvals);

	// 5. M = U * D⁻¹
	carma_dgmm<float>(cublas_handle(),CUBLAS_SIDE_RIGHT,this->ninflu,this->ninflu,
				d_geocov->getData(), this->ninflu, d_eigenvals_inv->getData(),1,
				d_tmp->getData(), this->ninflu);

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
  this->ndm = ndm;
}

sutra_dms::~sutra_dms() {
  for (std::map<type_screen, sutra_dm *>::iterator it = this->d_dms.begin();
      this->d_dms.end() != it; it++) {
    delete it->second;
  }
  this->d_dms.clear();
}

int sutra_dms::add_dm(carma_context *context, const char* type, float alt,
    long dim, long ninflu, long influsize, long ninflupos, long n_npoints,
    float push4imat, int device) {
  this->d_dms[make_pair(type, alt)] = new sutra_dm(context, type, dim, ninflu,
      influsize, ninflupos, n_npoints, push4imat, device);

  return EXIT_SUCCESS;
}

int sutra_dms::remove_dm(const char* type, float alt) {
  delete this->d_dms[make_pair(type, alt)];
  this->d_dms.erase(make_pair(type, alt));

  return EXIT_SUCCESS;
}

// Florian features
int sutra_dm::kl_floloadarrays(float *covmat, float *filter, float *evals,
    float *bas) {
  this->d_kl->d_covmat->host2device(covmat);
  this->d_kl->d_filter->host2device(filter);
  this->d_kl->d_bas->host2device(bas);
  this->d_kl->d_evals->host2device(evals);

  return EXIT_SUCCESS;
}

int sutra_dms::nact_total() {
  map<type_screen, sutra_dm *>::iterator p;
  p = d_dms.begin();
  int idx = 0;
  while (p != d_dms.end()) {
    idx += p->second->ninflu;
    p++;
  }
  return idx;
}
