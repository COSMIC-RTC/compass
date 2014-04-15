#include <sutra_controller_mv.h>
#include <string>

sutra_controller_mv::sutra_controller_mv(carma_context *context, long nvalid,
    long nactu, long delay) :
    sutra_controller(context, nvalid * 2, nactu) {

  this->nvalid = nvalid;
  this->nslope = 2 * nvalid;
  this->nactu = nactu;
  this->delay = delay;
  this->device = device;
  this->gain = 0.0f;

  this->nstreams = 1; //nvalid/10;
  while (nactu % this->nstreams != 0)
    nstreams--;
  cerr << "controller uses " << nstreams << " streams" << endl;
  streams = new carma_streams(nstreams);
  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = nslope;
  dims_data2[2] = nactu;
  this->d_imat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = nactu;
  dims_data2[2] = nslope;
  this->d_cmat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = dims_data2[2] = nactu;
  d_U = new carma_obj<float>(this->current_context, dims_data2);
  this->d_cenbuff = 0L;
  if (delay > 0) {
    dims_data2[1] = 2 * nvalid;
    dims_data2[2] = delay + 1;
    this->d_cenbuff = new carma_obj<float>(this->current_context, dims_data2);
  }

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nslope < nactu ? nslope : nactu;
  this->d_eigenvals = new carma_obj<float>(this->current_context, dims_data1);
  this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<float>(this->current_context, dims_data1);
  this->d_noisemat = new carma_obj<float>(this->current_context, dims_data1);
  dims_data1[1] = nactu;
  this->d_com = new carma_obj<float>(this->current_context, dims_data1);
  this->d_com1 = new carma_obj<float>(this->current_context, dims_data1);
  this->d_com2 = new carma_obj<float>(this->current_context, dims_data1);
  this->d_err = new carma_obj<float>(this->current_context, dims_data1);
  this->d_gain = new carma_obj<float>(this->current_context, dims_data1);

  // Florian features
  dims_data2[1] = nactu;//64564;
  dims_data2[2] = nactu;
  this->d_covmat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = nactu;
  dims_data2[2] = nactu;
  this->d_KLbasis = new carma_obj<float>(this->current_context, dims_data2);

  cublas_handle = current_context->get_cublasHandle();
  //carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

sutra_controller_mv::~sutra_controller_mv() {
  current_context->set_activeDevice(device);
  delete this->d_U;

  delete this->streams;
  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_gain;

  delete this->d_eigenvals;
  delete this->h_eigenvals;

  if (this->delay > 0)
    delete this->d_cenbuff;
  delete this->d_com1;
  delete this->d_com2;
  delete this->d_err;
// Florian features
  delete this->d_covmat;
  delete this->d_KLbasis;

  //carma_checkCublasStatus(cublasDestroy(this->cublas_handle));

  //delete this->current_context;
}

string sutra_controller_mv::get_type() {
  return "mv";
}

int sutra_controller_mv::set_gain(float gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

int sutra_controller_mv::load_mgain(float *mgain) {
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

int sutra_controller_mv::load_noisemat(float *noise) {
  this->d_noisemat->host2device(noise);
  return EXIT_SUCCESS;
}
// Florian features
int
sutra_controller_mv::do_covmat(sutra_dm *ydm, char *method, int *indx_pup,long dim, float *xpos, float *ypos, float norm){
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = dim;
	dims_data[2] = this->nactu;
	carma_obj<float> *d_IF = new carma_obj<float>(this->current_context, dims_data);
	dims_data[1] = this->nactu;
	carma_obj<float> *d_geocov = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_statcov = new carma_obj<float>(current_context, dims_data);

	long dims_data2[2];
	dims_data2[0] = 1;
	dims_data2[1] = this->nactu;
	carma_obj<float> *d_KLcov = new carma_obj<float>(current_context, dims_data2);

	dims_data2[1] = dim;
	carma_obj<int> *d_indx = new carma_obj<int>(this->current_context, dims_data2);

	// Get influence functions of the DM
	d_indx->host2device(indx_pup);
	ydm->get_IF(d_IF->getData(),d_indx->getData(),dim);

	// Compute geometric matrix (to be CUsparsed)
	this->do_geomat(d_geocov,d_IF,dim);

	delete d_IF;

	//Compute the statistic matrix from actuators positions & Kolomogorov statistic
	carma_obj<float> *d_xpos = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_ypos = new carma_obj<float>(current_context, dims_data);
	d_xpos->host2device(xpos);
	d_ypos->host2device(ypos);

	do_statmat(d_statcov->getData(), this->nactu,d_xpos->getData(), d_ypos->getData(), norm, device);
	// Compute and apply piston filter
	this->piston_filt(d_statcov);

	// Double diagonalisation to obtain KL basis on actuators
	this->DDiago(d_statcov,d_geocov);

	// Computation of covariance matrix
	// 1. Computation of covariance matrix in KL basis
	carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

	carma_syevd<float,1>('N', d_statcov, h_eigenvals);

	if(strcmp(method, "inv") == 0){
		// Inversion of the KL covariance matrix
		carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
		for (int i=0 ; i<this->nactu ; i++){
				// Valeurs propores négatives.... A voir et inverser ordre si valeurs propres positives
				h_eigenvals_inv->getData()[i] = -1./h_eigenvals->getData()[i];
			}

		h_eigenvals_inv->getData()[this->nactu - 1] = 0.;
		d_KLcov->host2device(*h_eigenvals_inv);

		// 2. Inversion of the KL basis
		dims_data2[0] = 1;
		dims_data2[1] = this->nactu;
		carma_obj<float> *d_eigen = new carma_obj<float>(current_context, dims_data2);
		dims_data[1] = this->nactu;
		dims_data[2] = this->nactu;
		carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
		carma_obj<float> *d_Ukl = new carma_obj<float>(current_context, dims_data);
		carma_obj<float> *d_Vkl = new carma_obj<float>(current_context, dims_data);
		carma_host_obj<float> *h_KL = new carma_host_obj<float>(dims_data, MA_PAGELOCK);
		carma_host_obj<float> *h_U = new carma_host_obj<float>(dims_data, MA_PAGELOCK);
		carma_host_obj<float> *h_Vt = new carma_host_obj<float>(dims_data, MA_PAGELOCK);

		h_KL->cpy_obj(this->d_KLbasis,cudaMemcpyDeviceToHost);

		carma_svd_cpu<float>(h_KL, h_eigenvals, h_U, h_Vt);

		d_Ukl->host2device(*h_Vt);
		d_Vkl->host2device(*h_U);

		for (int i=0 ; i<this->nactu ; i++){
				h_eigenvals_inv->getData()[i] = 1./h_eigenvals->getData()[i];
			}

		d_eigen->host2device(*h_eigenvals_inv);

		carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nactu,nactu,d_Vkl->getData(),nactu,d_eigen->getData(),1,d_tmp->getData(),nactu);
		carma_gemm(cublas_handle, 't', 't', nactu, nactu, nactu, 1.0f,
								  d_tmp->getData(), nactu, d_Ukl->getData(), nactu, 0.0f,
								  this->d_KLbasis->getData(), nactu);

		// 3. Projection of KL covariance matrix in the DM basis

		carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nactu,nactu,d_KLbasis->getData(),nactu,d_KLcov->getData(),1,d_tmp->getData(),nactu);
		carma_gemm(cublas_handle, 't', 'n', nactu, nactu, nactu, 1.0f,
					  d_tmp->getData(), nactu, this->d_KLbasis->getData(), nactu, 0.0f,
					  this->d_covmat->getData(), nactu);

		delete d_eigen;
		delete d_tmp;
		delete d_Ukl;
		delete d_Vkl;
		delete h_KL;
		delete h_U;
		delete h_Vt;
		delete h_eigenvals_inv;
	}

	else if(strcmp(method, "n") == 0){
		dims_data[1] = this->nactu;
		dims_data[2] = this->nactu;
		carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
		for (int i=0 ; i<this->nactu ; i++){
						// Valeurs propres négatives.... A voir et inverser ordre si valeurs propres positives
						h_eigenvals->getData()[i] = -h_eigenvals->getData()[i];
					}
		d_KLcov->host2device(*h_eigenvals);

		carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, nactu,nactu,d_KLbasis->getData(),nactu,d_KLcov->getData(),1,d_tmp->getData(),nactu);
		carma_gemm(cublas_handle, 'n', 't', nactu, nactu, nactu, 1.0f,
							  d_tmp->getData(), nactu, this->d_KLbasis->getData(), nactu, 0.0f,
							  this->d_covmat->getData(), nactu);

		delete d_tmp;
	}
	else {}// y_error("Specifiy the computation method for mv : inv or n \n");}

	delete d_geocov;
	delete d_statcov;
	delete h_eigenvals;
	delete d_KLcov;
	delete d_indx;
	delete d_xpos;
	delete d_ypos;

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::do_geomat(carma_obj<float> *d_geocov, carma_obj<float> *d_IF, long n_pts){
	carma_gemm(cublas_handle, 't', 'n', nactu, nactu, n_pts, 1.0f,
	      d_IF->getData(), n_pts, d_IF->getData(), n_pts, 0.0f,
	      d_geocov->getData(), nactu);

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::piston_filt(carma_obj<float> *d_statcov){
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = this->nactu;
	dims_data[2] = this->nactu;
	carma_obj<float> *d_F = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);

	int N = this->nactu * this->nactu;
	fill_filtmat(d_F->getData(),this->nactu, N, device);

	carma_gemm(cublas_handle, 'n', 'n', nactu, nactu, nactu, 1.0f,
		      d_F->getData(), nactu, d_statcov->getData(), nactu, 0.0f,
		      d_tmp->getData(), nactu);
	carma_gemm(cublas_handle, 'n', 't', nactu, nactu, nactu, 1.0f,
			      d_tmp->getData(), nactu, d_F->getData(), nactu, 0.0f,
			      d_statcov->getData(), nactu);

	delete d_tmp;
	delete d_F;

	return EXIT_SUCCESS;
}
/*
int sutra_controller_mv::do_statmat(float *statcov, float *xpos, float *ypos){
	int dim_x = sizeof(xpos)/sizeof(xpos[0]);
	int ind;
	for (i=0 ; i<dim_x ; i++){
		for(j=0 ; j<dim_x ; j++){
			ind = i*dim_x + j;
			statcov[ind] = 6.88 * pow(sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2)),5./3);
		}
	}

	return EXIT_SUCCESS;
}
*/
int sutra_controller_mv::DDiago(carma_obj<float> *d_statcov, carma_obj<float> *d_geocov){
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = this->nactu;
	dims_data[2] = this->nactu;
	carma_obj<float> *d_M1 = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context,dims_data);
	long dims_data2[2];
	dims_data2[0] = 1;
	dims_data2[1] = this->nactu;
	carma_obj<float> *d_eigenvals = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_sqrt = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,dims_data2);
	carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_sqrt = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

	// 1. SVdec(geocov,U) --> Ut * geocov * U = D²
	carma_syevd<float,1>('V', d_geocov, h_eigenvals);

	d_eigenvals->host2device(*h_eigenvals);
	for (int i=0 ; i<this->nactu ; i++){
		h_eigenvals_sqrt->getData()[i] = sqrt(h_eigenvals->getData()[i]); // D = sqrt(D²)
		h_eigenvals_inv->getData()[i] = 1./sqrt(h_eigenvals->getData()[i]);// D⁻¹ = 1/sqrt(D²)
	}
	d_eigenvals_sqrt->host2device(*h_eigenvals_sqrt);
	d_eigenvals_inv->host2device(*h_eigenvals_inv);

	// 2. M⁻¹ = sqrt(eigenvals) * Ut : here, we have transpose(M⁻¹)

	carma_dgmm<float>(cublas_handle,CUBLAS_SIDE_RIGHT,this->nactu,this->nactu,
					d_geocov->getData(), this->nactu, d_eigenvals_sqrt->getData(),1,
					d_M1->getData(), this->nactu);

	// 3. C' = M⁻¹ * statcov * M⁻¹t
	carma_gemm<float>(cublas_handle, 't', 'n', nactu, nactu, nactu, 1.0f,
		      d_M1->getData(), nactu, d_statcov->getData(), nactu, 0.0f,
		      d_tmp->getData(), nactu);
	carma_gemm<float>(cublas_handle, 'n', 'n', nactu, nactu, nactu, 1.0f,
			      d_tmp->getData(), nactu, d_M1->getData(), nactu, 0.0f,
			      d_tmp2->getData(), nactu);

	DEBUG_TRACE("here\n");
	// 4. SVdec(C',A)
	carma_syevd<float,1>('V', d_tmp2, h_eigenvals);

	// 5. M = U * D⁻¹
	carma_dgmm<float>(cublas_handle,CUBLAS_SIDE_RIGHT,this->nactu,this->nactu,
				d_geocov->getData(), this->nactu, d_eigenvals_inv->getData(),1,
				d_tmp->getData(), this->nactu);

	// 6. B = M * A;
	carma_gemm<float>(cublas_handle, 'n', 'n', nactu, nactu, nactu, 1.0f,
			      d_tmp->getData(), nactu, d_tmp2->getData(), nactu, 0.0f,
			      d_KLbasis->getData(), nactu);

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

int sutra_controller_mv::load_covmat(float *covmat) {
  this->d_covmat->host2device(covmat);
  return EXIT_SUCCESS;
}
int sutra_controller_mv::load_klbasis(float *klbasis) {
  this->d_KLbasis->host2device(klbasis);
  return EXIT_SUCCESS;
}
int sutra_controller_mv::set_delay(int delay) {
  this->delay = delay;
  return EXIT_SUCCESS;
}

// Florian features
int sutra_controller_mv::build_cmat(const char *dmtype, char *method) {
  float one = 1.;
  float zero = 0.;

  if(strcmp(method,"inv") == 0){
	  carma_obj<float> *d_tmp;
	  carma_obj<float> *d_tmp2;
	  carma_obj<float> *d_tmp3;
	  long *dims_data2 = new long[3];
	  dims_data2[0] = 2;
	  dims_data2[1] = nactu;
	  dims_data2[2] = nactu;
	  d_tmp = new carma_obj<float>(current_context, dims_data2);
	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
	  carma_obj<float> *d_U = new carma_obj<float>(current_context, dims_data2);
	  carma_obj<float> *d_Vt = new carma_obj<float>(current_context, dims_data2);
	  dims_data2[1] = nslope;
	  dims_data2[2] = nactu;
	  d_tmp3 = new carma_obj<float>(current_context,dims_data2);

	  carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nslope,nactu,d_imat->getData(),nslope,d_noisemat->getData(),1,d_tmp3->getData(),nslope);
	  carma_gemm(cublas_handle, 't', 'n', nactu, nactu, nslope, one,
			d_tmp3->getData(), nslope, d_imat->getData(), nslope, zero,
			d_tmp2->getData(), nactu);
	  carma_geam(cublas_handle, 'n', 'n', nactu, nactu, one, d_tmp2->getData(),
		  nactu, one, d_covmat->getData(), nactu, d_tmp->getData(), nactu);

	  long *dims_data = new long[3];
	  dims_data[0] = 1;
	  dims_data[1] = nactu;
	  dims_data2[1] = nactu;
	  dims_data2[2] = nactu;
	  carma_host_obj<float> *h_tmp = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	  carma_host_obj<float> *h_U = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	  carma_host_obj<float> *h_Vt = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	  carma_host_obj<float> *h_eigen = new carma_host_obj<float>(dims_data, MA_PAGELOCK);
	  carma_obj<float> *d_eigen = new carma_obj<float>(current_context, dims_data);
	  h_tmp->cpy_obj(d_tmp,cudaMemcpyDeviceToHost);

	  carma_svd_cpu<float>(h_tmp, h_eigen, h_U, h_Vt);

	  d_U->host2device(*h_Vt);
	  d_Vt->host2device(*h_U);

	  for (int i=0 ; i<this->nactu ; i++){
		h_eigen->getData()[i] = 1./h_eigen->getData()[i];
	  }
	  d_eigen->host2device(*h_eigen);

		carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nactu,nactu,d_Vt->getData(),nactu,d_eigen->getData(),1,d_tmp2->getData(),nactu);
		carma_gemm(cublas_handle, 't', 't', nactu, nactu, nactu, 1.0f,
								  d_tmp2->getData(), nactu, d_U->getData(), nactu, 0.0f,
								  d_tmp->getData(), nactu);

	  dims_data2[1] = nactu;
	  dims_data2[2] = nslope;
	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
	  carma_gemm(cublas_handle, 'n', 't', nactu, nslope, nactu, one,
		  d_tmp->getData(), nactu, d_imat->getData(), nslope, zero,
		  d_tmp2->getData(), nactu);


	  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, nactu,nslope,d_tmp2->getData(),nactu,d_noisemat->getData(),1,d_cmat->getData(),nactu);

	  delete d_tmp;
	  delete d_tmp2;
	  delete d_tmp3;
	  delete h_tmp;
	  delete h_U;
	  delete h_Vt;
	  delete h_eigen;
	  delete d_eigen;
  }

  else if(strcmp(method,"n") == 0){
  	  carma_obj<float> *d_tmp;
  	  carma_obj<float> *d_tmp2;
  	  carma_obj<float> *d_tmp3;
  	  carma_obj<float> *d_tmp4;
  	  long *dims_data2 = new long[3];
  	  dims_data2[0] = 2;
  	  dims_data2[1] = nslope;
  	  dims_data2[2] = nslope;
  	  carma_obj<float> *d_U = new carma_obj<float>(current_context, dims_data2);
  	  carma_obj<float> *d_Vt = new carma_obj<float>(current_context, dims_data2);
  	  d_tmp = new carma_obj<float>(current_context, dims_data2);
  	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
  	  d_tmp4 = new carma_obj<float>(current_context, dims_data2);
  	  dims_data2[1] = nslope;
  	  dims_data2[2] = nactu;
  	  d_tmp3 = new carma_obj<float>(current_context,dims_data2);

  	  carma_gemm(cublas_handle, 'n', 'n', nslope, nactu, nactu, one,
  			d_imat->getData(), nslope, d_covmat->getData(), nactu, zero,
  			d_tmp3->getData(), nslope);
 	  carma_gemm(cublas_handle, 'n', 't', nslope, nslope, nactu, one,
  			d_tmp3->getData(), nslope, d_imat->getData(), nslope, zero,
  			d_tmp2->getData(), nslope);
 	  add_md(d_tmp4->getData(),d_tmp4->getData(),d_noisemat->getData(), nslope, this->device);
  	  carma_geam(cublas_handle, 'n', 'n', nslope, nslope, one, d_tmp2->getData(),
  		  nslope, one, d_tmp4->getData(), nslope, d_tmp->getData(), nslope);

  	  long *dims_data = new long[3];
  	  dims_data[0] = 1;
  	  dims_data[1] = nslope;
  	  dims_data2[1] = nslope;
  	  dims_data2[2] = nslope;
  	  carma_host_obj<float> *h_tmp = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
  	  carma_host_obj<float> *h_U = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
  	  carma_host_obj<float> *h_Vt = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
  	  carma_host_obj<float> *h_eigen = new carma_host_obj<float>(dims_data, MA_PAGELOCK);
  	  carma_obj<float> *d_eigen = new carma_obj<float>(current_context, dims_data);
  	  h_tmp->cpy_obj(d_tmp,cudaMemcpyDeviceToHost);

  	  carma_svd_cpu<float>(h_tmp, h_eigen, h_U, h_Vt);

  	  d_U->host2device(*h_Vt);
  	  d_Vt->host2device(*h_U);

  	  for (int i=0 ; i<this->nactu ; i++){
  		h_eigen->getData()[i] = 1./h_eigen->getData()[i];
  	  }
  	  d_eigen->host2device(*h_eigen);

  		carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nslope,nslope,d_Vt->getData(),nslope,d_eigen->getData(),1,d_tmp2->getData(),nslope);
  		carma_gemm(cublas_handle, 't', 't', nslope, nslope, nslope, 1.0f,
  								  d_tmp2->getData(), nslope, d_U->getData(), nslope, 0.0f,
  								  d_tmp->getData(), nslope);

  	  dims_data2[1] = nactu;
  	  dims_data2[2] = nslope;
  	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
  	  carma_gemm(cublas_handle, 't', 'n', nactu, nslope, nslope, one,
  		  d_imat->getData(), nslope, d_tmp->getData(), nslope, zero,
  		  d_tmp2->getData(), nactu);
  	  carma_gemm(cublas_handle, 'n', 'n', nactu, nslope, nactu, one,
  		  d_covmat->getData(), nactu, d_tmp2->getData(), nactu, zero,
  		  d_cmat->getData(), nactu);


  	  delete d_tmp;
  	  delete d_tmp2;
  	  delete d_tmp3;
  	  delete d_tmp4;
  	  delete h_tmp;
  	  delete h_U;
  	  delete h_Vt;
  	  delete h_eigen;
  	  delete d_eigen;
    }
  else { }//y_error("Specify the computation method for mv : inv or n \n");}
  return EXIT_SUCCESS;
}

int sutra_controller_mv::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if (delay > 0) {
    for (int cc = 0; cc < delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope]), 1,
          this->nslope, this->device);

    cutilSafeCall(
        cudaMemcpy(&(this->d_cenbuff->getData()[delay * this->nslope]),
            this->d_centroids->getData(), sizeof(float) * this->nslope,
            cudaMemcpyDeviceToDevice));

    cutilSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
            sizeof(float) * this->nslope, cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

int sutra_controller_mv::comp_com() {

  long dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = nslope;
  carma_obj<float> d_tmp(current_context, dims_data2);
  carma_obj<float> d_tmp2(current_context, dims_data2);
  carma_obj<float> d_olmeas(current_context, dims_data2);
  dims_data2[1] = nactu;
  carma_obj<float> d_tmp3(current_context, dims_data2);
  carma_obj<float> d_tmp4(current_context, dims_data2);
  this->d_com2->copy(this->d_com1, 1, 1);
  this->d_com1->copy(this->d_com, 1, 1);
  // POLC equations
  carma_geam<float>(cublas_handle, 'n', 'n', nactu, 1, (float) delay, *d_com2,
      nactu, 1.0f - delay, *d_com1, nactu, d_tmp, nactu);
  //carma_gemv<float>(cublas_handle,'n',nslope,nactu,-1.0f,*d_imat,nslope,d_tmp,1,1.0f,*d_centroids,1);
  carma_gemv<float>(cublas_handle, 'n', nslope, nactu, 1.0f, *d_imat, nslope,
      d_tmp, 1, 0.0f, d_tmp2, 1);
  carma_geam<float>(cublas_handle, 'n', 'n', nslope, 1, 1.0f, *d_centroids,
      nslope, -1.0f, d_tmp2, nslope, d_olmeas, nslope);

  if (this->nstreams > 1) {
    int nstreams = this->nstreams;
    float alpha = -1.0f;
    float beta = 0.0f;
    //cout << this->d_cmat->getDims(1) << " x " << this->d_cmat->getDims(2) << endl;
    //cout << this->d_centroids->getNbElem()<< endl;
    //cout << this->d_err->getNbElem()<< endl;
    for (int i = 0; i < nstreams; i++) {
      int istart1 = i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1)
          / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      //cout << istart1 << " " << istart2 << endl;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      // carma_checkCublasStatus(
      //     cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
      //         this->d_cmat->getDims(2), &alpha,
      //         &((this->d_cmat->getData())[istart1]),
      //         this->d_cmat->getDims(1) / nstreams, this->d_centroids->getData(),
      //         1, &beta, &((this->d_err->getData())[istart2]), 1));
      carma_checkCublasStatus(
          cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
              this->d_cmat->getDims(2), &alpha,
              &((this->d_cmat->getData())[istart1]),
              this->d_cmat->getDims(1) / nstreams, d_olmeas, 1, &beta,
              &((this->d_err->getData())[istart2]), 1));
    }

    //mult_int(this->d_com->getData(), this->d_err->getData(),
    //   this->d_gain->getData(), this->gain, this->nactu, this->device,
    //   this->streams);

    this->streams->wait_all_streams();
    //this->d_com->copy(this->d_err,1,1);
  } else {
    // compute error
    //  this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
    //      this->d_centroids, 1, 0.0f, 1);
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
        &d_olmeas, 1, 0.0f, 1);
    // this->d_com->copy(this->d_err,1,1);
    // apply modal gain & loop gain
    //mult_int(this->d_com->getData(), this->d_err->getData(),
    //    this->d_gain->getData(), this->gain, this->nactu, this->device);
  }

  mult_int(d_tmp4, this->d_err->getData(), this->d_gain->getData(), this->gain,
      this->nactu, this->device);
  mult_int(d_tmp3, this->d_com1->getData(), this->d_gain->getData(),
      1 - (this->gain), this->nactu, this->device);
  carma_geam<float>(cublas_handle, 'n', 'n', nactu, 1, 1.0f, d_tmp4, nactu,
      1.0f, d_tmp3, nactu, *d_com, nactu);

  return EXIT_SUCCESS;
}
