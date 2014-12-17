#include <sutra_controller_mv.h>
#include <sutra_controller_utils.h>
#include <string>

sutra_controller_mv::sutra_controller_mv(carma_context *context, long nvalid,
    long nactu_, long delay) :
    sutra_controller(context, nvalid * 2, nactu_) {


  this->delay = delay;
  this->device = device;
  this->gain = 0.0f;

//  this->nstreams = 1; //nvalid/10;
//  while (nactu() % this->nstreams != 0)
//    nstreams--;
//  cerr << "controller uses " << nstreams << " streams" << endl;
//  streams = new carma_streams(nstreams);
  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = nslope();
  dims_data2[2] = nactu();
  this->d_imat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = nactu();
  dims_data2[2] = nslope();
  this->d_cmat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = nactu() - 2;
  dims_data2[2] = nslope();
  this->d_Cphim = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = dims_data2[2] = nactu();
  d_U = new carma_obj<float>(this->current_context, dims_data2);
  this->d_cenbuff = 0L;
  if (delay > 0) {
    dims_data2[1] = 2 * nvalid;
    dims_data2[2] = delay + 1;
    this->d_cenbuff = new carma_obj<float>(this->current_context, dims_data2);
  }
  dims_data2[1] = nslope();
  dims_data2[2] = nslope();
  this->d_Cmm = new carma_obj<float>(this->current_context, dims_data2);

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nslope() < nactu() ? nslope() : nactu();
  this->d_eigenvals = new carma_obj<float>(this->current_context, dims_data1);
  this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  dims_data1[1] = nslope();
  this->d_noisemat = new carma_obj<float>(this->current_context, dims_data1);
  this->d_olmeas = new carma_obj<float>(this->current_context, dims_data1);
  this->d_compbuff = new carma_obj<float>(this->current_context, dims_data1);
  this->d_compbuff2 = new carma_obj<float>(this->current_context, dims_data1);
  dims_data1[1] = nactu();
  this->d_com1 = new carma_obj<float>(this->current_context, dims_data1);
  this->d_com2 = new carma_obj<float>(this->current_context, dims_data1);
  this->d_err = new carma_obj<float>(this->current_context, dims_data1);
  this->d_gain = new carma_obj<float>(this->current_context, dims_data1);

  // Florian features
  dims_data2[1] = nactu();//64564;
  dims_data2[2] = nactu();
  this->d_covmat = new carma_obj<float>(this->current_context, dims_data2);
  dims_data2[1] = nactu();
  dims_data2[2] = nactu();
  this->d_KLbasis = new carma_obj<float>(this->current_context, dims_data2);

  cublas_handle = current_context->get_cublasHandle();
  //carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

sutra_controller_mv::~sutra_controller_mv() {
  current_context->set_activeDevice(device);
  delete this->d_U;

  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_gain;
  delete this->d_Cmm;
  delete this->d_Cphim;
  delete this->d_olmeas;
  delete this->d_compbuff;
  delete this->d_compbuff2;
  delete this->d_noisemat;

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
sutra_controller_mv::compute_Cmm(sutra_atmos *atmos, sutra_sensors *sensors, float *L0, float *cn2, float *alphaX, float *alphaY, float diamTel, float cobs){

	struct gtomo_struct g_tomo;
	init_tomo_gpu_gb(&g_tomo, atmos, sensors, diamTel, cobs);
	update_tomo_sys_gpu_gb(&g_tomo, sensors, alphaX, alphaY);
	update_tomo_atm_gpu_gb(&g_tomo, sensors, atmos, L0, cn2, alphaX, alphaY);
	matcov_gpu_4(this->d_Cmm->getData(), this->nslope(), this->nslope(), 0, 0, this->nslope(), &g_tomo, atmos, sensors, alphaX, alphaY);
	free_tomo_gpu_gb(&g_tomo);

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::compute_Cphim(sutra_atmos *atmos, sutra_sensors *sensors, sutra_dms *dms, float *L0, float *cn2, float *alphaX, float *alphaY, float *X, float *Y,float *xactu, float *yactu, float diamTel, float k2, float *Nact) {

	// Find number of actuators without TT DM
	struct cphim_struct cphim_struct;
	int Nactu = 0;
	map<type_screen, sutra_dm *>::iterator p;
	p = dms->d_dms.begin();
      while (p != dms->d_dms.end()) {
    	  if (p->second->type != "tt"){
	      	  Nactu += p->second->ninflu;
    	  }
	      p++;
      }

    // Compute Cphim matrix
	init_cphim_struct(&cphim_struct, atmos, sensors, dms, diamTel);
	update_cphim_sys(&cphim_struct, sensors, alphaX, alphaY, xactu, yactu,X,Y);
	tab_u831J0(&cphim_struct,-4.0f,10.0f,this->device);
	update_cphim_atm(&cphim_struct, sensors, atmos, L0, cn2, alphaX, alphaY);
	CPHIM(this->d_Cphim->getData(), Nactu, this->nslope(), 0, 0, Nactu, &cphim_struct, atmos, sensors, alphaX, alphaY, k2, this->device);
	free_cphim_struct(&cphim_struct);

	// Piston filter
	piston_filt_cphim(this->d_Cphim);
	//Init and inverse the coupling matrix
	long dims_data2[3];
	dims_data2[0] = 2;
	dims_data2[1] = nactu() - 2;
	dims_data2[2] = nactu() - 2;
	carma_obj<float> *d_Nact = new carma_obj<float>(current_context, dims_data2);
	dims_data2[2] = nslope();
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data2);
	d_Nact->host2device(Nact);
	carma_potri(d_Nact);
	carma_gemm(cublas_handle, 'n', 'n', nactu()-2, nslope(), nactu()-2, 1.0f,
				  d_Nact->getData(), nactu()-2, this->d_Cphim->getData(), nactu()-2, 0.0f,
				  d_tmp->getData(), nactu()-2);
	this->d_Cphim->copy(d_tmp,1,1);

	delete d_Nact;
	delete d_tmp;


	return EXIT_SUCCESS;
}


int
sutra_controller_mv::do_covmat(sutra_dm *ydm, char *method, int *indx_pup,long dim, float *xpos, float *ypos, long Nkl, float norm, float ampli){
		long dims_data[3];
		dims_data[0] = 2;
		dims_data[1] = Nkl;
		dims_data[2] = Nkl;
		carma_obj<float> *d_statcov = new carma_obj<float>(current_context, dims_data);
		long dims_data2[2];
		dims_data2[0] = 1;
		dims_data2[1] = this->nactu();
		carma_obj<float> *d_KLcov = new carma_obj<float>(current_context, dims_data2);

		dims_data2[1] = dim;
		carma_obj<int> *d_indx = new carma_obj<int>(this->current_context, dims_data2);

		//Compute the statistic matrix from actuators positions & Kolmogorov statistic
		carma_obj<float> *d_xpos = new carma_obj<float>(current_context, dims_data);
		carma_obj<float> *d_ypos = new carma_obj<float>(current_context, dims_data);
		d_xpos->host2device(xpos);
		d_ypos->host2device(ypos);
		do_statmat(d_statcov->getData(), Nkl,d_xpos->getData(), d_ypos->getData(), norm, device);
		// Compute and apply piston filter
		this->piston_filt(d_statcov);

		if (ydm->type == "pzt"){
			this->d_covmat->copy(d_statcov,1,1);
			delete d_statcov;
			delete d_KLcov;
			delete d_indx;
			delete d_xpos;
			delete d_ypos;

			return EXIT_SUCCESS;
			/*
			dims_data[1] = dim;
			dims_data[2] = this->nactu();
			carma_obj<float> *d_IF = new carma_obj<float>(this->current_context, dims_data);
			dims_data[1] = this->nactu();
			dims_data[2] = this->nactu();
			carma_obj<float> *d_geocov = new carma_obj<float>(current_context, dims_data);

			// Get influence functions of the DM (to be CUsparsed)
			d_indx->host2device(indx_pup);
			ydm->get_IF(d_IF->getData(),d_indx->getData(),dim);

			// Compute geometric matrix (to be CUsparsed)
			this->do_geomat(d_geocov,d_IF,dim,ampli);

			delete d_IF;

			// Double diagonalisation to obtain KL basis on actuators
			this->DDiago(d_statcov,d_geocov);

			delete d_geocov;
			*/
		}
		// Computation of covariance matrix
		// 1. Computation of covariance matrix in KL basis
		carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

		carma_syevd<float,1>('N', d_statcov, h_eigenvals);

		if (ydm->type == "kl"){
			dims_data2[1] = this->nactu();
			carma_host_obj<float> *h_KLcov = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
			if(strcmp(method, "inv") == 0){
				for (int i=0 ; i < this->nactu() ; i++){
					h_KLcov->getData()[i] = -1./h_eigenvals->getData()[i];
				}
				if (Nkl == this->nactu()){
					h_KLcov->getData()[this->nactu() - 1] = 0.;
				}
				d_KLcov->host2device(*h_KLcov);
				add_md(this->d_covmat->getData(),this->d_covmat->getData(),d_KLcov->getData(),this->nactu(),this->device);
			}
			if(strcmp(method, "n") == 0){
				for (int i=0 ; i < this->nactu() ; i++){
					h_KLcov->getData()[i] = -h_eigenvals->getData()[i];
				}
				d_KLcov->host2device(*h_KLcov);
				add_md(this->d_covmat->getData(),this->d_covmat->getData(),d_KLcov->getData(),this->nactu(),this->device);
			}
			delete h_KLcov;
		}
		if (ydm->type == "pzt"){
			if(strcmp(method, "inv") == 0){
				// Inversion of the KL covariance matrix
				carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
				for (int i=0 ; i<this->nactu() ; i++){
						// Valeurs propres négatives.... A voir et inverser ordre si valeurs propres positives
						h_eigenvals_inv->getData()[i] = -1./h_eigenvals->getData()[i];
					}

				h_eigenvals_inv->getData()[this->nactu() - 1] = 0.;
				d_KLcov->host2device(*h_eigenvals_inv);

				// 2. Inversion of the KL basis
				dims_data2[0] = 1;
				dims_data2[1] = this->nactu();
				carma_obj<float> *d_eigen = new carma_obj<float>(current_context, dims_data2);
				dims_data[1] = this->nactu();
				dims_data[2] = this->nactu();
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

				for (int i=0 ; i<this->nactu() ; i++){
						h_eigenvals_inv->getData()[i] = 1./h_eigenvals->getData()[i];
					}

				d_eigen->host2device(*h_eigenvals_inv);

				carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nactu(),nactu(),d_Vkl->getData(),nactu(),d_eigen->getData(),1,d_tmp->getData(),nactu());
				carma_gemm(cublas_handle, 't', 't', nactu(), nactu(), nactu(), 1.0f,
										  d_tmp->getData(), nactu(), d_Ukl->getData(), nactu(), 0.0f,
										  this->d_KLbasis->getData(), nactu());

				// 3. Projection of KL covariance matrix in the DM basis

				carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nactu(),nactu(),d_KLbasis->getData(),nactu(),d_KLcov->getData(),1,d_tmp->getData(),nactu());
				carma_gemm(cublas_handle, 't', 'n', nactu(), nactu(), nactu(), 1.0f,
							  d_tmp->getData(), nactu(), this->d_KLbasis->getData(), nactu(), 0.0f,
							  this->d_covmat->getData(), nactu());

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
				dims_data[1] = this->nactu();
				dims_data[2] = this->nactu();
				carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
				for (int i=0 ; i<this->nactu() ; i++){
								// Valeurs propres négatives.... A voir et inverser ordre si valeurs propres positives
								h_eigenvals->getData()[i] = -h_eigenvals->getData()[i];
								cout << h_eigenvals->getData()[i] << endl;
							}
				h_eigenvals->getData()[this->nactu() - 1] = 0.;
				d_KLcov->host2device(*h_eigenvals);

				carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, nactu(),nactu(),d_KLbasis->getData(),nactu(),d_KLcov->getData(),1,d_tmp->getData(),nactu());
				carma_gemm(cublas_handle, 'n', 't', nactu(), nactu(), nactu(), 1.0f,
									  d_tmp->getData(), nactu(), this->d_KLbasis->getData(), nactu(), 0.0f,
									  this->d_covmat->getData(), nactu());

				delete d_tmp;
			}
		}

		delete d_statcov;
		delete h_eigenvals;
		delete d_KLcov;
		delete d_indx;
		delete d_xpos;
		delete d_ypos;

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::do_geomat(carma_obj<float> *d_geocov, carma_obj<float> *d_IF, long n_pts, float ampli){
	carma_gemm(cublas_handle, 't', 'n', nactu(), nactu(), n_pts, 1.0f,
	      d_IF->getData(), n_pts, d_IF->getData(), n_pts, 0.0f,
	      d_geocov->getData(), nactu());
	mult_vect(d_geocov->getData(),ampli,this->nactu()*this->nactu(),this->device);

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::piston_filt(carma_obj<float> *d_statcov){
	long Nmod = d_statcov->getDims()[1];
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = Nmod;
	dims_data[2] = Nmod;
	carma_obj<float> *d_F = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);

	int N = d_statcov->getDims()[1] * d_statcov->getDims()[1];
	fill_filtmat(d_F->getData(),Nmod, N, device);

	carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
			  d_F->getData(), Nmod, d_statcov->getData(), Nmod, 0.0f,
			  d_tmp->getData(), Nmod);
	carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
			      d_tmp->getData(), Nmod, d_F->getData(), Nmod, 0.0f,
			      d_statcov->getData(), Nmod);

	delete d_tmp;
	delete d_F;

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::piston_filt_cphim(carma_obj<float> *d_cphim){
	long Nmod = d_cphim->getDims()[1];
	long dims_data[3];
	dims_data[0] = 2;
	dims_data[1] = Nmod;
	dims_data[2] = Nmod;
	carma_obj<float> *d_F = new carma_obj<float>(current_context, dims_data);
	dims_data[2] = nslope();
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);

	int N = Nmod * Nmod;
	fill_filtmat(d_F->getData(),Nmod, N, device);

	carma_gemm(cublas_handle, 'n', 'n', Nmod, nslope(), Nmod, 1.0f,
			  d_F->getData(), Nmod, d_cphim->getData(), Nmod, 0.0f,
			  d_tmp->getData(), Nmod);
	d_cphim->copy(d_tmp,1,1);

	delete d_tmp;
	delete d_F;

	return EXIT_SUCCESS;
}

int
sutra_controller_mv::invgen(carma_obj<float> *d_mat, float cond, int job){

	const long dims_data[3] = {2, d_mat->getDims()[1], d_mat->getDims()[2]};
	carma_obj<float> *d_U = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context,dims_data);
	int i;

	const long dims_data2[2]={1,d_mat->getDims()[1]};
	carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,dims_data2);
	carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

	d_U->copy(d_mat,1,1);
	carma_syevd<float,1>('V',d_U,h_eigenvals);
	//syevd_f('V',d_U,h_eigenvals);
	if (job == 1){ // Conditionnement
		float maxe = h_eigenvals->getData()[d_mat->getDims()[1] - 1];

		for (i=0 ; i<d_mat->getDims()[1] ; i++){
			if(h_eigenvals->getData()[i] < maxe/cond)
				h_eigenvals_inv->getData()[i] = 0.;
			else
				h_eigenvals_inv->getData()[i] = 1./h_eigenvals->getData()[i];
		}
	}
	if (job == 0){ // Filtre #cond modes
		for (i=0 ; i<d_mat->getDims()[1] ; i++){
					if(i < cond)
						h_eigenvals_inv->getData()[i] = 0.;
					else
						h_eigenvals_inv->getData()[i] = 1./h_eigenvals->getData()[i];
				}
	}

	d_eigenvals_inv->host2device(*h_eigenvals_inv);

	carma_dgmm(cublas_handle,CUBLAS_SIDE_RIGHT,d_mat->getDims()[1],d_mat->getDims()[2],
			d_U->getData(),d_mat->getDims()[1],d_eigenvals_inv->getData(),1,
			d_tmp->getData(),d_mat->getDims()[1]);
	carma_gemm<float>(cublas_handle, 'n', 't', d_mat->getDims()[1], d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
			      d_tmp->getData(), d_mat->getDims()[1], d_U->getData(), d_mat->getDims()[1], 0.0f,
			      d_mat->getData(), d_mat->getDims()[1]);

	delete d_U;
	delete d_tmp;
	delete d_eigenvals_inv;
	delete h_eigenvals;
	delete h_eigenvals_inv;


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
	const long dims_data[3] = {2, this->nactu(), this->nactu()};
	carma_obj<float> *d_M1 = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context,dims_data);
	carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context,dims_data);

	const long dims_data2[2]={1,this->nactu()};
	carma_obj<float> *d_eigenvals = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_sqrt = new carma_obj<float>(current_context,dims_data2);
	carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,dims_data2);
	carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
	carma_host_obj<float> *h_eigenvals_sqrt = new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

	// 1. SVdec(geocov,U) --> Ut * geocov * U = D²
	carma_syevd<float,1>('V', d_geocov, h_eigenvals);

	d_eigenvals->host2device(*h_eigenvals);
	for (int i=0 ; i<this->nactu() ; i++){
		h_eigenvals_sqrt->getData()[i] = sqrt(h_eigenvals->getData()[i]); // D = sqrt(D²)
		h_eigenvals_inv->getData()[i] = 1./sqrt(h_eigenvals->getData()[i]);// D⁻¹ = 1/sqrt(D²)
	}
	d_eigenvals_sqrt->host2device(*h_eigenvals_sqrt);
	d_eigenvals_inv->host2device(*h_eigenvals_inv);

	// 2. M⁻¹ = sqrt(eigenvals) * Ut : here, we have transpose(M⁻¹)

	carma_dgmm<float>(cublas_handle,CUBLAS_SIDE_RIGHT,this->nactu(),this->nactu(),
					d_geocov->getData(), this->nactu(), d_eigenvals_sqrt->getData(),1,
					d_M1->getData(), this->nactu());

	// 3. C' = M⁻¹ * statcov * M⁻¹t
	carma_gemm<float>(cublas_handle, 't', 'n', nactu(), nactu(), nactu(), 1.0f,
		      d_M1->getData(), nactu(), d_statcov->getData(), nactu(), 0.0f,
		      d_tmp->getData(), nactu());

	carma_gemm<float>(cublas_handle, 'n', 'n', nactu(), nactu(), nactu(), 1.0f,
			      d_tmp->getData(), nactu(), d_M1->getData(), nactu(), 0.0f,
			      d_tmp2->getData(), nactu());

	// 4. SVdec(C',A)
	carma_syevd<float,1>('V', d_tmp2, h_eigenvals);

	// 5. M = U * D⁻¹
	carma_dgmm<float>(cublas_handle,CUBLAS_SIDE_RIGHT,this->nactu(),this->nactu(),
				d_geocov->getData(), this->nactu(), d_eigenvals_inv->getData(),1,
				d_tmp->getData(), this->nactu());

	// 6. B = M * A;
	carma_gemm<float>(cublas_handle, 'n', 'n', nactu(), nactu(), nactu(), 1.0f,
			      d_tmp->getData(), nactu(), d_tmp2->getData(), nactu(), 0.0f,
			      d_KLbasis->getData(), nactu());

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
int sutra_controller_mv::build_cmat(float *Dm, float *Dtt, float cond){
	long *dims_data = new long[3];
	dims_data[0] = 2;
	dims_data[1] = nslope();
	dims_data[2] = nactu() - 2;
	carma_obj<float> *d_Dm = new carma_obj<float>(current_context, dims_data);
	dims_data[2] = 2;
	carma_obj<float> *d_Dtt = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_M = new carma_obj<float>(current_context, dims_data);
	dims_data[1] = nactu() - 2;
	dims_data[2] = nslope();
	carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_cmat_tt = new carma_obj<float>(current_context, dims_data);
	carma_obj<float> *d_Dm1 = new carma_obj<float>(current_context, dims_data);
	dims_data[2] = nactu() - 2;
	carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context, dims_data);
	dims_data[1] = nactu() - 2;
	dims_data[2] = 2;
	carma_obj<float> *d_TT2ho = new carma_obj<float>(current_context, dims_data);
	dims_data[1] = 2;
	dims_data[2] = nslope();
	carma_obj<float> *d_M1 = new carma_obj<float>(current_context, dims_data);
	dims_data[2] = 2;
	carma_obj<float> *d_tmp3 = new carma_obj<float>(current_context, dims_data);
	dims_data[1] = nslope();
	dims_data[2] = nslope();
	carma_obj<float> *d_Ftt = new carma_obj<float>(current_context, dims_data);

	// (Cmm + Cn)⁻¹
	add_md(this->d_Cmm->getData(),this->d_Cmm->getData(),this->d_noisemat->getData(), nslope(), this->device);
	invgen(this->d_Cmm,(float)(nslope()-nactu())/2,0);
	// Cphim * (Cmm + Cn)⁻¹
	carma_gemm(cublas_handle, 'n', 'n', nactu() - 2, nslope(), nslope(), 1.0f,
				this->d_Cphim->getData(), nactu() - 2, this->d_Cmm->getData(), nslope(), 0.0f,
				d_tmp->getData(), nactu() - 2);

	// Imat decomposition TT
	d_Dm->host2device(Dm);
	d_Dtt->host2device(Dtt);
	// Dm⁻¹
	carma_gemm(cublas_handle, 't', 'n', nactu() - 2, nactu() - 2, nslope(), 1.0f,
					d_Dm->getData(), nslope(), d_Dm->getData(), nslope(), 0.0f,
					d_tmp2->getData(), nactu() - 2);
	invgen(d_tmp2,cond,1);
	carma_gemm(cublas_handle, 'n', 't', nactu() - 2, nslope(), nactu() - 2, 1.0f,
						d_tmp2->getData(), nactu() - 2, d_Dm->getData(), nslope(), 0.0f,
						d_Dm1->getData(), nactu() - 2);

	// TT2ho = Dm⁻¹ * Dtt
	carma_gemm(cublas_handle, 'n', 'n', nactu() - 2, 2, nslope(), 1.0f,
							d_Dm1->getData(), nactu() - 2, d_Dtt->getData(), nslope(), 0.0f,
							d_TT2ho->getData(), nactu() - 2);
	// M = Dm * TT2ho
	carma_gemm(cublas_handle, 'n', 'n', nslope(), 2, nactu() - 2, 1.0f,
								d_Dm->getData(), nslope(), d_TT2ho->getData(), nactu() - 2, 0.0f,
								d_M->getData(), nslope());

	// M⁻¹
	carma_gemm(cublas_handle, 't', 'n', 2, 2, nslope(), 1.0f,
						d_M->getData(), nslope(), d_M->getData(), nslope(), 0.0f,
						d_tmp3->getData(), 2);
	invgen(d_tmp3,0.0f,0);
	carma_gemm(cublas_handle, 'n', 't', 2, nslope(), 2, 1.0f,
						d_tmp3->getData(), 2, d_M->getData(), nslope(), 0.0f,
						d_M1->getData(), 2);

	// M*M⁻¹
	carma_gemm(cublas_handle, 'n', 'n', nslope(), nslope(), 2, 1.0f,
								d_M->getData(), nslope(), d_M1->getData(), 2, 0.0f,
								d_Ftt->getData(), nslope());

	// TT filter
	TT_filt(d_Ftt->getData(),nslope(),this->device);

	//cmat without TT
	carma_gemm(cublas_handle, 'n', 'n', nactu()-2, nslope(), nslope(), 1.0f,
									d_tmp->getData(), nactu() - 2, d_Ftt->getData(), nslope(), 0.0f,
									d_cmat_tt->getData(), nactu() - 2);

	// Fill CMAT
	fill_cmat(this->d_cmat->getData(),d_cmat_tt->getData(),d_M1->getData(),nactu(),nslope(),this->device);

	delete d_Dm;
	delete d_Dtt;
	delete d_M;
	delete d_tmp;
	delete d_tmp2;
	delete d_tmp3;
	delete d_cmat_tt;
	delete d_Dm1;
	delete d_TT2ho;
	delete d_M1;
	delete d_Ftt;

	return EXIT_SUCCESS;
}
int sutra_controller_mv::build_cmat(const char *dmtype, char *method) {
  float one = 1.;
  float zero = 0.;

  if(strcmp(method,"inv") == 0){
	  //  R = (Dt*Cn⁻¹*D + Cphi⁻¹)⁻¹*Dt*Cn⁻¹

	  carma_obj<float> *d_tmp;
	  carma_obj<float> *d_tmp2;
	  carma_obj<float> *d_tmp3;
	  long *dims_data2 = new long[3];
	  dims_data2[0] = 2;
	  dims_data2[1] = nactu();
	  dims_data2[2] = nactu();
	  d_tmp = new carma_obj<float>(current_context, dims_data2);
	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
	  dims_data2[1] = nslope();
	  dims_data2[2] = nactu();
	  d_tmp3 = new carma_obj<float>(current_context,dims_data2);

	  carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, nslope(),nactu(),d_imat->getData(),nslope(),d_noisemat->getData(),1,d_tmp3->getData(),nslope());
	  carma_gemm(cublas_handle, 't', 'n', nactu(), nactu(), nslope(), one,
			d_tmp3->getData(), nslope(), d_imat->getData(), nslope(), zero,
			d_tmp2->getData(), nactu());
	  carma_geam(cublas_handle, 'n', 'n', nactu(), nactu(), one, d_tmp2->getData(),
		  nactu(), one, d_covmat->getData(), nactu(), d_tmp->getData(), nactu());

	  carma_potri<float>(d_tmp);

	  dims_data2[1] = nactu();
	  dims_data2[2] = nslope();
	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
	  carma_gemm(cublas_handle, 'n', 't', nactu(), nslope(), nactu(), one,
		  d_tmp->getData(), nactu(), d_imat->getData(), nslope(), zero,
		  d_tmp2->getData(), nactu());


	  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, nactu(),nslope(),d_tmp2->getData(),nactu(),d_noisemat->getData(),1,d_cmat->getData(),nactu());

	  delete d_tmp;
	  delete d_tmp2;
	  delete d_tmp3;
  }

  else if(strcmp(method,"n") == 0){
	  //  R = Cphi*Dt*(D*Cphi*Dt + Cn)⁻¹

  	  carma_obj<float> *d_tmp;
  	  carma_obj<float> *d_tmp2;
  	  carma_obj<float> *d_tmp3;
  	  carma_obj<float> *d_tmp4;
  	  long *dims_data2 = new long[3];
  	  dims_data2[0] = 2;
  	  dims_data2[1] = nslope();
  	  dims_data2[2] = nslope();
  	  d_tmp = new carma_obj<float>(current_context, dims_data2);
  	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
  	  d_tmp4 = new carma_obj<float>(current_context, dims_data2);
  	  d_U = new carma_obj<float>(current_context, dims_data2);
  	  dims_data2[1] = nslope();
  	  dims_data2[2] = nactu();
  	  d_tmp3 = new carma_obj<float>(current_context,dims_data2);
	  long *dims_data = new long[2];
	  dims_data[0] = 1;
	  dims_data[1] = this->nactu();
	  //carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data, MA_PAGELOCK);
	  //carma_obj<float> *d_eigenvals = new carma_obj<float>(current_context, dims_data);

  	  carma_gemm(cublas_handle, 'n', 'n', nslope(), nactu(), nactu(), one,
  			d_imat->getData(), nslope(), d_covmat->getData(), nactu(), zero,
  			d_tmp3->getData(), nslope());
 	  carma_gemm(cublas_handle, 'n', 't', nslope(), nslope(), nactu(), one,
  			d_tmp3->getData(), nslope(), d_imat->getData(), nslope(), zero,
  			d_tmp2->getData(), nslope());
 	  add_md(d_tmp4->getData(),d_tmp4->getData(),d_noisemat->getData(), nslope(), this->device);
  	  carma_geam(cublas_handle, 'n', 'n', nslope(), nslope(), one, d_tmp2->getData(),
  		  nslope(), one, d_tmp4->getData(), nslope(), d_tmp->getData(), nslope());

  	  carma_potri<float>(d_tmp);

  	  dims_data2[1] = nactu();
  	  dims_data2[2] = nslope();
  	  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
  	  carma_gemm(cublas_handle, 't', 'n', nactu(), nslope(), nslope(), one,
  		  d_imat->getData(), nslope(), d_tmp->getData(), nslope(), zero,
  		  d_tmp2->getData(), nactu());
  	  carma_gemm(cublas_handle, 'n', 'n', nactu(), nslope(), nactu(), one,
  		  d_covmat->getData(), nactu(), d_tmp2->getData(), nactu(), zero,
  		  d_cmat->getData(), nactu());

  	  delete d_tmp;
  	  delete d_tmp2;
  	  delete d_tmp3;
  	  delete d_tmp4;
    }
  else { }//y_error("Specify the computation method for mv : inv or n \n");}
  return EXIT_SUCCESS;
}

int sutra_controller_mv::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if (delay > 0) {
    for (int cc = 0; cc < delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope()]), 1,
          this->nslope(), this->device);

    cutilSafeCall(
        cudaMemcpy(&(this->d_cenbuff->getData()[delay * this->nslope()]),
            this->d_centroids->getData(), sizeof(float) * this->nslope(),
            cudaMemcpyDeviceToDevice));

    cutilSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
            sizeof(float) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

int sutra_controller_mv::comp_com() {

  this->frame_delay();
  this->d_com2->copy(this->d_com1, 1, 1);
  this->d_com1->copy(this->d_com, 1, 1);
  // POLC equations
  carma_geam<float>(cublas_handle, 'n', 'n', nactu(), 1, (float)(delay-1), *d_com2,
      nactu(), 1.0f - (delay-1), *d_com1, nactu(), *d_compbuff, nactu());
  carma_gemv<float>(cublas_handle, 'n', nslope(), nactu(), 1.0f, *d_imat, nslope(),
      *d_compbuff, 1, 0.0f, *d_compbuff2, 1);
  carma_geam<float>(cublas_handle, 'n', 'n', nslope(), 1, 1.0f, *d_centroids,
      nslope(), -1.0f, *d_compbuff2, nslope(), *d_olmeas, nslope());

  int nstreams = this->streams->get_nbStreams();
  if (nstreams > 1) {
    float alpha = -1.0f;
    float beta = 0.0f;

    for (int i = 0; i < nstreams; i++) {
      int istart1 = i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1)
          / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      carma_checkCublasStatus(
          cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
              this->d_cmat->getDims(2), &alpha,
              &((this->d_cmat->getData())[istart1]),
              this->d_cmat->getDims(1) / nstreams, *d_olmeas, 1, &beta,
              &((this->d_err->getData())[istart2]), 1));
    }

    this->streams->wait_all_streams();

  } else {
    // compute error
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
        d_olmeas, 1, 0.0f, 1); // POLC --> d_olmeas
  }
  	  /*
  mult_int(this->d_com->getData(), this->d_err->getData(),
          this->d_gain->getData(), this->gain, this->nactu(), this->device);*/

  carma_geam<float>(cublas_handle, 'n', 'n', nactu(), 1, this->gain, *d_err, nactu(),
        1.0f - this->gain, *d_com1, nactu(), *d_com, nactu());


  return EXIT_SUCCESS;
}
