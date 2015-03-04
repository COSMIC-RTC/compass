#include <sutra_controller_geo.h>

sutra_controller_geo::sutra_controller_geo(carma_context *context, long nactu, long Nphi,
      long delay):
      sutra_controller(context, 0, nactu){

	this->delay = delay;
	this->device = device;
	this->gain = 0.0f;
	this->Nphi = Nphi;

//	long dims_data2[3];
//	dims_data2[0] = 2;
//	dims_data2[1] = nactu;
//	dims_data2[2] = Nphi;
	//this->d_proj = new carma_obj<float>(this->current_context, dims_data2);
	this->d_proj = 0L;
	this->d_geocov = 0L;
	this->d_IFsparse = 0L;
	/*
	if (delay > 0) {
	    dims_data2[1] = Nphi;
	    dims_data2[2] = delay + 1;
	    this->d_cenbuff = new carma_obj<float>(this->current_context, dims_data2);
	  }
	  */
	long dims_data1[2];
	dims_data1[0] = 1;
	dims_data1[1] = nactu;
	this->d_gain = new carma_obj<float>(this->current_context, dims_data1);
	this->d_compfloat = new carma_obj<float>(this->current_context, dims_data1);
	this->d_compdouble = new carma_obj<double>(this->current_context, dims_data1);
	dims_data1[1] = Nphi;
	this->d_phi = new carma_obj<double>(this->current_context, dims_data1);
	this->d_indx_pup = new carma_obj<int>(this->current_context, dims_data1);

}

sutra_controller_geo::~sutra_controller_geo() {
  current_context->set_activeDevice(device,1);
  delete this->d_proj;
  delete this->d_gain;
  delete this->d_indx_pup;
  delete this->d_phi;
  delete this->d_compfloat;
  delete this->d_compdouble;

}

string sutra_controller_geo::get_type() {
  return "geo";
}

int sutra_controller_geo::set_gain(float gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

int sutra_controller_geo::set_delay(int delay) {
  this->delay = delay;
  return EXIT_SUCCESS;
}

int sutra_controller_geo::load_mgain(float *mgain) {
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

int
sutra_controller_geo::init_proj(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup){
  current_context->set_activeDevice(device,1);
	long dims_data[3] = {2,this->Nphi,nactu()};
	carma_obj<float> d_IF(current_context,dims_data);
	dims_data[1] = nactu();
	carma_obj<float> d_tmp(current_context,dims_data);
	long dims_data1[2] = {1,this->Nphi * dms->ndm};
	carma_obj<int> d_indx(current_context,dims_data1, indx_dm);

	this->d_indx_pup->host2device(indx_pup);

	// Get influence functions in d_IF
	int indx_start = 0;
	int ind = 0;
	map<type_screen, sutra_dm *>::iterator p;
	p = dms->d_dms.begin();
	while (p != dms->d_dms.end()) {
	  p->second->get_IF<float>(d_IF.getData(indx_start * this->Nphi), d_indx.getData(this->Nphi * ind), this->Nphi, 1.0f/*unitpervolt[ind]*/);
	  indx_start += p->second->ninflu;
	  ind++;
	  p++;
	}

	// d_tmp = (transpose(d_IF) * d_IF)⁻¹
	carma_gemm<float>(cublas_handle(), 't', 'n', nactu(), nactu(), this->Nphi, 1.0f,
	  			      d_IF.getData(), d_IF.getDims()[1], d_IF.getData(), d_IF.getDims()[1], 0.0f,
	  			      d_tmp.getData(), d_tmp.getDims()[1]);
	carma_potri(&d_tmp);
	//invgen(d_tmp,1000.0f,1);

	//d_proj = d_tmp * transpose(d_IF)
	carma_gemm<float>(cublas_handle(), 'n', 't', nactu(), this->Nphi, nactu(), 1.0f,
		  			      d_tmp.getData(), d_tmp.getDims()[1], d_IF.getData(), d_IF.getDims()[1], 0.0f,
		  			      this->d_proj->getData(), this->d_proj->getDims()[1]);

	return EXIT_SUCCESS;
}

int
sutra_controller_geo::init_proj_sparse(sutra_dms *dms, int *indx_dm, float *unitpervolt, int *indx_pup){

  current_context->set_activeDevice(device,1);
  carma_sparse_obj<double> *d_IFi[dms->ndm];
	long dims_data1[2] = {1,this->Nphi * dms->ndm};
	carma_obj<int> d_indx(current_context,dims_data1, indx_dm);

	this->d_indx_pup->host2device(indx_pup);

	// Get influence functions of the DM #ind in d_IFi
	int indx_start = 0;
	int ind = 0;
	int nnz = 0;
	int NNZ[dms->ndm];
	int Nact[dms->ndm];
	map<type_screen, sutra_dm *>::iterator p;
	p = dms->d_dms.begin();
	while (p != dms->d_dms.end()) {
	  p->second->get_IF_sparse<double>(d_IFi[ind], d_indx.getData(this->Nphi * ind), this->Nphi, 1.0f,1);
	  NNZ[ind] = d_IFi[ind]->nz_elem;
	  Nact[ind] = p->second->ninflu;
	  nnz += d_IFi[ind]->nz_elem;
	  indx_start += p->second->ninflu;
	  ind++;
	  p++;
	}
	//Create global d_IF_sparse from array of d_IFi
	long dims_data[2] = {1,nnz};
	carma_obj<double> d_val(current_context,dims_data);
	carma_obj<int> d_col(current_context,dims_data);
	dims_data[1] = this->nactu() + 1;
	carma_obj<int> d_row(current_context,dims_data);
	int cpt[dms->ndm];
	int nact = 0;
	cpt[0] = 0;
	p = dms->d_dms.begin();

	for (int i = 0; i < dms->ndm; i++) {
		cutilSafeCall(
				cudaMemcpyAsync(d_val.getData(cpt[i]), d_IFi[i]->d_data,
						sizeof(double) * d_IFi[i]->nz_elem, cudaMemcpyDeviceToDevice));
		cutilSafeCall(
				cudaMemcpyAsync(d_col.getData(cpt[i]), d_IFi[i]->d_colind,
						sizeof(int) * d_IFi[i]->nz_elem, cudaMemcpyDeviceToDevice));
		if(i == 0)
			cutilSafeCall(
						cudaMemcpyAsync(d_row.getData(), d_IFi[i]->d_rowind, sizeof(int) * (p->second->ninflu + 1),
							cudaMemcpyDeviceToDevice));
		else
			cutilSafeCall(
						cudaMemcpyAsync(d_row.getData(nact+1), &(d_IFi[i]->d_rowind[1]), sizeof(int) * (p->second->ninflu),
							cudaMemcpyDeviceToDevice));
		cpt[i+1] = cpt[i] + d_IFi[i]->nz_elem;
		nact += p->second->ninflu;
		p++;
		delete d_IFi[i];
	}

	if(dms->ndm > 1){
		dims_data[1] = dms->ndm;
		carma_obj<int>d_NNZ(current_context,dims_data);
		carma_obj<int>d_nact(current_context,dims_data);
		d_NNZ.host2device(NNZ);
		d_nact.host2device(Nact);

		adjust_csr_index(d_row.getData(1),d_NNZ.getData(),d_nact.getData(),this->nactu(),Nact[0],current_context->get_device(device));
	}

	long dims_data2[3] = {2,dms->nact_total(),Nphi};
	this->d_IFsparse = new carma_sparse_obj<double>(current_context, dims_data2,
				d_val.getData(), d_col.getData(), d_row.getData(),nnz,false);

	// d_geocov = (transpose(d_IF) * d_IF)⁻¹
	carma_sparse_obj<double> *d_tmp = new carma_sparse_obj<double>(this->current_context);
	dims_data2[2] = dms->nact_total();
	this->d_geocov = new carma_obj<float>(current_context,dims_data2);
	carma_obj<double> *d_tmp2 = new carma_obj<double>(this->current_context,this->d_geocov->getDims());

	carma_gemm<double>(cusparse_handle(),'n','t',this->d_IFsparse,this->d_IFsparse,d_tmp);
	carma_csr2dense<double>(d_tmp, d_tmp2->getData());
	doubletofloat(d_tmp2->getData(),this->d_geocov->getData(),this->d_geocov->getNbElem(),current_context->get_device(device));

	carma_potri(d_geocov);

	delete d_tmp;
	delete d_tmp2;

	return EXIT_SUCCESS;
}

int
sutra_controller_geo::comp_dphi(sutra_source *target){

  current_context->set_activeDevice(device,1);
	// Get the target phase in the pupil
	get_pupphase(this->d_phi->getData(),target->d_phase->d_screen->getData(), this->d_indx_pup->getData(), this->Nphi, this->current_context->get_device(device));
	remove_avg(this->d_phi->getData(), this->d_phi->getNbElem(), current_context->get_device(device));

	return EXIT_SUCCESS;
}

int
sutra_controller_geo::comp_com(){
	// Project the phase on the actuators
	/*
	//Dense version
	carma_gemv(cublas_handle(),'n', nactu(), this->Nphi, -1.0f, this->d_proj->getData(),this->d_proj->getDims()[1],
			this->d_phi->getData(),1, 0.0f, this->d_com->getData(),1);
*/
  current_context->set_activeDevice(device,1);
	// Sparse version
	carma_gemv(cusparse_handle(),'n',1.0,this->d_IFsparse,this->d_phi,0.0,this->d_compdouble);
	doubletofloat(this->d_compdouble->getData(),this->d_compfloat->getData(),this->nactu(),current_context->get_device(device));
	carma_gemv(cublas_handle(),'n', nactu(), nactu(), -1.0f, this->d_geocov->getData(),this->d_geocov->getDims()[1],
				this->d_compfloat->getData(),1, 0.0f, this->d_com->getData(),1);

	return EXIT_SUCCESS;
}
