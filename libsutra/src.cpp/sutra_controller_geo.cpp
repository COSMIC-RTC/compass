#include <sutra_controller_geo.h>

sutra_controller_geo::sutra_controller_geo(carma_context *context, long nactu, long Nphi,
      long delay):
      sutra_controller(context, 0, nactu){

	this->delay = delay;
	this->device = device;
	this->gain = 0.0f;
	this->Nphi = Nphi;

	long dims_data2[3];
	dims_data2[0] = 2;
	dims_data2[1] = nactu;
	dims_data2[2] = Nphi;
	this->d_proj = new carma_obj<float>(this->current_context, dims_data2);
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
	dims_data1[1] = Nphi;
	this->d_phi = new carma_obj<float>(this->current_context, dims_data1);
	this->d_indx_pup = new carma_obj<int>(this->current_context, dims_data1);

}

sutra_controller_geo::~sutra_controller_geo() {
  current_context->set_activeDevice(device);
  delete this->d_proj;
  delete this->d_gain;
  delete this->d_indx_pup;
  delete this->d_phi;

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
	  p->second->get_IF(d_IF.getData(indx_start * this->Nphi), d_indx.getData(this->Nphi * ind), this->Nphi, 1.0f/*unitpervolt[ind]*/);
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
sutra_controller_geo::comp_dphi(sutra_source *target){

	// Get the target phase in the pupil
	get_pupphase(this->d_phi->getData(),target->d_phase->d_screen->getData(), this->d_indx_pup->getData(), this->Nphi, this->current_context->get_device(device));
	remove_avg(this->d_phi->getData(), this->d_phi->getNbElem(), current_context->get_device(device));

	return EXIT_SUCCESS;
}

int
sutra_controller_geo::comp_com(){
	// Project the phase on the actuators
	carma_gemv(cublas_handle(),'n', nactu(), this->Nphi, -1.0f, this->d_proj->getData(),this->d_proj->getDims()[1],
			this->d_phi->getData(),1, 0.0f, this->d_com->getData(),1);

	return EXIT_SUCCESS;
}
