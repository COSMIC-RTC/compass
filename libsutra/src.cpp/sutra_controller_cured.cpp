#include <sutra_controller_cured.h>
#include <string>

sutra_controller_cured::sutra_controller_cured(carma_context *context,
		long nvalid, long nactu, long delay) :
		sutra_controller(context, nvalid * 2, nactu) {

	this->gain = 0;
	this->h_centroids = 0L;
	this->h_err = 0L;
	this->h_parcure = 0L;
	this->h_syscure = 0L;

	long dims_data1[2] = { 1, 0 };

	dims_data1[1] = nvalid * 2;
	this->d_centroids = new carma_obj<float>(context, dims_data1);

	dims_data1[1] = nactu;
	this->h_err = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

}

sutra_controller_cured::~sutra_controller_cured() {
	current_context->set_activeDevice(device);

	if (this->h_centroids != 0L)
		delete this->h_centroids;
	if (this->h_err != 0L)
		delete this->h_err;

}

string sutra_controller_cured::get_type() {
	return "cured";
}

int sutra_controller_cured::set_gain(float gain) {
	this->gain = gain;
	return EXIT_SUCCESS;
}

int sutra_controller_cured::comp_com() {

	h_centroids->cpy_obj(this->d_centroids, cudaMemcpyDeviceToHost);

	cured(this->h_syscure, this->h_parcure, this->h_centroids->getData(),
			this->h_err->getData(), this->gain);

	h_err->cpy_obj(this->d_com, cudaMemcpyHostToDevice);

	//mult_int(this->d_com->getData(),this->d_err->getData(),this->d_gain->getData(),this->gain,this->nactu,this->device);

	return EXIT_SUCCESS;
}

int sutra_controller_cured::init_cured(int nxsubs, int *isvalid) {
	this->h_syscure = cureSystem(nxsubs, this->nslope() / 2., this->nactu(),
			isvalid, 1);
	this->h_parcure = cureInit(this->h_syscure);
	return EXIT_SUCCESS;
}

