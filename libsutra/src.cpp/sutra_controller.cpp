#include <sutra_controller.h>
#include <string>

sutra_controller::sutra_controller(carma_context* context, int nslope,
    int nactu) {
  this->current_context = context;
  this->device = context->get_activeDevice();

  int nstreams = 1; //nvalid/10;
  while (nactu % nstreams != 0)
    nstreams--;

  cerr << "controller uses " << nstreams << " streams" << endl;
  streams = new carma_streams(nstreams);

  long dims_data1[2] = { 1, 0 };

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<float>(context, dims_data1);
  dims_data1[1] = nactu;
  this->d_com = new carma_obj<float>(context, dims_data1);

}

sutra_controller::~sutra_controller() {
  delete this->streams;

  delete this->d_centroids;
  delete this->d_com;
}
int
sutra_controller::syevd_f(char meth, carma_obj<float> *d_U, carma_host_obj<float> *h_eigenvals){
	// Init double arrays
	const long dims_data[3] = {2, d_U->getDims()[1], d_U->getDims()[2]};
	carma_obj<double> *d_Udouble = new carma_obj<double>(current_context,dims_data);
	const long dims_data2[2] = {1, h_eigenvals->getDims()[1]};
	carma_host_obj<double> *h_eigen_double = new carma_host_obj<double>(dims_data2,MA_PAGELOCK);

	// Copy float array in double array
	floattodouble(d_U->getData(),d_Udouble->getData(),d_U->getNbElem(), this->device);

	//Doing syevd<double>
	carma_syevd<double,1>(meth,d_Udouble,h_eigen_double);

	// Reverse copy
	doubletofloat(d_Udouble->getData(),d_U->getData(),d_U->getNbElem(), this->device);
	for(int cc=0 ; cc<h_eigenvals->getNbElem() ; cc++){
		h_eigenvals->getData()[cc] = (float)h_eigen_double->getData()[cc];
	}

	delete d_Udouble;
	delete h_eigen_double;

	return EXIT_SUCCESS;
}
