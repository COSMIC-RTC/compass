#include <sutra_centroider_pyr.h>
#include <iostream>

using namespace std;

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context,
		sutra_sensors *sensors, int nwfs,
		long nvalid, float offset,
		float scale, int device) {
	if ((sensors->d_wfs[nwfs]->type != "pyr") && (sensors->d_wfs[nwfs]->type
			!= "pyrhr"))
		throw "sutra_centroider_roof expect a sutra_wfs_pyr_pyr4 sensor";
	this->current_context = context;

	this->device = device;
	context->set_activeDevice(device,1);
	this->wfs = sensors->d_wfs[nwfs];
	this->nwfs = nwfs;
	this->nvalid = nvalid;
	this->offset = offset;
	this->scale = scale;
	this->pyr_type = sensors->d_wfs[nwfs]->type;

}

sutra_centroider_pyr::~sutra_centroider_pyr() {

}

string sutra_centroider_pyr::get_type() {
	return this->pyr_type;
}

int sutra_centroider_pyr::get_cog(carma_streams *streams, float *cube,
		float *subsum, float *centroids, int nvalid,
		int npix, int ntot) {
	//TODO: Implement sutra_centroider_pyr::get_cog
	cerr << "get_cog not implemented\n";

	return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_pyr(float *cube, float *subsum, float *centroids,
		int *subindx, int *subindy, int nvalid,
		int ns, int nim) {

	pyr_subsum(this->wfs->d_subsum->getData(), this->wfs->d_binimg->getData(),
			this->wfs->d_validsubsx->getData(), this->wfs->d_validsubsy->getData(),
			this->wfs->nfft / this->wfs->nrebin, this->wfs->nvalid,
			this->current_context->get_device(device));


	current_context->set_activeDevice(device,1);
	if (this->pyr_type == "pyr")
		pyr_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim,
				this->current_context->get_device(device));
	else if (this->pyr_type == "pyrhr")
		pyr2_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid,
				this->current_context->get_device(device));
	else
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(float *subsum, float *slopes) {
	if (this->pyr_type == "pyr")
		return this->get_pyr(*(wfs->d_bincube), *(wfs->d_subsum), slopes,
				*(wfs->d_validsubsx), *(wfs->d_validsubsy),
				wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
	else if (this->pyr_type == "pyrhr")
		return this->get_pyr(*(wfs->d_binimg), *(wfs->d_subsum), slopes,
				*(wfs->d_validsubsx), *(wfs->d_validsubsy),
				wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
	return EXIT_FAILURE;
}

int sutra_centroider_pyr::get_cog() {
	return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes));
}

