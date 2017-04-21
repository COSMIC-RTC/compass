#include <sutra_centroider_pyr.h>
#include <iostream>
#include <string>

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context,
        sutra_sensors *sensors, int nwfs, long nvalid, float offset,
        float scale, int device) {
    if ((sensors->d_wfs[nwfs]->type != "pyr")
            && (sensors->d_wfs[nwfs]->type != "pyrhr"))
        throw "sutra_centroider_roof expect a sutra_wfs_pyr_pyr4 sensor";
    this->current_context = context;

    this->device = device;
    context->set_activeDevice(device, 1);
    this->wfs = sensors->d_wfs[nwfs];
    this->nwfs = nwfs;
    this->nvalid = nvalid;
    this->offset = offset;
    this->scale = scale;
    this->valid_thresh = 1e-4;
    this->pyr_type = sensors->d_wfs[nwfs]->type;
    this->method = global;
}

sutra_centroider_pyr::~sutra_centroider_pyr() {
}

string sutra_centroider_pyr::get_type() {
    return this->pyr_type;
}

int sutra_centroider_pyr::set_valid_thresh(float valid_thresh) {
    this->valid_thresh = valid_thresh;
    return EXIT_SUCCESS;
}
float sutra_centroider_pyr::get_valid_thresh() {
    return this->valid_thresh;
}

int sutra_centroider_pyr::set_method(int type){
    switch (type) {
        case 0: this->method = local ; break;
        case 1: this->method = global; break;
        default: return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


int sutra_centroider_pyr::get_cog(carma_streams *streams, float *cube,
        float *subsum, float *centroids, int nvalid, int npix, int ntot) {
    // TODO(Implement sutra_centroider_pyr::get_cog)
    std::cerr << "get_cog not implemented" << std::endl;

    return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_pyr(float *cube, float *subsum, float *centroids,
        int *subindx, int *subindy, int nvalid, int ns, int nim) {
    current_context->set_activeDevice(device, 1);
    if (this->pyr_type == "pyr" || this->pyr_type == "roof") {
        DEBUG_TRACE("Not working anymore");
        throw "Not working anymore";
        // pyr_subsum(subsum, cube, subindx, subindy, ns, nvalid, nim,
        //         this->current_context->get_device(device));
        //
        // pyr_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim,
        //         this->current_context->get_device(device));
    } else if (this->pyr_type == "pyrhr") {
        pyr_subsum(subsum, cube, subindx, subindy, ns, nvalid,
                this->current_context->get_device(device));

        if(this->method == global){
            int blocks, threads;
        //  getNumBlocksAndThreads(current_context->get_device(device), this->d_binimg->getNbElem(),
        //      blocks, threads);
        	this->current_context->set_activeDevice(device,1);
        	sumGetNumBlocksAndThreads(nvalid,
        			this->current_context->get_device(device), blocks, threads);
            reduce(nvalid, threads, blocks, subsum, centroids);

            fillvalues(subsum, &centroids[0], nvalid,
                this->current_context->get_device(device));
        }

        pyr2_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, this->scale,
                this->valid_thresh, this->current_context->get_device(device));
    } else {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(float *subsum, float *slopes, bool noise) {
	if (this->pyr_type == "pyr" || this->pyr_type == "roof")
		return this->get_pyr(*(wfs->d_bincube), subsum, slopes,
				*(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
				wfs->nfft / wfs->nrebin, 4);
	else if (this->pyr_type == "pyrhr"){
		if(noise || wfs->error_budget == false){
			return this->get_pyr(*(wfs->d_binimg), subsum, slopes,
				*(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
				wfs->nfft / wfs->nrebin, 4);
		}
		else
			return this->get_pyr(*(wfs->d_binimg_notnoisy), subsum, slopes,
			*(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
			wfs->nfft / wfs->nrebin, 4);
	}
	return EXIT_FAILURE;
}

int sutra_centroider_pyr::get_cog() {
	return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes),true);
}
