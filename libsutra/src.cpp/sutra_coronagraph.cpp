// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_coronagraph.cpp
//! \ingroup   libsutra
//! \class     SutraCoronagraph
//! \brief     this class provides the coronagraph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <sutra_coronagraph.h>

SutraCoronagraph::SutraCoronagraph(CarmaContext *context, std::string type, SutraSource *d_source, 
                                    int dimx, int dimy, float *wavelength, int nWavelength, 
                                    int device):
    current_context(context),
    device(device),
    type(type),
    d_source(d_source),
    imageDimx(dimx),
    imageDimy(dimy),
    cntPsf(0),
    cntImg(0) {
    current_context->set_active_device(device, 1);
    long dims[3] = {2, dimx, dimy};
    d_image_se = new CarmaObj<float>(current_context, dims);
    d_image_le = new CarmaObj<float>(current_context, dims);
    d_psf_se = new CarmaObj<float>(current_context, dims);
    d_psf_le = new CarmaObj<float>(current_context, dims);
    d_pupil = d_source->d_pupil;
    pupDimx = d_source->d_phase->d_screen->get_dims(1);
    pupDimy = d_source->d_phase->d_screen->get_dims(2);
    d_electric_field = new CarmaObj<cuFloatComplex>(current_context, 
                                                    d_source->d_phase->d_screen->get_dims());
    d_complex_image = new CarmaObj<cuFloatComplex>(current_context, dims);
    for (int i = 0; i < nWavelength; i++) {
        this->wavelength.push_back(wavelength[i]);
        amplitude.push_back(new CarmaObj<float>(current_context, d_source->d_phase->d_screen->get_dims()));
        amplitude.back()->memset(1);
    }
}

int SutraCoronagraph::reset() {
    d_image_se->reset();
    d_image_le->reset();
    d_psf_se->reset();
    d_psf_le->reset();
    cntPsf = 0;
    cntImg = 0;
    return EXIT_SUCCESS;
}

int SutraCoronagraph::compute_electric_field(int wavelengthIndex) {
    float scale = 2 * CARMA_PI / wavelength[wavelengthIndex];
    ::compute_electric_field(d_electric_field->get_data(), d_source->d_phase->d_screen->get_data(), 
                            scale, amplitude[wavelengthIndex]->get_data(), d_pupil->get_data(), 
                            pupDimx, pupDimy, current_context->get_device(device));
    return EXIT_SUCCESS;
}

int SutraCoronagraph::mft(CarmaObj<cuFloatComplex> *A, CarmaObj<cuFloatComplex> *B, 
                CarmaObj<cuFloatComplex> *Ainput,
                CarmaObj<cuFloatComplex> *input, CarmaObj<cuFloatComplex> *output, float norm) {
    cuFloatComplex alpha(make_float2(norm, 0));
    cuFloatComplex beta(make_float2(0, 0));
    carma_gemm<cuFloatComplex>(current_context->get_cublas_handle(), 'n', 'n', A->get_dims(1), 
                input->get_dims(2), input->get_dims(1), alpha, A->get_data(), A->get_dims(1),
                input->get_data(), input->get_dims(1), 
                beta, Ainput->get_data(), Ainput->get_dims(1));
    alpha.x = 1.f;
    carma_gemm<cuFloatComplex>(current_context->get_cublas_handle(), 'n', 'n', Ainput->get_dims(1), 
                B->get_dims(2), B->get_dims(1), alpha, Ainput->get_data(), Ainput->get_dims(1),
                B->get_data(), B->get_dims(1), 
                beta, output->get_data(), output->get_dims(1));
    return EXIT_SUCCESS;
}

int SutraCoronagraph::set_amplitude(float *ampli) {
    long dims[3] = {2, pupDimx, pupDimy};
    for (int i = 0; i < wavelength.size() ; i++) {
        amplitude[i]->host2device(ampli + i * dims[1] * dims[2]);
    }

    return EXIT_SUCCESS;
}
