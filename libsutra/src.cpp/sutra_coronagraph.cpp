// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_coronagraph.cpp
//! \ingroup   libsutra
//! \class     SutraCoronagraph
//! \brief     this class provides the coronagraph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_coronagraph.hpp>

SutraCoronagraph::SutraCoronagraph(CarmaContext *context, std::string type, SutraSource *d_source,
                                    int32_t dimx, int32_t dimy, float *wavelength, int32_t nWavelength,
                                    int32_t device):
    current_context(context),
    device(device),
    type(type),
    d_source(d_source),
    imageDimx(dimx),
    imageDimy(dimy),
    cntPsf(0),
    cntImg(0) {
    current_context->set_active_device(device, 1);
    int64_t dims[3] = {2, dimx, dimy};
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
    for (int32_t i = 0; i < nWavelength; i++) {
        this->wavelength.push_back(wavelength[i]);
        amplitude.push_back(new CarmaObj<float>(current_context, d_source->d_phase->d_screen->get_dims()));
        amplitude.back()->memset(1);
    }
}

int32_t SutraCoronagraph::reset() {
    d_image_se->reset();
    d_image_le->reset();
    d_psf_se->reset();
    d_psf_le->reset();
    cntPsf = 0;
    cntImg = 0;
    return EXIT_SUCCESS;
}

int32_t SutraCoronagraph::compute_electric_field(int32_t wavelengthIndex) {
    float scale = 2 * CARMA_PI / wavelength[wavelengthIndex];
    ::compute_electric_field(d_electric_field->get_data(), d_source->d_phase->d_screen->get_data(),
                            scale, amplitude[wavelengthIndex]->get_data(), d_pupil->get_data(),
                            pupDimx, pupDimy, current_context->get_device(device));
    return EXIT_SUCCESS;
}

int32_t SutraCoronagraph::mft(CarmaObj<cuFloatComplex> *A, CarmaObj<cuFloatComplex> *B,
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

int32_t SutraCoronagraph::set_amplitude(float *ampli) {
    int64_t dims[3] = {2, pupDimx, pupDimy};
    for (int32_t i = 0; i < wavelength.size() ; i++) {
        amplitude[i]->host2device(ampli + i * dims[1] * dims[2]);
    }

    return EXIT_SUCCESS;
}
