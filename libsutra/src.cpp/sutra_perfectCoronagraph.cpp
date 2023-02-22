// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_perfectCoronograph.cpp
//! \ingroup   libsutra
//! \class     SutraPerfectCoronograph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <sutra_perfectCoronagraph.h>

SutraPerfectCoronagraph::SutraPerfectCoronagraph(CarmaContext *context, SutraSource *d_source, 
                                    int im_dimx, int im_dimy, float *wavelength, int nWavelength, 
                                    int device):
    SutraCoronagraph(context, "perfect", d_source, im_dimx, im_dimy, wavelength, nWavelength, device),
    tmp_mft(nullptr) {}

int SutraPerfectCoronagraph::_set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm0, bool psf) {
    std::vector<CarmaObj<cuFloatComplex>*> *vecA = &AA;
    std::vector<CarmaObj<cuFloatComplex>*> *vecB = &BB;
    std::vector<float> *vecNorm = &norm;
    if(psf) {
        vecA = &AA_psf;
        vecB = &BB_psf;
        vecNorm = &norm_psf;
    }

    if(!vecA->empty()) {
        vecA->clear();
    }
    if(!vecB->empty()) {
        vecB->clear();
    }
    if(!vecNorm->empty()) {
        vecNorm->clear();
    }

    long dims[3];
    dims[0] = 2;
    for (int i = 0; i < wavelength.size() ; i++) {
        dims[1] = imageDimx;
        dims[2] = pupDimx;
        vecA->push_back(new CarmaObj<cuFloatComplex>(current_context, dims, A + i * imageDimx * pupDimx));
        dims[1] = pupDimy;
        dims[2] = imageDimy;
        vecB->push_back(new CarmaObj<cuFloatComplex>(current_context, dims, B + i * imageDimy * pupDimy));
        vecNorm->push_back(norm0[i]);
    }
    if (tmp_mft == nullptr) {
        dims[1] = imageDimx;
        dims[2] = pupDimy;
        tmp_mft = new CarmaObj<cuFloatComplex>(current_context, dims);
    }
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm0) {
    return _set_mft(A, B, norm0, false);
}

int SutraPerfectCoronagraph::set_mft_psf(cuFloatComplex *A, cuFloatComplex *B, float* norm0) {
    return _set_mft(A, B, norm0, true);
}

int SutraPerfectCoronagraph::_compute_image(bool psf, bool accumulate, bool remove_coro) {
    CarmaObj<float> *img_se = d_image_se;
    CarmaObj<float> *img_le = d_image_le;
    std::vector<CarmaObj<cuFloatComplex>*> *mftA = &AA;
    std::vector<CarmaObj<cuFloatComplex>*> *mftB = &BB;
    std::vector<float> *mftNorm = &norm;
    if (psf) {
        img_se = d_psf_se;
        img_le = d_psf_le;
        mftA = &AA_psf;
        mftB = &BB_psf;
        mftNorm = &norm_psf;
    }

    img_se->reset();
    for (int i = 0 ; i < wavelength.size(); i++) {
        compute_electric_field(i);
        if(!remove_coro) {
            cuFloatComplex sum = d_electric_field->sum();
            remove_complex_avg(d_electric_field->get_data(), d_electric_field->sum(),
                            d_pupil->get_data(), d_source->d_wherephase->get_nb_elements(),
                            pupDimx, pupDimy, current_context->get_device(device));
        }
        mft((*mftA)[i], (*mftB)[i], tmp_mft, d_electric_field, d_complex_image, (*mftNorm)[i]);
        accumulate_abs2(d_complex_image->get_data(), img_se->get_data(), 
                        img_se->get_nb_elements(), current_context->get_device(device));
    }

    if(accumulate) {
        img_le->axpy(1.0f, img_se, 1, 1);
        cnt += 1;
    }
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::compute_psf(bool accumulate) {
    return _compute_image(true, accumulate, true);
}

int SutraPerfectCoronagraph::compute_image(bool accumulate) {
    return _compute_image(false, accumulate, false);
}