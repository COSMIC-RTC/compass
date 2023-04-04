// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_perfectCoronograph.cpp
//! \ingroup   libsutra
//! \class     SutraPerfectCoronograph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <sutra_perfectCoronagraph.h>

SutraPerfectCoronagraph::SutraPerfectCoronagraph(CarmaContext *context, SutraSource *d_source, 
                                    int im_dimx, int im_dimy, float *wavelength, int nWavelength, 
                                    int device):
    SutraCoronagraph(context, "perfect", d_source, im_dimx, im_dimy, wavelength, nWavelength, device),
    tmp_mft(nullptr) {

        AA = {
            {"img", std::vector<CarmaObj<cuFloatComplex>*>(nWavelength, nullptr)},
            {"psf", std::vector<CarmaObj<cuFloatComplex>*>(nWavelength, nullptr)}
        };
        BB = {
            {"img", std::vector<CarmaObj<cuFloatComplex>*>(nWavelength, nullptr)},
            {"psf", std::vector<CarmaObj<cuFloatComplex>*>(nWavelength, nullptr)}
        };
        norm = {
            {"img", std::vector<float>(nWavelength, 1)},
            {"psf", std::vector<float>(nWavelength, 1)}
        };
    }

int SutraPerfectCoronagraph::set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm0, 
                                        std::string mftType) {
    if(AA.count(mftType) < 1) {
        std::cout << "Invalid mftType. Must be img or psf" << std::endl;
        return EXIT_FAILURE;
    }

    long dims[3];
    dims[0] = 2;
    for (int i = 0; i < wavelength.size() ; i++) {
        dims[1] = imageDimx;
        dims[2] = pupDimx;
        AA[mftType][i] = new CarmaObj<cuFloatComplex>(current_context, dims, A + i * imageDimx * pupDimx);
        dims[1] = pupDimy;
        dims[2] = imageDimy;
        BB[mftType][i] = new CarmaObj<cuFloatComplex>(current_context, dims, B + i * imageDimy * pupDimy);
        norm[mftType][i] = norm0[i];
    }

    if (tmp_mft == nullptr) {
        dims[1] = imageDimx;
        dims[2] = pupDimy;
        tmp_mft = new CarmaObj<cuFloatComplex>(current_context, dims);
    }
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::_compute_image(bool psf, bool accumulate, bool remove_coro) {
    CarmaObj<float> *img_se = d_image_se;
    CarmaObj<float> *img_le = d_image_le;
    std::vector<CarmaObj<cuFloatComplex>*> mftA = AA["img"];
    std::vector<CarmaObj<cuFloatComplex>*> mftB = BB["img"];
    std::vector<float> mftNorm = norm["img"];
    if (psf) {
        img_se = d_psf_se;
        img_le = d_psf_le;
        mftA = AA["psf"];
        mftB = BB["psf"];
        mftNorm = norm["psf"];
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
        mft(mftA[i], mftB[i], tmp_mft, d_electric_field, d_complex_image, mftNorm[i]);
        accumulate_abs2(d_complex_image->get_data(), img_se->get_data(), 
                        img_se->get_nb_elements(), current_context->get_device(device));
    }

    if(accumulate) {
        img_le->axpy(1.0f, img_se, 1, 1);
        if(psf)
            cntPsf += 1;
        else 
            cntImg += 1;
    }
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::compute_psf(bool accumulate) {
    return _compute_image(true, accumulate, true);
}

int SutraPerfectCoronagraph::compute_image(bool accumulate) {
    return _compute_image(false, accumulate, false);
}