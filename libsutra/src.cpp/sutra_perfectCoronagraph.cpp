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
                                    int dimx, int dimy, float *wavelength, int nWavelength, 
                                    int device):
    SutraCoronagraph(context, "perfect", d_source, dimx, dimy, wavelength, nWavelength, device),
    focalPlaneDimx(0),
    focalPlaneDimy(0),
    remove_coro(false) {}

int SutraPerfectCoronagraph::set_remove_coro(bool remove) {
    remove_coro = remove;
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm0) {
    if(!AA.empty()) {
        AA.clear();
    }
    if(!BB.empty()) {
        BB.clear();
    }

    for (int i = 0; i < wavelength.size() ; i++) {
        long dims[3] = {2, imageDimx, pupDimx};
        AA.push_back(new CarmaObj<cuFloatComplex>(current_context, dims, A + i * imageDimx * pupDimx));
        dims[1] = pupDimy;
        dims[2] = imageDimy;
        BB.push_back(new CarmaObj<cuFloatComplex>(current_context, dims, B + i * imageDimy * pupDimy));
        norm.push_back(norm0[i]);
    }
    dims[1] = imageDimx;
    dims[2] = pupDimy;
    tmp_mft = new CarmaObj<cuFloatComplex>(current_context, dims);
    return EXIT_SUCCESS;
}

int SutraPerfectCoronagraph::compute_image() {
    d_image_se->reset();
    for (int i = 0 ; i < wavelength.size(); i++) {
        compute_electric_field(i);
        if(!remove_coro)
            remove_complex_avg(d_electric_field->get_data(), d_electric_field->sum(),
                            d_pupil->get_data(), d_source->d_wherephase->get_nb_elements(),
                            pupDimx, pupDimy, current_context->get_device(device));
        mft(AA[i], BB[i], tmp_mft, d_electric_field, d_complex_image, norm[i]);
        accumulate_abs2(d_complex_image->get_data(), d_image_se->get_data(), 
                        d_image_se->get_nb_elements(), current_context->get_device(device));
    }
}