// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_perfectCoronagraph.h
//! \ingroup   libsutra
//! \class     SutraPerfectCoronagraph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#ifndef _SUTRA_PERFECT_CORONAGRAPH_H_
#define _SUTRA_PERFECT_CORONAGRAPH_H_

#include <sutra_coronagraph.h>

class SutraPerfectCoronagraph : public SutraCoronagraph {
    public:
        std::vector<CarmaObj<cuFloatComplex>*> AA;
        std::vector<CarmaObj<cuFloatComplex>*> BB;
        std::vector<float> norm;
        std::vector<CarmaObj<cuFloatComplex>*> AA_psf;
        std::vector<CarmaObj<cuFloatComplex>*> BB_psf;
        std::vector<float> norm_psf;
        CarmaObj<cuFloatComplex> *tmp_mft;

    public:
        SutraPerfectCoronagraph(CarmaContext *context, SutraSource *d_source,int im_dimx, 
                                int im_dimy, float *wavelength, int nWavelength, int device);
        ~SutraPerfectCoronagraph()=default;
        int compute_image(bool accumulate);
        int compute_psf(bool accumulate);
        int set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm);
        int set_mft_psf(cuFloatComplex *A, cuFloatComplex *B, float* norm);

    private:
        int _set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm, bool psf);
        int _compute_image(bool psf, bool remove_coro, bool accumulate);
};

#endif //_SUTRA_PERFECT_CORONAGRAPH_H_