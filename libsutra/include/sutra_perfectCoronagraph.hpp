// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_perfectCoronagraph.hpp
//! \ingroup   libsutra
//! \class     SutraPerfectCoronagraph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_PERFECT_CORONAGRAPH_H_
#define _SUTRA_PERFECT_CORONAGRAPH_H_

#include <sutra_coronagraph.hpp>
#include <map>

class SutraPerfectCoronagraph : public SutraCoronagraph {
    public:
        std::map<std::string, std::vector<CarmaObj<cuFloatComplex>*>> AA;
        std::map<std::string, std::vector<CarmaObj<cuFloatComplex>*>> BB;
        std::map<std::string, std::vector<float>> norm;

        CarmaObj<cuFloatComplex> *tmp_mft;

    public:
        SutraPerfectCoronagraph(CarmaContext *context, SutraSource *d_source,int32_t im_dimx,
                                int32_t im_dimy, float *wavelength, int32_t nWavelength, int32_t device);
        ~SutraPerfectCoronagraph()=default;
        int32_t compute_image(bool accumulate);
        int32_t compute_psf(bool accumulate);
        int32_t set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm, std::string mftType);

    private:
        int32_t _compute_image(bool psf, bool remove_coro, bool accumulate);
};

#endif //_SUTRA_PERFECT_CORONAGRAPH_H_