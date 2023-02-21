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
        bool remove_coro;
        std::vector<CarmaObj<cuFloatComplex>*> AA;
        std::vector<CarmaObj<cuFloatComplex>*> BB;
        std::vector<float> norm;
        CarmaObj<cuFloatComplex> *tmp_mft;

    public:
        ~SutraPerfectCoronagraph()=default;
        int compute_image(bool accumulate);
        int set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm, int dimx, int dimy);
        int set_remove_coro(bool remove);
    protected:
        int propagate();
        SutraPerfectCoronagraph(CarmaContext *context, std::string type, SutraSource *d_source, 
                            int im_dimx, int im_dimy, float *wavelength, int nWavelength, 
                            int device);

}

#endif //_SUTRA_PERFECT_CORONAGRAPH_H_