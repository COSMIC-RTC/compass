// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_stellarCoronagraph.h
//! \ingroup   libsutra
//! \class     SutraStellarCoronagraph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_STELLAR_CORONAGRAPH_H_
#define _SUTRA_STELLAR_CORONAGRAPH_H_

#include <sutra_coronagraph.h>
#include <map>
#include <tuple>

using std::vector;
using std::map;
using std::string;
using std::tuple;

typedef vector<CarmaObj<cuFloatComplex>*> mftVec;
typedef tuple<mftVec, vector<int64_t>> mftTuple;

class SutraStellarCoronagraph : public SutraCoronagraph {
    public:
        bool babinet;
        int32_t fpmDimx;
        int32_t fpmDimy;
        map<string, mftTuple> AA;
        map<string, mftTuple> BB;
        map<string, vector<float>> norm;
        map<string, tuple<CarmaObj<cuFloatComplex>*, vector<int64_t>>> tmp_mft;

        CarmaObj<float> *d_lyot_stop;
        CarmaObj<float> *d_apodizer;
        vector<CarmaObj<float>*> focal_plane_mask;
        CarmaObj<cuFloatComplex> *d_electric_field_fpm;
        CarmaObj<cuFloatComplex> *d_electric_field_babinet;

    public:
        SutraStellarCoronagraph(CarmaContext *context, SutraSource *d_source,int32_t im_dimx,
                                int32_t im_dimy, int32_t fpm_dimx, int32_t fpm_dimy,
                                float *wavelength, int32_t nWavelength, bool babinet, int32_t device);
        ~SutraStellarCoronagraph()=default;
        int32_t compute_image(bool accumulate);
        int32_t compute_psf(bool accumulate);
        int32_t compute_image_normalization();
        int32_t set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm, std::string mftType);
        int32_t set_focal_plane_mask(float *mask);
        int32_t set_apodizer(float *apodizer);
        int32_t set_lyot_stop(float *lyot_stop);

    private:
        int32_t _compute_image(bool center_on_pixel, bool accumulate, bool no_fpm);
};

#endif //_SUTRA_STELLAR_CORONAGRAPH_H_