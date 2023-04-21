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
//! \version   5.4.3
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
typedef tuple<mftVec, vector<long>> mftTuple;

class SutraStellarCoronagraph : public SutraCoronagraph {
    public:
        bool babinet;
        int fpmDimx;
        int fpmDimy;
        map<string, mftTuple> AA;
        map<string, mftTuple> BB;
        map<string, vector<float>> norm;
        map<string, tuple<CarmaObj<cuFloatComplex>*, vector<long>>> tmp_mft;

        CarmaObj<float> *d_lyot_stop;
        CarmaObj<float> *d_apodizer;
        vector<CarmaObj<float>*> focal_plane_mask;
        CarmaObj<cuFloatComplex> *d_electric_field_fpm;
        CarmaObj<cuFloatComplex> *d_electric_field_babinet;

    public:
        SutraStellarCoronagraph(CarmaContext *context, SutraSource *d_source,int im_dimx, 
                                int im_dimy, int fpm_dimx, int fpm_dimy,
                                float *wavelength, int nWavelength, bool babinet, int device);
        ~SutraStellarCoronagraph()=default;
        int compute_image(bool accumulate);
        int compute_psf(bool accumulate);
        int compute_image_normalization();
        int set_mft(cuFloatComplex *A, cuFloatComplex *B, float* norm, std::string mftType);
        int set_focal_plane_mask(float *mask);
        int set_apodizer(float *apodizer);
        int set_lyot_stop(float *lyot_stop);

    private:
        int _compute_image(bool center_on_pixel, bool accumulate, bool no_fpm);
};

#endif //_SUTRA_STELLAR_CORONAGRAPH_H_