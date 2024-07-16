// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_stellarCoronagraph.hpp
//! \ingroup   libsutra
//! \class     SutraStellarCoronagraph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_STELLAR_CORONAGRAPH_H_
#define _SUTRA_STELLAR_CORONAGRAPH_H_

#include <sutra_coronagraph.hpp>
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