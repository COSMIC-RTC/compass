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

//! \file      sutra_coronagraph.hpp
//! \ingroup   libsutra
//! \class     SutraCoronagraph
//! \brief     this class provides the coronagraph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CORONAGRAPH_H_
#define _SUTRA_CORONAGRAPH_H_

#include <carma_cublas.hpp>
#include <sutra_utils.hpp>
#include <sutra_source.hpp>
#include <tuple>
#include <vector>

class SutraCoronagraph {
    public:
        int32_t device;
        std::string type;
        int64_t cntPsf;
        int64_t cntImg;
        int32_t imageDimx;
        int32_t imageDimy;
        int32_t pupDimx;
        int32_t pupDimy;
        std::vector<float> wavelength;
        CarmaContext *current_context;
        CarmaObj<float> *d_image_se;
        CarmaObj<float> *d_image_le;
        CarmaObj<float> *d_psf_se;
        CarmaObj<float> *d_psf_le;
        std::vector<CarmaObj<float>*> amplitude;

        CarmaObj<cuFloatComplex> *d_electric_field;
        CarmaObj<cuFloatComplex> *d_complex_image;
        CarmaObj<float> *d_pupil;

        SutraSource* d_source;

    public:
        virtual ~SutraCoronagraph()=default;
        virtual int32_t compute_image(bool accumulate) = 0;
        virtual int32_t compute_psf(bool accumulate) = 0;
        int32_t reset();
        int32_t compute_electric_field(int32_t wavelengthIndex);
        int32_t set_amplitude(float* amplitude);

    protected:
        SutraCoronagraph(CarmaContext *context, std::string type, SutraSource *d_source,
                            int32_t dimx, int32_t dimy,float *wavelength, int32_t nWavelength, int32_t device);
        int32_t mft(CarmaObj<cuFloatComplex> *A, CarmaObj<cuFloatComplex> *B,
                CarmaObj<cuFloatComplex> *Ainput,
                CarmaObj<cuFloatComplex> *input, CarmaObj<cuFloatComplex> *output, float norm);

};

int32_t compute_electric_field(cuFloatComplex *electric_field, float* phase_opd, float scale,
                            float* amplitude, float* mask, int32_t dimx, int32_t dimy, CarmaDevice *device);
int32_t remove_complex_avg(cuFloatComplex *electric_field, cuFloatComplex sum, float* mask, int32_t Nvalid,
                        int32_t dimx, int32_t dimy, CarmaDevice *device);
int32_t accumulate_abs2(cuFloatComplex *img, float* abs2img, int32_t N, CarmaDevice *device);
int32_t apply_mask(cuFloatComplex *electric_field, float* mask, int32_t N, CarmaDevice *device);
#endif //_SUTRA_CORONAGRAPH_H_
