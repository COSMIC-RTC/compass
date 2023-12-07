// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_coronagraph.h
//! \ingroup   libsutra
//! \class     SutraCoronagraph
//! \brief     this class provides the coronagraph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CORONAGRAPH_H_
#define _SUTRA_CORONAGRAPH_H_

#include <carma_cublas.h>
#include <sutra_utils.h>
#include <sutra_source.h>
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
