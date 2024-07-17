// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_perfectCoronagraph.hpp
//! \ingroup   libsutra
//! \class     SutraPerfectCoronagraph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
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