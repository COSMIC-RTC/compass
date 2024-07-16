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

//! \file      sutra_centroider_utils.cuh
//! \ingroup   libsutra
//! \class     SlopesIndex
//! \brief     this struct allows easy manipulation of x-y slopes indices
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2016/01/24

#ifndef SUTRA_CENTROIDER_UTILS_CUH
#define SUTRA_CENTROIDER_UTILS_CUH

#include <sutra_centroider.hpp>

namespace sutra
{

    struct SlopesIndex {
        int32_t factor, offset;

        __host__ __device__ SlopesIndex(int32_t nvalid, SlopeOrder so):
            factor((so == SlopeOrder::untied) ? 1 : 2),
            offset((so == SlopeOrder::untied) ? nvalid : 1)
        {}

        __host__ __device__
        constexpr int32_t x(std::size_t pos) const noexcept
        {
            return factor * pos;
        }

        __host__ __device__
        constexpr int32_t y(std::size_t pos) const noexcept
        {
            return factor * pos + offset;
        }

    };


} // namespace sutra

#endif //SUTRA_CENTROIDER_UTILS_CUH