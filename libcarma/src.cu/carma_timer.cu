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

//! \file      carma_timer.cu
//! \ingroup   libcarma
//! \brief     this file provides timer CUDA kernels
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#include <carma_timer.hpp>

__global__ void clockKernel(int64_t* clock_counter) {
  *clock_counter = clock64();
}

__global__ void clockKernel(double* time_buffer, int64_t* clock_counter,
                            double gpu_freq) {
  time_buffer[0] = (clock64() - *clock_counter) / gpu_freq;
}

void get_clock_count(int64_t* clock_counter) {
  clockKernel<<<1, 1>>>(clock_counter);
}

void get_clock_count(double* time_buffer, int64_t* clock_counter,
                   double gpu_freq) {
  clockKernel<<<1, 1>>>(time_buffer, clock_counter, gpu_freq);
}
