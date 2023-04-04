// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_timer.cu
//! \ingroup   libcarma
//! \brief     this file provides timer CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24


#include <carma_timer.h>

__global__ void clockKernel(long long int* clock_counter) {
  *clock_counter = clock64();
}

__global__ void clockKernel(double* time_buffer, long long int* clock_counter,
                            double gpu_freq) {
  time_buffer[0] = (clock64() - *clock_counter) / gpu_freq;
}

void get_clock_count(long long int* clock_counter) {
  clockKernel<<<1, 1>>>(clock_counter);
}

void get_clock_count(double* time_buffer, long long int* clock_counter,
                   double gpu_freq) {
  clockKernel<<<1, 1>>>(time_buffer, clock_counter, gpu_freq);
}
