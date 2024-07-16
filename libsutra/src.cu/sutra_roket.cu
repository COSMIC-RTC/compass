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

//! \file      sutra_roket.cu
//! \ingroup   libsutra
//! \class     SutraRoket
//! \brief     this class provides the roket features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_roket.hpp>

__global__ void separate_modes_krnl(float *modes, float *filtmodes, int32_t nmodes,
                                    int32_t nfilt) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t bornemin = nmodes - nfilt - 2;
  int32_t bornemax = nmodes - 2;
  while (tid < nmodes) {
    if (tid >= bornemin && tid < bornemax) {
      filtmodes[tid] = modes[tid];
      modes[tid] = 0.0f;
    } else
      filtmodes[tid] = 0.0f;

    tid += blockDim.x * gridDim.x;
  }
}

int32_t separate_modes(float *modes, float *filtmodes, int32_t nmodes, int32_t nfilt,
                   CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nmodes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  separate_modes_krnl<<<grid, threads>>>(modes, filtmodes, nmodes, nfilt);
  carma_check_msg("separate_modes_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
