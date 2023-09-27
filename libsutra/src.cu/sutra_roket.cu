// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_roket.cu
//! \ingroup   libsutra
//! \class     SutraRoket
//! \brief     this class provides the roket features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_roket.h>

__global__ void separate_modes_krnl(float *modes, float *filtmodes, int nmodes,
                                    int nfilt) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int bornemin = nmodes - nfilt - 2;
  int bornemax = nmodes - 2;
  while (tid < nmodes) {
    if (tid >= bornemin && tid < bornemax) {
      filtmodes[tid] = modes[tid];
      modes[tid] = 0.0f;
    } else
      filtmodes[tid] = 0.0f;

    tid += blockDim.x * gridDim.x;
  }
}

int separate_modes(float *modes, float *filtmodes, int nmodes, int nfilt,
                   CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nmodes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  separate_modes_krnl<<<grid, threads>>>(modes, filtmodes, nmodes, nfilt);
  carma_check_msg("separate_modes_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
