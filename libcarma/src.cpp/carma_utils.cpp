// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_utils.cpp
//! \ingroup   libcarma
//! \brief     this file provides tools to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_context.hpp>
#include <carma_utils.hpp>
#include <cuda_profiler_api.h>
#include <cmath>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <sstream>
#include <string>

#include <sys/ioctl.h>


void get_num_blocks_and_threads(CarmaDevice *device, int32_t n, int32_t &blocks,
                            int32_t &threads) {
  // int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  // blocks = device->get_properties().multiProcessorCount *
  //          8;  // device->get_cores_per_sm();
  // threads = (n + blocks - 1) / blocks;

  // if (threads > maxThreads) {
  //   threads = maxThreads;
  //   blocks = MIN(device->get_properties().maxGridSize[0],
  //                (n + threads - 1) / threads);
  // }

  threads = device->get_properties().maxThreadsPerBlock;
  blocks = (n + threads - 1) / threads;
  if (blocks > device->get_properties().maxGridSize[0]) {
    throw std::runtime_error("problem is too big!");
  }

}

void sum_get_num_blocks_and_threads(int32_t n, CarmaDevice *device, int32_t &blocks,
                               int32_t &threads) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  int32_t maxBlocks = device->get_properties().multiProcessorCount;

  threads = (n < maxThreads * 2) ? next_pow2((n + 1) / 2) : maxThreads;
  blocks = (n + (threads * 2 - 1)) / (threads * 2);
  blocks = MIN(maxBlocks, blocks);
}

void carma_start_profile() {
  printf("CUDA Profiling started\n");
  cudaProfilerStart();
}
void carma_stop_profile() {
  printf("CUDA Profiling stoped\n");
  cudaProfilerStop();
}

float ran1() {
  float norm;
  norm = 2147483647.f;

  // srand(time(NULL)); // initialisation de rand
  // return rand()/norm;
  return random() / norm;
}

extern "C" {
int32_t _dist(float *d, int64_t dimx, int64_t dimy, float xc, float yc) {
  /* Declarations */
  int64_t i, j;

  /* Loop and fill d with distance values */
  for (i = 0; i < dimx; ++i) {
    for (j = 0; j < dimy; ++j) {
      d[i + j * dimx] = (float)sqrt((xc - (float)i) * (xc - (float)i) +
                                    (yc - (float)j) * (yc - (float)j));
    }
  }
  return EXIT_SUCCESS;
}

void _poidev(float *xmv, int64_t n)
/* */
{
  float gammln(float xx);
  /*  float ran1(int64_t *idum);*/
  static float sq, alxm, g, oldm = (-1.0);
  float xm, em, t, y;
  int64_t i;

  for (i = 0; i < n; i++) {
    xm = (float)xmv[i];
    if (xm == 0.0f) continue;
    if (xm < 20.0) { /* Use direct method. */
      if (xm != oldm) {
        oldm = xm;
        g = exp(-xm); /* If xm is new, compute the exponential. */
      }
      em = -1;
      t = 1.0;
      do {
        ++em;
        t *= ran1();
      } while (t > g);
    } else { /* Use rejection method. */
      if (xm != oldm) {
        oldm = xm;
        sq = sqrt(2.0 * xm);
        alxm = log(xm);
        g = xm * alxm - gammln(xm + 1.0);
      }
      do {
        do {
          y = tan(CARMA_PI * ran1());
          em = sq * y + xm;
        } while (em < 0.0);
        em = floor(em);
        t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
      } while (ran1() > t);
    }
    xmv[i] = (float)em;
  }
}

float gammln(float xx) {
  /* Returns the value ln[?(xx)] for xx>0. */
  float x, y, tmp, ser;
  static float cof[6] = {76.18009172947146,     -86.50532032941677,
                         24.01409824083091,     -1.231739572450155,
                         0.1208650973866179e-2, -0.5395239384953e-5};
  int32_t j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}
}

int32_t print_mem_info() {
  size_t free_mem;
  size_t total_mem;
  float free_float;
  float total_float;
  float used_mem;

  carma_safe_call(cudaMemGetInfo(&free_mem, &total_mem));
  free_float = (float)free_mem / 1000000.;
  total_float = (float)total_mem / 1000000.;
  used_mem = total_float - free_float;
  printf("GPU Memory usage : used memory = %f MB, free memory = %f MB\n",
         used_mem, free_float);

  return EXIT_SUCCESS;
}
