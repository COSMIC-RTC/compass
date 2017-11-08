#include <carma_utils.h>
#include <cuda_profiler_api.h>
#include <carma_context.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include <sys/ioctl.h>

int quick_pow10(int n) {
  static int pow10[10] = {
    1, 10, 100, 1000, 10000,
    100000, 1000000, 10000000, 100000000, 1000000000
  };

  return pow10[n];
}

std::string disp(int idx, int width, char align='r') {
  std::ostringstream stm;
  if(idx >= quick_pow10(width)) {
    for( int i(0); i<width-1; i++) {
      stm << "#";
    }
    return stm.str();
  }
  if(align == 'l') {
    stm << idx;
  }
  while(idx < quick_pow10(width)-1) {
    stm << " ";
    width--;
  }
  if(align != 'l') {
    stm << idx;
  }
  return stm.str();
}

carma_utils::ProgressBar::ProgressBar(int max_): max(max_), start(std::chrono::system_clock::now()) {
  ndigits = 0;
  while(max_ > quick_pow10(ndigits)) ndigits++;
  struct winsize w;
  ioctl(0, TIOCGWINSZ, &w);
  barWidth = w.ws_col - 2*(ndigits) - 24 ;

};

void carma_utils::ProgressBar::update() {
  count++;
  progress = (double)(count) / max;
  if(2*barWidth*progress > prev) {
    prev++;
    int n = (int)(progress*100.0);
    std::cout << disp(n, 3) << "% [";
    int pos = barWidth * progress;
    for(int i=0; i<barWidth; i++) {
      if(i<pos) std::cout << "=";
      else if(i==pos) std::cout << ">";
      else std::cout << " ";
    }
    std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = time - start;
    int eta = elapsed.count()>0?(int)(elapsed.count() / progress - elapsed.count()):0;
    std::cout << "] " << disp(count, ndigits) << "/" << disp(max, ndigits, 'l') << " ETA " << disp(eta,4) << "s \r" << std::flush;
  }
};

void carma_utils::ProgressBar::finish() {
  std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = time - start;
  std::cout << std::endl << "done in " << elapsed.count() << " s" << std::endl;
};

void getNumBlocksAndThreads(carma_device *device, int n, int &blocks, int &threads) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  blocks = device->get_properties().multiProcessorCount * 8; //device->get_cores_per_sm();
  threads = (n + blocks - 1) / blocks;

  if (threads > maxThreads) {
    threads = maxThreads;
    blocks = MIN( device->get_properties().maxGridSize[0],(n + threads - 1) / threads);
  }
  /*
  threads = device->get_cores_per_sm();
  blocks = (n + threads - 1) / threads;
  if (blocks > device->get_properties().maxGridSize[0]) {
    blocks=device->get_properties().maxGridSize[0];
    threads = (n + blocks - 1) / blocks;
  }
  */
}

void sumGetNumBlocksAndThreads(int n, carma_device *device, int &blocks, int &threads) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int maxBlocks = device->get_properties().multiProcessorCount;

  threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
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

  //srand(time(NULL)); // initialisation de rand
  //return rand()/norm;
  return random() / norm;
}

extern "C" {
  int _dist(float *d, long dimx, long dimy, float xc, float yc) {
    /* Declarations */
    long i, j;

    /* Loop and fill d with distance values */
    for (i = 0; i < dimx; ++i) {
      for (j = 0; j < dimy; ++j) {
        d[i + j * dimx] = (float) sqrt(
                            (xc - (float) i) * (xc - (float) i)
                            + (yc - (float) j) * (yc - (float) j));
      }
    }
    return EXIT_SUCCESS;
  }

  void _poidev(float *xmv, long n)
  /* */
  {
    float
    gammln(float xx);
    /*  float ran1(long *idum);*/
    static float sq, alxm, g, oldm = (-1.0);
    float xm, em, t, y;
    long i;

    for (i = 0; i < n; i++) {
      xm = (float) xmv[i];
      if (xm == 0.0f)
        continue;
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
            y = tan(3.1415926535897932384626433832 * ran1());
            em = sq * y + xm;
          } while (em < 0.0);
          em = floor(em);
          t = 0.9 * (1.0 + y * y) * exp(em * alxm - gammln(em + 1.0) - g);
        } while (ran1() > t);
      }
      xmv[i] = (float) em;
    }
  }

  float gammln(float xx) {
    /* Returns the value ln[?(xx)] for xx>0. */
    float x, y, tmp, ser;
    static float cof[6] = { 76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
                            -0.5395239384953e-5
                          };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
      ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
  }

}

int printMemInfo() {
  size_t free_mem;
  size_t total_mem;
  float free_float;
  float total_float;
  float used_mem;

  carmaSafeCall(cudaMemGetInfo(&free_mem, &total_mem));
  free_float = (float) free_mem / 1000000.;
  total_float = (float) total_mem / 1000000.;
  used_mem = total_float - free_float;
  printf("GPU Memory usage : used memory = %f MB, free memory = %f MB\n",
         used_mem, free_float);

  return EXIT_SUCCESS;
}
