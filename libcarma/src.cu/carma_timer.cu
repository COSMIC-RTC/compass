#include <carma_timer.h>

__global__ void clockKernel(long long int* clockCounter) {
  *clockCounter = clock64();
}

__global__ void clockKernel(double* timeBuffer, long long int* clockCounter,
                            double GPUfreq) {
  timeBuffer[0] = (clock64() - *clockCounter) / GPUfreq;
}

void getClockCount(long long int* clockCounter) {
  clockKernel<<<1, 1>>>(clockCounter);
}

void getClockCount(double* timeBuffer, long long int* clockCounter,
                   double GPUfreq) {
  clockKernel<<<1, 1>>>(timeBuffer, clockCounter, GPUfreq);
}
