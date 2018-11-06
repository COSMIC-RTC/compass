#!/bin/bash
#
# Prints the compute capability of the first CUDA device installed
# on the system, or alternatively the device whose index is the
# first command-line argument

timestamp=$(date +%s.%N)
gcc_binary=$(which g++)
gcc_binary=${gcc_binary:-g++}
cuda_root=${CUDA_ROOT:-/usr/local/cuda}
CUDA_INCLUDE_DIRS=${CUDA_INCLUDE_DIRS:-${cuda_root}/include}
CUDA_CUDART_LIBRARY=${CUDA_CUDART_LIBRARY:-${cuda_root}/lib64/libcudart.so}
generated_binary="/tmp/cuda-compute-version-helper-$$-$timestamp"
# create a 'here document' that is code we compile and use to probe the card
source_code="$(cat << EOF
#include <stdio.h>
#include <set>
#include <cuda_runtime_api.h>

int main()
{
    cudaDeviceProp prop;
    cudaError_t status;
    int device_count;
    status = cudaGetDeviceCount(&device_count);
    if (status != cudaSuccess) {
            fprintf(stderr,"cudaGetDeviceCount() failed: %s\n", cudaGetErrorString(status));
            return -1;
    }
    if (device_count == 0) {
            fprintf(stderr, "No cuda device available\n");
            return -1;
    }
    std::set<int> me_set;
    for(int device_index = 0; device_index < device_count; ++device_index)
    {
        status = cudaGetDeviceProperties(&prop, device_index);
        if (status != cudaSuccess) {
                fprintf(stderr,"cudaGetDeviceProperties() for device ${device_index} failed: %s\n", cudaGetErrorString(status));
                return -1;
        }
        me_set.emplace(prop.major * 10 + prop.minor);
    }
    for(int e : me_set)
        printf("%d ", e);
    return 0;
}
EOF
)"
echo "$source_code" | $gcc_binary -x c++ -I"$CUDA_INCLUDE_DIRS" -o "$generated_binary" - -x none "$CUDA_CUDART_LIBRARY"

# probe the card and cleanup

CUDA_SM=$($generated_binary)
rm $generated_binary
echo -n $CUDA_SM
exit 0
