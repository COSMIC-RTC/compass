#!/bin/env bash

CC=/opt/cuda/bin/gcc CXX=/opt/cuda/bin/g++ cmake ..

# CC=/opt/cuda/bin/gcc CXX=/opt/cuda/bin/g++ cmake -GNinja ..
# time ninja
# real	0m31,448s
# user	3m19,223s
# sys	0m14,071s

# CC=/opt/cuda/bin/gcc CXX=/opt/cuda/bin/g++ cmake ..
# time make -j8
# real	0m30,622s
# user	2m54,184s
# sys	0m12,232s
