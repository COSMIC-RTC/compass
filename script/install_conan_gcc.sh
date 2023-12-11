#!/bin/bash

rm -rf ~/.conan2

# install dependencies
pip install -U -r script/requirements.txt
pip install -U -r script/requirements-dev.txt

# Adds obspm conan repository if it is not already the case.
conan remote list | grep obspm || conan remote add obspm https://conan.obspm.fr/conan

conan profile detect > /dev/null

# for Arch: 10
# for Ubuntu 18.04: 7.5
if [ ! -z $1 ]
then
    conan profile update settings.compiler.version=$1 default

    # for Arch: CUDA_ROOT=/opt/cuda
    # for Ubuntu 18.04: CUDA_ROOT=/usr/local/cuda
    if [ ! -z $CUDA_ROOT ]
    then
        conan profile update env.CC=$CUDA_ROOT/bin/gcc default
        conan profile update env.CXX=$CUDA_ROOT/bin/g++ default
    fi
fi
echo "tools.build:skip_test=True" >> ~/.conan2/global.conf 
