#!/bin/bash

rm -rf ~/.conan

# install dependencies
conda install --file script/requirements-conda.txt -y
pip install -r script/requirements-dev.txt --upgrade

# Adds obspm conan repository if it is not already the case.
conan remote list | grep obspm || conan remote add obspm https://conan.obspm.fr/conan False

# OPTIONAL: Adds hippo6 conan repository if it is not already the case.
#conan remote list | grep hippo6 || conan remote add hippo6 https://hippo6.obspm.fr/conan False

conan profile new default --detect --force > /dev/null
conan profile update settings.compiler=clang default
conan profile update settings.compiler.version=13 default
conan profile update settings.compiler.libcxx=libstdc++11 default
conan profile update env.CC=/usr/bin/clang default
conan profile update env.CXX=/usr/bin/clang++ default

