#!/bin/bash

rm -rf ~/.conan

# install dependencies
pip install "conan<1.48" cmake --upgrade

# Adds cosmic conan repository if it is not already the case.
conan remote list | grep cosmic || conan remote add cosmic https://odp2.jfrog.io/artifactory/api/conan/cosmic
conan remote list | grep obspm || conan remote add obspm https://conan.obspm.fr/conan

# OPTIONAL: Adds hippo6 conan repository if it is not already the case.
#conan remote list | grep hippo6 || conan remote add hippo6 https://hippo6.obspm.fr/conan False

conan profile new default --detect --force > /dev/null
conan profile update settings.compiler=clang default
conan profile update settings.compiler.version=13 default
conan profile update settings.compiler.libcxx=libstdc++11 default
conan profile update env.CC=/usr/bin/clang default
conan profile update env.CXX=/usr/bin/clang++ default

