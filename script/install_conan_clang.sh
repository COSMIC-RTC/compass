#!/bin/bash

rm -rf ~/.conan2

# install dependencies
pip install -U -r script/requirements.txt
pip install -U -r script/requirements-dev.txt

# Adds obspm conan repository if it is not already the case.
conan remote list | grep obspm || conan remote add obspm https://conan.obspm.fr/conan

CC=/usr/bin/clang CXX=/usr/bin/clang++ conan profile detect > /dev/null

echo "[buildenv]" >> ~/.conan2/profiles/default
echo "CC=/usr/bin/clang" >> ~/.conan2/profiles/default
echo "CXX=/usr/bin/clang++" >> ~/.conan2/profiles/default

echo "tools.build:skip_test=True" >> ~/.conan2/global.conf 