#!/bin/bash
LOCAL_DIR="$(realpath --relative-to=$(pwd) $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd ))"

# install dependencies
pip install -U pip
pip install -r $LOCAL_DIR/requirements.txt || exit 0

# Adds cosmic conan repository if it is not already the case.
conan remote list | grep cosmic || conan remote add cosmic https://api.bintray.com/conan/odp/cosmic

# OPTIONAL: Adds hippo6 conan repository if it is not already the case.
#conan remote list | grep hippo6 || conan remote add hippo6 https://hippo6.obspm.fr/conan False

conan install $LOCAL_DIR --build=missing
conan profile update settings.compiler.libcxx=libstdc++11 default
