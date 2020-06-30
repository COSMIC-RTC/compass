#!/bin/bash
LOCAL_DIR="$(realpath --relative-to=$(pwd) $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd ))"

# install dependencies
pip install -r $LOCAL_DIR/requirements.txt || exit 0

# Adds cosmic conan repository if it is not already the case.
conan remote list | grep cosmic || conan remote add cosmic https://api.bintray.com/conan/odp/cosmic