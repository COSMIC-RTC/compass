#!/bin/bash

if [ -z $1 ]
then
    VCPKG_INSTALL_DIR=$HOME/vcpkg
else
    VCPKG_INSTALL_DIR=$1
fi

git clone https://github.com/Microsoft/vcpkg.git $VCPKG_INSTALL_DIR
$VCPKG_INSTALL_DIR/bootstrap-vcpkg.sh

$VCPKG_INSTALL_DIR/vcpkg integrate install

# Add the vcpkg variables to your .bashrc file
echo "export VCPKG_ROOT=$VCPKG_INSTALL_DIR" >> ~/.bashrc
echo "export PATH=\$VCPKG_ROOT:\$PATH" >> ~/.bashrc
source ~/.bashrc
