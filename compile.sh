
if [ -z $1 ]
then
    PYTHON_INSTALL_PATH=$COMPASS_ROOT/local
else
    PYTHON_INSTALL_PATH=$1
fi

if [ ! -d build ]
then
        echo "Create build directory"
        mkdir -p build
fi

if [ -z $COMPASS_DO_HALF ]
then
    COMPASS_DO_HALF="OFF"
fi

cd build

# BUILD_TOOL="ninja" # make or ninja
if [ "$BUILD_TOOL" = "ninja" ]
then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PYTHON_INSTALL_PATH -Ddo_half=$COMPASS_DO_HALF -GNinja
    ninja install
else
    cmake .. -DCMAKE_INSTALL_PREFIX=$PYTHON_INSTALL_PATH -Ddo_half=$COMPASS_DO_HALF
    make -j8
    make install
fi

if [ ! -z $OCTOPUS_ROOT ]
then
    cp -r $OCTOPUS_ROOT/Octopus $PYTHON_INSTALL_PATH/python
fi
