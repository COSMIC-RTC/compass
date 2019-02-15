
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

cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PYTHON_INSTALL_PATH
make -j8
make install
