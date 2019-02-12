
if [ -z $1 ]
then
    PYTHON_INSTALL_PATH=python
else
    PYTHON_INSTALL_PATH=$1
fi

if [ ! -d $PYTHON_INSTALL_PATH ]
then
        echo "Create $PYTHON_INSTALL_PATH directory"
        mkdir -p $PYTHON_INSTALL_PATH
fi

easy_install -d $PYTHON_INSTALL_PATH . $OCTOPUS_ROOT
