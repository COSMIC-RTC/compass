#!/bin/bash

YORICK_PATH=$(which yorick)
if [[ $YORICK_PATH=="yorick not found" ]]
then
    if [[ $1 == "" ]]
    then
	echo "yorick is not in the path, use $0 yorick_path"
	exit
    else
	YORICK_PATH=$1
    fi
fi

echo "using $YORICK_PATH to update Makefile"
$YORICK_PATH -batch make.i || exit
(cd yoga_ao && $YORICK_PATH -batch make.i) || exit

(cd libyoga && make clean && make -j) || exit
(make clean install) || exit
(cd yoga_ao/libyoga_ao && make clean && make -j) || exit
(cd yoga_ao && make clean install) || exit
