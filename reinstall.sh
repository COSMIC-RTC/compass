#!/bin/sh

YORICK_PATH=`which yorick`

if [ $# -gt 0 ]; then
    YORICK_PATH=$1
fi

if [ -z "$YORICK_PATH" ]; then
    echo "yorick is not in the path, use $0 full_yorick_path"
    exit
else
    if [ -x $YORICK_PATH ]; then
	echo "using $YORICK_PATH to update Makefile"
	(cd yoga && $YORICK_PATH -batch make.i) || exit
	(cd yoga_ao && $YORICK_PATH -batch make.i) || exit
	
	(cd libcarma && make clean && make -j) || exit
	(cd yoga && make clean install) || exit
	(cd libsutra && make clean && make -j) || exit
	(cd yoga_ao && make clean install) || exit
    else
	echo "yorick is not in the path, use $0 full_yorick_path"
	exit
    fi
fi

