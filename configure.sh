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
	(cd yoga && cp Makefile.in Makefile && $YORICK_PATH -batch make.i) || exit
	(cd yoga_ao && cp Makefile.in Makefile && $YORICK_PATH -batch make.i) || exit
    else
	echo "yorick is not in the path, use $0 full_yorick_path"
	exit
    fi
fi

