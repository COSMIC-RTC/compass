#!/bin/bash

FILES="scao_16x16_8pix.par
scao_16x16_16pix.par 
scao_40x40_8pix.par 
scao_40x40_16pix.par 
scao_64x64_8pix.par
scao_64x64_16pix.par
scao_80x80_8pix.par 
scao_80x80_16pix.par"

YORICK_PATH=`which yorick`

if [ $# -gt 0 ]; then
    YORICK_PATH=$1
fi

if [ -z "$YORICK_PATH" ]; then
    echo "yorick is not in the path, use $0 full_yorick_path"
    exit
else
    for f in $FILES
    do
        for CTR in "ls" "modopti" "mv"
        do
	        for COG in "cog" "tcog" "bpcog" "geom"
            do
	            CMD=$YORICK_PATH" -batch benchmark_script.i "$f" "$COG" "$CTR
	            echo "execute $CMD"
	            $CMD
            done
        done
    done
fi
