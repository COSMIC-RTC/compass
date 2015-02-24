#!/bin/bash

FILES="scao_16x16_8pix.par
scao_16x16_16pix.par 
scao_40x40_8pix.par 
scao_40x40_16pix.par 
scao_64x64_8pix.par
scao_64x64_16pix.par
scao_80x80_8pix.par 
scao_80x80_16pix.par"
#scao_16x16_8pix_noisy.par
#scao_16x16_16pix_noisy.par 
#scao_40x40_8pix_noisy.par 
#scao_40x40_16pix_noisy.par 
#scao_64x64_8pix_noisy.par
#scao_64x64_16pix_noisy.par
#scao_80x80_8pix_noisy.par 
#scao_80x80_16pix_noisy.par

YORICK_PATH=`which yorick`
#YORICK_PATH="/home/ferreira/yorick-2.2/relocate/bin/yorick"


if [ $# -gt 0 ]; then
    YORICK_PATH=$1
fi

DATE=`date +%F_%Hh%M`

if [ -z "$YORICK_PATH" ]; then
    echo "yorick is not in the path, use $0 full_yorick_path"
    exit
else
    for f in $FILES
    do
        for CTR in "ls" "modopti" "mv" "geo"
        do
	        for COG in "cog" "tcog" "bpcog" "geom"
            do
	            CMD=$YORICK_PATH" -batch benchmark_script.i "$f" "$COG" "$CTR
	            echo "execute $CMD" >> $DATE\_$HOSTNAME\_outputfile
	            $CMD 2>> $DATE\_$HOSTNAME\_outputfile >> $DATE\_$HOSTNAME\_outputfile
            done
        done
    done
fi

FILES="scao_16x16_8pix_lgs.par
scao_40x40_10pix_lgs.par  
scao_64x64_16pix_lgs.par
scao_80x80_20pix_lgs.par"

if [ $# -gt 0 ]; then
    YORICK_PATH=$1
fi

if [ -z "$YORICK_PATH" ]; then
    echo "yorick is not in the path, use $0 full_yorick_path"
    exit
else
    for f in $FILES
    do
        for CTR in "ls" "modopti" "mv" "geo"
        do
	        for COG in "wcog" "corr"
            do
	            CMD=$YORICK_PATH" -batch benchmark_script.i "$f" "$COG" "$CTR
	            echo "execute $CMD" >> $DATE\_$HOSTNAME\_outputfile
	            $CMD 2>> $DATE\_$HOSTNAME\_outputfile >> $DATE\_$HOSTNAME\_outputfile
            done
        done
    done
fi
