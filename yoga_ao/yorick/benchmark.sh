#!/bin/bash

FILES="scao_16x16_8pix.par
scao_16x16_16pix.par 
scao_40x40_8pix.par 
scao_40x40_16pix.par 
scao_64x64_8pix.par
scao_64x64_16pix.par
scao_80x80_8pix.par 
scao_80x80_16pix.par"

for f in $FILES
do
    for CTR in "ls" "modopti" "mv"
    do
	for COG in "cog" "tcog" "bpcog" "geom"
        do
	        CMD="/home/ferreira/yorick-2.2/relocate/bin/yorick -batch benchmark_script.i "$f" "$COG" "$CTR
	        echo "execute $CMD"
	        $CMD
        done
    done
done
