#!/bin/bash

FILES="scao_16x16_8pix.py  scao_40x40_8pix.py  scao_64x64_8pix.py  scao_80x80_8pix.py "
FILES+="scao_16x16_16pix.py scao_40x40_16pix.py scao_64x64_16pix.py scao_80x80_16pix.py"
#FILES+="scao_16x16_8pix_noisy.par scao_40x40_8pix_noisy.par scao_64x64_8pix_noisy.par scao_80x80_8pix_noisy.par 
#FILES+="scao_16x16_16pix_noisy.par scao_40x40_16pix_noisy.par scao_64x64_16pix_noisy.par scao_80x80_16pix_noisy.par"

YORICK_PATH=`which yorick`

if [ $# -gt 0 ]; then
    YORICK_PATH=$1
fi

DATE=`date +%F_%Hh%M`
SVN=`svnversion`
OUTPUT="$CHAKRA_AO/data/bench-results/outputfile_$DATE\_$HOSTNAME\_r$SVN"

script="$CHAKRA_AO/test/benchmark_script.py"

#if [ -z "$YORICK_PATH" ]; then
#    echo "yorick is not in the path, use $0 full_yorick_path"
#    exit
#else
    for f in $FILES
    do
        for CTR in "ls" "modopti" "mv" "geo"
        do
	        for COG in "cog" "tcog" "bpcog" "geom"
            do
	            CMD="python $script $f $COG $CTR $1"
	            echo "execute $CMD" >> $OUTPUT
	            $CMD 2>> $OUTPUT >> $OUTPUT
            done
        done
    done
#fi

FILES_LGS="scao_16x16_8pix_lgs.py"
FILES_LGS+="scao_40x40_10pix_lgs.par" 
FILES_LGS+="scao_64x64_16pix_lgs.par"
FILES_LGS+="scao_80x80_20pix_lgs.par"

for f in $FILES_LGS
do
    for CTR in "ls" "modopti" "mv" "geo"
    do
        for COG in "wcog" "corr"
        do
            CMD="python $script $f $COG $CTR $1"
            echo "execute $CMD" >> $OUTPUT
            $CMD 2>> $OUTPUT >> $OUTPUT
        done
    done
done
