#!/bin/bash

PARFILE="$SHESHA_ROOT/data/par/MICADO/micado_39m_PYR_ELTPupil.py"
OUTPUT="$SHESHA_ROOT/test/scripts/resultatsScripts/script39mPYRLog3.txt"
rm $OUTPUT
echo "writing output in "$OUTPUT
echo "To monitor process: tail -f" $OUTPUT
#script="$SHESHA_ROOT/test/scripts/script_PYR39m.py"
script="$SHESHA_ROOT/test/scripts/script_PYR39m_optimGain.py"
# Relevant parameters for pyramid:
# REQ, MODU, GAIN, MAG, NKL
SIMULNAME="PYR_39m_ELTPupil_TestNbKL"
GPU="5"
for FREQ in "500"
do
    for RONS in "0.1"
    do
        for MODU in "5"
        do
           for MAG in "11" "12" "13" "14" "15" "16" "17" "18"
           #for MAG in "18"
            do
               for GAIN in "0.6"
               #for GAIN in "1.0"
                do
                    for KLFILT in "2000" "1600" "1000" "500" "2"
                    do
                        CMD="python $script $PARFILE $FREQ $RONS $MODU $GAIN $MAG $KLFILT $SIMULNAME $GPU"
                        echo "execute $CMD" >> $OUTPUT
                        $CMD 2>> $OUTPUT >> $OUTPUT
                    done
                done
            done
        done
    done
done

# To monitor the script log:
# tail -f resultatsScripts/script39mPYRLog.txt


#for f in $FILES
#do
#    for CTR in "ls" "modopti" "mv" "geo"
#    do
#        for COG in "cog" "tcog" "bpcog" "geom" #"pyr" #
#        do
#            CMD="python -i $script $f $COG $CTR $DEVICE"
#            echo "execute $CMD" >> $OUTPUT
#            $CMD 2>> $OUTPUT >> $OUTPUT
#        done
#    done
#done
echo "Script 39mPYR Done"

#FILES_LGS="scao_16x16_8pix_lgs.py"
#FILES_LGS+="scao_40x40_10pix_lgs.par"
#FILES_LGS+="scao_64x64_16pix_lgs.par"
#FILES_LGS+="scao_80x80_20pix_lgs.par"
