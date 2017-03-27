#!/bin/bash

#PARFILE="$SHESHA_ROOT/data/par/MICADO/micado_39m_PYR.py" # scao_40x40_8pix.py  scao_64x64_8pix.py  scao_80x80_8pix.py "
PARFILE="$SHESHA_ROOT/data/par/MICADO/micado_39m_PYR.py" # scao_40x40_8pix.py  scao_64x64_8pix.py  scao_80x80_8pix.py "

DATE=`date +%F_%Hh%M`
OUTPUT="$SHESHA_ROOT/test/scripts/resultatsScripts/script39mPYRLog.txt"
rm $OUTPUT
echo "writing output in "$OUTPUT
echo "To monitor process: tail -f" $OUTPUT
#script="$SHESHA_ROOT/test/scripts/script_PYR39m.py"
script="$SHESHA_ROOT/test/scripts/script_PYR39m_optimGain.py"
# Relevant parameters for pyramid:
# REQ, MODU, GAIN, MAG, NKL
SIMULNAME="PYR_39m_RoundPupil500"

for FREQ in "500"
do
    for RONS in "0.1"
    do
        for MODU in "5" "10" "20"
        do
            for MAG in "13" "14" "15" "16" "17" "18" "19" "20"
            do
               for GAIN in "0.2" "0.4" "0.6" "0.8" "1.0"
                do
                    for KLFILT in "450"
                    do
                        CMD="python $script $PARFILE $FREQ $RONS $MODU $GAIN $MAG $KLFILT $SIMULNAME"
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
