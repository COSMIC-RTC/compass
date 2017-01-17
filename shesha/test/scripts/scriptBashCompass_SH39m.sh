#!/bin/bash

PARFILE="/home/fvidal/compass/shesha/data/par/MICADO/micado_39m_SH.py" # scao_40x40_8pix.py  scao_64x64_8pix.py  scao_80x80_8pix.py "

DATE=`date +%F_%Hh%M`
SVN=`svnversion`
OUTPUT="$SHESHA_ROOT/test/scripts/resultatsScripts/scriptLog.txt"
rm $OUTPUT
echo "writing output in "$OUTPUT
echo "To monitor process: tail -f" $OUTPUT
script="$SHESHA_ROOT/test/scripts/script_SH39m.py"
for FREQ in "500" #"1000"
do
    for NPIX in "8" #"8"
    do
        for PIXSIZE in "1." #"1.5"
        do
          #for GAIN in "0.1" "0.2" "0.3" "0.4" "0.5" #"pyr" #
          for GAIN in "0.1" "0.2" "0.3" "0.4" "0.5" #"pyr" #
            #for GAIN in "0.5" #"pyr" #
            do
                for BP in "10" # nb of Brightest pixels
                #for TH in "1" #"pyr" #
                do
                    for RON in "3" # RON
                    do
                        #for MAG in "11" "12" "13" "14" "15" "16" "17" #"pyr" #
                        for MAG in "11" "12" "13" "13.5" "14" "14.5" "15" "15.5" "16" #"pyr" #
                        #for MAG in "11" "12" #"pyr" #
                        do
                          #for KLFILT in "500" "1000" "1500"
                          for KLFILT in "1000"
                            do
                                CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE $GAIN $BP $MAG $RON $KLFILT"
                                echo "execute $CMD" >> $OUTPUT
                                $CMD 2>> $OUTPUT >> $OUTPUT
                            done
                        done
                    done
                done
            done
        done
    done
done

#FREQ="500"
#NPIX="6"
#PIXSIZE="1"
#CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE 0.5 0 16 1000"
#echo "execute $CMD" >> $OUTPUT
#$CMD 2>> $OUTPUT >> $OUTPUT
#sleep 5

#CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE 0.5 0 16 1500"
#echo "execute $CMD" >> $OUTPUT
#$CMD 2>> $OUTPUT >> $OUTPUT
#sleep 5

#CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE 0.5 0 17 500"
#echo "execute $CMD" >> $OUTPUT
#$CMD 2>> $OUTPUT >> $OUTPUT
#sleep 5

#CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE 0.5 0 17 1000"
#echo "execute $CMD" >> $OUTPUT
#$CMD 2>> $OUTPUT >> $OUTPUT
#sleep 5

#CMD="python $script $PARFILE $FREQ $NPIX $PIXSIZE 0.5 0 17 1500"
#echo "execute $CMD" >> $OUTPUT
#$CMD 2>> $OUTPUT >> $OUTPUT
#sleep 5



# To monitor the script log:
# tail -f resultatsScripts/scriptLog.txt


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
echo "Script Done"

#FILES_LGS="scao_16x16_8pix_lgs.py"
#FILES_LGS+="scao_40x40_10pix_lgs.par"
#FILES_LGS+="scao_64x64_16pix_lgs.par"
#FILES_LGS+="scao_80x80_20pix_lgs.par"
