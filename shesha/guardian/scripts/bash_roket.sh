#!/bin/bash
DIAM="8"
NSSP="40"
NFILT="20"
NPIX="2 4 6"
PIXSIZE="0.2 0.4 0.6 0.8"
DATE=`date +%F_%Hh%M`
OUTPUT="$SHESHA_ROOT/guardian/scripts/outputfile_$DATE"

echo "writing output in "$OUTPUT

script="$SHESHA_ROOT/guardian/scripts/script_roket.py"

for n in $NPIX
do
    for p in $PIXSIZE
    do
        CMD="ipython $script $SHESHA_ROOT/data/par/par4roket/roket_8m_1layer.py -- --d $DIAM --npix $n --pixsize $p --nfilt $NFILT --nssp $NSSP -s roket_"$DIAM"m_nssp"$NSSP"_npix"$n"_psize"$p".h5"
        echo "execute $CMD" >> $OUTPUT
        $CMD 2>> $OUTPUT >> $OUTPUT
    done
done
