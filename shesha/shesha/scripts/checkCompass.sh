script="$SHESHA_ROOT/shesha/scripts/check.py"
for file in $SHESHA_ROOT/data/par/par4tests/*.py
do
    CMD="python $script $file"
    echo "execute $CMD"
    $CMD &> /dev/null
done
CMD="python $script osef --displayResult"
$CMD
