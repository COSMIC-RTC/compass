#!/bin/bash

coverage run -m pytest $COMPASS_ROOT/libcarma/python_module/test $SHESHA_ROOT/tests/pytest/rtc #|| echo "It's just a scratch"

# script="$SHESHA_ROOT/shesha/tests/check.py"
rm -f $COMPASS_ROOT/check.h5
script="tests.check"
conf_path="$SHESHA_ROOT/data/par/par4tests"
nb_test=$(ls -1 $conf_path/*.py | wc -l)
current_test=1
for file in $conf_path/*.py
do
    name=$(basename $file ".py")
    CMD="coverage run --append -m $script $file"
    echo "[$current_test/$nb_test] running $name"
    $CMD > /dev/null
    # let "current_test++"
    current_test=$(expr $current_test + 1)
done
rm -f $COMPASS_ROOT/check.h5
