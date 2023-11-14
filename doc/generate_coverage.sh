#!/bin/bash
# Exit on first error.
set -e

coverage run -m pytest $COMPASS_ROOT/python_module/carmaWrap/test 
coverage run -a -m pytest $SHESHA_ROOT/tests/pytest/rtc $SHESHA_ROOT/tests/pytest/rtc_standalone $SHESHA_ROOT/tests/pytest/supervisor #|| echo "It's just a scratch"

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
