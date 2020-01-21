#!/bin/bash

# pytest --html=doc/html/report_unit_test.html --self-contained-html libcarma/python_module/test shesha/tests/pytest/rtc

# pytest --cov-report html:doc/html/carma_cov_html --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report html:doc/html/sutra_cov_html --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report html:doc/html/shesha_cov_html --cov=shesha shesha/tests/pytest

# pytest --cov-report xml:carma_cov.xml --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report xml:sutra_cov.xml --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report xml:shesha_cov.xml --cov=shesha shesha/tests/pytest

# pytest --cov-report html:doc/html/coverage --cov-report xml:compass_cov.xml --cov=carmaWrap --cov=sutraWrap --cov=shesha libcarma/python_module/test shesha/tests/pytest/rtc

# ./shesha/tests/checkCompass.sh

coverage run -m pytest --verbose --junit-xml compass_tests.xml $COMPASS_ROOT/libcarma/python_module/test $SHESHA_ROOT/tests/pytest/rtc || echo "It's just a scratch"

# script="$SHESHA_ROOT/shesha/tests/check.py"
rm -f check.h5
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
rm -f check.h5

coverage xml -o compass_cov.xml
coverage html -d public
