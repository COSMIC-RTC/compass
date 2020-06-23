#!/bin/bash

# pytest --html=doc/html/report_unit_test.html --self-contained-html libcarma/python_module/test shesha/tests/pytest/rtc

# pytest --cov-report html:doc/html/carma_cov_html --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report html:doc/html/sutra_cov_html --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report html:doc/html/shesha_cov_html --cov=shesha shesha/tests/pytest

# pytest --cov-report xml:carma_cov.xml --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report xml:sutra_cov.xml --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report xml:shesha_cov.xml --cov=shesha shesha/tests/pytest

rm $COMPASS_ROOT/doc/html/report_unit_test.html
# TODO: generate large file with errors... Need to be fixed
# pytest --html=$COMPASS_ROOT/doc/html/report_unit_test.html --self-contained-html --cov-report html:$COMPASS_ROOT/doc/html/coverage --cov=carmaWrap --cov=sutraWrap --cov=shesha $COMPASS_ROOT/libcarma/python_module/test $COMPASS_ROOT/shesha/tests/pytest

# script="$SHESHA_ROOT/shesha/tests/check.py"
rm -f $COMPASS_ROOT/check.h5
script="tests.check"
conf_path="$SHESHA_ROOT/data/par/par4tests"
nb_test=$(ls -1 $conf_path/*.py | wc -l)
current_test=1
for file in $conf_path/*.py
do
    name=$(basename $file ".py")
    CMD="python -m $script $file"
    echo "[$current_test/$nb_test] running $name"
    $CMD > /dev/null
    # let "current_test++"
    current_test=$(expr $current_test + 1)
done

python -m $script osef --displayResult --repportResult=$COMPASS_ROOT/report_E2E.md

rm -rf $COMPASS_ROOT/public

doxygen $COMPASS_ROOT/doc/Doxyfile

mv doc/doxygen-doc/public .
$COMPASS_ROOT/doc/correctDoxygen.sh

mkdir -p public/coverage
coverage html -d public/coverage

echo 'Documentation generated in public/. To Publish it:'
echo 'rsync -PaW --inplace --del public/* lesia:compass-doc/html'