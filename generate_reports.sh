#!/bin/bash

pytest --html=doc/html/report_unit_test.html --self-contained-html libcarma/python_module/test shesha/tests/pytest/rtc

# pytest --cov-report html:doc/html/carma_cov_html --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report html:doc/html/sutra_cov_html --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report html:doc/html/shesha_cov_html --cov=shesha shesha/tests/pytest

# pytest --cov-report xml:carma_cov.xml --cov=carmaWrap libcarma/python_module/test
# pytest --cov-report xml:sutra_cov.xml --cov=sutraWrap shesha/tests/pytest
# pytest --cov-report xml:shesha_cov.xml --cov=shesha shesha/tests/pytest

# pytest --cov-report html:doc/html/coverage --cov-report xml:compass_cov.xml --cov=carmaWrap --cov=sutraWrap --cov=shesha libcarma/python_module/test shesha/tests/pytest/rtc

./shesha/tests/checkCompass.sh
