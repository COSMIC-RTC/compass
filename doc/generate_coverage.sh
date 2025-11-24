#!/bin/bash
# Exit on first error.
set -e

# Initialize coverage database
coverage erase

# Run unit tests with coverage
coverage run -m pytest --junitxml=report.xml \
  $SHESHA_ROOT/tests/unit \
  $SHESHA_ROOT/tests/rtc \
  $SHESHA_ROOT/tests/supervisor

# Run CARMA tests with coverage (append to existing)
if [ -d "$COMPASS_ROOT/python_module/carma/test" ]; then
  coverage run -a -m pytest $COMPASS_ROOT/python_module/carma/test || true
fi

# Run additional tests with coverage (append to existing)
if [ -d "$SHESHA_ROOT/tests/rtc_standalone" ]; then
  coverage run -a -m pytest $SHESHA_ROOT/tests/rtc_standalone || true
fi

# Run end-to-end tests if available
rm -f check.h5
if [ -f "$SHESHA_ROOT/tests/check.py" ]; then
  script="$SHESHA_ROOT/tests/check.py"
  conf_path="$SHESHA_ROOT/data/par/par4tests"
  
  if [ -d "$conf_path" ]; then
    nb_test=$(ls -1 $conf_path/*.py 2>/dev/null | wc -l)
    
    if [ $nb_test -gt 0 ]; then
      current_test=1
      for file in $conf_path/*.py
      do
        name=$(basename $file ".py")
        echo "[$current_test/$nb_test] running $name"
        coverage run --append $script $file > /dev/null 2>&1 || true
        current_test=$(expr $current_test + 1)
      done
      
      # Generate E2E report if script supports it
      python $script osef --displayResult --repportResult=report_E2E.md || true
    fi
  fi
fi

# Generate coverage reports
echo ""
echo "=== Coverage Summary ==="
coverage report

# Generate HTML coverage report for GitLab pages
coverage html -d htmlcov

# Generate XML coverage report (for potential CI tools)
coverage xml -o coverage.xml

echo "Coverage reports generated:"
echo "  - Text report (above)"
echo "  - HTML report: htmlcov/index.html"
echo "  - XML report: coverage.xml"
