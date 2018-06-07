#!/bin/bash

#echo "Checking CPP/CU files..."

# astyle options
# EXE="astyle"
# OPTIONS="--style=google --indent=tab"

# find . -type f -a \( -type d -path .git -prune -o -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" -o -name "*.cuh" \) -exec clang-format -i -style=Google {} \;

# clang-format options
EXE="clang-format"
OPTIONS="-style=Google"

FORMATTER=$(which $EXE)
if [ $? -ne 0 ];
then
  echo "[!] $EXE not installed. Unable to check source file format policy." >&2
  exit 1
fi

RETURN=0
git diff --cached --name-status --diff-filter=ACMR | {
  # Command grouping to workaround subshell issues. When the while loop is
  # finished, the subshell copy is discarded, and the original variable
  # RETURN of the parent hasn't changed properly.
  while read STATUS FILE;
  do
    if [[ "$FILE" =~ \.+(c|cpp|h|cu|cuh)$ ]];
    then
      $FORMATTER $OPTIONS $FILE > $FILE.beautified
      md5sum -b $FILE | { read stdin; echo $stdin.beautified; } | md5sum -c --status -
      if [ $? -ne 0 ];
      then
        # echo "[!] $FILE does not respect the agreed coding standards." >&2
        mv $FILE.beautified $FILE
        RETURN=1
      else
        # echo "[!] $FILE respects the agreed coding standards." >&2
        rm $FILE.beautified
      fi
    # else
    #   echo "[!] $FILE not checked" >&2
    fi
  done

  exit $RETURN
}
