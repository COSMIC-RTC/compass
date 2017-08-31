#!/bin/bash
# Installation:
#   cd my_gitproject
#   wget -O pre-commit.sh http://tinyurl.com/mkovs45
#   ln -s ../../pre-commit.sh .git/hooks/pre-commit.legacy
#   chmod +x pre-commit.sh

OPTIONS="-A8 -t8 --lineend=linux"

RETURN=0
ASTYLE=$(which astyle)
if [ $? -ne 0 ]; then
	echo "[!] astyle not installed. Unable to check source file format policy." >&2
	exit 1
fi

FILES=`git diff --cached --name-only --diff-filter=ACMR | grep -E "\.(c|cpp|h)$"`
FILES=`find . -name "*.cu" -or -name "*.hpp" -or -name "*.cuh" -or -name "*.h"`
for FILE in $FILES; do
	$ASTYLE $OPTIONS < $FILE | cmp -s $FILE -
	if [ $? -ne 0 ]; then
		echo "[!] $FILE does not respect the agreed coding style." >&2
		RETURN=1
	fi
done

if [ $RETURN -eq 1 ]; then
	echo "" >&2
	echo "Make sure you have run astyle with the following options:" >&2
	echo $OPTIONS >&2
fi

exit $RETURN


CONFIG_REGMAP_I2C=y
CONFIG_I2C=y
CONFIG_I2C_BOARDINFO=y
CONFIG_I2C_ALGOBIT=y
CONFIG_SND_SOC_I2C_AND_SPI=y
