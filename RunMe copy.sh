#!/usr/bin/env sh

HOME=$(pwd)

BLADE=$HOME'/bin/'

echo "****************************************\n"
echo "********* 3D BLADE PARAMETRIZER ********\n"
echo "****************************************\n\n"
echo "Add Following in your .bashrc\n"
echo "\texport BLADE_HOME="$HOME
echo "\texport BIN="$BLADE
echo "\tPATH=\$PATH:\$BIN\n\n"
echo "*************************************\n\n"

PYTHON3=$(which python3)

if [ "$PYTHON3" != "/usr/bin/python3" ]; then
  echo "IMPORTANT !!! \n Python3 directory is different than that used in the bin file. \nPlease change the header in bin/MakeBlade.py file to:\n"
  echo $PYTHON3
  echo "\n"
fi
