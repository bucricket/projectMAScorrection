#!/bin/sh
#
# MW emissivity atlas test control script
#
# 16/07/2010   Created (P. Brunel, J. Hocking)
#

# Set BIN directory if supplied
BIN=`perl -e 'for(@ARGV){m/^BIN=(\S+)$/o && print "$1";}' $*`
if [ "x$BIN" = "x" ]
then
  BIN=bin
fi

######## Edit this section for your pathnames #######

# Relative to the rttov_test/test_emisatlas.1 directory:
DATART=../../rtcoef_rttov11/rttov7pred54L
DATAATLAS=../../emis_data
EXEC=../../$BIN/rttov_cnrmmwatlas_test.exe
REF_TEST=../test_emisatlas.2

# Relative to the rttov_test directory:
TEST=./test_emisatlas.1

TEST_NUMBER=01
TEST_INPUT_FILE="rttov_cnrmmwatlas_test_input"
MONTH=06
COEF_FILENAME="rtcoef_noaa_15_amsua.dat"

###################################################

ARG_ARCH=`perl -e 'for(@ARGV){m/^ARCH=(\S+)$/o && print "$1";}' $*`
if [ ! "x$ARG_ARCH" = "x" ]; then
  ARCH=$ARG_ARCH
fi
if [ "x$ARCH" = "x" ];
then
  echo 'Please supply ARCH'
  exit 1
fi

if [ ! -d $TEST ]; then
  echo 'Test directory required, containing profiles_mw'
  echo 'Also, make sure to run directly from the rttov_test directory'
  exit
fi

cwd=`pwd`
cd $TEST

echo " "
echo " "
echo " Test CNRM MW emissivity atlas "
echo " "

cat  > $TEST_INPUT_FILE  <<EOF
$COEF_FILENAME
$DATAATLAS
$MONTH          # Month for which to load emissivity data
650             # number of profiles
85, -0.20       # latitude start, latitude step
20,  0.00       # longitude start, longitude step
 0,  0.08       # Zenith angle (degrees)
6               # number of channels
1,2,3,4,7,15    # list of channels
EOF

rm -f $COEF_FILENAME
if [ -s  $DATART/$COEF_FILENAME ] ; then
  ln -s $DATART/$COEF_FILENAME $COEF_FILENAME
else
  echo "rtcoef file not found"
  cd $cwd
  exit
fi

if [ ! -d $DATAATLAS ] ; then
  echo "Atlas emissivity directory not found"
  cd $cwd
  exit
fi


$EXEC < $TEST_INPUT_FILE

if [ $? -ne 0 ]; then
  echo " "
  echo "TEST FAILED"
  echo " "
  exit 1
fi

mv output_cnrmmwatlas.ascii  output_cnrmmwatlas.$TEST_NUMBER.$ARCH

if [ $? -ne 0 ]; then
  echo "Expected output file not found"
  exit 1
fi

echo
echo "Output is in the file ${TEST}/output_cnrmmwatlas.$TEST_NUMBER.$ARCH"

DIFF_FILE=diff_cnrmmwatlas.$TEST_NUMBER.$ARCH

if [ -f $REF_TEST/output_cnrmmwatlas.$TEST_NUMBER ] ; then
  diff -biw output_cnrmmwatlas.$TEST_NUMBER.$ARCH \
    $REF_TEST/output_cnrmmwatlas.$TEST_NUMBER > $DIFF_FILE

  if [ -s $DIFF_FILE ]; then
    echo "--- Diff file contents: ---"
    cat $DIFF_FILE
    echo "---------------------------"
  else
    echo " "
    echo "Diff file has zero size: TEST SUCCESSFUL"
    echo " "
  fi
else
  echo "Test reference output not found"
fi
echo

rm -f $TEST_INPUT_FILE
rm -f $COEF_FILENAME


cd $cwd

exit
