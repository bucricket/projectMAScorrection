#!/bin/sh
#
# MW emissivity atlas test control script
#
# 14/06/2010   Created (J. Hocking)
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
EXEC=../../$BIN/rttov_mwatlas_test.exe
REF_TEST=../test_emisatlas.2

# Relative to the rttov_test directory:
TEST=./test_emisatlas.1

TEST_INPUT_FILE="rttov_mwatlas_test_input"
MONTH=08
# Coefficient filenames are specified below

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

ATLAS_EMIS_FILE="ssmi_mean_emis_climato_${MONTH}_cov_interpol"
ATLAS_COR_FILE="correlations"

if [ ! -d $TEST ]; then
  echo 'Test directory required, containing profiles_mw'
  echo 'Also, make sure to run directly from the rttov_test directory'
  exit
fi

cwd=`pwd`
cd $TEST

echo " "
echo " "
echo " Test TELSEM MW emissivity atlas "
echo " "

for TEST_NUMBER in 01 02 
# for TEST_NUMBER in 01 
do

  case $TEST_NUMBER in
  
    '01') COEF_FILENAME=rtcoef_noaa_15_amsua.dat
   
    cat  > $TEST_INPUT_FILE  <<EOF
$MONTH          # Month for which to load emissivity data
$COEF_FILENAME
23              # Zenith angle (degrees)
$DATAATLAS
EOF

    ;;

    '02') COEF_FILENAME=rtcoef_dmsp_13_ssmi.dat

    cat  > $TEST_INPUT_FILE  <<EOF
$MONTH          # Month for which to load emissivity data
$COEF_FILENAME
57              # Zenith angle (degrees)
$DATAATLAS
EOF

    ;;
  esac

  rm -f $COEF_FILENAME
  if [ -s  $DATART/$COEF_FILENAME ] ; then
    ln -s $DATART/$COEF_FILENAME $COEF_FILENAME
  else
    echo "rtcoef file not found"
    cd $cwd
    exit
  fi

  if [ ! -s $DATAATLAS/$ATLAS_EMIS_FILE ] ; then
    echo "Atlas emissivity file not found"
    cd $cwd
    exit
  fi
  
  if [ ! -s $DATAATLAS/$ATLAS_COR_FILE ] ; then
    echo "Atlas correlations file not found"
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

  mv output_mwatlas.ascii  output_mwatlas.$TEST_NUMBER.$ARCH

  if [ $? -ne 0 ]; then
    echo "Expected output file not found"
    exit 1
  fi

  echo
  echo "Output is in the file ${TEST}/output_mwatlas.$TEST_NUMBER.$ARCH"

  DIFF_FILE=diff_mwatlas.$TEST_NUMBER.$ARCH

  if [ -f $REF_TEST/output_mwatlas.$TEST_NUMBER ] ; then
    diff -biw output_mwatlas.$TEST_NUMBER.$ARCH \
      $REF_TEST/output_mwatlas.$TEST_NUMBER > $DIFF_FILE

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

done

cd $cwd

exit
