#!/bin/sh
#
# RTTOV_SCATT tests control script
#
# This test requires the (~1GB) Mie Tables to be available
# and placed in a directory pointed to by DATAMIE
# 
# 13/12/2007   Alan Geer   First version
# 24/02/2010   Alan Geer   Temporary version for RTTOV-10 - only N15 AMSU-A
# 30/07/2010   J Hocking   Add SSM/I test
#

# Set BIN directory if supplied
BIN=`perl -e 'for(@ARGV){m/^BIN=(\S+)$/o && print "$1";}' $*`
if [ "x$BIN" = "x" ]
then
  BIN=bin
fi

######## Edit this section for your pathnames #######

# Relative to the rttov_test/test_rttovscatt.1 directory:
DATART=../../rtcoef_rttov11/rttov7pred54L
DATAMIE=../../rtcoef_rttov11/mietable
EXEC=../../$BIN/rttovscatt_test.exe
REF_TEST=../test_rttovscatt.2

# Relative to the rttov_test directory:
TEST=./test_rttovscatt.1

TEST_INPUT_FILE=rttovscatt_test_input

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

if [ ! -d ${TEST} ]; then
  echo 'Test directory required, containing profiles2_fmt'
  echo 'Also, make sure to run directly from the rttov_test directory'
  exit
fi

cwd=`pwd`
cd $TEST

echo " "
echo " "
echo " Test RTTOV-SCATT "
echo " "

for TEST_NUMBER in 01 02 
do
   
#  export PROGINF=DETAIL
#  export DR_HOOK=1
#  export DR_HOOK_OPT='hpmprof'
#  export OMP_NUM_THREADS=1

  case ${TEST_NUMBER} in
  
    '01') COEF_FILENAME=rtcoef_dmsp_13_ssmi.dat
          MIETABLE_FILENAME=mietable_dmsp_ssmi.dat
   
    cat  > ${TEST_INPUT_FILE}  <<EOF
$COEF_FILENAME
57  # zenith angle (degrees)
EOF

    ;;

    '02') COEF_FILENAME=rtcoef_noaa_15_amsua.dat
          MIETABLE_FILENAME=mietable_noaa_amsua.dat

    cat  > ${TEST_INPUT_FILE}  <<EOF
$COEF_FILENAME
23  # zenith angle (degrees)
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

  rm -f $MIETABLE_FILENAME
  if [ -s  $DATAMIE/$MIETABLE_FILENAME ] ; then
    ln -s $DATAMIE/$MIETABLE_FILENAME $MIETABLE_FILENAME
  else
    echo "Mie tables not found"
    cd $cwd
    exit
  fi

  $EXEC < ${TEST_INPUT_FILE}

  if [ $? -ne 0 ]; then
    echo " "
    echo "TEST FAILED"
    echo " "
    exit 1
  fi

#  mv drhook.prof.1      drhook.prof.${TEST_NUMBER}.rttov_scatt
  mv outputscatt_full.ascii  output_full.${TEST_NUMBER}.rttov_scatt.${ARCH}
  mv outputscatt.ascii       output.${TEST_NUMBER}.rttov_scatt.${ARCH}

  if [ $? -ne 0 ]; then
    echo "Expected output file not found"
    exit 1
  fi

  echo
  echo "Output is in the file ${TEST}/output.${TEST_NUMBER}.rttov_scatt.${ARCH}"

  DIFF_FILE=diff.${TEST_NUMBER}.${ARCH}

  if [ -f $REF_TEST/output.${TEST_NUMBER}.rttov_scatt ] ; then
    diff -biw output.${TEST_NUMBER}.rttov_scatt.${ARCH} \
      $REF_TEST/output.${TEST_NUMBER}.rttov_scatt > $DIFF_FILE

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

  rm -f ${TEST_INPUT_FILE}
done

cd $cwd

exit
