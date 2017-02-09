#!/bin/sh

# User test

# Tests the full RTTOV-11 code for various MW and IR instruments.

ARG_ARCH=`perl -e 'for(@ARGV){m/^ARCH=(\S+)$/o && print "$1";}' $*`
if [ ! "x$ARG_ARCH" = "x" ]; then
  ARCH=$ARG_ARCH
fi
if [ "x$ARCH" = "x" ];
then
  echo 'Please supply ARCH'
  exit 1
fi

set -x

SESSION=test_rttov11
OPTS="IGNORETINY=1 $*"
WHAT="DIRECT=1 TL=1 AD=1 K=1"
CHECK="CHECK=1 TEST_REF=$SESSION.2"

./rttov_test.pl SESSION=$SESSION $WHAT $CHECK ARCH=$ARCH $OPTS -- << EOF
  TEST_LIST=amsre/001
  TEST_LIST=amsua/001,amsua/021clw
  TEST_LIST=amsub/001
  TEST_LIST=msu/001                 REFRACTION=1
  TEST_LIST=ssmis/001,ssmis/021
  TEST_LIST=windsat/001
  TEST_LIST=hirs/001                REFRACTION=1 APPLY_REG_LIMITS=1
  TEST_LIST=modis/021
  TEST_LIST=seviri/222              SOLAR=1
EOF

