#!/bin/sh

# User test

# Tests PC-RTTOV.

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

SESSION=test_pc
OPTS="IGNORETINY=1 TINYABS=1.E-8 $*"
WHAT="DIRECT=1 TL=1 AD=1 K=1"
CHECK="CHECK=1 TEST_REF=$SESSION.2"

./rttov_test.pl SESSION=$SESSION $WHAT $CHECK ARCH=$ARCH $OPTS -- << EOF
  TEST_LIST=airs/pc241 SWITCHRAD=0
  TEST_LIST=iasi/pc261
EOF
