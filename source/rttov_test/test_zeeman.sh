#!/bin/sh

# User test

# Tests with Zeeman coefficient files using profiles selected for this purpose.

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

SESSION=test_zeeman
OPTS="IGNORETINY=1 $*"
WHAT="DIRECT=1 TL=1 AD=1 K=1"
CHECK="CHECK=1 TEST_REF=$SESSION.2"

./rttov_test.pl SESSION=$SESSION $WHAT $CHECK ARCH=$ARCH $OPTS -- << EOF
  TEST_LIST=amsua/zm01
  TEST_LIST=ssmis/zm01,ssmis/zm02,ssmis/zm13,ssmis/zm24
EOF

