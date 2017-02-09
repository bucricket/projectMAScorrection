#!/bin/sh

# User test

# Tests the full RTTOV-11 code for hi-res IR instruments including cloud/aerosol profiles.

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

SESSION=test_rttov11_hires
OPTS="IGNORETINY=1 TINYABS=1.E-8 $*"
WHAT="DIRECT=1 TL=1 AD=1 K=1"
CHECK="CHECK=1 TEST_REF=$SESSION.2"

./rttov_test.pl SESSION=$SESSION $WHAT $CHECK ARCH=$ARCH $OPTS -- << EOF
  TEST_LIST=airs/001,airs/241
  TEST_LIST=airs/281,airs/282,airs/283
  TEST_LIST=iasi/241,iasi/261           REFRACTION=1 SOLAR=1
  TEST_LIST=iasi/281,iasi/282,iasi/283  REFRACTION=1
EOF
