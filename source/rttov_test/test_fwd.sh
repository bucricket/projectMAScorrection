#!/bin/sh

# User test

# Tests forward model for a wide range of instruments.

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

SESSION=test_fwd
OPTS="IGNORETINY=1 $*"
WHAT="DIRECT=1"
CHECK="CHECK=1 TEST_REF=$SESSION.2"

./rttov_test.pl SESSION=$SESSION $WHAT $CHECK ARCH=$ARCH $OPTS -- << EOF
  TEST_LIST=airs/001,airs/241
  TEST_LIST=amsre/001
  TEST_LIST=amsua/001,amsua/002,amsua/003,amsua/004,amsua/005
  TEST_LIST=amsua/021clw,amsua/099
  TEST_LIST=amsub/001
  TEST_LIST=atsr/001,atsr/002,atsr/003
  TEST_LIST=atsr/201,atsr/202,atsr/203 SOLAR=1
  TEST_LIST=avhrr/011,avhrr/012,avhrr/013
  TEST_LIST=avhrr/211,avhrr/212,avhrr/213 SOLAR=1
  TEST_LIST=hirs/015,hirs/016,hirs/017
  TEST_LIST=iasi/001,iasi/241
  TEST_LIST=imager/001,imager/002,imager/009,imager/010
  TEST_LIST=mhs/001,mhs/002,mhs/003
  TEST_LIST=modis/001,modis/002
  TEST_LIST=modis/201,modis/202 SOLAR=1
  TEST_LIST=msu/001,msu/002
  TEST_LIST=mviri/001,mviri/002,mviri/003
  TEST_LIST=seviri/001,seviri/002,seviri/099
  TEST_LIST=seviri/201,seviri/202,seviri/299 SOLAR=1
  TEST_LIST=sounder/007,sounder/008
  TEST_LIST=ssmi/007,ssmi/008
  TEST_LIST=ssmis/001,ssmis/002,ssmis/003,ssmis/021
  TEST_LIST=tmi/001
  TEST_LIST=vissr/001,vissr/002
  TEST_LIST=windsat/001,windsat/021
EOF
