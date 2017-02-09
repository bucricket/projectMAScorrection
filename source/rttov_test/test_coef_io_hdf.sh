#!/bin/sh

# Tests coefficient I/O routines for HDF5

# Does not need reference (.2 directory) for this test.

# grok ARCH environment variable

ARG_ARCH=`perl -e 'for(@ARGV){m/^ARCH=(\S+)$/o && print "$1";}' $*`
if [ ! "x$ARG_ARCH" = "x" ]; then
  ARCH=$ARG_ARCH
fi
if [ "x$ARCH" = "x" ];
then
  echo 'Please supply ARCH'
  exit 1
fi

TEST_LIST="amsua/001,seviri/081,hirs/001,ssu/101pmc,seviri/201,airs/061_ascii,iasi/261varichan_ascii,iasi/282varichan_ascii,iasi/283varichan_ascii,iasi/pc241_ascii"

if [ -d "coefs.1.$ARCH" ]
then

for t in $(perl -e '@x = split (m/,/o, $ARGV[0]); print "@x"' $TEST_LIST)
do
  for f in $(find coefs.1.$ARCH/$t -type l)
  do
    l=$(perl -e 'use File::Basename;my $f = shift; print dirname ($f) . "/" . readlink ($f) ' "$f")
    \rm -f "$l"
  done
  \rm -rf coefs.1.$ARCH/$t
done

fi


set -x


\rm -rf test_coef_io_hdf.*.1.$ARCH


# usual way to proceed (read only channels of interest from coefficient file), create a reference

./rttov_test.pl DIRECT=1 PRINT=1 TEST_LIST=$TEST_LIST ARCH=$ARCH \
                SESSION=test_coef_io_hdf $*

# extract required channels from coefficient files, save them to disk and read them again
# compare to reference

./rttov_test.pl DIRECT=1 PRINT=1 TEST_LIST=$TEST_LIST ARCH=$ARCH \
                SESSION=test_coef_io_hdf.hdf5 COEF_EXTRACT=1 \
                COEF_FORMAT=hdf5 TEST_REF=test_coef_io_hdf.1.$ARCH $*

