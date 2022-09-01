#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./xsample7* ./xsample04
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./xsample07 ./xsample07_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./xsample07 ./xsample07_avx2

make rebuild
cp ./xsample07 ./xsample07_release

