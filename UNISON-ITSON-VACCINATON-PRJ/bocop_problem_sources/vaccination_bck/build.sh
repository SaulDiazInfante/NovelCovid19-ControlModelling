#!/bin/bash

## shell script for Bocop build
## Pierre Martinon, Inria
## 2019

## set build folder
if ! [ -d build ]; then
  mkdir build
fi
cd build

## default cmake options
NLPPATH=${BOCOPPATH:-"/home/saul/Desktop/Bocop-2.2.1-linux"}
platform=""
buildtype="Release" 

## set specific cmake options
#if [[ "$(uname -s)" =~ "MSYS" ]]; then
#  platform="-G \"MSYS Makefiles\" "
#fi
KERNEL=`uname -s`
case "$KERNEL" in
*"MSYS"*) platform="-G \"MSYS Makefiles\" ";;
esac

while getopts d option
do
    case "${option}" in d) buildtype="Debug";;
    esac
done

## launch cmake, make and go back to problem folder
cmake ${platform} -DCMAKE_BUILD_TYPE=${buildtype} ${NLPPATH}
make -j
cd -
