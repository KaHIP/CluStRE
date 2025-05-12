#!/bin/bash

NCORES=12
unamestr=`uname`

rm -rf deploy
rm -rf build
mkdir build
cd build 
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_CXX_FLAGS="-w" ../
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg ../
make -j $NCORES
cd ..

mkdir deploy
cp ./build/clustre deploy/

#rm -rf build
