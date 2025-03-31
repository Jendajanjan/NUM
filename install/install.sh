#!/bin/bash

cd ..
mkdir -p Build/Release

cd Build/Release
cmake -D CMAKE_BUILD_TYPE=Release ../../src

make -j 8
