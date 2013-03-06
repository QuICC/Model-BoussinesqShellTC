#!/usr/bin/bash

sed -e "s/++TESTNAME++/$1/g" TemplateTest.cpp > $2/$1Test.cpp
mkdir -p $2/cmake.d
cp TemplateTestSourcesList.cmake $2/cmake.d/$1TestSourcesList.cmake
