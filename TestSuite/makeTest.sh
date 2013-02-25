#!/usr/bin/bash

sed -e "s/++TESTNAME++/$1/g" TemplateTest.cpp > $2/$1Test.cpp
cp TemplateTestSourcesList.cmake $2/$1TestSourcesList.cmake
