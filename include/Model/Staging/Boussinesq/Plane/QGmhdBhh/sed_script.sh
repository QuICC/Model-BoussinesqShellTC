#!/bin/bash

for file in *
do

sed s/QGMHDBHH/QGMHDBHH/g < $file > new; mv new $file
sed s/QGmhdBhh/QGmhdBhh/g < $file > new; mv new $file

echo $file
 
done
