#!/bin/sh

gzip -r9 ./*.dat ./*.f90

if [ ! -e $1 ]
then
   mkdir ${1}
fi
scp -r ./*.dat.gz ./*.f90.gz giga20:${1}

