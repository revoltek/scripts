#!/bin/bash

# to copy COLUMN2 in COLUMN1
#use: ./copy-col.sh file.ms COLUMN1 COLUMN2

if [ x$1 == 'x' -o x$2 == 'x' -o x$3 == 'x' ]
then
	echo 'use: ./copy-col.sh file.ms TO_COLUMN FROM_COLUMN'
else
	taql "using style python update $1 set $2 = $3"
fi
