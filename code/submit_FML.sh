#!/bin/bash -e


RUNPATH=$HOME/scripts/Fortran_Matrix_Lab/code
cd $RUNPATH

if [ -n "$3" ]; then
	./FML $1 $2 $3
else
	./FML $1 $2
fi
