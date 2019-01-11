#!/bin/bash -e


RUNPATH=$HOME/scripts/TMM_model/code
cd $RUNPATH


if [ -n "$3" ]; then
	./FML $1 $2 $3
else
	./FML $1 $2
fi