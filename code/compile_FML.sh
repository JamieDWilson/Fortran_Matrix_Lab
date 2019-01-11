#!/bin/bash -e

gfortran -c -O3 fml_lib.f90 
gfortran -c -O3 -I/share/apps/netcdf-4.0/gcc-4.4.7/include  tm_module.f90 
gfortran -c -O3 bg_module.f90 
gfortran -c -O3  fml.f90

gfortran -O3 fml.o fml_lib.o tm_module.o bg_module.o -L/share/apps/netcdf-4.0/gcc-4.4.7/lib -lnetcdf -o FML




