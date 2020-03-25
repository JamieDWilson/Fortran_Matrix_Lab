# compile Fortran_Matrix_Lab

# load user options
include user.mak

# define variables
objects = fml.o fml_lib.o io_module.o tm_module.o bg_module.o
switch = -O3 -fdefault-real-8

# Makefile

FML: $(objects)
	$(compiler) $(switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

fml_lib.mod: fml_lib.o fml_lib.f90
	$(compiler) -c $(switch) fml_lib.f90
fml_lib.o: fml_lib.f90
	$(compiler) -c $(switch) fml_lib.f90

io_module.mod: io_module.o io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include io_module.f90
io_module.o: fml_lib.mod io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include io_module.f90

tm_module.mod: tm_module.o tm_module.f90
	$(compiler) -c $(switch) tm_module.f90
tm_module.o: fml_lib.mod tm_module.f90
	$(compiler) -c $(switch) tm_module.f90

bg_module.mod: bg_module.o bg_module.f90
	$(compiler) -c $(switch) bg_module.f90
bg_module.o: fml_lib.mod bg_module.f90
	$(compiler) -c $(switch) bg_module.f90

fml.o: fml_lib.mod io_module.mod tm_module.mod bg_module.mod fml.f90
	$(compiler) -c $(switch) fml.f90

# Profiling
profile: switch = -O3 -fdefault-real-8 -pg -g
profile: $(objects)
	$(compiler) $(switch) $(profile_switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

# Cleaning
clean:
	rm fml_lib.mod fml_lib.o io_module.mod io_module.o tm_module.mod tm_module.o bg_module.mod bg_module.o fml.o

# End of makefile
