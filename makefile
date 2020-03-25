# compile Fortran_Matrix_Lab

# load user options
include user.mak

# define variables
objects = code/fml.o code/fml_lib.o code/io_module.o code/tm_module.o code/bg_module.o
switch = -O3 -fdefault-real-8
src_dir = code

# Makefile

FML: $(objects)
	$(compiler) $(switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

code/fml_lib.mod: code/fml_lib.o code/fml_lib.f90
	$(compiler) -c $(switch) code/fml_lib.f90 -o $@
code/fml_lib.o: code/fml_lib.f90
	$(compiler) -c $(switch) code/fml_lib.f90 -o $@

code/io_module.mod: code/io_module.o code/io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include code/io_module.f90 -o $@
code/io_module.o: code/fml_lib.mod code/io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include code/io_module.f90 -o $@

code/tm_module.mod: code/tm_module.o code/tm_module.f90
	$(compiler) -c $(switch) code/tm_module.f90 -o $@
code/tm_module.o: code/fml_lib.mod code/tm_module.f90
	$(compiler) -c $(switch) code/tm_module.f90 -o $@

code/bg_module.mod: code/bg_module.o code/bg_module.f90
	$(compiler) -c $(switch) code/bg_module.f90 -o $@
code/bg_module.o: code/fml_lib.mod code/bg_module.f90
	$(compiler) -c $(switch) code/bg_module.f90 -o $@

code/fml.o: code/fml_lib.mod code/io_module.mod code/tm_module.mod code/bg_module.mod code/fml.f90
	$(compiler) -c $(switch) code/fml.f90 -o $@

# Profiling
profile: switch = -O3 -fdefault-real-8 -pg -g
profile: $(objects)
	$(compiler) $(switch) $(profile_switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

# Cleaning
clean:
	rm code/*.mod code/*.o
	#rm fml_lib.mod fml_lib.o io_module.mod io_module.o tm_module.mod tm_module.o bg_module.mod bg_module.o fml.o

# End of makefile
