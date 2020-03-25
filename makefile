# compile Fortran_Matrix_Lab

# load user options
include user.mak

# define variables
objects = src/fml.o src/fml_lib.o src/io_module.o src/tm_module.o src/bg_module.o
switch = -O3 -fdefault-real-8
src_dir = src

# Makefile

FML: $(objects)
	$(compiler) $(switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

src/fml_lib.mod: src/fml_lib.o src/fml_lib.f90
	$(compiler) -c $(switch) src/fml_lib.f90 -o $@
src/fml_lib.o: src/fml_lib.f90
	$(compiler) -c $(switch) src/fml_lib.f90 -o $@

src/io_module.mod: src/io_module.o src/io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include src/io_module.f90 -o $@
src/io_module.o: src/fml_lib.mod src/io_module.f90
	$(compiler) -c $(switch) -I$(netcdf_switch)/include src/io_module.f90 -o $@

src/tm_module.mod: src/tm_module.o src/tm_module.f90
	$(compiler) -c $(switch) src/tm_module.f90 -o $@
src/tm_module.o: src/fml_lib.mod src/tm_module.f90
	$(compiler) -c $(switch) src/tm_module.f90 -o $@

src/bg_module.mod: src/bg_module.o src/bg_module.f90
	$(compiler) -c $(switch) src/bg_module.f90 -o $@
src/bg_module.o: src/fml_lib.mod src/bg_module.f90
	$(compiler) -c $(switch) src/bg_module.f90 -o $@

src/fml.o: src/fml_lib.mod src/io_module.mod src/tm_module.mod src/bg_module.mod src/fml.f90
	$(compiler) -c $(switch) src/fml.f90 -o $@

# Profiling
profile: switch = -O3 -fdefault-real-8 -pg -g
profile: $(objects)
	$(compiler) $(switch) $(profile_switch) $(objects) -L$(netcdf_switch)/lib $(netcdf_lib) -o FML

# Cleaning
clean:
	rm src/*.mod src/*.o

# End of makefile
