# compile Fortran_Matrix_Lab
# defining variables
objects = fml.o fml_lib.o tm_module.o bg_module.o
compiler = gfortran
switch = -O3
netcdf_switch = -L/share/apps/netcdf-4.0/gcc-4.4.7/lib -lnetcdf

# Makefile

FML: $(objects)
	$(compiler) $(switch) $(objects) $(netcdf_switch) -o FML

fml_lib.mod: fml_lib.o fml_lib.f90
	$(compiler) -c $(switch) fml_lib.f90
fml_lib.o: fml_lib.f90
	$(compiler) -c $(switch) fml_lib.f90

tm_module.mod: tm_module.o tm_module.f90
	$(compiler) -c $(switch) -I/share/apps/netcdf-4.0/gcc-4.4.7/include tm_module.f90
tm_module.o: fml_lib.mod tm_module.f90
	$(compiler) -c $(switch) -I/share/apps/netcdf-4.0/gcc-4.4.7/include tm_module.f90

bg_module.mod: bg_module.o bg_module.f90
	$(compiler) -c $(switch) bg_module.f90
bg_module.o: fml_lib.mod bg_module.f90
	$(compiler) -c $(switch) bg_module.f90

fml.o: fml_lib.mod tm_module.mod bg_module.mod fml.f90
	$(compiler) -c $(switch) fml.f90

# Cleaning 
clean:
	rm fml_lib.mod fml_lib.o tm_module.mod tm_module.o bg_module.mod bg_module.o fml.o

# End of makefile
