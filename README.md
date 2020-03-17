# Fortran_Matrix_Lab

Fortran (Transport) Matrix Lab (FML) is a OCMIP-style carbon cycle model based on the 2.8 deg transport matrix (Khatiwala et al., 2005). The model is written in Fortran.

FML was developed and used in Wilson et al., (2019) Biogeosciences.

Current status (2020): I am updating this repository so that the code is available to use.

# Setting up FML
1) download the code
2) edit the makefile with your compiler of choice and netcdf library location
3) download the transport matrix files

# Running FML
1) create an experiment file in the experiment directory
2) within the code directory, run the makefile to compile the code
3) to run the experiment from cold: ./FML <experiment_name> <number of years>
4) to run from a restart: ./FML <experiment_name> <number of years> <restart_experiment_name>
5) output is written to netcdf and text files in the output directory under the experiment name
  
# Documentation
Documentation of the code is contained in the doc directory.

To view, either open the index.html file or click on the link below.
http://htmlpreview.github.io/?https://github.com/JamieDWilson/Fortran_Matrix_Lab/blob/master/doc/index.html
