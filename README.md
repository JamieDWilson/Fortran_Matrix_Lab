# Fortran_Matrix_Lab

Fortran (Transport) Matrix Lab (FML) is a OCMIP-style carbon cycle model based on the 2.8 deg transport matrix (Khatiwala et al., 2005). The model is written in Fortran.

FML was developed and used in Wilson et al., (2019) Biogeosciences.

Current status (2020): I am updating this repository so that the code is available to use.

# Setting up FML
1) download or clone the code here
2) edit `user.mak` with your compiler, netcdf library location and library name
3) download the transport matrix files

# Running FML
1) create an experiment file in the `experiment` directory (see `expeiment/example_experiment`)
2) run `make` to compile the code
3) to run the experiment, run the FML executible with 2 or 3 arguments:
  - the experiment name
  - number of years to run
  - (optional) name of an experiment to restart from
4) e.g., `./FML initial_experiment 10 restart_experiment`
5) output is written to netcdf and text files in the `output` directory under the experiment name
  
# Documentation
Documentation of the code is contained in the doc directory.

To view, either open the index.html file in a browser or click on the link:
http://htmlpreview.github.io/?https://github.com/JamieDWilson/Fortran_Matrix_Lab/blob/master/doc/index.html
