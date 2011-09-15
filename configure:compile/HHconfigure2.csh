#!/bin/csh

make distclean
module unload hdf5
module load hdf5-parallel
setenv CC  cc
setenv CXX "CC"
setenv CXXFLAGS -DMPICH_IGNORE_CXX_SEEK
./configure 
module unload hdf5-parallel
module load hdf5
make
