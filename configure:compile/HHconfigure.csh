#!/bin/csh

make distclean
module unload hdf5
module load hdf5-parallel/1.8.3.1
setenv CC  cc
setenv CXX "CC $HDF5_POST_LINK_OPTS $HDF5_INCLUDE_OPTS"
setenv CXXFLAGS -DMPICH_IGNORE_CXX_SEEK
./configure 
module unload hdf5-parallel
module load hdf5
make
