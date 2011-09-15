#!/bin/csh

make distclean
./configure CXX=CC CC=cc CPPFLAGS="-I/opt/cray/hdf5-parallel/1.8.3.1/hdf5-parallel-pgi/include" LIBS="-L/opt/cray/hdf5-parallel/1.8.3.1/hdf5-parallel-pgi/lib -lhdf5_hl -lhdf5 -lz" 
make
mv autodock4 autodock4.lga.MPI
