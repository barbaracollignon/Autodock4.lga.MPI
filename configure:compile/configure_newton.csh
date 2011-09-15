#!/bin/csh

make distclean
./configure CXX=mpicxx CC=cc CPPFLAGS="-I/data/apps/hdf5/1.6.10/include/" LIBS="-L/data/apps/hdf5/1.6.10/lib -lhdf5_hl -lhdf5 -lz" 
make
mv autodock4 autodock4.lga.MPI
