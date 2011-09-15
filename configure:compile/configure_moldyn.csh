#!/bin/csh

make distclean
./configure CC=mpicc  CXX="mpicxx -I/share/apps/hdf5/include -I/share/apps/mpi/gcc/openmpi-1.2.8/include -L/share/apps/hdf5/lib -lhdf5_hl -lhdf5 -lz -L/share/apps/mpi/gcc/openmpi-1.2.8/lib"
make
mv autodock4 autodock4.lga.MPI
