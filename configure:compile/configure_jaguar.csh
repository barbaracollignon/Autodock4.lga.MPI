#!/bin/csh


##me clean
##module load hdf5-parallel/
##module load perftools
##make clean
cp main.cc main.c
./configure CXX=CC CC=CC CPPFLAGS="-I/opt/cray/hdf5-parallel/1.8.4.1/hdf5-parallel-pgi/include" LIBS="-L/opt/cray/hdf5-parallel/1.8.4.1/hdf5-parallel-pgi/lib -lhdf5_hl -lhdf5 -lz" 
make
