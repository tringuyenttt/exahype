#!/bin/bash

function verbose { echo $@; $@; }

set -e

TEST="TEST_NEW_PDE_AUTONOMOUSLY"

# instead of softlinks
FOR="../../GRMHD/" # fortran includes
Fortran="../../GRMHD/Fortran"
FortranInitialData="../../GRMHD/InitialData"

# WAS Dim3
DIMSTRING="Dim2"

verbose g++ -c -g3 --std=c++11 -D$DIMSTRING -D$TEST -I$FOR -Wall -I../ -I../../Peano/ \
	indi-test.cc ../InitialData.cpp \
	../PDE/Cons2Prim.cpp ../PDE/Prim2Cons.cpp ../PDE/PDE.cpp 
	
verbose g++ -g3 -oinditest.out \
	$FortranInitialData/InitialDataFort.o \
	$Fortran/{abort,C2P-GRMHD,C2PRoutines,Metric,Parameters,PDE}.o \
	{indi-test,InitialData,Cons2Prim,PDE,Prim2Cons}.o -lgfortran
