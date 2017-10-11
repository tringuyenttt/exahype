SHELL=/bin/bash

TEST=TEST_NEW_PDE_AUTONOMOUSLY

# path to GRMHD base directory
GRMHD=../../..

# GRMHD Fortran formulation
GRMHD_Finc=$(GRMHD)/GRMHD
FortranObjects=$(GRMHD)/GRMHD/Fortran
FortranInitialData=$(GRMHD)/GRMHD/InitialData

# GRMHD cpp formulation
GRMHD_Cinc=$(GRMHD)/GRMHD_cpp

# WAS Dim3
DIMSTRING=Dim2

all:
	g++ -c -g3 --std=c++11 -D$(DIMSTRING) -D$(TEST) -I$(GRMHD_Finc) -I$(GRMHD_Cinc) -I../ -I../../Peano/ -Wall \
		indi-test.cc $(GRMHD_Cinc)/InitialData.cpp \
		../PDE/Cons2Prim.cpp ../PDE/Prim2Cons.cpp ../PDE/PDE.cpp 
	
	g++ -g3 -oinditest.out \
		$(FortranInitialData)/InitialDataFort.o \
		$(FortranObjects)/{abort,C2P-GRMHD,C2PRoutines,Metric,Parameters,PDE}.o \
		{indi-test,InitialData,Cons2Prim,PDE,Prim2Cons}.o -lgfortran

clean:
	rm -f *.o *.out
