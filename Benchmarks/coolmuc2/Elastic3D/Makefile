# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           None              OMP, TBB                     Shared-memory parallelisation
# DISTRIBUTEDMEM      None              MPI                          Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On                           IPO during compilation (intel only)
# ******************************************************************************************************

PEANO_KERNEL_PEANO_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Peano/peano
PEANO_KERNEL_TARCH_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Peano/mpibalancing
EXAHYPE_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/ExaHyPE
PROJECT_PATH=/home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/Benchmarks/coolmuc2/Elastic3D
EXECUTABLE=ExaHyPE-ElasticWaveEquation3D 

# Several MPI variants face problems with multithreaded MPI. As we run into 
# such issues multiple times, we disable by default multithreaded MPI in ExaHyE. 
# However, feel free to give it a try in your code on your system by disabling 
# this flag. 
PROJECT_CFLAGS+= -DnoMultipleThreadsMayTriggerMPICalls

PROJECT_CFLAGS+= -DDim3

ARCHITECTURE=snb
PROJECT_CFLAGS+= -DALIGNMENT=32

# Default options required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=TBB
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=MPI
endif



# Main Makefile include
-include /home/dominic/dev/codes/c/upstream/ExaHyPE-Engine/ExaHyPE/Makefile


all: 
	@echo ElasticWaveEquation3D 
	@echo =================
	@echo An ExaHyPE solver
