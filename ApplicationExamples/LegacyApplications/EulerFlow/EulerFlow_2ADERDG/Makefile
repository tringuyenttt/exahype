# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#--------------------------------------------------------------------------------------------------
# ARCHITECTURE        CPU               Phi                       Hardware-plattfrom
# COMPILER            Intel             GNU                       Used compiler (and linker)
# MODE                Release           Debug, Profile, Asserts   Verbosity and Debug level
# SHAREDMEM           None              OMP, TBB                  Shared-memory parallelisation
# DISTRIBUTEDMEM      None              MPI                       Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                  Type of boundary conditions
# *********************************************************************************************
PEANO_KERNEL_PEANO_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Peano/mpibalancing
EXAHYPE_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./ExaHyPE
PROJECT_PATH=/home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./Applications/EulerFlow_2ADERDG
EXECUTABLE=ExaHyPE-Euler 


PROJECT_CFLAGS+=-DALIGNMENT=16
# Several MPI variants face problems with multithreaded MPI. As we run into 
# such issues multiple times, we disable by default multithreaded MPI in ExaHyE. 
# However, feel free to give it a try in your code on your system by disabling 
# this flag. 
PROJECT_CFLAGS+=-DnoMultipleThreadsMayTriggerMPICalls
PROJECT_CFLAGS+=-DDim2




-include /home/hpc/pr63so/ga96nuv3/exahype/ExaHyPE/Code/./ExaHyPE/Makefile




all: 
	@echo Euler 
	@echo =================
	@echo An ExaHyPE solver
