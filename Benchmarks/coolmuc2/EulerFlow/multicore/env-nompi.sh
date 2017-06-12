source /etc/profile.d/modules.sh

module unload python
module load python/3.5_intel
module unload java
module load java
module unload intel
module load intel/16.0
module unload mpi.intel
module load mpi.intel/5.1
module unload gcc
module load gcc
module unload tbb
module load tbb

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=None
