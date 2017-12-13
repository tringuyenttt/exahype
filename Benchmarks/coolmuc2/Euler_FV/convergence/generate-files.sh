#!/bin/bash
#
# Perform convergence tests on coolmuc2.
#
# coolmuc2 uses SLURM. SLURM supports array jobs.
#

project=Euler_FV

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
hMax=( 0.11112 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
T=0.03 # SIMULATION END TIME # chosen the same for all patch and mesh sizes

kernels=gengodunov # this is just an identifier; actual kernels must be chosen before building the executables # gengodunov or genmusclhancock

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do 
  for patchSize in 7 11 15 19 # corresponds to orders=3 5 7 9
  do 
    # Create script
    prefix=$project-$kernels
    if [ "$fuseAlgorithmicSteps" == "on" ]; then
      prefix+="-fused"
    else
      prefix+="-nonfused"
    fi
    script=convergence/coolmuc2.slurm-script
    newScript=convergence/coolmuc2-$prefix-p${patchSize}.slurm-script
    cp $script $newScript
   
    sed -i 's,prefix='$project',prefix='$prefix',g' $newScript
    sed -i 's,kernels=gen,kernels='$kernels',g' $newScript
    sed -i 's,p3,p'$patchSize',g' $newScript
    sed -i 's,regular-0,'$mesh',g' $newScript
    sed -i 's,script=convergence/coolmuc2.slurm-script,script='$newScript',g' $newScript
  
    for i in 0 1 2
    do
      mesh=regular-$i
      h=${hMax[i]}
      t=$T
       
      specPrefix=$prefix-$mesh
    
      # Create spec files
      coresPerTask=24
      spec=convergence/$project.exahype
      filename=convergence/$specPrefix-p$patchSize
      newSpec=$filename'.exahype'

      cp $spec $newSpec
      
      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
      sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
     
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
      sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
    
      sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$patchSize',' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.|-|\+|e|E)*),maximum-mesh-size\1=\2'$h',' $newSpec
    done
  done
done
