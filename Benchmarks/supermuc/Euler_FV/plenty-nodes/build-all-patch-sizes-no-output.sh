directory=plenty-nodes

exe=ExaHyPE-Euler
spec=$directory/Euler_FV-no-output.exahype

# save original file
cp $spec ${spec}_tmp

for m in 1 2
do
  if (( m == 1 )); then
    make clean
    export SHAREDMEM=TBB
  else
    make clean
    export SHAREDMEM=None
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 7 11 15 19
  do 
    rm *.o
    sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$p',' $spec
    cat $spec
    module unload gcc/4.9
    module load python/3.3_anaconda_nompi
    $directory/configure-no-output.sh
    module unload python/3.3_anacobda_nompi
    module load gcc/4.9
    make -j28 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${spec}_tmp $spec
