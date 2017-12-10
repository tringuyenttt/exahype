exe=ExaHyPE-Euler
spec=convergence/Euler_FV.exahype

# save original file
cp $spec ${spec}_file

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

  for p in 7 11 15 19 # 3 5 7 9
  do 
    rm *.o
    sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$p',' $spec
    cat $spec
    convergence/configure.sh
    make -j24 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${spec}_file $spec
