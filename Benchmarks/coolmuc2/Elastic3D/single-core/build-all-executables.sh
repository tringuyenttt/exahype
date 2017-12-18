export MODE="Release"
export COMPILER="Intel"
export SHAREDMEM="None"
export DISTRIBUTEDMEM="None"
export USE_IPO="off"

CORES=28

FOLDER=single-core
SPEC=${FOLDER}/Elastic3D-no-output.exahype
APP=ExaHyPE-ElasticWaveEquation3D
#ARCH=knl
ARCH=hsw

# A1
# pre: remove -O from the Makefile

declare -a CFLAGS=(" -O2 -no-vec -no-fma -no-ip" " -O2 -vec -fma -ip " " -O3 -vec -fma -ip ")
declare -a suffixes=("O2-novec" "O2-vec" "O3-vec")

for arch in noarch $ARCH
do 
  sed -i -r 's,architecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
  
  for i in 0 1 2
  do  
    export COMPILER_CFLAGS="${CFLAGS[i]}"
    make clean
 
    for order in 3 6 9
    do
      rm *.o cipofiles.mk ffiles.mk cfiles.mk
      sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
      
      ${FOLDER}/configure-no-output.sh
      cat $SPEC

      # O2 no-vec no-fma
      make -j$cores
      mv $APP ExaHyPE-Elastic3d-${arch}-${suffixes[i]}-p${order}

    done   
  done
done 
