export MODE="Release"
export COMPILER="Intel"
export SHAREDMEM="None"
export DISTRIBUTEDMEM="None"

CORES=4

FOLDER=single-core
SPEC=${FOLDER}/Elastic3D-no-output.exahype
APP=ExaHyPE-Elastic3D
#ARCH=knl
ARCH=hsw

# A1
# pre: remove -O from the Makefile
for order in 3 6 9
do
  for arch in noarch $ARCH
  do
   
    sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
    sed -i -r 's,archictecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
    ${FOLDER}/configure-no-output.sh

    # O2 no-vec no-fma
    make clean
    export COMPILER_CFLAGS=" -O2 -no-vec -no-fma -no-ip"
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}-O2-novec--p${order}
   
    # O2 
    make clean
    export COMPILER_CFLAGS=" -O2 -vec -fma -ip "
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}-O2-vec-p${order}

    # O3 
    make clean
    export COMPILER_CFLAGS=" -O3 -vec -fma -ip "
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}-O3-vec-p${order}
    
  done
done 
