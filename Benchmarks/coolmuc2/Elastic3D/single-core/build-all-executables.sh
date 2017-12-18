export MODE="Release"
export COMPILER="Intel"
export SHAREDMEM="None"
export DISTRIBUTEDMEM="None"

SPEC=Elastic3D.exahype
APP=ExaHyPE-Elastic3D
CORES=4
#ARCH=knl
ARCH=hsw

# A1
# pre: remove -O from the Makefile
for order in 3 6 9
do
  for arch in noarch $ARCH
  do
   
    sed -i -r 's,archictecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
    sed -i -r 's,order(\s*)const(\s*)=(\s*).+,order\1const\2=\3'$order',' $SPEC
    ./configure.sh

    # O0 no-vec
    make clean
    export COMPILER_CFLAGS=" -O0 -no-vec "
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}_O0-novec-p${order}
    # O2 no-vec
    make clean
    export COMPILER_CFLAGS=" -O2 -no-vec "
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}-O2-novec-p${order}

    # O3 no-vec
    make clean
    export COMPILER_CFLAGS=" -O3 -no-vec "
    make -j$cores
    mv $APP ExaHyPE-Elastic3d-${arch}-O3-novec-p${order}

    # O3 
    make clean
    export COMPILER_CFLAGS=" -O3 "
    mv $APP ExaHyPE-Elastic3d-${arch}_O3-p${order}
    make -j$cores
    
  done
done 
