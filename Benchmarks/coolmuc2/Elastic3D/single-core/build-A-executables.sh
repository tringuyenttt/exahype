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
# pre: FORTRAN: opt-matmul,fma,
# pre: remove -O from the Makefile

for arch in noarch $ARCH
do

sed -i -r 's,archictecture(\s*)const(\s*)=(\s*).+,architecture\1const\2=\3'$arch',' $SPEC
./configure.sh

# O0 no-vec
make clean
export COMPILER_CFLAGS=" -O0 -no-vec "
make -j$cores
mv $APP ExaHyPE-Elastic3d_${arch}_O0_no-vec
# O2 no-vec
make clean
export COMPILER_CFLAGS=" -O2 -no-vec "
make -j$cores
mv $APP ExaHyPE-Elastic3d_${arch}_O2_no-vec

# O3 no-vec
make clean
export COMPILER_CFLAGS=" -O3 -no-vec "
make -j$cores
mv $APP ExaHyPE-Elastic3d_${arch}_O3_no-vec

# O3 
make clean
export COMPILER_CFLAGS=" -O3 "
mv $APP ExaHyPE-Elastic3d_${arch}_O3
make -j$cores

done
