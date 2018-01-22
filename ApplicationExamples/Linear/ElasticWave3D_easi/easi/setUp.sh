#!/bin/bash

module load intel/16.0
module load git
module load cmake
module load netcdf

mkdir libs

cd libs

c_compiler=icc
cpp_compiler=icpc
num_lib_dir="$HOME/software/numa/lib"


cmake_set_compiler="-DCMAKE_CXX_COMPILER=${cpp_compiler} -DCMAKE_C_COMPILER=${c_compiler}"


function clone_or_pull {
    if [ ! -d ${1} ]; then
	git clone --recursive ${2} ${1}
    else
	cd ${1}
	git pull
	cd ..
    fi
}

#clone and build impalaJIT
impala_repo="https://github.com/Manuel1605/ImpalaJIT.git"
impala_build="ImpalaJIT_build"
impala_src="ImpalaJIT"
if [ ! -d ${impala_build} ]; then
    echo "Building Impala"
    clone_or_pull ${impala_src} ${impala_repo}
    mkdir -p ${impala_build}
    cd ${impala_build}
    cmake ${cmake_set_compiler} ../${impala_src}
    make -j8
    cd ..
    echo "IMPALA Build"
fi

#clone and build yaml
yaml_repo="https://github.com/jbeder/yaml-cpp.git"
yaml_build="yaml_build"
yaml_src="yaml-cpp"
if [ ! -d ${yaml_build} ]; then
    echo "Building YAML"
    clone_or_pull ${yaml_src} ${yaml_repo}
    mkdir -p ${yaml_build}
    cd ${yaml_build}
    cmake_flags="${cmake_set_compiler} -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_FLAGS=--std=c++11"
    cmake ${cmake_flags} ../${yaml_src}
    make -j8
    cd ..
    echo "YAML Build"
fi

#clone and build ASAGI
asagi_repo="https://github.com/TUM-I5/ASAGI.git"
asagi_build="asagi_build"
asagi_src="ASAGI"
if [ ! -d ${asagi_build} ]; then
    echo "Building ASAGI"
    clone_or_pull ${asagi_src} ${asagi_repo}
    mkdir -p ${asagi_build}
    cd ${asagi_build}
    cmake_flags="-DCMAKE_BUILD_TYPE=Debug -DNOMPI=ON -DNONUMA=ON ${cmake_set_compiler} -DCMAKE_PREFIX_PATH=${num_lib_dir}"
    cmake ${cmake_flags}  ../${asagi_src}
    make -j8
    cd ..
    echo "ASAGI Build"
fi

#clone easi
echo "Cloning/Pulling easi"
easi_repo="https://github.com/SeisSol/easi.git"
easi_src="easi"
clone_or_pull $easi_src $easi_repo
cd ..

echo "Done."
cd ..