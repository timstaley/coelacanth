#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"
source ${SCRIPT_DIR}/setup_env_vars.sh

mkdir -p ${COELA_INSTALL}/lib
mkdir -p ${COELA_INSTALL}/include


cd ${SCRIPT_DIR}/deps/src/
unzip -o unittest-cpp-1.4.zip
tar -xvf Minuit2-5.28.00.tar.gz
tar -xvf rngstreams-1.0.1.tar.gz
tar -xvf unuran-1.8.1.tar.gz

#Unittest++
cd ${SCRIPT_DIR}/deps/src/UnitTest++ 
make -j$NCPU
cp -v libUnitTest++.a ${COELA_INSTALL}/lib
cp -r src ${COELA_INSTALL}/include/UnitTest++

#RNGstreams
cd ${SCRIPT_DIR}/deps/src/rngstreams-1.0.1
./configure --enable-shared=yes --with-pic --prefix=${COELA_INSTALL}
make install -j$NCPU

#UNURan
cd ${SCRIPT_DIR}/deps/src/unuran-1.8.1
./configure --with-urng-rngstream --enable-shared=yes --prefix=${COELA_INSTALL} --with-pic --with-urng-default=rngstream
make install -j$NCPU

#Minuit2
cd ${SCRIPT_DIR}/deps/src/Minuit2-5.28.00
./configure --with-pic --disable-openmp --prefix=${COELA_INSTALL}
make install -j$NCPU

#TBB
cd ${SCRIPT_DIR}/deps/src
tar -xvf tbb41_20130314oss_src.tgz
cd ${SCRIPT_DIR}/deps/src/tbb41_20130314oss
make -j$NCPU
mv include/tbb ${COELA_INSTALL}/include
mv -v build/linux_*/*.so* ${COELA_INSTALL}/lib

cd ${SCRIPT_DIR}
